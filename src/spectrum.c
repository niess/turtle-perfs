/* The standard library */
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* The PUMAS engine */
#include "pumas.h"
/* The TURTLE library */
#include "turtle.h"
/* For time measurements */
#include <sys/time.h>
/* For external ray tracers */
#include <dlfcn.h>

/* Prototypes for any external ray tracer library */
void (*raytracer_create)(const char *) = NULL;
void (*raytracer_destroy)(void) = NULL;
void (*raytracer_initialise)(const double *) = NULL;
double (*raytracer_medium)(const double *, const double *, double) = NULL;

/* This is a chrono */
struct tictac {
        struct timeval t0, t1;
};

static void tictac_start(struct tictac * tictac)
{
        gettimeofday(&tictac->t0, NULL);
}

static double tictac_stop(struct tictac * tictac)
{
        gettimeofday(&tictac->t1, NULL);
        return (double)tictac->t1.tv_sec - (double)tictac->t0.tv_sec +
            ((double)tictac->t1.tv_usec - (double)tictac->t0.tv_usec) * 1E-06;
}

/* Global handles */
static struct turtle_map * map = NULL;
static struct turtle_map * geoid = NULL;
static struct turtle_stack * stack = NULL;
static struct turtle_stepper * stepper = NULL;
static struct pumas_context * context = NULL;

static void gracefully_exit(int code)
{
        turtle_stepper_destroy(&stepper);
        turtle_map_destroy(&map);
        turtle_stack_destroy(&stack);
        pumas_recorder_destroy(&context->recorder);
        pumas_context_destroy(&context);

        exit(code);
}

static void handle_pumas(enum pumas_return rc, pumas_function_t * function,
    const char * message)
{
        fprintf(stderr, "A PUMAS library error occured:\n%s\n", message);
        gracefully_exit(EXIT_FAILURE);
}

/* The geometry */
static struct {
        int geometry;
        double z_min;
        double z_max;
        int steps;
} geometry = { 1, -1E+03, 2E+03, 0 };

#define MAX_LAYERS 100
static struct {
        int index;
        struct pumas_medium * medium[MAX_LAYERS];
        double grammage[MAX_LAYERS];
} layer;

static double locals_air(struct pumas_medium * medium,
    struct pumas_state * state, struct pumas_locals * locals_)
{
        locals_->density = 1.205;
        return 0.;
}

static double locals_rock(struct pumas_medium * medium,
    struct pumas_state * state, struct pumas_locals * locals_)
{
        locals_->density = 2.65E+03;
        return 0.;
}

static struct pumas_medium media[2] = {{ 1, &locals_air }, { 0, &locals_rock }};

static double medium(struct pumas_context * context, struct pumas_state * state,
    struct pumas_medium ** medium_p)
{
        if (geometry.geometry == 1) {
                int depth;
                double altitude, ground, step;
                turtle_stepper_step(stepper, state->position, NULL, NULL, NULL,
                    &altitude, &ground, &step, &depth);
                if ((depth < 0) || (altitude < geometry.z_min) ||
                    (altitude > geometry.z_max)) {
                        if (medium_p != NULL) *medium_p = NULL;
                        return 0.;
                } else {
                        if (medium_p != NULL)
                                *medium_p = (altitude > ground) ?
                                    media : media + 1;
                        return step;
                }
        } else if (geometry.geometry == 0) {
                const double a = 6378137. + geometry.z_max;
                const double b = 6356752.3142 + geometry.z_max;
                const double x = state->position[0] / a;
                const double y = state->position[1] / a;
                const double z = state->position[2] / b;
                if (x * x + y * y + z * z > 1.) {
                        if (medium_p != NULL) *medium_p = NULL;
                        return 0.;
                } else {
                        if (medium_p != NULL) *medium_p = media;
                        return geometry.z_max;
                }
        } else {
                if (medium_p != NULL)
                        *medium_p = layer.medium[layer.index];
                if (geometry.geometry == 2)
                        return 0.;

                double step = raytracer_medium(
                    state->position, state->direction, state->distance);
                if (step < 0) step = 0.;
                else if (step < 1E-05) step = 1E-05;
                return step;
        }
}

void record_steps(struct pumas_context * context, struct pumas_state * state,
    struct pumas_medium * medium, enum pumas_event event)
{
        if (!(event & PUMAS_EVENT_START))
                geometry.steps++;
}

struct trex_data {
        int index;
#define MT_PERIOD 624
        unsigned long buffer[MT_PERIOD];
};

/* Set the MT initial state */
static void trex_initialise(struct trex_data * data)
{
        /* Initialise the random seed from /dev/urandom */
        unsigned long seed;
        static const char * urandom = "/dev/urandom";
        FILE * fp = fopen(urandom, "rb");
        fread(&seed, sizeof(long), 1, fp);
        fclose(fp);

        /* Set the MT initial state */
        data->buffer[0] = seed & 0xffffffffUL;
        int j;
        for (j = 1; j < MT_PERIOD; j++) {
                data->buffer[j] = (1812433253UL *
                        (data->buffer[j - 1] ^ (data->buffer[j - 1] >> 30)) +
                    j);
                data->buffer[j] &= 0xffffffffUL;
        }
        data->index = MT_PERIOD;
}

/* Uniform pseudo random distribution from a Mersenne Twister */
static double trex_uniform01(struct pumas_context * context)
{
        /* Check the buffer. */
        struct trex_data * data = (struct trex_data *)context->user_data;
        if (data->index < MT_PERIOD - 1) {
                data->index++;
        } else {
                /* Update the MT state. */
                const int M = 397;
                const unsigned long UPPER_MASK = 0x80000000UL;
                const unsigned long LOWER_MASK = 0x7fffffffUL;
                static unsigned long mag01[2] = { 0x0UL, 0x9908b0dfUL };
                unsigned long y;
                int kk;
                for (kk = 0; kk < MT_PERIOD - M; kk++) {
                        y = (data->buffer[kk] & UPPER_MASK) |
                            (data->buffer[kk + 1] & LOWER_MASK);
                        data->buffer[kk] =
                            data->buffer[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
                }
                for (; kk < MT_PERIOD - 1; kk++) {
                        y = (data->buffer[kk] & UPPER_MASK) |
                            (data->buffer[kk + 1] & LOWER_MASK);
                        data->buffer[kk] = data->buffer[kk + (M - MT_PERIOD)] ^
                            (y >> 1) ^ mag01[y & 0x1UL];
                }
                y = (data->buffer[MT_PERIOD - 1] & UPPER_MASK) |
                    (data->buffer[0] & LOWER_MASK);
                data->buffer[MT_PERIOD - 1] =
                    data->buffer[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
                data->index = 0;
        }

        /* Tempering. */
        unsigned long y = data->buffer[data->index];
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);

        /* Convert to a floating point and return. */
        return y * (1.0 / 4294967295.0);
}

int main(int argc, char * argv[])
{
        /* Parse the arguments */
        const char * site = "CDC", * mode = "DETAILED";
        int events = 1000;
        if (argc && --argc) site = *++argv;
        if (argc && --argc) mode = *++argv;
        if (argc && --argc) events = atoi(*++argv); 
        if (argc && --argc) geometry.geometry = atoi(*++argv);

        /* Load any external ray tracer */
        void * lib = NULL;
        if (geometry.geometry > 2) {
                const char * path = (geometry.geometry == 3) ?
                        "lib/libraytracer-tree.so" :
                        "lib/libraytracer-mesh.so";
                lib = dlopen(path, RTLD_NOW);
                if (!lib) {
                        fprintf(stderr, "%s\n", dlerror());
                        exit(EXIT_FAILURE);
                }
                dlerror();

                raytracer_create = dlsym(lib, "raytracer_create");
                raytracer_destroy = dlsym(lib, "raytracer_destroy");
                raytracer_initialise = dlsym(lib, "raytracer_initialise");
                raytracer_medium = dlsym(lib, "raytracer_medium");

                if (strcmp(site, "ULS") == 0) {
                        if (geometry.geometry == 3)
                                path = "share/data/ulastai-2.ply";
                        else
                                path = "share/data/ulastai-1.png";
                } else {
                        if (geometry.geometry == 3)
                                path = "share/data/pdd-1.ply";
                        else
                                path = "share/data/pdd-1.png";
                }
                raytracer_create(path);
        }

        /* Set the site */
        double latitude, longitude, height, azimuth0, daz, del;
        int nx, ny;
        if (strcmp(site, "CDC") == 0) {
                turtle_map_load(&map, "share/data/pdd-2.5m.png");
                latitude = 45.76415653;
                longitude = 2.95536402;
                height = 0.5;
                geometry.z_max = 2E+03;
                azimuth0 = -4.;
                daz = 0.1;
                del = 0.1;
                nx = 601;
                ny = 301;
        } else if (strcmp(site, "TDF") == 0) {
                turtle_map_load(&map, "share/data/pdd-2.5m.png");
                latitude = 45.769961;
                longitude = 2.980687;
                height = 0.5;
                geometry.z_max = 2E+03;
                azimuth0 = 250.;
                daz = 0.1;
                del = 0.1;
                nx = 601;
                ny = 301;
        } else if (strcmp(site, "ULS") == 0) {
                turtle_stack_create(&stack, "share/SRTMGL1", 0, NULL, NULL);
                latitude = 42.9242110;
                longitude = 86.6982730;
                height = 100.;
                geometry.z_max = 7.0E+03;
                azimuth0 = 0.;
                daz = 0.2;
                del = 0.05;
                nx = 1801;
                ny = 241;
        } else {
                fprintf(stderr, "invalid site\n");
                exit(EXIT_FAILURE);
        }

        /* Initialise the stepper */
        turtle_stepper_create(&stepper);
        turtle_map_load(&geoid, "share/data/egm96.png");
        turtle_stepper_geoid_set(stepper, geoid);
        turtle_stepper_add_flat(stepper, 0.);
        if (stack != NULL)
                turtle_stepper_add_stack(stepper, stack, 0.);
        if (map != NULL)
                turtle_stepper_add_map(stepper, map, 0.);

        /* Initialise PUMAS */
        pumas_error_handler_set(&handle_pumas);
        const char * dump = ".materials.dump";
        FILE * fid = fopen(dump, "r");
        if (fid != NULL) {
                pumas_load(fid);
                fclose(fid);
        } else {
                pumas_initialise(PUMAS_PARTICLE_MUON,
                    "share/materials/materials.xml", "share/materials/dedx");
                fid = fopen(dump, "w+");
                pumas_dump(fid);
                fclose(fid);
        }

        pumas_context_create(&context, sizeof(struct trex_data));
        context->forward = 0;
        context->medium = &medium;
        context->random = &trex_uniform01;
        context->event |= PUMAS_EVENT_LIMIT_KINETIC | PUMAS_EVENT_MEDIUM;
        if (strcmp(mode, "DETAILED") == 0) {
                context->longitudinal = 0;
                context->scheme = PUMAS_SCHEME_DETAILED;
        } else if (strcmp(mode, "HYBRID") == 0) {
                context->longitudinal = 1;
                context->scheme = PUMAS_SCHEME_HYBRID;
        } else {
                fprintf(stderr, "invalid mode: %s\n", mode);
                gracefully_exit(EXIT_FAILURE);
        }
        context->kinetic_limit = 1E+08;

        pumas_material_index("Rock", &media[0].material);
        pumas_material_index("Air", &media[1].material);

        /* Set the recorder */
        pumas_recorder_create(0, &context->recorder);
        context->recorder->period = 1;
        context->recorder->record = &record_steps;

        /* Initialise the random engine */
        trex_initialise((struct trex_data *)context->user_data);

        /* Scan the line of sights */
        char filename[128];
        sprintf(filename, "share/results/%s_%.3lE_%.3lE_%.3lE_%s_%d.txt",
            site, 1E-02, 1., 0.4, mode, geometry.geometry);

        fid = fopen(filename, "w+");
        fclose(fid);

        struct tictac tictac;
        const double dle = log(1E+07);
        fputs("scanning depths ...", stdout);
        fflush(stdout);
        double total_time = 0.;
        long long total_iterations = 0;
        for (int iy = 0; iy < ny; iy++) {
                const double elevation = del * iy;
                for (int ix = 0; ix < nx; ix++) {
                        const double azimuth = daz * ix + azimuth0;

                        enum pumas_event default_event = context->event;
                        if (geometry.geometry >= 2) {
                                const int geogeo = geometry.geometry;
                                enum pumas_scheme scheme = context->scheme;
                                context->scheme = PUMAS_SCHEME_NO_LOSS;

                                struct pumas_state state = {
                                    .charge = -1., .kinetic = 1E+06,
                                    .weight = 1. };
                                turtle_stepper_position(stepper, latitude,
                                    longitude, height, 0, state.position, NULL);
                                turtle_ecef_from_horizontal(latitude, longitude,
                                    azimuth + 180, -elevation, state.direction);

                                geometry.geometry = 1;
                                layer.index = 0;
                                enum pumas_event event = PUMAS_EVENT_NONE;
                                for (;;) {
                                        if (layer.index >= MAX_LAYERS) {
                                                fprintf(stderr,
                                                    "max layers exceeded\n");
                                                gracefully_exit(EXIT_FAILURE);
                                        }
                                        struct pumas_medium * med[2];
                                        pumas_transport(context, &state,
                                            &event, med);
                                        layer.medium[layer.index] = med[0];
                                        layer.grammage[layer.index++] =
                                            state.grammage;
                                        if (med[1] == NULL)
                                                break;
                                }
                                layer.medium[layer.index++] = NULL;

                                geometry.geometry = geogeo;
                                context->event |= PUMAS_EVENT_LIMIT_GRAMMAGE;
                                context->scheme = scheme;
                        }

                        tictac_start(&tictac);
                        geometry.steps = 0;
                        double rock_depth = 0.;
                        for (int j = 0; j < events; j++) {
                                const double energy =
                                    1E-03 * exp(trex_uniform01(context) * dle);
                                struct pumas_state state = {
                                    .charge = -1., .kinetic = energy,
                                    .weight = energy * dle };
                                turtle_stepper_position(stepper, latitude,
                                    longitude, height, 0, state.position, NULL);
                                turtle_ecef_from_horizontal(latitude, longitude,
                                    azimuth + 180, -elevation, state.direction);

                                if (geometry.geometry > 2)
                                        raytracer_initialise(state.position);
                                layer.index = 0;
                                enum pumas_event event = PUMAS_EVENT_NONE;
                                while (event != PUMAS_EVENT_LIMIT_KINETIC) {
                                        if (geometry.geometry >= 2) {
                                                context->grammage_max =
                                                    layer.grammage[layer.index];
                                        }
                                        struct pumas_medium * med[2];
                                        const double d0 = state.distance;
                                        pumas_transport(context, &state,
                                            &event, med);
                                        if (med[0] == media + 1) {
                                                rock_depth +=
                                                    state.distance - d0;
                                        }
                                        if (med[1] == NULL)
                                                break;
                                        if (geometry.geometry >= 2) {
                                                if (layer.medium[++layer.index]
                                                    == NULL)
                                                        break;
                                        }
                                }
                        }

                        if (geometry.geometry >= 2) {
                                context->event = default_event;
                        }

                        const double dt = tictac_stop(&tictac) / events;
                        int steps = geometry.steps / events;
                        total_time += dt;
                        total_iterations += steps;

                        fid = fopen(filename, "a");
                        fprintf(fid, "%d %d %.9lE %.9lE %d\n",
                            iy, ix, rock_depth / events, dt, steps);
                        fclose(fid);
                }
                printf("\rscanning depths ... %d", iy + 1);
                fflush(stdout);
        }
        puts("\rall depths scanned!       ");
        printf("time       = %.5lE s\n", total_time);
        printf("iterations = %.5lld\n", total_iterations);

        /* Unload any ray tracer */
        if (geometry.geometry > 2) {
                raytracer_destroy();
                dlclose(lib);
        }

        /* Exit to the OS */
        exit(EXIT_SUCCESS);
}
