/* The standard library */
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* The PUMAS engine */
#include "pumas.h"
/* The interfaces */
#include "geometry.h"
#include "simulator.h"
#include "tictac.h"


#define WGS84_RADIUS_A 6378137.0
#define WGS84_RADIUS_B 6356752.314245


static struct pumas_context * context = NULL;
static struct geometry * geometry;
static struct simulator * simulator;
static int events = 100;
static int precompute = 0;
static int steps = 0;
static double zmax = 0.;


static void gracefully_exit(int code)
{
        pumas_context_destroy(&context);
        exit(code);
}


static void handle_pumas(enum pumas_return rc, pumas_function_t * function,
    const char * message)
{
        fprintf(stderr, "A PUMAS library error occured:\n%s\n", message);
        gracefully_exit(EXIT_FAILURE);
}


/* The pre-computed geometry */
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


static struct {
        int volume;
        int next;
        double position[3];
        double direction[3];
        double step;
        int event;
        int first;
        int start_step;
        double start_distance;
        int outside;
} stepping = { 0 };

static void stepping_initialise(int event, const struct pumas_state * state)
{
        int i;

        stepping.volume = geometry->volume_at(state->position);
        stepping.next = stepping.volume;
        stepping.step = 0.;
        stepping.event = event;
        stepping.first = 1;
        stepping.start_step = steps;
        stepping.start_distance = state->distance;
        memcpy(stepping.position, state->position, sizeof stepping.position);
        for (i = 0; i < 3; i++) {
                stepping.direction[i] = -state->direction[i];
        }
        stepping.outside = 0;
}


static struct pumas_medium * volume_to_medium(int volume)
{
        if (volume >= 0) {
                return media + geometry->medium_at(volume);
        } else {
                return NULL;
        }
}


static void medium(struct pumas_context * context, struct pumas_state * state,
    struct pumas_medium ** medium_p, double * step_p)
{
        if ((medium_p == NULL) && (step_p == NULL)) /* This case is not used */
                return;

        if (precompute) {
                if (medium_p != NULL)
                        *medium_p = layer.medium[layer.index];
                if (step_p != NULL)
                        *step_p = 0.;
        } else if (stepping.outside) {
                /* The particle has escaped the DEM. No further intersection
                 * with the topography is expected. Let us transport it up to
                 * the primary altitude
                 */
                if (medium_p != NULL) *medium_p = media;
                if (step_p != NULL) *step_p = 0.;
        } else {
                /* The particle is within the topography described by
                 * the DEM. Let us compute its local coordinate w.r.t. the
                 * last computed distance
                 */
                int i;
                double step = 0.;
                for (i = 0; i < 3; i++) {
                        step += (state->position[i] - stepping.position[i]) *
                            stepping.direction[i];
                }

                if (step_p == NULL) {
                        /* This is a volume request subsequent to a distance
                         * computation. Let us return the corresponding volume.
                         */
                        int volume = (step >= stepping.step) ? stepping.next :
                                                               stepping.volume;
                        *medium_p = volume_to_medium(volume);
                        return;
                } else if (step >= stepping.step) {
                        /* The last geometric step was accepted. Let us update
                         * the volume accordingly.
                         */
                        stepping.volume = stepping.next;
                }

                /* Fetch the next medium along the line of sight */
                const double direction[3] = {
                    -state->direction[0], -state->direction[1],
                    -state->direction[2] };
                stepping.next = geometry->distance_to(stepping.volume,
                    state->position, direction, &step);

                if (step < 0) {
                        /* Void case */
                        if (medium_p != NULL) *medium_p = NULL;
                        *step_p =  0.;
                        return;
                } else {
                        if (medium_p != NULL)
                                *medium_p = volume_to_medium(stepping.volume);
                        *step_p = (step < 1E-06) ? 1E-06 : step;
                }

                /* Update the geometry stepping data */
                stepping.step = step;
                memcpy(stepping.position, state->position,
                    sizeof stepping.position);
                memcpy(stepping.direction, direction,
                    sizeof stepping.direction);

        }
}


void record_steps(struct pumas_context * context, struct pumas_state * state,
    struct pumas_medium * medium, enum pumas_event event)
{
        if (simulator->record != NULL) {
                int step_index;
                if (stepping.first) {
                        stepping.first = 0;
                        step_index = 0;
                } else {
                        step_index = steps - stepping.start_step + 1;
                }

                int medium_index = (medium == NULL) ? -1 : medium - media;
                simulator->record(stepping.event, step_index, event, medium_index,
                    state->kinetic, state->distance - stepping.start_distance,
                    state->position);
        }

        if (!(event & PUMAS_EVENT_START))
                steps++;
}


struct prng_data {
        int index;
#define MT_PERIOD 624
        unsigned long buffer[MT_PERIOD];
};


/* Set the MT initial state */
static void prng_initialise(struct prng_data * data)
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
static double prng_uniform01(struct pumas_context * context)
{
        /* Check the buffer. */
        struct prng_data * data = (struct prng_data *)context->user_data;
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


static int primary_set(struct pumas_state * state)
{
        /* Limit the state propagation distance to
         * the intersection with primary altitude
         */

        const double * const r = state->position;
        const double * const u = state->direction;
        const double ai = 1. / (WGS84_RADIUS_A + zmax);
        const double bi = 1. / (WGS84_RADIUS_B + zmax);

        const double rx = ai * r[0];
        const double ry = ai * r[1];
        const double rz = bi * r[2];
        const double r2 = rx * rx + ry * ry + rz * rz;
        if (r2 >= 1) return 0;

        const double ux = -ai * u[0];
        const double uy = -ai * u[1];
        const double uz = -bi * u[2];
        const double a = ux * ux + uy * uy + uz * uz;
        const double b = rx * ux + ry * uy + rz * uz;

        const double d2 = b * b + a * (1. - r2);
        const double d = (d2 <= 0.) ? 0. : sqrt(d2);
        double step = (d - b) / a;
        if (step < 1E-03) step = 1E-03;

        context->event |= PUMAS_EVENT_LIMIT_DISTANCE;
        context->distance_max = state->distance + step;
        stepping.outside = 1;

        return 1;
}

static void primary_unset(void)
{
        if (context->event & PUMAS_EVENT_LIMIT_DISTANCE)
                context->event -= PUMAS_EVENT_LIMIT_DISTANCE;
        context->distance_max = 0.;
        stepping.outside = 0;
}


static int simulator_simulate(struct simulator * simulator_,
    struct geometry * geometry_, const double position[3],
    const double direction[3], double * rock_depth_p, double * flux_p,
    double * time_p)
{
        geometry = geometry_;
        simulator = simulator_;

        zmax = (geometry->location == GEOMETRY_LOCATION_PDD) ?
            GEOMETRY_ZMAX_PDD : GEOMETRY_ZMAX_TIANSHAN;

        enum pumas_event default_event = context->event;
        if (precompute) {
                precompute = 0;
                context->scheme = PUMAS_SCHEME_NO_LOSS;

                struct pumas_state state = {
                    .charge = -1.,
                    .kinetic = 1E+06,
                    .weight = 1.,
                    .position = { position[0], position[1], position[2] },
                    .direction = { -direction[0], -direction[1], -direction[2] }
                };

                stepping_initialise(-1, &state);
                layer.index = 0;
                enum pumas_event event = PUMAS_EVENT_NONE;
                for (;;) {
                        if (layer.index >= MAX_LAYERS) {
                                fprintf(stderr, "max layers exceeded\n");
                                gracefully_exit(EXIT_FAILURE);
                        }
                        struct pumas_medium * med[2];
                        pumas_transport(context, &state, &event, med);
                        layer.medium[layer.index] = med[0];
                        layer.grammage[layer.index++] = state.grammage;
                        if ((med[1] == NULL) && !primary_set(&state))
                                break;
                        else if (event & PUMAS_EVENT_LIMIT_DISTANCE)
                                break;
                }
                primary_unset();
                layer.medium[layer.index++] = NULL;

                context->event |= PUMAS_EVENT_LIMIT_GRAMMAGE;
                context->scheme = PUMAS_SCHEME_DETAILED;
                precompute = 1;
        }

        struct tictac tictac;
        tictac_initialise(&tictac);
        tictac.start(&tictac);

        int i;
        const double dle = log(1E+07);
        steps = 0;
        double rock_depth = 0., flux = 0.;

        for (i = 0; i < events; i++) {
                const double energy =
                    1E-03 * exp(prng_uniform01(context) * dle);
                struct pumas_state state = {
                    .charge = -1.,
                    .kinetic = energy,
                    .weight = energy * dle,
                    .position = { position[0], position[1], position[2] },
                    .direction = { -direction[0], -direction[1], -direction[2] }
                };

                stepping_initialise(i, &state);
                layer.index = 0;
                enum pumas_event event = PUMAS_EVENT_NONE;
                while (event != PUMAS_EVENT_LIMIT_KINETIC) {
                        if (precompute) {
                                context->grammage_max =
                                    layer.grammage[layer.index];
                        }
                        struct pumas_medium * med[2];
                        const double d0 = state.distance;
                        pumas_transport(context, &state, &event, med);
                        if (med[0] == media + 1) {
                                rock_depth += state.distance - d0;
                        }
                        if ((med[1] == NULL) && !primary_set(&state))
                                break;
                        else if (event & PUMAS_EVENT_LIMIT_DISTANCE)
                                break;
                        if ((precompute) &&
                            (layer.medium[++layer.index] == NULL)) break;
                }
                primary_unset();

                if (event != PUMAS_EVENT_LIMIT_KINETIC)
                        flux += state.weight * pow(state.kinetic, -2.7);
        }
        context->event = default_event;

        *time_p = tictac.stop(&tictac) / events;
        *rock_depth_p = rock_depth / events;
        *flux_p = flux / events;
        return steps / events;
}


extern void simulator_initialise(struct simulator * simulator)
{
        pumas_error_handler_set(&handle_pumas);
        const char * dump = SIMULATOR_CACHE_FILE;
        FILE * fid = fopen(dump, "r");
        if (fid != NULL) {
                pumas_load(fid);
                fclose(fid);
        } else {
                pumas_initialise(PUMAS_PARTICLE_MUON,
                    SIMULATOR_MDF_FILE, SIMULATOR_DEDX_DIR);
                fid = fopen(dump, "w+");
                pumas_dump(fid);
                fclose(fid);
        }

        pumas_context_create(&context, sizeof(struct prng_data));
        context->forward = 0;
        context->medium = &medium;
        context->random = &prng_uniform01;
        context->event |= PUMAS_EVENT_LIMIT_KINETIC | PUMAS_EVENT_MEDIUM;
        context->longitudinal = 0;
        context->scheme = PUMAS_SCHEME_DETAILED;
        context->kinetic_limit = 1E+08;

        pumas_material_index("Rock", &media[0].material);
        pumas_material_index("Air", &media[1].material);

        /* Set the recorder */
        pumas_recorder_create(&context->recorder, 0);
        context->recorder->period = 1;
        context->recorder->record = &record_steps;

        /* Initialise the random engine */
        prng_initialise((struct prng_data *)context->user_data);

        /* Initialise the simulator */
        simulator->simulate = &simulator_simulate;
        simulator->record = NULL;
}


extern void simulator_configure(
    struct simulator * simulator, int events_, int precompute_)
{
        if (events_ > 0) events = events_;
        if (precompute_ >= 0) precompute = precompute_;
}
