/* C89 library */
#include <float.h>
#include <stdlib.h>
#include <string.h>
/* Interfaces */
#include "downsampler.h"
#include "geometry.h"
#include "turtle.h"


static struct turtle_map * geoid = NULL;
static struct turtle_map * map = NULL;
static struct turtle_stepper * stepper = NULL;
static double zmax = 0.;


static int opti_volume_at(const double position[3])
{
        double r[3] = { position[0], position[1], position[2] };
        double altitude, ground;
        int layer[2];
        double latitude, longitude;
        turtle_stepper_step(stepper, r, NULL, &latitude, &longitude,
            &altitude, &ground, NULL, layer);

        if ((layer[0] < 0) || (altitude <= 0.) || (altitude > zmax))
                return GEOMETRY_MEDIUM_VOID;
        else return (altitude > ground) ? GEOMETRY_MEDIUM_AIR :
                                          GEOMETRY_MEDIUM_ROCK;
}


static int opti_distance_to(int volume, const double position[3],
    const double direction[3], double * distance)
{
        double r[3] = { position[0], position[1], position[2] };
        double altitude, ground;
        int layer;
        turtle_stepper_step(stepper, r, direction, NULL, NULL, &altitude,
            &ground, distance, &layer);

        if ((altitude <= 0.) || (altitude > zmax)) {
                *distance = -1;
                return GEOMETRY_MEDIUM_VOID;
        } else {
                if ((*distance >= 0.) && (*distance < FLT_EPSILON))
                        *distance = FLT_EPSILON;
                if (layer < 0) {
                        return -1;
                } else {
                        return (altitude >= ground) ? GEOMETRY_MEDIUM_AIR :
                                                      GEOMETRY_MEDIUM_ROCK;
                }
        }
}


static enum geometry_medium opti_medium_at(int volume)
{
        return (enum geometry_medium)volume;
}


static void opti_clear(void)
{
        turtle_stepper_destroy(&stepper);
        turtle_map_destroy(&map);
        turtle_map_destroy(&geoid);
}


extern void geometry_initialise_opti(
    struct geometry * geometry, const char * path, int period)
{
        if (stepper == NULL) {
                /* Create the stepper */
                downsampler_map_load(&map, path, period);
                turtle_map_load(&geoid, GEOMETRY_GEOID);

                turtle_stepper_create(&stepper);
                turtle_stepper_geoid_set(stepper, geoid);
                turtle_stepper_add_map(stepper, map, 0.);

                if (strstr(path, "tianshan") != NULL)
                        zmax = GEOMETRY_ZMAX_TIANSHAN;
                else
                        zmax = GEOMETRY_ZMAX_PDD;
        }

        /* Initialise the corresponding geometry interface */
        geometry->algorithm = "opti";
        if (strstr(path, "tianshan") != NULL)
                geometry->location = GEOMETRY_LOCATION_TIANSHAN;
        else
                geometry->location = GEOMETRY_LOCATION_PDD;
        geometry->volume_at = &opti_volume_at;
        geometry->distance_to = &opti_distance_to;
        geometry->medium_at = &opti_medium_at;
        geometry->clear = &opti_clear;
}


extern void geometry_configure_opti(struct geometry * geometry,
    double resolution, double range, double slope)
{
        turtle_stepper_resolution_set(stepper, resolution);
        turtle_stepper_range_set(stepper, range);
        turtle_stepper_slope_set(stepper, slope);
}
