/* C89 includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* The TURTLE library */
#include "turtle.h"
/* The view object */
#include "view.h"


static void get_direction(const struct view * view, double azimuth,
    double elevation, double direction[3])
{
        turtle_ecef_from_horizontal(view->latitude, view->longitude,
            azimuth, elevation, direction);
}


extern void view_initialise(struct view * view, const char * site)
{
        if (strcmp(site, "CDC") == 0) {
                view->site = "CDC";
                view->LOS = 601 * 301;
                view->position[0] = 4451978.8500717236;
                view->position[1] = 229840.64926429678;
                view->position[2] = 4547788.6190320663;
        } else if (strcmp(site, "ULS") == 0) {
                view->site = "ULS";
                view->LOS = 1801 * 241;
                view->position[0] = 269523.22314760089 ;
                view->position[1] = 4671932.518579076;
                view->position[2] = 4323208.8337448305;
        } else {
                fprintf(stderr, "invalid site name: %s\n", site);
                exit(EXIT_FAILURE);
        }

        double altitude;
        turtle_ecef_to_geodetic(view->position, &view->latitude,
            &view->longitude, &altitude);

        view->direction = &get_direction;
}
