/* C89 standard library */
#include <stdlib.h>
/* Interfaces */
#include "downsampler.h"
#include "turtle.h"


extern int downsampler_map_load(
    struct turtle_map ** map, const char * path, int period)
{
        if (period <= 1) {
                return turtle_map_load(map, path);
        }

        /* Load the full map */
        struct turtle_map * full;
        struct turtle_map_info info;
        const char * projection;

        turtle_map_load(&full, path);
        turtle_map_meta(full, &info, &projection);

        /* Create the downsampled one */
        int nx = info.nx, ny = info.ny;
        info.encoding = NULL;
        const double dx = (info.x[1] - info.x[0]) / (info.nx - 1);
        info.nx /= period;
        info.x[1] = info.x[0] + (info.nx - 1) * period * dx;
        const double dy = (info.y[1] - info.y[0]) / (info.ny - 1);
        info.ny /= period;
        info.y[1] = info.y[0] + (info.ny - 1) * period * dy;

        turtle_map_create(map, &info, projection);
        int ix, iy;
        for (iy = 0; iy < info.ny; iy++) for (ix = 0; ix < info.nx; ix++) {
                double z;
                turtle_map_node(full, ix * period, iy * period, NULL, NULL, &z);
                turtle_map_fill(*map, ix, iy, z);
        }

        turtle_map_destroy(&full);
        return EXIT_SUCCESS;
}
