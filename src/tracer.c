/* C89 standard library */
#include <stdlib.h>
/* Interfaces */
#include "geometry.h"
#include "tictac.h"
#include "tracer.h"


static struct tictac tictac;


static int trace(const struct tracer * tracer, struct geometry * geometry,
    const double position[3], const double direction[3], double * depth,
    double * time)
{
        tictac.start(&tictac);
        int volume = geometry->volume_at(position);

        int n = 0;
        *depth = 0;
        double r[3] = { position[0], position[1], position[2] };
        double s = 0.;
        for (;;) {
                if (tracer->record != NULL)
                        tracer->record(n, geometry->medium_at(volume), s, r);
                if (volume == -1) break;

                double d;
                const int next = geometry->distance_to(
                                     volume, r, direction, &d);
                if (d < 0) break;
                n++;
                if (geometry->medium_at(volume) == GEOMETRY_MEDIUM_ROCK)
                    *depth += d;
                volume = next;

                int i;
                for (i = 0; i < 3; i++)
                        r[i] += d * direction[i];

                if (tracer->record != NULL) s += d;
        }
        *time = tictac.stop(&tictac);

        return n;
}


extern void tracer_initialise(struct tracer * tracer)
{
        tictac_initialise(&tictac);

        tracer->trace = &trace;
        tracer->record = NULL;
}
