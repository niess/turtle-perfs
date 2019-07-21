/* C89 standard library */
#include <stdlib.h>
/* The clock object */
#include "tictac.h"


static void start(struct tictac * tictac)
{
        gettimeofday(&tictac->base, NULL);
}


static double stop(struct tictac * tictac)
{
        struct timeval t1, * t0 = &tictac->base;
        gettimeofday(&t1, NULL);

        return (double)t1.tv_sec - (double)t0->tv_sec +
                    ((double)t1.tv_usec - (double)t0->tv_usec) * 1E-06;
}


void tictac_initialise(struct tictac * tictac)
{
        tictac->start = &start;
        tictac->stop = &stop;
}
