#pragma once


struct tracer {
        int (*trace)(const struct tracer * trace, struct geometry * geometry,
            const double position[3], const double direction[3],
            double * depth, double * time);
        void (*record)(
            int iter, int medium, double step, const double position[3]);
};


extern void tracer_initialise(struct tracer * tracer);
