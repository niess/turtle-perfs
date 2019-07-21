#pragma once


#ifndef INSTALL_PREFIX
#define INSTALL_PREFIX "."
#endif

#define SIMULATOR_CACHE_FILE INSTALL_PREFIX "/share/materials/.cache"
#define SIMULATOR_DEDX_DIR INSTALL_PREFIX "/share/materials/dedx"
#define SIMULATOR_MDF_FILE INSTALL_PREFIX "/share/materials/materials.xml"


struct simulator {
        int (*simulate)(struct simulator * simulator,
            struct geometry * geometry, const double position[3],
            const double direction[3], double * depth, double * flux,
            double * time);
        void (*record)(
            int event, int step, int flag, int medium, double energy,
            double distance, const double position[3]);
};


extern void simulator_initialise(struct simulator * simulator);

extern void simulator_configure(
    struct simulator * simulator, int events, int precompute);
