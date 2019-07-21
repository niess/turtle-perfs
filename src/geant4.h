#pragma once


enum geant4_particle {
        GEANT4_PARTICLE_GEANTINO = 0,
        GEANT4_PARTICLE_MUON_MINUS,
        GEANT4_PARTICLE_MUON_PLUS
};


enum geant4_mode {
        GEANT4_MODE_G4TURTLE = 0,
        GEANT4_MODE_TESSELLATE,
        GEANT4_MODE_STL
};


struct geant4 {
        void (*clear)(void);
        int (*run)(enum geant4_particle, double energy, double position[3],
                    double direction[3], double * depth, double * time);
};


void geant4_initialise(struct geant4 * geant4, const char * topography,
    int period, enum geant4_mode mode, int kill_secondaries, double cut);

void geant4_configure(struct geant4 * geant4, int events, int verbosity);
