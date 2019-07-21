#pragma once


#ifndef INSTALL_PREFIX
#define INSTALL_PREFIX "."
#endif

#define GEOMETRY_GEOID INSTALL_PREFIX "/share/topography/egm96.png"

#define GEOMETRY_ZMAX_PDD 2E+03
#define GEOMETRY_ZMAX_TIANSHAN 7.5E+03


enum geometry_medium {
        GEOMETRY_MEDIUM_VOID = -1,
        GEOMETRY_MEDIUM_AIR = 0,
        GEOMETRY_MEDIUM_ROCK = 1
};


enum geometry_location {
        GEOMETRY_LOCATION_PDD = 0,
        GEOMETRY_LOCATION_TIANSHAN
};


struct geometry {
        enum geometry_location location;
        const char * algorithm;
        int (*volume_at)(const double position[3]);
        int (*distance_to)(int volume, const double position[3],
            const double direction[3], double * distance);
        enum geometry_medium (*medium_at)(int volume);
        void (*clear)(void);
};


extern void geometry_initialise_bvh(
    struct geometry * geometry, const char * path, int period);

extern void geometry_configure_bvh(
    struct geometry * geometry, int mode);


extern void geometry_initialise_embree(
    struct geometry * geometry, const char * path, int period);

extern void geometry_configure_embree(
    struct geometry * geometry, int mode);


extern void geometry_initialise_mesh(
    struct geometry * geometry, const char * path, int period);


extern void geometry_initialise_opti(
    struct geometry * geometry, const char * path, int period);

extern void geometry_configure_opti(struct geometry * geometry,
    double resolution, double range, double slope);
