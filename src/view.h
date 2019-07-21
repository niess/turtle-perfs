#pragma once


struct view {
        const char * site;
        int LOS;
        double position[3];
        double latitude;
        double longitude;
        void (*direction)(const struct view * view, double azimuth,
            double elevation, double direction[3]);
};


extern void view_initialise(struct view * view, const char * site);
