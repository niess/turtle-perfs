#pragma once


struct turtle_map;
extern int downsampler_map_load(
    struct turtle_map ** map, const char * path, int period);
