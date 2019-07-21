/* C89 standard library */
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* For the ray tracing */
#include "embree3/rtcore.h"
/* The TURTLE library */
#include "turtle.h"
/* The interfaces */
#include "downsampler.h"
#include "geometry.h"
#include "tictac.h"


struct Vertex { float x, y, z; };
struct Triangle { unsigned int v0, v1, v2; };


static struct turtle_map * map = NULL;
static struct turtle_map * geoid = NULL;
static const struct turtle_projection * projection = NULL;

static RTCDevice device = NULL;
static RTCScene scene = NULL;
static RTCGeometry topography = NULL;
static unsigned int topographyId = 0;
static struct Vertex origin = { 0 };


static void handleError(void * userPtr, enum RTCError code, const char * str)
{
        fprintf(stderr, "error: %s\n", str);
        exit(EXIT_FAILURE);
}


static void get_node(int ix, int iy, struct Vertex * node)
{
        double x, y, z;
        turtle_map_node(map, ix, iy, &x, &y, &z);

        if (projection != NULL) {
                turtle_projection_unproject(
                    projection, x, y, &y, &x);
        }

        double undulation;
        turtle_map_elevation(geoid, x, y, &undulation, NULL);
        z += undulation;

        double ecef[3];
        turtle_ecef_from_geodetic(y, x, z, ecef);
        node->x = ecef[0];
        node->y = ecef[1];
        node->z = ecef[2];
}


/* Create the scene from a topography map */
static void embree_create(const char * path, int period_)
{
        /* Configure the device */
        device = rtcNewDevice("threads=1");
        rtcSetDeviceErrorFunction(device, &handleError, NULL);

        /* Load the map data */
        downsampler_map_load(&map, path, period_);
        projection = turtle_map_projection(map);
        turtle_map_load(&geoid, GEOMETRY_GEOID);

        struct turtle_map_info info;
        turtle_map_meta(map, &info, NULL);
        const int nx = info.nx, ny = info.ny;

        /* Set the origin */
        get_node(nx / 2, ny / 2, &origin);

        /* Build the nodes */
        const int n_nodes = nx * ny;
        int period = n_nodes / 100;
        if (period <= 0) period = 1;
        struct tictac tictac;
        tictac_initialise(&tictac);
        printf("# building nodes ...");
        tictac.start(&tictac);

        /* Instanciate the topography */
        topography = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
        struct Vertex * nodes = rtcSetNewGeometryBuffer(
            topography, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
            sizeof(struct Vertex), n_nodes);
        struct Vertex * node;
        int ix, iy, it;
        for (iy = 0, it = 0, node = nodes; iy < ny; iy++)
            for (ix = 0; ix < nx; ix++, it++, node++) {
                get_node(ix, iy, node);
                node->x -= origin.x;
                node->y -= origin.y;
                node->z -= origin.z;

                if ((it + 1) % period == 0) {
                        printf("\r# building nodes ... %d %%",
                            (it + 1) / period);
                        fflush(stdout);
                }
        }

        double dt = tictac.stop(&tictac);
        printf("\r# %d nodes loaded in %.3f s\n", n_nodes, dt);
        fflush(stdout);

        /* Build the facets */
        printf("# building facets ...");
        tictac.start(&tictac);

        const int n_iter = (nx - 1) * (ny - 1);
        period = n_iter / 100;
        if (period <= 0) period = 1;

        struct Triangle * triangles = rtcSetNewGeometryBuffer(
            topography, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
            sizeof(struct Triangle), (long long int)n_iter * 2);
        struct Triangle * triangle;
        for (iy = 0, it = 0, triangle = triangles; iy < ny - 1; iy++)
            for (ix = 0; ix < nx - 1; ix++, it++, triangle++) {
                triangle->v0 = iy * nx + ix;
                triangle->v1 = iy * nx + ix + 1;
                triangle->v2 = (iy + 1) * nx + ix + 1;
                triangle++;
                triangle->v0 = (iy + 1) * nx + ix + 1;
                triangle->v1 = (iy + 1) * nx + ix;
                triangle->v2 = iy * nx + ix;

                if ((it + 1) % period == 0) {
                        printf("\r# building facets ... %d %%",
                            (it + 1) / period);
                        fflush(stdout);
                }
        }
        dt = tictac.stop(&tictac);
        printf("\r# %d facets built in %.3f s\n", 2 * n_iter, dt);

        /* Free the maps */
        turtle_map_destroy(&map);
        turtle_map_destroy(&geoid);
        projection = NULL;

        /* Finalise the geometry */
        tictac.start(&tictac);
        rtcCommitGeometry(topography);

        /* Set the scene */
        scene = rtcNewScene(device);
        rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);
        topographyId = rtcAttachGeometry(scene, topography);
        rtcReleaseGeometry(topography);
        rtcCommitScene(scene);

        dt = tictac.stop(&tictac);
        printf("\r# geometry commited in %.3f s\n", dt);
}


static void embree_clear(void)
{
        /* Clean and exit */
        rtcReleaseScene(scene);
        rtcReleaseDevice(device);
}


static int embree_distance_to(int volume, const double position[3],
    const double direction[3], double * distance)
{
        /* Intersect the topography */
        struct RTCIntersectContext context = {
            .flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT };
        struct RTCRayHit rayHit = {
            .ray = {
                .org_x = position[0] - origin.x,
                .org_y = position[1] - origin.y,
                .org_z = position[2] - origin.z,
                .dir_x = direction[0],
                .dir_y = direction[1],
                .dir_z = direction[2],
                .tnear = 0.0f,
                .tfar = FLT_MAX },
            .hit = {
                .instID = RTC_INVALID_GEOMETRY_ID,
                .geomID = RTC_INVALID_GEOMETRY_ID }};

        rtcInitIntersectContext(&context);
        rtcIntersect1(scene, &context, &rayHit);

        if (rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
                if (distance != NULL) *distance = -1.;
                return GEOMETRY_MEDIUM_VOID;
        } else {
                if (distance != NULL)
                    *distance = rayHit.ray.tfar + 1E-02; /* Protect against rounding */
                return (volume == GEOMETRY_MEDIUM_ROCK) ? GEOMETRY_MEDIUM_AIR :
                                                          GEOMETRY_MEDIUM_ROCK;
        }
}


static int embree_volume_at(const double position_[3])
{
        /* Translate to the local frame */
        double position[3] = {
            position_[0] - origin.x,
            position_[1] - origin.y,
            position_[2] - origin.z };

        /* 1st check the bounding box of the tree */
        struct RTCBounds box;
        rtcGetSceneBounds(scene, &box);
        if ((position[0] < box.lower_x) || (position[0] > box.upper_x) ||
            (position[1] < box.lower_y) || (position[1] > box.upper_y) ||
            (position[2] < box.lower_z) || (position[2] > box.upper_z)) {
                return GEOMETRY_MEDIUM_AIR;
        }

        /* Draw a ray along the local vertical */
        double latitude, longitude, altitude, vertical[3];
        turtle_ecef_to_geodetic(position_, &latitude, &longitude, &altitude);
        turtle_ecef_from_horizontal(latitude, longitude, 0., 90., vertical);

        double d;
        embree_distance_to(GEOMETRY_MEDIUM_AIR, position_, vertical, &d);

        /* Count the number of unique intersections */
        const int n = (d > 0.);

        /* An odd number of intersections implies that the position is
         * inside the volume
         */
        return (n % 2) ? GEOMETRY_MEDIUM_ROCK : GEOMETRY_MEDIUM_AIR;
}


static struct {
        int n;
#define EMBREE_MAX_INTERSECTIONS 1024
        float distance[EMBREE_MAX_INTERSECTIONS];
} intersectionData;


static void topographyFilter(const struct RTCFilterFunctionNArguments * args)
{
        if (intersectionData.n < EMBREE_MAX_INTERSECTIONS - 2) {
                struct RTCRay * ray = (struct RTCRay *)args->ray;
                if (ray->tfar < FLT_MAX)
                        intersectionData.distance[intersectionData.n++] =
                            ray->tfar;
        }

        int * valid = (int *)args->valid;
        *valid = 0;
}


static int embree_distance_to_straight(int volume, const double position[3],
    const double direction[3], double * distance)
{
        static double dlast = 0.;
        static double dirlast[3] = { 0., 0., 0. };
        static int it = 0;

        if ((direction[0] != dirlast[0]) ||
            (direction[1] != dirlast[1]) ||
            (direction[2] != dirlast[2]))  {
                /* Set the intersection callback */
                rtcSetGeometryIntersectFilterFunction(
                    topography, &topographyFilter);
                intersectionData.n = 0;

                /* Do the ray tracing */
                embree_distance_to(volume, position, direction, NULL);

                /* Remove the intersection callback */
                rtcSetGeometryIntersectFilterFunction(topography, NULL);

                /* Reset the iteration count */
                it = 0;
                dlast = 0.;
                memcpy(dirlast, direction, sizeof dirlast);
        }

        /* Yield back the next result */
        if (it < intersectionData.n) {
                const double d = intersectionData.distance[it] - dlast;
                dlast = intersectionData.distance[it++];
                *distance = d;
                return (volume == GEOMETRY_MEDIUM_ROCK) ? GEOMETRY_MEDIUM_AIR :
                                                          GEOMETRY_MEDIUM_ROCK;
        } else {
                *distance = -1.;
                 return GEOMETRY_MEDIUM_VOID;
        }
}


static enum geometry_medium embree_medium_at(int volume)
{
        return (enum geometry_medium)volume;
}


enum embree_mode {
        EMBREE_MODE_FIRST = 0,
        EMBREE_MODE_STRAIGHT
};


extern void geometry_initialise_embree(
    struct geometry * geometry, const char * path, int period)
{
        /* Create the BVH */
        embree_create(path, period);

        /* Initialise the corresponding geometry interface */
        geometry->algorithm = "embree";
        if (strstr(path, "tianshan") != NULL)
                geometry->location = GEOMETRY_LOCATION_TIANSHAN;
        else
                geometry->location = GEOMETRY_LOCATION_PDD;
        geometry->volume_at = &embree_volume_at;
        geometry->medium_at = &embree_medium_at;
        geometry->clear = &embree_clear;

        geometry_configure_embree(geometry, EMBREE_MODE_FIRST);
}


extern void geometry_configure_embree(
    struct geometry * geometry, int mode)
{
        if (mode == EMBREE_MODE_FIRST) {
                geometry->distance_to = &embree_distance_to;
        } else if (mode == EMBREE_MODE_STRAIGHT) {
                geometry->distance_to = &embree_distance_to_straight;
        }
}
