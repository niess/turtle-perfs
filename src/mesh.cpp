/* Standard library */
#include <iostream>
#include <list>
/* The TURTLE library */
#include "turtle.h"
/* The CGAL library */
#include <CGAL/Simple_cartesian.h>
/* The C interfaces */
extern "C" {
#include "downsampler.h"
#include "geometry.h"
}


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT Scalar;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Plane_3 Plane;


struct Cell {
        int volume;
        Plane faces[9];

        void set_face(int index, const Point & a, const Point & b,
            const Point & c);
        void build(int volume);
        bool inside(const Point & point);
        int distance(const Point & point, const Vector & direction,
            double * distance);
};


static struct turtle_map * gMap = NULL;
static struct turtle_map * gGeoid = NULL;
struct turtle_map_info gInfo;
static double gZBottom = -1E+03;
static double gZTop = 1E+04;
static struct Cell gCell;
static Vector gVertical;


static void get_nodes(int ix, int iy, int iz, Point & r_top, Point & r_bottom)
{
        if (ix < 0) ix = 0;
        else if (ix >= gInfo.nx) ix = gInfo.nx - 1;
        if (iy < 0) iy = 0;
        else if (iy >= gInfo.ny) iy = gInfo.ny - 1;

        double x, y, z;
        turtle_map_node(gMap, ix, iy, &x, &y, &z);
        const struct turtle_projection * p = turtle_map_projection(gMap);
        if (p != NULL) {
                const double xmap = x, ymap = y;
                turtle_projection_unproject(p, xmap, ymap, &y, &x);
        }
        double undulation;
        turtle_map_elevation(gGeoid, x, y, &undulation, NULL);
        z += undulation;

        double ecef[3];
        turtle_ecef_from_geodetic(y, x, z, ecef);
        if (iz == 0) {
                r_top = Point(ecef[0], ecef[1], ecef[2]);
                r_bottom = r_top + (gZBottom - z) * gVertical;
        } else {
                r_bottom = Point(ecef[0], ecef[1], ecef[2]);
                r_top = r_bottom + (gZTop - z) * gVertical;
        }
}


static int index_to_volume(int ix, int iy, int iz)
{
        if (iz < 0) return -1;
        else return gInfo.nx * (iz * gInfo.ny + iy) + ix;
}


static void volume_to_index(int volume, int * ix, int * iy, int * iz)
{
        if (volume < 0) {
                *ix = *iy = *iz = -1;
        } else {
                *iz = volume / (gInfo.nx * gInfo.ny);
                *iy = volume / gInfo.nx - (*iz) * gInfo.ny;
                *ix = volume % gInfo.nx;
        }
}


void Cell::set_face(
    int index, const Point & a, const Point & b, const Point & c)
{
       faces[index] = Plane(a, b, c);
}


void Cell::build(int volume_)
{
        if (volume_ == volume) return;
        volume = volume_;

        int ix, iy, iz;
        volume_to_index(volume, &ix, &iy, &iz);

        /* Get the node data in ECEF */
        Point nodes[8];
        get_nodes(ix, iy, iz, nodes[0], nodes[4]);
        get_nodes(ix + 1, iy , iz, nodes[1], nodes[5]);
        get_nodes(ix + 1, iy + 1, iz, nodes[2], nodes[6]);
        get_nodes(ix, iy + 1, iz, nodes[3], nodes[7]);

        /* Model the cell */
        int index = 0;
        set_face(index++, nodes[0], nodes[1], nodes[2]);
        set_face(index++, nodes[2], nodes[3], nodes[0]);

        for (int i = 0; i < 4; i++) {
                const int i00 = i;
                const int i10 = (i + 1) % 4;
                const int i01 = i00 + 4;
                set_face(index++, nodes[i10], nodes[i00],
                    nodes[i01]);
        }
        set_face(index++, nodes[4], nodes[6], nodes[5]);
        set_face(index++, nodes[6], nodes[4], nodes[7]);
        set_face(index++, nodes[2], nodes[0], nodes[4]);
}


bool Cell::inside(const Point & point)
{
        for (int i = 0; i < 8; i++) {
                if (faces[i].has_on_positive_side(point))
                        return false;
        }
        return true;
}


static double distance_to_face(const Point & point, const Vector & direction,
    const Plane & face)
{
        const double a = face.a();
        const double b = face.b();
        const double c = face.c();
        const double d = face.d();
        const double dp = a * direction[0] + b * direction[1] +
            c * direction[2];
        if (fabs(dp) < DBL_EPSILON)
                return DBL_MAX;
        else
                return -(a * point[0] + b * point[1] +
                    c * point[2] + d) / dp;
}


int Cell::distance(
    const Point & point, const Vector & direction, double * distance)
{
        if (volume < 0) {
                *distance = -1.;
                return -1;
        }

        /* Locate the closest face */
        const int positive[] = { 0, 2, 3, 6, 8 };
        const int negative[] = { 1, 4, 5, 7, 8 };
        const int * iface = faces[8].has_on_positive_side(point) ?
            positive : negative;
        double distance_ = DBL_MAX;
        int index = -1;
        for (int ii = 0; ii < 5; ii++) {
                const int i = iface[ii];
                const double di = distance_to_face(
                    point, direction, faces[i]);
                if ((di < DBL_EPSILON) || (di == DBL_MAX))
                        continue;
                else if (di < distance_) {
                        distance_ = di;
                        index = i;
                }
        }
        if (index < 0) {
                *distance = -1.;
                return -1;
        } else {
                *distance = distance_;
        }

        /* Locate the next volume */
        if (index == 8) return volume;

        int ix, iy, iz;
        volume_to_index(volume, &ix, &iy, &iz);
        if (index < 2) {
                if (iz == 0) iz  = 1;
                else iz = -1;
        } else if (index == 2) {
                if (--iy < 0) iz = -1;
        } else if (index == 3) {
                if (++ix >= gInfo.nx - 1) iz = -1;
        } else if (index == 4) {
                if (++iy >= gInfo.ny - 1) iz = -1;
        } else if (index == 5) {
                if (--ix < 0) iz = -1;
        } else {
                if (iz == 1) iz = 0;
                else iz = -1;
        }

        return index_to_volume(ix, iy, iz);
}


/* Build the mesh from a topography file */
static void mesh_create(const char * filename, int period)
{
        /* Load the topography file */
        downsampler_map_load(&gMap, filename, period);
        turtle_map_load(&gGeoid, GEOMETRY_GEOID);

        /* Set the top altitude */
        if (strstr(filename, "tianshan") != NULL) {
                gZTop = GEOMETRY_ZMAX_TIANSHAN;
        } else {
                gZTop = GEOMETRY_ZMAX_PDD;
        }

        /* Get the vertical */
        turtle_map_meta(gMap, &gInfo, NULL);
        double x0 = 0.5 * (gInfo.x[0] + gInfo.x[1]);
        double y0 = 0.5 * (gInfo.y[0] + gInfo.y[1]);
        const struct turtle_projection * p = turtle_map_projection(gMap);
        if (p) {
                const double xx = x0;
                const double yy = y0;
                turtle_projection_unproject(p, xx, yy, &y0, &x0);
        }

        double n[3];
        turtle_ecef_from_horizontal(y0, x0, 0., 90., n);
        gVertical = Vector(n[0], n[1], n[2]);
}


static int mesh_volume_at(const double position[3])
{
        double latitude, longitude, altitude;
        turtle_ecef_to_geodetic(position, &latitude, &longitude, &altitude);
        double x, y, z;
        const struct turtle_projection * p = turtle_map_projection(gMap);
        if (p != NULL) {
                turtle_projection_project(p, latitude, longitude, &x, &y);
        } else {
                x = longitude;
                y = latitude;
        }
        int checkMap;
        turtle_map_elevation(gMap, x, y, &z, &checkMap);
        if (checkMap == 0) {
                return GEOMETRY_MEDIUM_VOID;
        }

        double undulation;
        turtle_map_elevation(gGeoid, longitude, latitude, &undulation, NULL);
        z += undulation;

        int ix, iy, iz;
        const double dx = (gInfo.x[1] - gInfo.x[0]) / (gInfo.nx - 1);
        const double dy = (gInfo.y[1] - gInfo.y[0]) / (gInfo.ny - 1);
        ix = (int)((x - gInfo.x[0]) / dx);
        if ((ix < 0) || (ix >= gInfo.nx)) return -1;
        iy = (int)((y - gInfo.y[0]) / dy);
        if ((iy < 0) || (iy >= gInfo.ny)) return -1;
        iz = (altitude > z) ? 1 : 0;

        Point r(position[0], position[1], position[2]);
        int volume = index_to_volume(ix, iy, iz);
        gCell.build(volume);

        double d;
        if (iz == 0) {
                gCell.distance(r, gVertical, &d);
        } else {
                gCell.distance(r, -gVertical, &d);
        }
        if (d < 0) {
                iz = !iz;
                volume = index_to_volume(ix, iy, iz);
        }

        return volume;
}


static int mesh_distance_to(int volume, const double position[3],
    const double direction[3], double * distance)
{
        Point r(position[0], position[1], position[2]);
        Vector u(direction[0], direction[1], direction[2]);

        gCell.build(volume);
        const int next = gCell.distance(r, u, distance);
        if (*distance > 0)
                *distance += 1E-06; /* Push out of the current cell */
        return next;
}


static enum geometry_medium mesh_medium_at(int volume)
{
        if (volume < 0) {
                return GEOMETRY_MEDIUM_VOID;
        } else {
                int ix, iy, iz;
                volume_to_index(volume, &ix, &iy, &iz);
                return iz ? GEOMETRY_MEDIUM_AIR : GEOMETRY_MEDIUM_ROCK;
        }
}


static void mesh_clear()
{
        turtle_map_destroy(&gMap);
        turtle_map_destroy(&gGeoid);
}


extern "C" void geometry_initialise_mesh(
    struct geometry * geometry, const char * path, int period)
{
        /* Create the Mesh */
        mesh_create(path, period);

        /* Initialise the corresponding geometry interface */
        geometry->algorithm = "mesh";
        if (strstr(path, "tianshan") != NULL)
                geometry->location = GEOMETRY_LOCATION_TIANSHAN;
        else
                geometry->location = GEOMETRY_LOCATION_PDD;
        geometry->volume_at = &mesh_volume_at;
        geometry->distance_to = &mesh_distance_to;
        geometry->medium_at = &mesh_medium_at;
        geometry->clear = &mesh_clear;
}
