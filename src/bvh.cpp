#include <iostream>
#include <list>
/* The TURTLE library */
#include "turtle.h"
/* The CGAL library */
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
/* The C interfaces */
extern "C" {
#include "downsampler.h"
#include "geometry.h"
#include "tictac.h"
}


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT Scalar;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;

typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

typedef Kernel::Ray_3 Ray;
typedef CGAL::AABB_traits<Kernel, Primitive>::Bounding_box Bounding_box;
typedef boost::optional<
    Tree::Intersection_and_primitive_id<Ray>::Type > Ray_intersection;


static Point gOrigin;
static std::list<Triangle> gTriangles;
static Tree gTree;

static struct turtle_map * gMap, * gGeoid;
static const struct turtle_projection * gProjection = NULL;


static Point get_node(int ix, int iy)
{
        double x, y, z;
        turtle_map_node(gMap, ix, iy, &x, &y, &z);

        if (gProjection != NULL) {
                turtle_projection_unproject(
                    gProjection, x, y, &y, &x);
        }

        double undulation;
        turtle_map_elevation(gGeoid, x, y, &undulation, NULL);
        z += undulation;

        double ecef[3];
        turtle_ecef_from_geodetic(y, x, z, ecef);
        return Point(ecef[0], ecef[1], ecef[2]);
}


/* Create the BVH from a topography map */
static void bvh_create(const char * path, int period_)
{
        /* Load the map data */
        downsampler_map_load(&gMap, path, period_);
        gProjection = turtle_map_projection(gMap);
        turtle_map_load(&gGeoid, GEOMETRY_GEOID);

        struct turtle_map_info info;
        turtle_map_meta(gMap, &info, NULL);
        const int nx = info.nx, ny = info.ny;

        /* Build the facets */
        std::cout << "# building facets ...";
        struct tictac tictac;
        tictac_initialise(&tictac);
        tictac.start(&tictac);

        gOrigin = get_node(nx / 2, ny / 2);
        Vector offset = gOrigin - CGAL::ORIGIN;

        const int n_iter = (nx - 1) * (ny - 1);
        int period = n_iter / 100;
        if (period <= 0) period = 1;

        int ix, iy, it = 0;
        for (iy = 0; iy < ny - 1; iy++) for (ix = 0; ix < nx - 1; ix++, it++) {
                Point n00 = get_node(ix, iy);
                Point n10 = get_node(ix + 1, iy);
                Point n01 = get_node(ix, iy + 1);
                Point n11 = get_node(ix + 1, iy + 1);

                Triangle triangle0(n00 - offset, n10 - offset, n11 - offset);
                gTriangles.push_back(triangle0);
                Triangle triangle1(n11 - offset, n01 - offset, n00 - offset);
                gTriangles.push_back(triangle1);

                if ((it + 1) % period == 0) {
                        std::cout << "\r# building facets ... "
                                  << (it + 1) / period << " %";
                        fflush(stdout);
                }
        }
        double dt = tictac.stop(&tictac);
        std::cout << "\r# " <<  2 * n_iter << " facets built in "
                  << dt << " s" << std::endl;

        /* Free the maps */
        turtle_map_destroy(&gMap);
        turtle_map_destroy(&gGeoid);
        gProjection = NULL;

        /* Construct the AABB tree */
        std::cout << "# building the AABB tree ...";
        tictac.start(&tictac);
        gTree.insert(gTriangles.begin(), gTriangles.end());
        gTree.build();
        dt = tictac.stop(&tictac);
        std::cout << "\r# AABB tree has been built in "
                  << dt << " s" << std::endl;
}


static void bvh_clear(void)
{
        gTree.clear();
        gTriangles.clear();
}


static double bounding_distance(
    const Bounding_box & box, const Point & position)
{
        Scalar r = 0;
        for (int i = 0; i < 3; i++) {
                const Scalar s0 = fabs(box.max(i) - position[i]);
                const Scalar s1 = fabs(position[i] - box.min(i));
                const Scalar s = s0 > s1 ? s0 : s1;
                if (s > r) r = s;
        }

        const Scalar safety = 1.1;
        return r * safety;
}


inline static Point to_local(const Point & global_position)
{
        Vector o = gOrigin - CGAL::ORIGIN;
        return global_position - o;
}


static int bvh_volume_at(const double position_[3])
{
        /* Translate to the local frame */
        const Point global_position(position_[0], position_[1], position_[2]);
        Point position = to_local(global_position);

        /* 1st check the bounding box of the tree */
        const Bounding_box box = gTree.bbox();
        if ((position.x() < box.xmin()) || (position.x() > box.xmax()) ||
            (position.y() < box.ymin()) || (position.y() > box.ymax()) ||
            (position.z() < box.zmin()) || (position.z() > box.zmax())) {
                return GEOMETRY_MEDIUM_VOID;
        }

        /* Draw a ray along the local vertical */
        const Scalar r = bounding_distance(box, position);
        double ecef[3] = { global_position.x(), global_position.y(),
            global_position.z() };
        double latitude, longitude, altitude;
        turtle_ecef_to_geodetic(ecef, &latitude, &longitude, &altitude);
        double vertical[3];
        turtle_ecef_from_horizontal(latitude, longitude, 0., 90., vertical);
        Vector u(vertical[0], vertical[1], vertical[2]);
        Point end = position + r * u;
        Ray ray_query(position, end);

        /* Count the number of unique intersections */
        std::list<Ray_intersection> intersections;
        std::set<Point> uniques;
        gTree.all_intersections(
            ray_query, std::back_inserter(intersections));
        for (std::list<Ray_intersection>::iterator it =
            intersections.begin(); it != intersections.end(); it++) {
                const Point * p = boost::get<Point>(&((*it)->first));
                if (p) uniques.insert(*p);
        }

        /* An odd number of intersections implies that the position is
         * inside the volume
         */
        return (uniques.size() % 2) ? GEOMETRY_MEDIUM_ROCK :
                                      GEOMETRY_MEDIUM_AIR;
}


static Ray get_ray(const Point & position, const Vector & direction)
{
        const Bounding_box box = gTree.bbox();
        const Scalar r = bounding_distance(box, position);
        Point end = position + r * direction;
        return Ray(position, end);
}


static int bvh_distance_to(int volume, const double position_[3],
    const double direction_[3], double * distance)
{
        /* Translate to the local frame */
        const Point global_position(position_[0], position_[1], position_[2]);
        Point position = to_local(global_position);

        /* Get a long enough ray */
        const Vector direction(direction_[0], direction_[1], direction_[2]);
        Ray ray_query = get_ray(position, direction);

        /* Get all intersections and return the closest one */
        Scalar closest = -1;
        std::list<Ray_intersection> intersections;
        gTree.all_intersections(ray_query, std::back_inserter(intersections));
        for (std::list<Ray_intersection>::iterator it =
            intersections.begin(); it != intersections.end(); it++) {
                const Point * point = boost::get<Point>(&((*it)->first));
                if (point) {
                        Vector v = position - *point;
                        Scalar d = sqrt(v.squared_length());
                        if ((closest < 0) || (d < closest))
                                closest = d;
                }
        }
        if (closest < 0) {
                *distance = closest;
                return GEOMETRY_MEDIUM_VOID;
        } else {
                *distance = closest + 1E-06; /* Protect against rounding errors */
                return (volume == GEOMETRY_MEDIUM_ROCK) ? GEOMETRY_MEDIUM_AIR :
                                                          GEOMETRY_MEDIUM_ROCK;
        }
}


static int bvh_distance_to_first(int volume, const double position_[3],
    const double direction_[3], double * distance)
{
        /* Translate to the local frame */
        const Point global_position(position_[0], position_[1], position_[2]);
        Point position = to_local(global_position);

        /* Get a long enough ray */
        const Vector direction(direction_[0], direction_[1], direction_[2]);
        Ray ray_query = get_ray(position, direction);

        /* Get the closest intersection */
        Scalar closest = -1;
        Ray_intersection intersection = gTree.first_intersection(ray_query);
        if (intersection) {
                const Point * point = boost::get<Point>(&intersection->first);
                Vector v = position - *point;
                closest = sqrt(v.squared_length());
        }

        if (closest < 0) {
                *distance = closest;
                return GEOMETRY_MEDIUM_VOID;
        } else {
                *distance = closest + 1E-06; /* Protect against rounding errors */
                return (volume == GEOMETRY_MEDIUM_ROCK) ? GEOMETRY_MEDIUM_AIR :
                                                          GEOMETRY_MEDIUM_ROCK;
        }
}


static int bvh_distance_to_straight(int volume, const double position_[3],
    const double direction_[3], double * distance)
{
        static std::set<Scalar> distances;
        static std::set<Scalar>::iterator it=distances.begin();
        static double dlast = 0.;
        static double dirlast[3] = { 0., 0., 0. };

        if ((direction_[0] != dirlast[0]) ||
            (direction_[1] != dirlast[1]) ||
            (direction_[2] != dirlast[2]))  {
                /* Translate to the local frame */
                const Point global_position(
                    position_[0], position_[1], position_[2]);
                Point position = to_local(global_position);

                /* Get a long enough ray */
                const Vector direction(
                    direction_[0], direction_[1], direction_[2]);
                Ray ray_query = get_ray(position, direction);

                /* Get all intersections and order them */
                std::list<Ray_intersection> intersections;
                gTree.all_intersections(
                    ray_query, std::back_inserter(intersections));
                distances.clear();
                for (std::list<Ray_intersection>::iterator i =
                    intersections.begin(); i != intersections.end(); i++) {
                        const Point * point =
                            boost::get<Point>(&((*i)->first));
                        if (point) {
                                Vector v = position - *point;
                                Scalar d = sqrt(v.squared_length());
                                distances.insert(d);
                        }
                }

                it = distances.begin();
                dlast = 0.;
                memcpy(dirlast, direction_, sizeof dirlast);
        }

        /* Yield back the next result */
        if (it != distances.end()) {
                const double d = *it - dlast;
                dlast = *it++;
                *distance = d;
                return (volume == GEOMETRY_MEDIUM_ROCK) ? GEOMETRY_MEDIUM_AIR :
                                                          GEOMETRY_MEDIUM_ROCK;
        } else {
                *distance = -1.;
                 return GEOMETRY_MEDIUM_VOID;
        }
}


static enum geometry_medium bvh_medium_at(int volume)
{
        return (enum geometry_medium)volume;
}


enum bvh_mode {
        BVH_MODE_ALL = 0,
        BVH_MODE_STRAIGHT,
        BVH_MODE_FIRST
};


extern "C" void geometry_initialise_bvh(
    struct geometry * geometry, const char * path, int period)
{
        /* Create the BVH */
        bvh_create(path, period);

        /* Initialise the corresponding geometry interface */
        geometry->algorithm = "BVH";
        if (strstr(path, "tianshan") != NULL)
                geometry->location = GEOMETRY_LOCATION_TIANSHAN;
        else
                geometry->location = GEOMETRY_LOCATION_PDD;
        geometry->volume_at = &bvh_volume_at;
        geometry->medium_at = &bvh_medium_at;
        geometry->clear = &bvh_clear;

        geometry_configure_bvh(geometry, BVH_MODE_ALL);
}


extern "C" void geometry_configure_bvh(
    struct geometry * geometry, int mode)
{
        if (mode == BVH_MODE_ALL) {
                geometry->distance_to = &bvh_distance_to;
        } else if (mode == BVH_MODE_STRAIGHT) {
                geometry->distance_to = &bvh_distance_to_straight;
        } else if (mode == BVH_MODE_FIRST)  {
                geometry->distance_to = &bvh_distance_to_first;
        }
}
