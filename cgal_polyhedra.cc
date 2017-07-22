#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef K::Ray_3 Ray;
typedef K::Vector_3 Vector;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron, CGAL::Default, CGAL::Tag_false> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef CGAL::Bbox_3 Bbox_3;

typedef struct {Polyhedron *poly; Tree *tree;} Polytree;
typedef struct {double x,y,z;} d3;

extern "C" int debuglevel;

extern "C" {

  void polyhedron_from_file (Polytree **pptree, const char *fname, int verbose, int * const ierr) {

    Point p(1.0, 0.0, 0.0);
    Point q(0.0, 1.0, 0.0);
    Point r(0.0, 0.0, 1.0);
    Point s(0.0, 0.0, 0.0);
    Polyhedron *polyhedron1 = new Polyhedron;
    polyhedron1->make_tetrahedron(p, q, r, s);

    Tree *tree = new Tree(polyhedron1->facets_begin(), polyhedron1->facets_end(), *polyhedron1);

    tree->accelerate_distance_queries();

    *pptree = new Polytree;
    (*pptree)->poly = polyhedron1;
    (*pptree)->tree = tree;

    return;
  }

  void polyhedron_closest (const Polytree *ptree, const d3 *query, d3 *near) {
    Point query_point(query->x,query->y,query->z);

    Point closest = ptree->tree->closest_point(query_point);

    near->x = closest.x();
    near->y = closest.y();
    near->z = closest.z();
    return;
  }

  bool polyhedron_intersects_ray(const Polytree *ptree, const d3 *origin, const d3 *vec){
    Ray ray(Point(origin->x,origin->y,origin->z),
            Vector(vec->x,vec->y,vec->z));
    try{
      return ptree->tree->do_intersect(ray);
    }
    catch (...) {
      std::cout << origin->x <<" "<< origin->y <<" "<< origin->z << std::endl;
      std::cout << vec->x <<" "<< vec->y <<" "<< vec->z << std::endl;
      return false;
    }
    return false;
  }

  void polyhedron_bbox(const Polytree *ptree, d3 *const min, d3 *const max){
    Bbox_3 bbox = ptree->tree->bbox();
    *min = {bbox.xmin(), bbox.ymin(), bbox.zmin()};
    *max = {bbox.xmax(), bbox.ymax(), bbox.zmax()};
  }

  void polyhedron_finalize(Polytree **pptree){
    delete (*pptree)->tree; (*pptree)->tree = NULL;
    delete (*pptree)->poly; (*pptree)->poly = NULL;
    delete *pptree; *pptree = NULL;
  }

}





/*   /1* bool polyhedron_inside(const Polytree *ptree, const d3 *query, const d3 *outside_ref){ *1/ */
/*   /1*   Segment seg(Point(query->x,query->y,query->z), *1/ */
/*   /1*               Point(outside_ref->x,outside_ref->y,outside_ref->z)); *1/ */

/*   /1*   std::list<Object_and_primitive_id> intersections; *1/ */

/*   /1*   ptree->tree->all_intersections(seg, std::back_inserter(intersections)); *1/ */


/*   /1*   std::vector<Point> points; *1/ */

/*   /1*   int i = 0; *1/ */
/*   /1*   for (auto iter = intersections.begin(); iter != intersections.end(); ++iter){ *1/ */
/*   /1*     i += 1; *1/ */
/*   /1*     // gets intersection object *1/ */
/*   /1*     Object_and_primitive_id op = *iter; *1/ */
/*   /1*     CGAL::Object object = op.first; *1/ */
/*   /1*     Point point; *1/ */
/*   /1*     if(CGAL::assign(point,object)) { *1/ */
/*   /1*       points.push_back(point); *1/ */
/*   /1*     } *1/ */

/*   /1*   } *1/ */

/*   /1*   int n_dist = 0; *1/ */
/*   /1*   // find how many of the points are distinct *1/ */
/*   /1*   for (std::vector<Point>::size_type i = 0; i < points.size(); ++i){ *1/ */
/*   /1*     bool distinct = true; *1/ */

/*   /1*     for (std::vector<Point>::size_type j = 0; j < i; ++j){ *1/ */
/*   /1*       Vector v = points[i] - points[j]; *1/ */
/*   /1*       distinct = ( v.squared_length() > 1e-10 ); *1/ */

/*   /1*       if (!distinct)  break; *1/ */
/*   /1*     } *1/ */
/*   /1*     if (distinct) n_dist += 1; *1/ */
/*   /1*   } *1/ */

/*   /1*   return n_dist%2 == 1; *1/ */
/*   /1* } *1/ */

/*   /1* bool polyhedron_intersects_ray(const Polytree *ptree, const d3 *origin, const d3 *vec){ *1/ */
/*   /1*   Ray ray(Point(origin->x,origin->y,origin->z), *1/ */
/*   /1*           Vector(vec->x,vec->y,vec->z)); *1/ */
/*   /1*   try{ *1/ */
/*   /1*     return ptree->tree->do_intersect(ray); *1/ */
/*   /1*   } *1/ */
/*   /1*   catch (...) { *1/ */
/*   /1*     cout << origin->x <<" "<< origin->y <<" "<< origin->z << endl; *1/ */
/*   /1*     cout << vec->x <<" "<< vec->y <<" "<< vec->z << endl; *1/ */
/*   /1*     return false; *1/ */
/*   /1*   } *1/ */
/*   /1* } *1/ */

/*   /1* void polyhedron_bbox(const Polytree *ptree, d3 *const min, d3 *const max){ *1/ */
/*   /1*   Bbox_3 bbox = ptree->tree->bbox(); *1/ */
/*   /1*   *min = {bbox.xmin(), bbox.ymin(), bbox.zmin()}; *1/ */
/*   /1*   *max = {bbox.xmax(), bbox.ymax(), bbox.zmax()}; *1/ */
/*   /1* } *1/ */

/*   /1* void polyhedron_finalize(Polytree **pptree){ *1/ */
/*   /1*   delete (*pptree)->tree; (*pptree)->tree = NULL; *1/ */
/*   /1*   delete (*pptree)->poly; (*pptree)->poly = NULL; *1/ */
/*   /1*   delete *pptree; *pptree = NULL; *1/ */
/*   /1* } *1/ */

/* } */
