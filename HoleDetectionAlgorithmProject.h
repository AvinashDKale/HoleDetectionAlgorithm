//// HoleDetectionAlgorithmProject.h : Include file for standard system include files,
//// or project specific include files.
//
//#pragma once
//#include <iostream>
//#include <chrono>
//#include <ctime>
//#include <vector>
//#include <boost/filesystem.hpp>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/Triangulation_hierarchy_2.h>
//#include <CGAL/Constrained_triangulation_plus_2.h>
//#include <CGAL/point_generators_2.h>
//#include<CGAL/convex_hull_2.h>
//#include <CGAL/enum.h>
//#include<CGAL/Polygon_2.h>
//#include <CGAL/Projection_traits_xy_3.h>
//#include <CGAL/Triangle_2.h>
//#include <pcl/point_cloud.h>
//#include <pcl/io/pcd_io.h>
//#include <pcl/io/obj_io.h>
//#include <pcl/common/io.h>
//#include <pcl/common/common.h>
//#include <pcl/point_types.h>
//#include <pcl/kdtree/kdtree.h>
//#include <pcl/kdtree/kdtree_flann.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/statistical_outlier_removal.h>
//#include <pcl/filters/extract_indices.h>
//#include <pcl/filters/uniform_sampling.h>
//#include <pcl/filters/voxel_grid_occlusion_estimation.h>
//#include <pcl/filters/radius_outlier_removal.h>
//#include <pcl/segmentation/progressive_morphological_filter.h>
//#include <pcl/segmentation/sac_segmentation.h>
//#include <pcl/segmentation/extract_clusters.h>
//#include <pcl/features/normal_3d.h>
//
//#include<stdio.h>
//#include<string>
//#include <exception>
//#include <algorithm>
//#include <boost/functional/hash.hpp> 
//#include <tuple>
//#include <unordered_map>
//#include <queue>
//#include <CGAL/nearest_neighbor_delaunay_2.h>
//#include <fstream>
//#include <sstream>
//#include <limits>
//#include <cstdlib>
////#include <GL/glut.h>
////#include <GL/glut.h>
//
//#include <CGAL/Alpha_shape_2.h>
//#include <CGAL/Alpha_shape_vertex_base_2.h>
//#include <CGAL/Alpha_shape_face_base_2.h>
//#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/algorithm.h>
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/centroid.h>
//#include <list>
//
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Projection_traits_xy_3<K>  Gt;
////typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
////typedef K::Point_3   Point;
//typedef K::Point_2 Point;
////typedef K::Polygon_2 Polygon;
//typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
////typedef CGAL::Gabriel_graph_2<DelaunayTriangulation> GabrielGraph;
////typedef CGAL::Polygon_2<K> Polygon;
////typedef Polygon::Vertex_iterator VertexIterator;
////typedef Polygon::Edge_const_iterator EdgeIterator;
//
//
//typedef CGAL::Triangle_3<K> Triangle;
//typedef std::list<Triangle> TriangleList;
//typedef std::vector<Point> Points;
//
//
//typedef K::FT                                                FT;
//typedef K::Segment_2                                         Segment;
//typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb1;
//typedef CGAL::Alpha_shape_face_base_2<K>                     Fb1;
//typedef CGAL::Triangulation_data_structure_2<Vb1, Fb1>          Tds;
//typedef CGAL::Delaunay_triangulation_2<K, Tds>                Triangulation_2;
//typedef CGAL::Alpha_shape_2<Triangulation_2>                 Alpha_shape_2;
//typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;
//
//
//
////typedef CGAL::Exact_predicates_tag Itag;
////typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
//typedef CGAL::Triangulation_vertex_base_2<K>Vbb;
//typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>Vb;
//typedef CGAL::Constrained_triangulation_face_base_2<K>Fb;
//typedef CGAL::Triangulation_data_structure_2<Vb, Fb>TDS;
//typedef CGAL::Exact_predicates_tag Itag;
//typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
//typedef CGAL::Triangulation_hierarchy_2<CDT> CDTH;
//typedef CGAL::Constrained_triangulation_plus_2<CDTH>Triangulation;
//typedef Triangulation::Point Point_T;
////typedef CDT::Point Point;
////typedef CDT::Edge  Edge;
//int cellularizationLevel = 48, reductionLevel = 28;
//
////Delaunay dt;
////float minX = std::numeric_limits<int>::max(), maxX = std::numeric_limits<int>::min();
////float minY = std::numeric_limits<int>::max(), maxY = std::numeric_limits<int>::min();
////Points inputPoints, exists;
////int numberOfHoles = 0;
//
//float v = 1; //parameter for outer bdry detection
//
//int window;
////float diagonalDistance;
////GLdouble width, height;
//
//class Edge
//{
//public:
//    Point source;
//    Point target;
//};
//
//typedef std::vector<Edge> Edges;
//
////Edges shape;
////Edges holeEdges;
//
//////Priority Queue
//struct node
//{
//    float length{ 0.0 };
//    Delaunay::All_faces_iterator ffi;
//    struct node* link;
//};
//
//class PriorityQueue
//{
//private:
//    node* front;
//public:
//    PriorityQueue()
//    {
//        front = NULL;
//    }
//
//    void pqInsert(float length, Delaunay::All_faces_iterator ffi)
//    {
//        node* tmp, * q;
//        tmp = new node;
//        tmp->ffi = ffi;
//        tmp->length = length;
//        if (front == NULL || length > front->length)
//        {
//            tmp->link = front;
//            front = tmp;
//        }
//        else
//        {
//            q = front;
//            while (q->link != NULL && q->link->length >= length)
//            {
//                q = q->link;
//            }
//            tmp->link = q->link;
//            q->link = tmp;
//        }
//    }
//
//    Delaunay::All_faces_iterator pqDelete()
//    {
//        node* tmp;
//        if (front == NULL)
//        {
//            std::cout << "Queue Underflow" << std::endl;
//        }
//        else
//        {
//            tmp = front;
//            Delaunay::All_faces_iterator frontFace;
//            frontFace = tmp->ffi;
//            front = front->link;
//            free(tmp);
//            return frontFace;
//        }
//    }
//
//    bool isEmpty()
//    {
//        if (front == NULL)
//            return 1;
//        return 0;
//    }
//};
//
////PriorityQueue pqOuterBoundary; //priority queue for outer boundary detection
////PriorityQueue pqInnerBoundary; //priority queue for inner boundary detection
//
//
//float area(Point a, Point b, Point c) /*compute the area of the triangle formed by the points a, b, and c*/
//{
//    return (float)std::abs(0.5 * (a.x() * (b.y() - c.y()) + b.x() * (c.y() - a.y()) + c.x() * (a.y() - b.y())));
//}
//
//double distance(Point a, Point b) /*compute the distance between two points*/
//{
//    return (sqrt(std::abs(((a.x() - b.x()) * (a.x() - b.x()))) + std::abs(((a.y() - b.y()) * (a.y() - b.y())))));
//}
//
//struct trianglePointIndices {
//    int i1, i2, i3;
//};
//
//// a useful general-purpose accessor
//auto as_tuple(trianglePointIndices const& v) -> decltype(auto)
//{
//    return std::tie(v.i1, v.i2, v.i3);
//}
//
//// equality implemented in terms of tuple, for simplicity
//bool operator==(trianglePointIndices const& l, trianglePointIndices const& r)
//{
//    return as_tuple(l) == as_tuple(r);
//}
//
//// hash_value implemented in terms of tuple, for consistency and simplicity
//std::size_t hash_value(trianglePointIndices const& v)
//{
//    using boost::hash_value;
//    return hash_value(as_tuple(v));
//}
//
//// the boring bit  injecting a hash specialisation into the std:: namespace
//// but let's derive from boost's hash class, which is much better
//// in that it allows easy hashing using free functions
//namespace std {
//    template<> struct hash<::trianglePointIndices> : boost::hash<::trianglePointIndices> {};
//}
//
//using trianglePointIndices_map = std::unordered_map<int, trianglePointIndices>;
//
//
//double circumradii(Point p1, Point p2, Point p3) {
//    return  std::sqrt(CGAL::squared_distance(p1, CGAL::circumcenter(p1, p2, p3)));
//}
//
//
//int  maximumOfThree(double a, double b, double c) {
//    auto m = std::max(a, std::max(b, c));
//    if (m == a) {
//        return 0;
//    }
//    else if (m == b) {
//        return 1;
//    }
//    else if (m == c) {
//        return 2;
//    }
//    else {
//        return -1;
//    }
//}
//
//
//
//
//
