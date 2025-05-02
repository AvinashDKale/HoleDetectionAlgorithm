#pragma once
#include <queue>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/point_generators_2.h>
#include<CGAL/convex_hull_2.h>
#include <CGAL/enum.h>
#include<CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/nearest_neighbor_delaunay_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/centroid.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/AABB_tree.h>


#define DEBUG_MODE 1

//CGAL typedefinations
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;
typedef CGAL::Delaunay_triangulation_2<Kernel> Delaunay;
typedef CGAL::Triangle_3<Kernel> Triangle;
typedef std::list<Triangle> TriangleList;
typedef std::vector<Point> Points;

//face iterator
typedef Delaunay::All_faces_iterator allFaceItr;
typedef Delaunay::Finite_faces_iterator finitFaceItr;
typedef Delaunay::Face_handle faceHandle;

//k-d tree
typedef CGAL::Search_traits_2<Kernel> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;


class Edge
{
public:
    Point source;
    Point target;
};

typedef std::vector<Edge> Edges;

//Priority Queue Node
class node
{
public:
    float length;
    allFaceItr ffi;
    node() :length(0.0), ffi(nullptr) {};
    node(float nlength, allFaceItr nffi) :length(nlength), ffi(nffi) {};

    //Overloading > operator.
    friend bool operator> (const node& node1, const node& node2)
    {
        return node1.length < node2.length;
    }
};

//Priority Queue
typedef std::priority_queue<node, std::vector<node>, std::greater<std::vector<node>::value_type> > PriorityQueue;

//FindHole class
class FindHole {

private:
    std::vector<Point> m_points;
    Delaunay m_delaunayTriangulation;

protected:
    // to find the index of point of hole from input point set
    void indicesOfHoles(std::vector<std::vector<Point>>& holes, std::vector<std::vector<int>>& holeIndices);

    // Is triangle is on boundary or not 
    // Checking whether the face is infinite or not
    bool isInfinite( allFaceItr afi);

    // Is triangle inside and not on boundary
    // Checking whether the face is finite or not
    bool isFinite(allFaceItr afi);

    // lies inside scalen
    // To check whether p3 is inside scaled diametric circle(p1,p2)
    bool liesInsideScaled( Point a, Point b, Point c);

    // lines insied unscalen
   // To check whether p3 is inside unscaled diametric circle(p1, p2)
    bool liesInsideUnscaled( Point a, Point b, Point c);

    // inside chord/mid circle
    // Checking whether mid point and chord circles are empty
    bool liesWithinException( Point a, Point b, Point c, double radius);

    // Reularity constrant
    // Checking whether a point is already in shape
    bool isInShape( Points& exists, Point a);

    // to update the outerboudary by adding hole as boundary
    // Inserting a new edge into boundary
    void updateNeighbors( PriorityQueue& pqOuterBoundary, allFaceItr currentFace, allFaceItr neighbor1, allFaceItr neighbor2, float d1, float d2);

    // inserting into shape/outer boundaary
    void insertToShape( Edges& shape, Point a, Point b);

    // to get maximum area triangle in the delaunay
    faceHandle getMaxAreaTriangleHandle( Edges& shape, finitFaceItr ffi);

    // condtion to check part of hole 
    bool liesInsideHole( Point p1, Point p2, Point p3);

    //Checking whether mid point and chord circles are empty
    bool liesWithinExceptionHole( Point a, Point b, Point c, double radius);

    //Inserting a new edge into (empty circle)ec-shape/boundary
    void insertToHoleBoundary( Edges& holeEdges, Point a, Point b);

    // updates the hole boundary
    void updateHole( PriorityQueue& pqInnerBoundary, allFaceItr currentFace, allFaceItr neighbor1, allFaceItr neighbor2);

    // Checking whether a point is already in shape
    bool inBoundary( Points exists, Point a);

    //Checking the status of mid point circle and chord circles for outerboundary
    bool outerBoundaryCondition(Point a, Point b, Point c, faceHandle fh);

    //Checking the status of mid point circle and chord circles for hole boundary
    bool holeDetectionCondition( Point a, Point b, Point c, faceHandle fh);


public:
    //scaling factor
    static constexpr float scaledFactor = 1;

    //default constructor not  used
    FindHole() = delete;   

    //parametric constructor
    FindHole(const std::vector<Point>&);

    //destructor
    virtual ~FindHole();

    // copy constructor
    FindHole(const FindHole&) = delete;

    // move constructor
    FindHole(FindHole&&) = delete;

    // assignment operator for copy constructor
    FindHole& operator=(const FindHole&) = delete;

    // assignment operator for move constructor
    FindHole& operator=(FindHole&&) = delete;

    // outer boundary
    void outerBoundary( Points& exists, Edges& shape, PriorityQueue& pqOuterBoundary);

    // holes detection
    void holeDetection( Points exists, Edges& shape, Edges& holeEdges, PriorityQueue& pqInnerBoundary, std::vector< std::vector<Point>>& holes, int& numberOfHoles);

    // single hole
    void singleHoleDetection( Edges& shape, Edges& holeEdges, std::vector< std::vector<Point>>& holes, int& numberOfHoles);

    //detecting and getting holes
    int getholes(std::vector<std::vector<Point>>& holes, std::vector<std::vector<int>>& holeIndices);

};
