#include "FindHole.h"

void saveDelaunayTriangulationTxtFile(const Delaunay& dt, std::vector<Point> points) {
    // Iterate through triangles and print them
    std::string str = "";
    int trnls = 0;
    for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) {
        if (dt.is_infinite(fit)) continue;  // Skip infinite faces
        Delaunay::Face_handle face = fit;
        ++trnls;
        Delaunay::Point p1 = fit->vertex(0)->point();
        Delaunay::Point p2 = fit->vertex(1)->point();
        Delaunay::Point p3 = fit->vertex(2)->point();
        std::string temp = std::to_string(p1.x()) + " " + std::to_string(p1.y()) + "\n" +
            std::to_string(p2.x()) + " " + std::to_string(p2.y()) + "\n" +
            std::to_string(p3.x()) + " " + std::to_string(p3.y()) + "\n" +
            std::to_string(p1.x()) + " " + std::to_string(p1.y()) + "\n";
        str = str + temp;
        str = str + "\n\n";

    }
    try {
        std::string filename = "E:\\Avinash C drive\\Downloads\\gp610-20230826-win64-mingw\\gnuplot\\bin\\sample2DPtsToTriangle.txt";
        std::ofstream file(filename);

        if (!file.is_open())
            throw "Error to Open File !!!" + filename;
       
        if (str.empty())
            throw "Write the file atleast one char !!!";
        file.write(str.c_str(), str.length());
        file.close();

        std::cout << "\n Delaunay triangulation file saved ! " << std::endl;
    }
    catch (...) {
        std::cout << "ERROR!!! 40 SampleCodeGeneration1.h  " << std::endl;
    }
}


FindHole::FindHole(const std::vector<Point>&points)
    :m_points(points)
{
    this->m_delaunayTriangulation.insert(points.begin(), points.end());
}


FindHole::~FindHole()
{
}


bool FindHole::isInfinite( allFaceItr afi)
{
    for (int i = 0; i < 3; ++i)
    {
        if (afi->vertex(i) == m_delaunayTriangulation.infinite_vertex())
            return true;
    }
    return false;
}


bool FindHole::isFinite( allFaceItr afi) 
{
    for (int i = 0; i < 3; ++i)
    {
        if (afi->vertex(i) != m_delaunayTriangulation.infinite_vertex())
            return true;
    }
    return false;
}


inline bool FindHole::liesInsideScaled(Point a, Point b, Point c)
{
    if (sqrt(CGAL::squared_distanceC2((a.x() + b.x()) / 2, (a.y() + b.y()) / 2, c.x(), c.y())) 
        <=
        scaledFactor * sqrt(CGAL::squared_distanceC2((a.x() + b.x()) / 2, (a.y() + b.y()) / 2, a.x(), a.y())))
        return true;
    return false;
}


inline bool FindHole::liesInsideUnscaled( Point a, Point b, Point c)
{
    if (sqrt(CGAL::squared_distanceC2((a.x() + b.x()) / 2, (a.y() + b.y()) / 2, c.x(), c.y()))
        <=
        sqrt(CGAL::squared_distanceC2((a.x() + b.x()) / 2, (a.y() + b.y()) / 2, a.x(), a.y())))
        return true;
    return false;
}


bool FindHole::liesWithinException(Point a, Point b, Point c, double radius)
{
    auto distanceAB = sqrt(CGAL::squared_distanceC2(a.x(), a.y(), b.x(), b.y()));
    if (distanceAB >= radius)
    {
        //(p1,p2) is the midpoint of edge (a,b)
        if (
            liesInsideScaled
            (
                Point(((a.x() + b.x()) / 2 - (radius / 2)), ((a.y() + b.y()) / 2), 1),
                Point(((a.x() + b.x()) / 2 + (radius / 2)), ((a.y() + b.y()) / 2), 1),
                c
            )
           )
            return true;
        return false;
    }
    radius = (radius / 2) * scaledFactor;
    double d = sqrt(std::abs((radius * radius) - ((distanceAB / 2) * (distanceAB / 2))));
    double a1New = 
        ((a.x() + b.x()) / 2) +
        (d / 
            (
                sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x())))
            )
        ) *
        (a.y() - b.y());
    double b1New = 
        ((a.y() + b.y()) / 2) + 
        (d /
            (
                sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x())))
                )
            ) *
        (b.x() - a.x());
    double a2New = 
        ((a.x() + b.x()) / 2) -
        (d /
            (
                sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x())))
                )
            ) *
        (a.y() - b.y());
    double b2New = 
        ((a.y() + b.y()) / 2) - 
        (d /
            (
                sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x())))
                )
            ) *
        (b.x() - a.x());

    if (liesInsideUnscaled( Point(a1New + radius, b1New, 1), Point(a1New - radius, b1New, 1), c))
        return true;
    if (liesInsideUnscaled( Point(a2New + radius, b2New, 1), Point(a2New - radius, b2New, 1), c))
        return true;
    return false;
}


inline bool FindHole::isInShape( Points& exists, Point a)
{
    for (int i = 0; i < exists.size(); i++)
    {   
        if (std::abs(exists[i].x() - a.x())<= DBL_EPSILON && 
            std::abs(exists[i].y() - a.y())<= DBL_EPSILON)
            return true;
    }
    return false;
}


void FindHole::updateNeighbors
( 
    PriorityQueue& pqOuterBoundary,
    allFaceItr currentFace,
    allFaceItr neighbor1,
    allFaceItr neighbor2,
    float area1, 
    float area2
)
{
    for (int i = 0; i < 3; ++i) {
        if (neighbor1->neighbor(i) == currentFace)
        {
            neighbor1->set_neighbor(i, this->m_delaunayTriangulation.infinite_face());
        }
        if (neighbor2->neighbor(i) == currentFace)
        {
            neighbor2->set_neighbor(i, this->m_delaunayTriangulation.infinite_face());
        }
    }
    this->m_delaunayTriangulation.delete_face(currentFace);
    pqOuterBoundary.push(node(area1, neighbor1));
    pqOuterBoundary.push(node(area2, neighbor2));

}


inline void FindHole::insertToShape( Edges& shape, Point a, Point b)
{
    Edge e;
    e.source = a;
    e.target = b;
    shape.push_back(e);
}


faceHandle FindHole::getMaxAreaTriangleHandle( Edges& shape, finitFaceItr ffi)
{
    double biggestArea = 0.0;
    bool presentInShape = false;
    Delaunay::Face_handle highestAreaTriangleHandle;
    ffi = this->m_delaunayTriangulation.finite_faces_begin();//output argument
    int flag = 0;

    for (Delaunay::Finite_faces_iterator ffi = this->m_delaunayTriangulation.finite_faces_begin();
        ffi != this->m_delaunayTriangulation.finite_faces_end();
        ffi++
        )
    {
        Point a = ffi->vertex(0)->point();
        Point b = ffi->vertex(1)->point();
        Point c = ffi->vertex(2)->point();

        
        for (int i = 0; i < shape.size(); i++)
        {
            for (int v = 0; v < 3; ++v) {

                if (
                    (
                        std::abs(ffi->vertex(v)->point().x() - shape[i].source.x())<=DBL_EPSILON && 
                        std::abs(ffi->vertex(v)->point().y() - shape[i].source.y()) <= DBL_EPSILON
                        ) 
                    ||
                    (std::abs(ffi->vertex(v)->point().x() - shape[i].target.x()) <=DBL_EPSILON && 
                        std::abs(ffi->vertex(v)->point().y() - shape[i].target.y()) <= DBL_EPSILON
                        )
                   )
                    presentInShape = true;
            }
        }

        if (!presentInShape)
        {
            if (isFinite( ffi->neighbor(0)) && 
                isFinite( ffi->neighbor(1)) && 
                isFinite( ffi->neighbor(2))
                )
            {
                flag = 1;
                
                if (
                    CGAL::area(ffi->vertex(0)->point(), ffi->vertex(1)->point(), ffi->vertex(2)->point()) >= biggestArea)
                {
                    biggestArea = CGAL::area(ffi->vertex(0)->point(), ffi->vertex(1)->point(), ffi->vertex(2)->point());
                    highestAreaTriangleHandle = ffi;
                }
            }
        }
        presentInShape = false;
    }
    if (flag == 0)
    {
        highestAreaTriangleHandle = this->m_delaunayTriangulation.infinite_face();
    }
    return highestAreaTriangleHandle;
}


inline bool FindHole::liesInsideHole( Point p1, Point p2, Point p3)
{    
    if (sqrt(CGAL::squared_distanceC2((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2, p3.x(), p3.y()))
        <= CGAL::squared_distanceC2((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2, p1.x(),p1.y()))
    {
        return true;
    }
    return 0;
}


bool FindHole::liesWithinExceptionHole( Point a, Point b, Point c, double radius)
{
    auto distanceAB= sqrt(CGAL::squared_distanceC2(a.x(), a.y(), b.x(), b.y()));

    if (distanceAB >= radius)
    {
        //(p1,p2) is the midpoint of edge (a,b)
        if (
            liesInsideScaled
            (
                Point(((a.x() + b.x()) / 2 - (radius / 2)), ((a.y() + b.y()) / 2), 1),
                Point(((a.x() + b.x()) / 2 + (radius / 2)), ((a.y() + b.y()) / 2), 1), 
                c
            )
           )
        {
            return true;
        }
        return false;
    }
    radius = (radius / 2);
    double d = sqrt(std::abs((radius * radius) - ((distanceAB / 2) * (distanceAB / 2))));
    double a1New =
        ((a.x() + b.x()) / 2) +
        (d /
            (
                sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x())))
                )
            ) *
        (a.y() - b.y());
    double b1New = 
        ((a.y() + b.y()) / 2) +
        (d / 
            (
                sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x())))
                )
            ) * 
        (b.x() - a.x());
    double a2New =
        ((a.x() + b.x()) / 2) - 
        (d /
            (
                sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x())))
                )
            ) *
        (a.y() - b.y());
    double b2New = 
        ((a.y() + b.y()) / 2) -
        (d /
            (
                sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x())))
                )
            ) *
        (b.x() - a.x());
    if (liesInsideUnscaled( Point(a1New + radius, b1New, 1), Point(a1New - radius, b1New, 1), c))
    {
        return true;
    }
    if (liesInsideUnscaled( Point(a2New + radius, b2New, 1), Point(a2New - radius, b2New, 1), c))
    {
        return true;
    }
    return false;
}


inline void FindHole::insertToHoleBoundary( Edges& holeEdges, Point a, Point b)
{
    Edge e;
    e.source = a;
    e.target = b;
    holeEdges.push_back(e);
}


void FindHole::updateHole( PriorityQueue& pqInnerBoundary, allFaceItr currentFace, allFaceItr neighbor1, allFaceItr neighbor2) 
{
    for (int i = 0; i < 3; ++i) {
        if (neighbor1->neighbor(i) == currentFace)
        {
            neighbor1->set_neighbor(i, this->m_delaunayTriangulation.infinite_face());
        }

        if (neighbor2->neighbor(i) == currentFace)
        {
            neighbor2->set_neighbor(i, this->m_delaunayTriangulation.infinite_face());
        }
    }
    
    this->m_delaunayTriangulation.delete_face(currentFace);
    pqInnerBoundary.push
    (
        node
        (
            CGAL::area(neighbor1->vertex(0)->point(), neighbor1->vertex(1)->point(), neighbor1->vertex(2)->point()),
            neighbor1
        )
    );
    pqInnerBoundary.push
    (
        node
        (
            CGAL::area(neighbor2->vertex(0)->point(), neighbor2->vertex(1)->point(), neighbor2->vertex(2)->point()),
            neighbor2
        )
    );
}


inline bool FindHole::inBoundary( Points exists, Point a)
{
    for (int i = 0; i < exists.size(); i++)
    {
        if (std::abs(exists[i].x() - a.x())<= DBL_EPSILON && 
            std::abs(exists[i].y() - a.y())<= DBL_EPSILON)
        {
            return true;
        }
    }
    return false;
}


bool FindHole::outerBoundaryCondition( Point a, Point b, Point c, faceHandle fh)
{  
    //Diametric circle is non-empty
    if (liesInsideScaled( a, b, c))
    {
        return true;
    }

    //Checking the status of mid point circle and chord circles
    auto distanceAB = sqrt(CGAL::squared_distanceC2(a.x(), a.y(), b.x(), b.y()));
    for (int n = 0; n < 3; ++n) {
        if (isFinite(fh->neighbor(n)))
        {
            if (fh->neighbor(n)->vertex(0)->point() != a &&
                fh->neighbor(n)->vertex(0)->point() != b && 
                fh->neighbor(n)->vertex(0)->point() != c)
            {
                if (liesWithinException(fh->neighbor(n)->vertex(1)->point(), 
                    fh->neighbor(n)->vertex(2)->point(),
                    fh->neighbor(n)->vertex(0)->point(),
                    scaledFactor * distanceAB))
                    return true;
            }
            else
            {
                if (fh->neighbor(n)->vertex(1)->point() != a && 
                    fh->neighbor(n)->vertex(1)->point() != b &&
                    fh->neighbor(n)->vertex(1)->point() != c)
                {
                    if (liesWithinException(fh->neighbor(n)->vertex(0)->point(),
                        fh->neighbor(n)->vertex(2)->point(),
                        fh->neighbor(n)->vertex(1)->point(), 
                        scaledFactor * distanceAB))
                        return true;
                }
                else
                {
                    if (fh->neighbor(n)->vertex(2)->point() != a && 
                        fh->neighbor(n)->vertex(2)->point() != b && 
                        fh->neighbor(n)->vertex(2)->point() != c)
                    {
                        if (liesWithinException(fh->neighbor(n)->vertex(1)->point(), 
                            fh->neighbor(n)->vertex(0)->point(),
                            fh->neighbor(n)->vertex(2)->point(),
                            scaledFactor * distanceAB))
                            return true;
                    }
                }
            }
        }

    }

    return false;
}


bool FindHole::holeDetectionCondition(Point a, Point b, Point c, faceHandle fh)
{

    auto distance = (sqrt(CGAL::squared_distanceC2(a.x(), a.y(), b.x(), b.y())));

    // Diametric circle is non-empty
    if (liesInsideHole(a, b, c))
    {
        return true;
    }

    //Checking the status of mid point circle and chord circles
    for (int n = 0; n < 3; ++n) {
        if (isFinite(fh->neighbor(n)))
        {
            if (fh->neighbor(n)->vertex(0)->point() != a && 
                fh->neighbor(n)->vertex(0)->point() != b &&
                fh->neighbor(n)->vertex(0)->point() != c)
            {
                if (
                    liesWithinExceptionHole
                    (
                        fh->neighbor(n)->vertex(1)->point(),
                        fh->neighbor(n)->vertex(2)->point(),
                        fh->neighbor(n)->vertex(0)->point(),
                        distance
                    )
                   )
                {
                    return true;
                }
            }
            else
            {
                if (fh->neighbor(n)->vertex(1)->point() != a &&
                    fh->neighbor(n)->vertex(1)->point() != b &&
                    fh->neighbor(n)->vertex(1)->point() != c)
                {
                    if (
                        liesWithinExceptionHole
                        (
                            fh->neighbor(n)->vertex(0)->point(),
                            fh->neighbor(n)->vertex(2)->point(),
                            fh->neighbor(n)->vertex(1)->point(),
                            distance
                        )
                       )
                    {
                        return true;
                    }
                }
                else
                {
                    if (fh->neighbor(n)->vertex(2)->point() != a && 
                        fh->neighbor(n)->vertex(2)->point() != b &&
                        fh->neighbor(n)->vertex(2)->point() != c)
                    {
                        if (
                            liesWithinExceptionHole
                            (
                                fh->neighbor(n)->vertex(1)->point(), 
                                fh->neighbor(n)->vertex(0)->point(),
                                fh->neighbor(n)->vertex(2)->point(),
                                distance
                            )
                           )
                        {
                            return true;
                        }
                    }
                }
            }
        }

    }

    return false;
}


void FindHole::outerBoundary(Points& exists, Edges& shape, PriorityQueue& pqOuterBoundary)
{
    for (
        Delaunay::Finite_faces_iterator ffi = this->m_delaunayTriangulation.finite_faces_begin(); 
        ffi != this->m_delaunayTriangulation.finite_faces_end(); 
        ffi++
        )
    {
        if (isInfinite(ffi->neighbor(0)) || 
            isInfinite(ffi->neighbor(1)) || 
            isInfinite(ffi->neighbor(2)))
        {
            if (isInfinite(ffi->neighbor(0)))
            {

                ffi->neighbor(0) = this->m_delaunayTriangulation.infinite_face(); //set the 0th neighbor as infinite face
                exists.push_back(ffi->vertex(1)->point());
                exists.push_back(ffi->vertex(2)->point());
                pqOuterBoundary.push
                (
                    node
                    (
                        sqrt(CGAL::squared_distanceC2(ffi->vertex(1)->point().x(),
                            ffi->vertex(1)->point().y(),
                            ffi->vertex(2)->point().x(),
                            ffi->vertex(2)->point().y())
                        ), 
                        ffi
                    )
                );
            }
            else
            {
                if (isInfinite(ffi->neighbor(1)))
                {
                    ffi->neighbor(1) = this->m_delaunayTriangulation.infinite_face();
                    exists.push_back(ffi->vertex(0)->point());
                    exists.push_back(ffi->vertex(2)->point());
                    pqOuterBoundary.push
                    (
                        node
                        (
                            sqrt(CGAL::squared_distanceC2(ffi->vertex(0)->point().x(),
                                ffi->vertex(0)->point().y(),
                                ffi->vertex(2)->point().x(),
                                ffi->vertex(2)->point().y())),
                            ffi
                        )
                    );
                }
                else
                {
                    if (isInfinite(ffi->neighbor(2)))
                    {
                        ffi->neighbor(2) = this->m_delaunayTriangulation.infinite_face();
                        exists.push_back(ffi->vertex(1)->point());
                        exists.push_back(ffi->vertex(0)->point());
                        pqOuterBoundary.push
                        (
                            node
                            (
                                sqrt(CGAL::squared_distanceC2(ffi->vertex(0)->point().x(),
                                    ffi->vertex(0)->point().y(),
                                    ffi->vertex(1)->point().x(),
                                    ffi->vertex(1)->point().y())),
                                ffi
                            )
                        );
                    }
                }
            }
        }
    }

    while (!pqOuterBoundary.empty())
    {
        allFaceItr currentFace = pqOuterBoundary.top().ffi;
        pqOuterBoundary.pop();

        if (isInfinite(currentFace->neighbor(0)) && 
            isFinite(currentFace->neighbor(1)) && 
            isFinite(currentFace->neighbor(2)))
        {
            if (
                outerBoundaryCondition
                (
                    currentFace->vertex(1)->point(),
                    currentFace->vertex(2)->point(),
                    currentFace->vertex(0)->point(),
                    currentFace
                )
               )/*If any of the circle is non-empty*/
            {
                if (!isInShape(exists, currentFace->vertex(0)->point()))
                {
                    exists.push_back(currentFace->vertex(0)->point());
                    updateNeighbors
                    (
                        pqOuterBoundary,
                        currentFace, currentFace->neighbor(1),
                        currentFace->neighbor(2),
                        sqrt
                        (
                            CGAL::squared_distanceC2
                            (
                                currentFace->vertex(0)->point().x(),
                                currentFace->vertex(0)->point().y(),
                                currentFace->vertex(2)->point().x(),
                                currentFace->vertex(2)->point().y()
                            )
                        ),
                        sqrt
                        (
                            CGAL::squared_distanceC2
                            (
                                currentFace->vertex(1)->point().x(),
                                currentFace->vertex(1)->point().y(),
                                currentFace->vertex(0)->point().x(),
                                currentFace->vertex(0)->point().y()
                            )
                        )
                    );
                }
                else
                {
                    insertToShape(shape, currentFace->vertex(1)->point(), currentFace->vertex(2)->point());
                }
            }
            else
            {
                insertToShape(shape, currentFace->vertex(1)->point(), currentFace->vertex(2)->point());
            }
        }

        else
        {
            if (isInfinite(currentFace->neighbor(1)) && 
                isFinite(currentFace->neighbor(0)) &&
                isFinite(currentFace->neighbor(2)))
            {
                if (
                    outerBoundaryCondition
                    (
                        currentFace->vertex(0)->point(),
                        currentFace->vertex(2)->point(),
                        currentFace->vertex(1)->point(),
                        currentFace
                    )
                   )
                {
                    if (!isInShape(exists, currentFace->vertex(1)->point()))
                    {
                        exists.push_back(currentFace->vertex(1)->point());
                        updateNeighbors
                        (
                            pqOuterBoundary,
                            currentFace, 
                            currentFace->neighbor(2),
                            currentFace->neighbor(0),
                            sqrt
                            (
                                CGAL::squared_distanceC2
                                (
                                    currentFace->vertex(0)->point().x(),
                                    currentFace->vertex(0)->point().y(),
                                    currentFace->vertex(1)->point().x(),
                                    currentFace->vertex(1)->point().y()
                                )
                            ),
                            sqrt(
                                CGAL::squared_distanceC2
                                (
                                    currentFace->vertex(1)->point().x(),
                                    currentFace->vertex(1)->point().y(),
                                    currentFace->vertex(2)->point().x(),
                                    currentFace->vertex(2)->point().y()
                                )
                            )
                        );
                    }
                    else
                    {
                        insertToShape(shape, currentFace->vertex(0)->point(), currentFace->vertex(2)->point());
                    }
                }
                else
                {
                    insertToShape(shape, currentFace->vertex(0)->point(), currentFace->vertex(2)->point());
                }
            }
            else
            {
                if (isInfinite(currentFace->neighbor(2)) &&
                    isFinite(currentFace->neighbor(1)) && 
                    isFinite(currentFace->neighbor(0)))
                {
                    if (
                        outerBoundaryCondition
                        (
                            currentFace->vertex(1)->point(),
                            currentFace->vertex(0)->point(),
                            currentFace->vertex(2)->point(),
                            currentFace
                        )
                       )
                    {
                        if (!isInShape(exists, currentFace->vertex(2)->point()))
                        {
                            exists.push_back(currentFace->vertex(2)->point());
                            updateNeighbors
                            (
                                pqOuterBoundary,
                                currentFace,
                                currentFace->neighbor(1),
                                currentFace->neighbor(0), 
                                sqrt
                                (
                                    CGAL::squared_distanceC2
                                    (
                                        currentFace->vertex(0)->point().x(),
                                        currentFace->vertex(0)->point().y(),
                                        currentFace->vertex(2)->point().x(),
                                        currentFace->vertex(2)->point().y()
                                    )
                                ),
                                sqrt
                                (
                                    CGAL::squared_distanceC2
                                    (
                                        currentFace->vertex(1)->point().x(),
                                        currentFace->vertex(1)->point().y(),
                                        currentFace->vertex(2)->point().x(),
                                        currentFace->vertex(2)->point().y()
                                    )
                                )
                            );
                        }
                        else
                        {
                            insertToShape(shape, currentFace->vertex(0)->point(), currentFace->vertex(1)->point());
                        }
                    }
                    else
                    {
                        insertToShape(shape, currentFace->vertex(0)->point(), currentFace->vertex(1)->point());
                    }
                }
                else
                {
                    if (isInfinite(currentFace->neighbor(2)) &&
                        isInfinite(currentFace->neighbor(1)) && 
                        isFinite(currentFace->neighbor(0)))
                    {
                        insertToShape(shape, currentFace->vertex(2)->point(), currentFace->vertex(0)->point());
                        insertToShape(shape, currentFace->vertex(0)->point(), currentFace->vertex(1)->point());
                    }
                    else
                    {
                        if (isInfinite(currentFace->neighbor(2)) && 
                            isFinite(currentFace->neighbor(1))   && 
                            isInfinite(currentFace->neighbor(0)))
                        {
                            insertToShape(shape, currentFace->vertex(2)->point(), currentFace->vertex(1)->point());
                            insertToShape(shape, currentFace->vertex(1)->point(), currentFace->vertex(0)->point());
                        }
                        else
                        {
                            if (isFinite(currentFace->neighbor(2))   && 
                                isInfinite(currentFace->neighbor(1)) && 
                                isInfinite(currentFace->neighbor(0)))
                            {
                                insertToShape(shape, currentFace->vertex(0)->point(), currentFace->vertex(2)->point());
                                insertToShape(shape, currentFace->vertex(2)->point(), currentFace->vertex(1)->point());
                            }
                        }
                    }
                }
            }
        }
    }

    if (DEBUG_MODE) {
        std::string str = "";

        std::cout << "Boundary Edge:" << std::endl;
        for (const auto& ed : shape) {

            std::string temp = std::to_string(ed.source.x()) + " " + std::to_string(ed.source.y()) + "\n" +
                std::to_string(ed.target.x()) + " " + std::to_string(ed.target.y()) + "\n";
            str = str + temp;
            str = str + "\n\n";

        }
        str = str + std::to_string(shape.at(0).source.x()) + " " + std::to_string(shape.at(0).source.y()) + "\n" +
            std::to_string(shape.at(0).target.x()) + " " + std::to_string(shape.at(0).target.y()) + "\n";
        str = str + "\n\n";
        std::string filename = "E:\\Avinash C drive\\Downloads\\gp610-20230826-win64-mingw\\gnuplot\\bin\\sampleBoundary.txt";
        std::ofstream file(filename);

        try {
            if (!file.is_open())
                throw "Error to Open File !!!" + filename;
            if (str.empty())
                throw "Write the file atleast one char !!!";
            file.write(str.c_str(), str.length());
            file.close();

            std::cout << "\n outerboundary file saved ! " << std::endl;

            std::cout << "No of boundary vertices:" << shape.size() << std::endl;
        }
        catch (...) {

            std::cout << "ERROR!!! 620 SampleCodeGeneration1.cpp" << std::endl;
        }
    }
}


void FindHole::holeDetection(
    Points exists, 
    Edges& shape,
    Edges& holeEdges, 
    PriorityQueue& pqInnerBoundary,
    std::vector< std::vector<Point>>& holes,
    int& numberOfHoles
)
{
    Delaunay::Finite_faces_iterator ffi = this->m_delaunayTriangulation.finite_faces_begin();

    Delaunay::Face_handle highestAreaTriangleHandle;//the face handle for the highest area triangle
    Delaunay::Face_handle faceHandle0, faceHandle1, faceHandle2;//facehandles for neighbouring triangles
    bool inHole;
    int tempVar = 0;//variable for for multiple hole detection
    bool zerothPointFaceHandle0, firstPointFaceHandle0, secondPointFaceHandle0;
    bool zerothPointFaceHandle1, firstPointFaceHandle1, secondPointFaceHandle1;
    bool zerothPointFaceHandle2, firstPointFaceHandle2, secondPointFaceHandle2;
    bool inOuterBdry = false;
    Delaunay::Face_handle faceHandle[3];//facehandles for neighbouring triangles
    bool PointFaceHandle[3][3]{ false };//PointFaceHandle[nth face][nth point]
    std::vector<Point> hole;

    //do while loop for multiple hole detection
    do //do of while ((count1-tempVar)>4);
    {
        std::vector<Point> hole;


        tempVar = holeEdges.size(); //variable for multiple hole detection
        inHole = false;//variable to check whether already the edge is in OHL
        zerothPointFaceHandle0 = false, firstPointFaceHandle0 = false, secondPointFaceHandle0 = false;
        zerothPointFaceHandle1 = false, firstPointFaceHandle1 = false, secondPointFaceHandle1 = false;
        zerothPointFaceHandle2 = false, firstPointFaceHandle2 = false, secondPointFaceHandle2 = false;
        inOuterBdry = false;

        highestAreaTriangleHandle = getMaxAreaTriangleHandle( shape, ffi);
        if (highestAreaTriangleHandle == this->m_delaunayTriangulation.infinite_face())
        {
            if (holeEdges.empty())
            {
                exit(0);
            }
            else
            {
                break;
            }
        }

        //checking any of the vertices of Highest Area Triangle are in the hole, if so do not put that it to Q
        for (int i = 0; i < holeEdges.size(); i++)
        {
            if (
                (
                    std::abs(highestAreaTriangleHandle->vertex(0)->point().x() - holeEdges[i].source.x()) <= DBL_EPSILON &&
                    std::abs(highestAreaTriangleHandle->vertex(0)->point().y() - holeEdges[i].source.y()) <= DBL_EPSILON
                )
                ||
                (
                    std::abs (highestAreaTriangleHandle->vertex(0)->point().x() - holeEdges[i].target.x()) <= DBL_EPSILON &&
                    std::abs (highestAreaTriangleHandle->vertex(0)->point().y() - holeEdges[i].target.y()) <= DBL_EPSILON
                )
                ||
                (
                    std::abs(highestAreaTriangleHandle->vertex(1)->point().x() - holeEdges[i].source.x()) <= DBL_EPSILON &&
                    std::abs(highestAreaTriangleHandle->vertex(1)->point().y() - holeEdges[i].source.y()) <= DBL_EPSILON
                )
                ||
                (
                    std::abs(highestAreaTriangleHandle->vertex(1)->point().x() - holeEdges[i].target.x()) <= DBL_EPSILON &&
                    std::abs(highestAreaTriangleHandle->vertex(1)->point().y() - holeEdges[i].target.y()) <= DBL_EPSILON
                )
                ||
                (
                    std::abs(highestAreaTriangleHandle->vertex(2)->point().x() - holeEdges[i].source.x()) <= DBL_EPSILON &&
                    std::abs(highestAreaTriangleHandle->vertex(2)->point().y() - holeEdges[i].source.y()) <= DBL_EPSILON
                 )
                ||
                (
                    std::abs(highestAreaTriangleHandle->vertex(2)->point().x() - holeEdges[i].target.x()) <= DBL_EPSILON &&
                    std::abs(highestAreaTriangleHandle->vertex(2)->point().y() - holeEdges[i].target.y()) <= DBL_EPSILON
                    )
                )
            {
                inHole = true;
                if (DEBUG_MODE)
                    std::cout << "line 1131 break" << std::endl;
                break;
            }
        }

        //checking any of the vertices of Highest Area Triangle are in the outerbdry, if so do not put that it to Q
        for (int i = 0; i < shape.size(); i++)
        {
            if (
                (
                    std::abs(highestAreaTriangleHandle->vertex(0)->point().x() - shape[i].source.x()) <= DBL_EPSILON && 
                    std::abs(highestAreaTriangleHandle->vertex(0)->point().y() - shape[i].source.y()) <= DBL_EPSILON
                    )
                ||
                (
                    std::abs(highestAreaTriangleHandle->vertex(0)->point().x() - shape[i].target.x()) <= DBL_EPSILON &&
                    std::abs(highestAreaTriangleHandle->vertex(0)->point().y() - shape[i].target.y()) <= DBL_EPSILON
                    )
                ||
                (
                    std::abs(highestAreaTriangleHandle->vertex(1)->point().x() - shape[i].source.x()) <= DBL_EPSILON && 
                    std::abs(highestAreaTriangleHandle->vertex(1)->point().y() - shape[i].source.y()) <= DBL_EPSILON
                    ) 
                ||
                (
                    std::abs(highestAreaTriangleHandle->vertex(1)->point().x() - shape[i].target.x()) <= DBL_EPSILON && 
                    std::abs(highestAreaTriangleHandle->vertex(1)->point().y() - shape[i].target.y()) <= DBL_EPSILON
                    )
                ||
                (
                    std::abs(highestAreaTriangleHandle->vertex(2)->point().x() - shape[i].source.x()) <= DBL_EPSILON &&
                    std::abs(highestAreaTriangleHandle->vertex(2)->point().y() - shape[i].source.y()) <= DBL_EPSILON
                    ) 
                ||
                (
                    std::abs(highestAreaTriangleHandle->vertex(2)->point().x() - shape[i].target.x()) <= DBL_EPSILON &&
                    std::abs(highestAreaTriangleHandle->vertex(2)->point().y() - shape[i].target.y()) <= DBL_EPSILON
                    )
                )
            {
                inOuterBdry = true;
                if (DEBUG_MODE)
                    std::cout << "line 964 break" << std::endl;
                break;
            }
        }
        
        if (inOuterBdry)
        {
            if (DEBUG_MODE)
                std::cout << "highest area triangle is outer bdry triangle and it is an invalid input" << std::endl;
            if (holeEdges.empty())
            {
                if (DEBUG_MODE)
                    std::cout << "line 629 exit" << std::endl;
                exit(0);
            }
            else
            {
                //break;
            }
        }

        if ((!inOuterBdry) && (!inHole))
        {
            // if the hole is not along the outer bdry, initialize the queue1 with the highest area traingle
            // the next three if loops are to find out the vertex number of the third vertex w.r.t 
            // facehandle0 when v1, v2, 
            // facehandle1 when v2, v0 and
            // facehandle2 when v0, v1 is processed. 
            // this is to use it in lies inside condition checking for dimeter-circle
            // to get the unvisited vertex of the neighboring triangle of 2nd vertex for cheking lies inside condition

            for (int f = 0; f < 3; ++f) {

                faceHandle[f] = highestAreaTriangleHandle->neighbor(f);
                pqInnerBoundary.push
                (
                    node
                    (
                        area(faceHandle[f]->vertex(0)->point(), faceHandle[f]->vertex(1)->point(), faceHandle[f]->vertex(2)->point()),
                        faceHandle[f]
                    )
                );

                for (int v = 0; v < 3; ++v) {
                    if (
                        (
                            std::abs(faceHandle[f]->vertex(v)->point().x() - highestAreaTriangleHandle->vertex((f + 1) % 3)->point().x()) <= DBL_EPSILON &&
                            std::abs(faceHandle[f]->vertex(v)->point().y() - highestAreaTriangleHandle->vertex((f + 1) % 3)->point().y()) <= DBL_EPSILON
                            ) 
                        ||
                        (
                            std::abs(faceHandle[f]->vertex(v)->point().x() - highestAreaTriangleHandle->vertex((f + 2) % 3)->point().x()) <= DBL_EPSILON &&
                            std::abs(faceHandle[f]->vertex(v)->point().y() - highestAreaTriangleHandle->vertex((f + 2) % 3)->point().y()) <= DBL_EPSILON
                            )
                        )
                    {
                        PointFaceHandle[f][v] = false;
                    }
                    else
                    {
                        PointFaceHandle[f][v] = true;
                    }

                    if (PointFaceHandle[f][v])
                    {
                        faceHandle[f]->set_neighbor(v, this->m_delaunayTriangulation.infinite_face());
                    }
                }
            }




            exists.push_back(highestAreaTriangleHandle->vertex(0)->point());
            exists.push_back(highestAreaTriangleHandle->vertex(1)->point());
            exists.push_back(highestAreaTriangleHandle->vertex(2)->point());
            this->m_delaunayTriangulation.delete_face(highestAreaTriangleHandle);
        }

        while (!pqInnerBoundary.empty())
        {
            highestAreaTriangleHandle = pqInnerBoundary.top().ffi;
            pqInnerBoundary.pop();

            for (int n = 0; n < 3; ++n) {
                if (isInfinite(highestAreaTriangleHandle->neighbor(n)) && 
                    isFinite(highestAreaTriangleHandle->neighbor((n + 1) % 3)) &&
                    isFinite(highestAreaTriangleHandle->neighbor((n + 2) % 3)))
                {
                    if (
                        holeDetectionCondition
                        (
                            highestAreaTriangleHandle->vertex(n)->point(),
                            highestAreaTriangleHandle->vertex((n + 1) % 3)->point(), 
                            highestAreaTriangleHandle->vertex((n + 2) % 3)->point(), 
                            highestAreaTriangleHandle
                        )
                       )
                    {
                        if (!inBoundary(exists, highestAreaTriangleHandle->vertex(n)->point()))
                        {
                            exists.push_back(highestAreaTriangleHandle->vertex(n)->point());
                            updateHole
                            (
                                pqInnerBoundary,
                                highestAreaTriangleHandle,
                                highestAreaTriangleHandle->neighbor((n + 1) % 3),
                                highestAreaTriangleHandle->neighbor((n + 2) % 3)
                            );
                        }
                        else
                        {
                            insertToHoleBoundary
                            (
                                holeEdges,
                                highestAreaTriangleHandle->vertex((n + 1) % 3)->point(),
                                highestAreaTriangleHandle->vertex((n + 2) % 3)->point()
                            );
                            hole.push_back(highestAreaTriangleHandle->vertex((n + 1) % 3)->point());
                            hole.push_back(highestAreaTriangleHandle->vertex((n + 2) % 3)->point());
                        }
                    }
                    else
                    {
                        insertToHoleBoundary
                        (
                            holeEdges,
                            highestAreaTriangleHandle->vertex((n + 1) % 3)->point(),
                            highestAreaTriangleHandle->vertex((n + 2) % 3)->point()
                        );
                        hole.push_back(highestAreaTriangleHandle->vertex((n + 1) % 3)->point());
                        hole.push_back(highestAreaTriangleHandle->vertex((n + 2) % 3)->point());
                    }
                }
                else
                {
                    if (isInfinite(highestAreaTriangleHandle->neighbor(2)) &&
                        isInfinite(highestAreaTriangleHandle->neighbor(1)) &&
                        isFinite(highestAreaTriangleHandle->neighbor(0)))
                    {
                        insertToHoleBoundary
                        (
                            holeEdges, 
                            highestAreaTriangleHandle->vertex(2)->point(),
                            highestAreaTriangleHandle->vertex(0)->point()
                        );
                        hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                        insertToHoleBoundary
                        (
                            holeEdges, 
                            highestAreaTriangleHandle->vertex(0)->point(),
                            highestAreaTriangleHandle->vertex(1)->point()
                        );
                        hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                        hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                    }
                    else if (isInfinite(highestAreaTriangleHandle->neighbor(2)) &&
                        isFinite(highestAreaTriangleHandle->neighbor(1)) &&
                        isInfinite(highestAreaTriangleHandle->neighbor(0)))
                    {
                        insertToHoleBoundary
                        (
                            holeEdges,
                            highestAreaTriangleHandle->vertex(2)->point(),
                            highestAreaTriangleHandle->vertex(1)->point()
                        );
                        hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                        hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                        insertToHoleBoundary
                        (
                            holeEdges,
                            highestAreaTriangleHandle->vertex(1)->point(),
                            highestAreaTriangleHandle->vertex(0)->point()
                        );
                        hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    }
                    else if (isFinite(highestAreaTriangleHandle->neighbor(2)) && 
                        isInfinite(highestAreaTriangleHandle->neighbor(1)) && 
                        isInfinite(highestAreaTriangleHandle->neighbor(0)))
                    {
                        insertToHoleBoundary
                        (
                            holeEdges,
                            highestAreaTriangleHandle->vertex(0)->point(),
                            highestAreaTriangleHandle->vertex(2)->point()
                        );
                        hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                        insertToHoleBoundary
                        (
                            holeEdges, 
                            highestAreaTriangleHandle->vertex(2)->point(),
                            highestAreaTriangleHandle->vertex(1)->point()
                        );
                        hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                        hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                    }
                }
            }
        }
        
        
        
        
        if (hole.size() > 4) {
            holes.push_back(hole);
            ++numberOfHoles;
        }


    } while ((holeEdges.size() - 1 - tempVar) > 4);

    if (DEBUG_MODE)
    std::cout << "holedetection :No of Holes:" << numberOfHoles << "," << holeEdges.size() << std::endl;
}


void FindHole::singleHoleDetection( Edges& shape, Edges& holeEdges, std::vector< std::vector<Point>>& holes, int& numberOfHoles)
{
    Points  exists;
    PriorityQueue pqOuterBoundary;
    PriorityQueue pqInnerBoundary;

    Delaunay::Finite_faces_iterator ffi = this->m_delaunayTriangulation.finite_faces_begin();
    Delaunay::Face_handle highestAreaTriangleHandle;//the face handle for the highest area triangle
    bool inHole = false;//variable to check whether already the edge is in Outer Hole List
    int tempVar = holeEdges.size();//variable for for multiple hole detection
    bool inOuterBdry = false;
    Delaunay::Face_handle faceHandle[3];//facehandles for neighbouring triangles
    bool PointFaceHandle[3][3]{false};//PointFaceHandle[nth face][nth point]
    std::vector<Point> hole;


    highestAreaTriangleHandle = getMaxAreaTriangleHandle( shape, ffi);

    if (highestAreaTriangleHandle == this->m_delaunayTriangulation.infinite_face())
    {
        if (holeEdges.size() == 0)
        {
            if (DEBUG_MODE)
            std::cout << "line 1110 exit" << std::endl;
            exit(0);
        }
        else
        {
            //break;
        }
    }

    //checking any of the vertices of Highest Area Triangle are in the hole, if so do not put that it to Q
    for (int i = 0; i < holeEdges.size(); i++)
    {
        if (
            (
                std::abs(highestAreaTriangleHandle->vertex(0)->point().x() - holeEdges[i].source.x())<=DBL_EPSILON &&
                std::abs(highestAreaTriangleHandle->vertex(0)->point().y() - holeEdges[i].source.y())<=DBL_EPSILON
                ) 
            ||
            (
                std::abs(highestAreaTriangleHandle->vertex(0)->point().x() - holeEdges[i].target.x())<=DBL_EPSILON &&
                std::abs(highestAreaTriangleHandle->vertex(0)->point().y() - holeEdges[i].target.y())<=DBL_EPSILON
                ) 
            ||
            (
                std::abs(highestAreaTriangleHandle->vertex(1)->point().x() - holeEdges[i].source.x())<=DBL_EPSILON && 
                std::abs(highestAreaTriangleHandle->vertex(1)->point().y() - holeEdges[i].source.y())<=DBL_EPSILON
                )
            ||
            (
                std::abs(highestAreaTriangleHandle->vertex(1)->point().x() - holeEdges[i].target.x())<=DBL_EPSILON &&
                std::abs(highestAreaTriangleHandle->vertex(1)->point().y() - holeEdges[i].target.y())<=DBL_EPSILON
                ) 
            ||
            (
                std::abs(highestAreaTriangleHandle->vertex(2)->point().x() - holeEdges[i].source.x())<=DBL_EPSILON &&
                std::abs(highestAreaTriangleHandle->vertex(2)->point().y() - holeEdges[i].source.y())<=DBL_EPSILON
                ) 
            ||
            (
                std::abs(highestAreaTriangleHandle->vertex(2)->point().x() - holeEdges[i].target.x())<=DBL_EPSILON && 
                std::abs(highestAreaTriangleHandle->vertex(2)->point().y() - holeEdges[i].target.y())<=DBL_EPSILON
                )
            )
        {
            inHole = true;
            if (DEBUG_MODE)
            std::cout << "line 1131 break" << std::endl;
            break;
        }
    }

    //checking any of the vertices of Highest Area Triangle are in the outerbdry, if so do not put that it to Q
    for (int i = 0; i < shape.size(); i++)
    {
        if (
            (
                std::abs(highestAreaTriangleHandle->vertex(0)->point().x() - shape[i].source.x())<=DBL_EPSILON &&
                std::abs(highestAreaTriangleHandle->vertex(0)->point().y() - shape[i].source.y())<=DBL_EPSILON
                ) 
            ||
            (
                std::abs(highestAreaTriangleHandle->vertex(0)->point().x() - shape[i].target.x())<=DBL_EPSILON &&
                std::abs(highestAreaTriangleHandle->vertex(0)->point().y() - shape[i].target.y())<=DBL_EPSILON
                ) 
            ||
            (
                std::abs(highestAreaTriangleHandle->vertex(1)->point().x() - shape[i].source.x())<=DBL_EPSILON && 
                std::abs(highestAreaTriangleHandle->vertex(1)->point().y() - shape[i].source.y())<=DBL_EPSILON
                ) 
            ||
            (
                std::abs(highestAreaTriangleHandle->vertex(1)->point().x() - shape[i].target.x())<=DBL_EPSILON &&
                std::abs(highestAreaTriangleHandle->vertex(1)->point().y() - shape[i].target.y())<=DBL_EPSILON
                ) 
            ||
            (
                std::abs(highestAreaTriangleHandle->vertex(2)->point().x() - shape[i].source.x())<=DBL_EPSILON &&
                std::abs(highestAreaTriangleHandle->vertex(2)->point().y() - shape[i].source.y())<=DBL_EPSILON
                )
            ||
            (
                std::abs(highestAreaTriangleHandle->vertex(2)->point().x() - shape[i].target.x())<=DBL_EPSILON &&
                std::abs(highestAreaTriangleHandle->vertex(2)->point().y() - shape[i].target.y())<=DBL_EPSILON
                )
            )
        {
            inOuterBdry = true;
            if (DEBUG_MODE)
            std::cout << "line 964 break" << std::endl;
            break;
        }
    }

    if (inOuterBdry)
    {
        if(DEBUG_MODE)
        std::cout << "highest area triangle is outer bdry triangle and it is an invalid input" << std::endl;
        if (holeEdges.empty())
        {
            if (DEBUG_MODE)
            std::cout << "line 1160 exit" << std::endl;
            exit(0);
        }
        else
        {
            //break;
        }
    }
    
    if ((!inOuterBdry) && (!inHole))
    {
        // if the hole is not along the outer bdry, initialize the queue1 with the highest area traingle
        // the next three if loops are to find out the vertex number of the third vertex w.r.t 
        // facehandle0 when v1, v2, 
        // facehandle1 when v2, v0 and
        // facehandle2 when v0, v1 is processed. 
        // this is to use it in lies inside condition checking for dimeter-circle
        // to get the unvisited vertex of the neighboring triangle of 2nd vertex for cheking lies inside condition

        for (int f = 0; f < 3; ++f) {

            faceHandle[f] = highestAreaTriangleHandle->neighbor(f);
            pqInnerBoundary.push(node(area(faceHandle[f]->vertex(0)->point(), faceHandle[f]->vertex(1)->point(), faceHandle[f]->vertex(2)->point()), faceHandle[f]));

            for (int v = 0; v < 3; ++v) {
                if ((
                    (std::abs(faceHandle[f]->vertex(v)->point().x() - highestAreaTriangleHandle->vertex((f + 1) % 3)->point().x()) <= DBL_EPSILON) &&
                    ( std::abs(faceHandle[f]->vertex(v)->point().y() - highestAreaTriangleHandle->vertex((f + 1) % 3)->point().y()) <= DBL_EPSILON)
                    )
                    ||
                    (
                     (std::abs(faceHandle[f]->vertex(v)->point().x() - highestAreaTriangleHandle->vertex((f + 2) % 3)->point().x()) <= DBL_EPSILON) &&
                     (std::abs(faceHandle[f]->vertex(v)->point().y() - highestAreaTriangleHandle->vertex((f + 2) % 3)->point().y()) <= DBL_EPSILON)
                    )
                    )
                {
                    PointFaceHandle[f][v] = false;
                }
                else
                {
                    PointFaceHandle[f][v] = true;
                }

                if (PointFaceHandle[f][v])
                {
                    faceHandle[f]->set_neighbor(v, this->m_delaunayTriangulation.infinite_face());
                }
            }
        }


        

        exists.push_back(highestAreaTriangleHandle->vertex(0)->point());
        exists.push_back(highestAreaTriangleHandle->vertex(1)->point());
        exists.push_back(highestAreaTriangleHandle->vertex(2)->point());
        this->m_delaunayTriangulation.delete_face(highestAreaTriangleHandle);
    }

    while (!pqInnerBoundary.empty())
    {
        highestAreaTriangleHandle = pqInnerBoundary.top().ffi;
        pqInnerBoundary.pop();

        for (int n = 0; n < 3; ++n) {
            if (isInfinite(highestAreaTriangleHandle->neighbor(n)) &&
                isFinite(highestAreaTriangleHandle->neighbor((n + 1) % 3)) && 
                isFinite(highestAreaTriangleHandle->neighbor((n + 2) % 3)))
            {
                if (
                    holeDetectionCondition
                    (
                        highestAreaTriangleHandle->vertex(n)->point(), 
                        highestAreaTriangleHandle->vertex((n + 1) % 3)->point(),
                        highestAreaTriangleHandle->vertex((n + 2) % 3)->point(),
                        highestAreaTriangleHandle
                    )
                   )
                {
                    if (!inBoundary(exists, highestAreaTriangleHandle->vertex(n)->point()))
                    {
                        exists.push_back(highestAreaTriangleHandle->vertex(n)->point());
                        updateHole
                        (
                            pqInnerBoundary,
                            highestAreaTriangleHandle,
                            highestAreaTriangleHandle->neighbor((n + 1) % 3), 
                            highestAreaTriangleHandle->neighbor((n + 2) % 3)
                        );
                    }
                    else
                    {
                        insertToHoleBoundary
                        (
                            holeEdges, 
                            highestAreaTriangleHandle->vertex((n + 1) % 3)->point(), 
                            highestAreaTriangleHandle->vertex((n + 2) % 3)->point()
                        );
                        hole.push_back(highestAreaTriangleHandle->vertex((n + 1) % 3)->point());
                        hole.push_back(highestAreaTriangleHandle->vertex((n + 2) % 3)->point());
                    }
                }
                else
                {
                    insertToHoleBoundary
                    (
                        holeEdges,
                        highestAreaTriangleHandle->vertex((n + 1) % 3)->point(),
                        highestAreaTriangleHandle->vertex((n + 2) % 3)->point()
                    );
                    hole.push_back(highestAreaTriangleHandle->vertex((n + 1) % 3)->point());
                     hole.push_back(highestAreaTriangleHandle->vertex((n + 2) % 3)->point());
                }
            }
            else
            {
                if (isInfinite(highestAreaTriangleHandle->neighbor(2)) && 
                    isInfinite(highestAreaTriangleHandle->neighbor(1)) &&
                    isFinite(highestAreaTriangleHandle->neighbor(0)))
                {
                    insertToHoleBoundary
                    (
                        holeEdges,
                        highestAreaTriangleHandle->vertex(2)->point(), 
                        highestAreaTriangleHandle->vertex(0)->point()
                    );
                    hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                    insertToHoleBoundary
                    (
                        holeEdges,
                        highestAreaTriangleHandle->vertex(0)->point(),
                        highestAreaTriangleHandle->vertex(1)->point()
                    );
                    hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                }
                else if (
                    isInfinite(highestAreaTriangleHandle->neighbor(2)) &&
                    isFinite(highestAreaTriangleHandle->neighbor(1)) &&
                    isInfinite(highestAreaTriangleHandle->neighbor(0)))
                {
                    insertToHoleBoundary
                    (
                        holeEdges,
                        highestAreaTriangleHandle->vertex(2)->point(),
                        highestAreaTriangleHandle->vertex(1)->point()
                    );
                    hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                    insertToHoleBoundary
                    (
                        holeEdges,
                        highestAreaTriangleHandle->vertex(1)->point(),
                        highestAreaTriangleHandle->vertex(0)->point()
                    );
                    hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                }
                else if
                    (
                        isFinite(highestAreaTriangleHandle->neighbor(2)) && 
                        isInfinite(highestAreaTriangleHandle->neighbor(1)) &&
                        isInfinite(highestAreaTriangleHandle->neighbor(0))
                    )
                {
                    insertToHoleBoundary
                    (
                        holeEdges,
                        highestAreaTriangleHandle->vertex(0)->point(),
                        highestAreaTriangleHandle->vertex(2)->point()
                    );
                    hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    insertToHoleBoundary
                    (
                        holeEdges,
                        highestAreaTriangleHandle->vertex(2)->point(),
                        highestAreaTriangleHandle->vertex(1)->point()
                    );
                    hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                }
            }
        }
    }

    holes.push_back(hole);
    numberOfHoles++;

    if(DEBUG_MODE)
    std::cout << "No of Holes:" << numberOfHoles << std::endl;
}


int indxPoint(const std::vector<Point>& points, Point p) {
    int inx = 0;
    int cnt = 0;
    for (auto a : points) {
        auto d = sqrt(CGAL::squared_distanceC2(p.x(), p.y(), a.x(), a.y()));
        if (d <= DBL_EPSILON) {
            inx = cnt;
            break;
        }
        ++cnt;
    }
    return inx;
}


void FindHole::indicesOfHoles(std::vector<std::vector<Point>>& holes,std::vector<std::vector<int>>& holeIndices) {
    
    // k-d tree
    Tree tree(this->m_points.begin(), this->m_points.end());

    for (auto nthHole : holes) {
        std::vector<int> indicesOfHole;
        for (auto pointOfHole : nthHole) {
            auto index = indxPoint(this->m_points,pointOfHole);
            indicesOfHole.emplace_back(index);

        }
        holeIndices.emplace_back(indicesOfHole);
    }

}


int FindHole::getholes(std::vector<std::vector<Point>>& holes, std::vector<std::vector<int>>& holeIndices) {

    
    PriorityQueue pqOuterBoundary;
    PriorityQueue pqInnerBoundary;
    Edges shape;
    Edges holeEdges;


    double minX = std::numeric_limits<int>::max(), maxX = std::numeric_limits<int>::min();
    double minY = std::numeric_limits<int>::max(), maxY = std::numeric_limits<int>::min();
    Points inputPoints, exists;

    int numberOfHoles = 0;

    if(DEBUG_MODE)
        saveDelaunayTriangulationTxtFile(this->m_delaunayTriangulation, this->m_points);

    for (
        Delaunay::Vertex_iterator vi = this->m_delaunayTriangulation.vertices_begin(); 
        vi != this->m_delaunayTriangulation.vertices_end();
        vi++
        )
    {
        if (vi->point().x() > maxX)
        {
            maxX = vi->point().x();
        }
        if (vi->point().y() > maxY)
        {
            maxY = vi->point().y();
        }
        if (vi->point().x() < minX)
        {
            minX = vi->point().x();
        }
        if (vi->point().y() < minY)
        {
            minY = vi->point().y();
        }

        inputPoints.push_back(vi->point());
    }

    //outer boundary
    FindHole::outerBoundary(exists, shape, pqOuterBoundary);

    //multiple holes
    //FindHole::holeDetection(exists, shape, holeEdges, pqInnerBoundary, holes, numberOfHoles);

    //single hole
    FindHole::singleHoleDetection(shape, holeEdges, holes, numberOfHoles);
    indicesOfHoles(holes, holeIndices);
    
    if (DEBUG_MODE)
        saveDelaunayTriangulationTxtFile(this->m_delaunayTriangulation, this->m_points);

    return numberOfHoles;

}