#include "FindHole.h"
float v = 1;

FindHole::FindHole(const std::vector<Point>&points):m_points(points)
{
    this->m_delaunayTriangulation.insert(points.begin(), points.end());
}


FindHole::~FindHole()
{
}


float FindHole::area(Point a, Point b, Point c)
{
    //compute the area of the triangle formed by the points a, b, and c
    return (float)std::abs(0.5 * (a.x() * (b.y() - c.y()) + b.x() * (c.y() - a.y()) + c.x() * (a.y() - b.y())));
}


double FindHole::distance(Point a, Point b)
{
    //compute the distance between two points
    return (sqrt(std::abs(((a.x() - b.x()) * (a.x() - b.x()))) + std::abs(((a.y() - b.y()) * (a.y() - b.y())))));
}


bool FindHole::isInfinite( allFaceItr afi)
{
    //Checking whether the face is infinite or not
    if (afi->vertex(0) == m_delaunayTriangulation.infinite_vertex() || afi->vertex(1) == m_delaunayTriangulation.infinite_vertex() || afi->vertex(2) == m_delaunayTriangulation.infinite_vertex())
        return true;
    return false;
}


bool FindHole::isFinite( allFaceItr afi) 
{
    //Checking whether the face is finite or not
    if (afi->vertex(0) != m_delaunayTriangulation.infinite_vertex() && afi->vertex(1) != m_delaunayTriangulation.infinite_vertex() && afi->vertex(2) != m_delaunayTriangulation.infinite_vertex())
        return true;
    return false;
}


bool FindHole::liesInsideScaled(Point a, Point b, Point c) /*To check whether p3 is inside scaled diametric circle(p1,p2)*/
{
    if (distance(Point((a.x() + b.x()) / 2, (a.y() + b.y()) / 2, 1), c) <= v * distance(Point((a.x() + b.x()) / 2, (a.y() + b.y()) / 2, 1), a))
        return true;
    return false;
}


bool FindHole::liesInsideUnscaled( Point a, Point b, Point c) /*To check whether p3 is inside unscaled diametric circle(p1,p2)*/
{
    if (distance(Point((a.x() + b.x()) / 2, (a.y() + b.y()) / 2, 1), c) <= distance(Point((a.x() + b.x()) / 2, (a.y() + b.y()) / 2, 1), a))
        return true;
    return false;
}


bool FindHole::liesWithinException(Point a, Point b, Point c, double radius)/*Checking whether mid point and chord circles are empty*/
{
    if (distance(a, b) >= radius)
    {
        Point p1 = Point(((a.x() + b.x()) / 2 - (radius / 2)), ((a.y() + b.y()) / 2), 1);
        Point p2 = Point(((a.x() + b.x()) / 2 + (radius / 2)), ((a.y() + b.y()) / 2), 1); /*(p1,p2) is the midpoint of edge (a,b)*/
        if (liesInsideScaled( p1, p2, c))
            return true;
        return false;
    }
    radius = (radius / 2) * v;
    double d = sqrt(std::abs((radius * radius) - ((distance(a, b) / 2) * (distance(a, b) / 2))));
    double a1New = ((a.x() + b.x()) / 2) + (d / (sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x()))))) * (a.y() - b.y());
    double b1New = ((a.y() + b.y()) / 2) + (d / (sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x()))))) * (b.x() - a.x());
    double a2New = ((a.x() + b.x()) / 2) - (d / (sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x()))))) * (a.y() - b.y());
    double b2New = ((a.y() + b.y()) / 2) - (d / (sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x()))))) * (b.x() - a.x());

    if (liesInsideUnscaled( Point(a1New + radius, b1New, 1), Point(a1New - radius, b1New, 1), c))
        return true;
    if (liesInsideUnscaled( Point(a2New + radius, b2New, 1), Point(a2New - radius, b2New, 1), c))
        return true;
    return false;
}


bool FindHole::isInShape( Points& exists, Point a)/*Checking whether a point is already in shape*/
{
    for (int i = 0; i < exists.size(); i++)
    {

        if (exists[i].x() == a.x() && exists[i].y() == a.y())
            return true;
    }
    return false;
}


void FindHole::updateNeighbors( PriorityQueue& pqOuterBoundary, allFaceItr currentFace, allFaceItr neighbor1, allFaceItr neighbor2, float d1, float d2)
{
    if (neighbor1->neighbor(0) == currentFace)
    {
        neighbor1->set_neighbor(0, this->m_delaunayTriangulation.infinite_face());
    }
    else
    {
        if (neighbor1->neighbor(1) == currentFace)
        {
            neighbor1->set_neighbor(1, this->m_delaunayTriangulation.infinite_face());
        }
        else
        {
            if (neighbor1->neighbor(2) == currentFace)
            {
                neighbor1->set_neighbor(2, this->m_delaunayTriangulation.infinite_face());
            }
        }
    }

    if (neighbor2->neighbor(0) == currentFace)
    {
        neighbor2->set_neighbor(0, this->m_delaunayTriangulation.infinite_face());
    }
    else
    {
        if (neighbor2->neighbor(1) == currentFace)
        {
            neighbor2->set_neighbor(1, this->m_delaunayTriangulation.infinite_face());
        }
        else
        {
            if (neighbor2->neighbor(2) == currentFace)
            {
                neighbor2->set_neighbor(2, this->m_delaunayTriangulation.infinite_face());
            }
        }
    }


    this->m_delaunayTriangulation.delete_face(currentFace);
    pqOuterBoundary.pqInsert(d1, neighbor1);
    pqOuterBoundary.pqInsert(d2, neighbor2);
}


void FindHole::insertToShape( Edges& shape, Point a, Point b)/*Inserting a new edge into boundary */
{
    Edge e;
    e.source = a;
    e.target = b;
    shape.push_back(e);
}


faceHandle FindHole::getMaxAreaTriangleHandle( Edges& shape, finitFaceItr ffi) /*returns the largest triangle*/
{
    double biggestArea = 0.0;
    bool presentInShape = false;
    Delaunay::Face_handle highestAreaTriangleHandle;
    ffi = this->m_delaunayTriangulation.finite_faces_begin();
    int fl = 0;

    for (Delaunay::Finite_faces_iterator ffi = this->m_delaunayTriangulation.finite_faces_begin(); ffi != this->m_delaunayTriangulation.finite_faces_end(); ffi++)
    {
        Point a = ffi->vertex(0)->point();
        Point b = ffi->vertex(1)->point();
        Point c = ffi->vertex(2)->point();

        for (int i = 0; i < shape.size(); i++)
        {
            if ((a.x() == shape[i].source.x() && a.y() == shape[i].source.y()) || (a.x() == shape[i].target.x() && a.y() == shape[i].target.y()) ||
                (b.x() == shape[i].source.x() && b.y() == shape[i].source.y()) || (b.x() == shape[i].target.x() && b.y() == shape[i].target.y()) ||
                (c.x() == shape[i].source.x() && c.y() == shape[i].source.y()) || (c.x() == shape[i].target.x() && c.y() == shape[i].target.y()))
                presentInShape = true;
        }

        if (!presentInShape)
        {
            if (isFinite( ffi->neighbor(0)) && isFinite( ffi->neighbor(1)) && isFinite( ffi->neighbor(2)))
            {
                fl = 1;
                if (area(ffi->vertex(0)->point(), ffi->vertex(1)->point(), ffi->vertex(2)->point()) >= biggestArea)
                {
                    biggestArea = area(ffi->vertex(0)->point(), ffi->vertex(1)->point(), ffi->vertex(2)->point());
                    highestAreaTriangleHandle = ffi;
                }
            }
        }
        presentInShape = false;
    }
    if (fl == 0)
    {
        highestAreaTriangleHandle = this->m_delaunayTriangulation.infinite_face();
    }
    return highestAreaTriangleHandle;
}


bool FindHole::liesInsideHole( Point p1, Point p2, Point p3)
{
    if (distance(Point((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2, 1), p3) <= distance(Point((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2, 1), p1))
    {
        return true;
    }
    return 0;
}


bool FindHole::liesWithinExceptionHole( Point a, Point b, Point c, double radius)/*Checking whether mid point and chord circles are empty*/
{
    if (distance(a, b) >= radius)
    {
        Point p1 = Point(((a.x() + b.x()) / 2 - (radius / 2)), ((a.y() + b.y()) / 2), 1);
        Point p2 = Point(((a.x() + b.x()) / 2 + (radius / 2)), ((a.y() + b.y()) / 2), 1);/*(p1,p2) is the midpoint of edge (a,b)*/
        if (liesInsideScaled( p1, p2, c))
        {
            return true;
        }
        return false;
    }
    radius = (radius / 2);
    double d = sqrt(std::abs((radius * radius) - ((distance(a, b) / 2) * (distance(a, b) / 2))));
    double a1New = ((a.x() + b.x()) / 2) + (d / (sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x()))))) * (a.y() - b.y());
    double b1New = ((a.y() + b.y()) / 2) + (d / (sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x()))))) * (b.x() - a.x());
    double a2New = ((a.x() + b.x()) / 2) - (d / (sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x()))))) * (a.y() - b.y());
    double b2New = ((a.y() + b.y()) / 2) - (d / (sqrt(((a.y() - b.y()) * (a.y() - b.y())) + ((b.x() - a.x()) * (b.x() - a.x()))))) * (b.x() - a.x());
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


void FindHole::insertToHoleBoundary( Edges& holeEdges, Point a, Point b)/*Inserting a new edge into ec-shape/boundary */
{
    Edge e;
    e.source = a;
    e.target = b;
    holeEdges.push_back(e);
}


void FindHole::updateHole( PriorityQueue& pqInnerBoundary, allFaceItr currentFace, allFaceItr neighbor1, allFaceItr neighbor2) /*updates the hole boundary*/
{
    if (neighbor1->neighbor(0) == currentFace)
    {
        neighbor1->set_neighbor(0, this->m_delaunayTriangulation.infinite_face());
    }
    else
    {
        if (neighbor1->neighbor(1) == currentFace)
        {
            neighbor1->set_neighbor(1, this->m_delaunayTriangulation.infinite_face());
        }
        else
        {
            if (neighbor1->neighbor(2) == currentFace)
            {
                neighbor1->set_neighbor(2, this->m_delaunayTriangulation.infinite_face());
            }
        }
    }

    if (neighbor2->neighbor(0) == currentFace)
    {
        neighbor2->set_neighbor(0, this->m_delaunayTriangulation.infinite_face());
    }
    else
    {
        if (neighbor2->neighbor(1) == currentFace)
        {
            neighbor2->set_neighbor(1, this->m_delaunayTriangulation.infinite_face());
        }
        else
        {
            if (neighbor2->neighbor(2) == currentFace)
            {
                neighbor2->set_neighbor(2, this->m_delaunayTriangulation.infinite_face());
            }
        }
    }

    this->m_delaunayTriangulation.delete_face(currentFace);

    pqInnerBoundary.pqInsert(area(neighbor1->vertex(0)->point(), neighbor1->vertex(1)->point(), neighbor1->vertex(2)->point()), neighbor1);
    pqInnerBoundary.pqInsert(area(neighbor2->vertex(0)->point(), neighbor2->vertex(1)->point(), neighbor2->vertex(2)->point()), neighbor2);
}


bool FindHole::inBoundary( Points exists, Point a)/*Checking whether a point is already in shape*/
{
    for (int i = 0; i < exists.size(); i++)
    {
        if (exists[i].x() == a.x() && exists[i].y() == a.y())
        {
            return true;
        }
    }
    return false;
}


bool FindHole::outerBoundaryCondition( Point a, Point b, Point c, faceHandle fh)
{
    if (liesInsideScaled( a, b, c)) /*Diametric circle is non-empty*/
    {
        return true;
    }

    /*Checking the status of mid point circle and chord circles*/
    if (isFinite( fh->neighbor(0)))
    {
        if (fh->neighbor(0)->vertex(0)->point() != a && fh->neighbor(0)->vertex(0)->point() != b && fh->neighbor(0)->vertex(0)->point() != c)
        {
            if (liesWithinException( fh->neighbor(0)->vertex(1)->point(), fh->neighbor(0)->vertex(2)->point(), fh->neighbor(0)->vertex(0)->point(), v * distance(a, b)))
                return true;
        }
        else
        {
            if (fh->neighbor(0)->vertex(1)->point() != a && fh->neighbor(0)->vertex(1)->point() != b && fh->neighbor(0)->vertex(1)->point() != c)
            {
                if (liesWithinException( fh->neighbor(0)->vertex(0)->point(), fh->neighbor(0)->vertex(2)->point(), fh->neighbor(0)->vertex(1)->point(), v * distance(a, b)))
                    return true;
            }
            else
            {
                if (fh->neighbor(0)->vertex(2)->point() != a && fh->neighbor(0)->vertex(2)->point() != b && fh->neighbor(0)->vertex(2)->point() != c)
                {
                    if (liesWithinException( fh->neighbor(0)->vertex(1)->point(), fh->neighbor(0)->vertex(0)->point(), fh->neighbor(0)->vertex(2)->point(), v * distance(a, b)))
                        return true;
                }
            }
        }
    }

    if (isFinite( fh->neighbor(1)))
    {
        if (fh->neighbor(1)->vertex(0)->point() != a && fh->neighbor(1)->vertex(0)->point() != b && fh->neighbor(1)->vertex(0)->point() != c)
        {
            if (liesWithinException(fh->neighbor(1)->vertex(1)->point(), fh->neighbor(1)->vertex(2)->point(), fh->neighbor(1)->vertex(0)->point(), v * distance(a, b)))
                return true;
        }
        else
        {
            if (fh->neighbor(1)->vertex(1)->point() != a && fh->neighbor(1)->vertex(1)->point() != b && fh->neighbor(1)->vertex(1)->point() != c)
            {
                if (liesWithinException(fh->neighbor(1)->vertex(2)->point(), fh->neighbor(1)->vertex(0)->point(), fh->neighbor(1)->vertex(1)->point(), v * distance(a, b)))
                    return true;
            }
            else
            {
                if (fh->neighbor(1)->vertex(2)->point() != a && fh->neighbor(1)->vertex(2)->point() != b && fh->neighbor(1)->vertex(2)->point() != c)
                {
                    if (liesWithinException( fh->neighbor(1)->vertex(0)->point(), fh->neighbor(1)->vertex(1)->point(), fh->neighbor(1)->vertex(2)->point(), v * distance(a, b)))
                        return true;
                }
            }
        }
    }

    if (isFinite( fh->neighbor(2)))
    {
        if (fh->neighbor(2)->vertex(0)->point() != a && fh->neighbor(2)->vertex(0)->point() != b && fh->neighbor(2)->vertex(0)->point() != c)
        {
            if (liesWithinException( fh->neighbor(2)->vertex(1)->point(), fh->neighbor(2)->vertex(2)->point(), fh->neighbor(2)->vertex(0)->point(), v * distance(a, b)))
                return true;
        }
        else
        {
            if (fh->neighbor(2)->vertex(1)->point() != a && fh->neighbor(2)->vertex(1)->point() != b && fh->neighbor(2)->vertex(1)->point() != c)
            {
                if (liesWithinException( fh->neighbor(2)->vertex(2)->point(), fh->neighbor(2)->vertex(0)->point(), fh->neighbor(2)->vertex(1)->point(), v * distance(a, b)))
                    return true;
            }
            else
            {
                if (fh->neighbor(2)->vertex(2)->point() != a && fh->neighbor(2)->vertex(2)->point() != b && fh->neighbor(2)->vertex(2)->point() != c)
                {
                    if (liesWithinException( fh->neighbor(2)->vertex(0)->point(), fh->neighbor(2)->vertex(1)->point(), fh->neighbor(2)->vertex(2)->point(), v * distance(a, b)))
                        return true;
                }
            }
        }
    }

    return false;
}


bool FindHole::holeDetectionCondition( Point a, Point b, Point c, faceHandle fh)
{
    if (liesInsideHole( a, b, c))/*Diametric circle is non-empty*/
    {
        return true;
    }

    /*Checking the status of mid point circle and chord circles*/
    if (isFinite( fh->neighbor(0)))
    {
        if (fh->neighbor(0)->vertex(0)->point() != a && fh->neighbor(0)->vertex(0)->point() != b && fh->neighbor(0)->vertex(0)->point() != c)
        {
            if (liesWithinExceptionHole( fh->neighbor(0)->vertex(1)->point(), fh->neighbor(0)->vertex(2)->point(), fh->neighbor(0)->vertex(0)->point(), distance(a, b)))
            {
                return true;
            }
        }
        else
        {
            if (fh->neighbor(0)->vertex(1)->point() != a && fh->neighbor(0)->vertex(1)->point() != b && fh->neighbor(0)->vertex(1)->point() != c)
            {
                if (liesWithinExceptionHole( fh->neighbor(0)->vertex(0)->point(), fh->neighbor(0)->vertex(2)->point(), fh->neighbor(0)->vertex(1)->point(), distance(a, b)))
                {
                    return true;
                }
            }
            else
            {
                if (fh->neighbor(0)->vertex(2)->point() != a && fh->neighbor(0)->vertex(2)->point() != b && fh->neighbor(0)->vertex(2)->point() != c)
                {
                    if (liesWithinExceptionHole( fh->neighbor(0)->vertex(1)->point(), fh->neighbor(0)->vertex(0)->point(), fh->neighbor(0)->vertex(2)->point(), distance(a, b)))
                    {
                        return true;
                    }
                }
            }
        }
    }

    if (isFinite( fh->neighbor(1)))
    {
        if (fh->neighbor(1)->vertex(0)->point() != a && fh->neighbor(1)->vertex(0)->point() != b && fh->neighbor(1)->vertex(0)->point() != c)
        {
            if (liesWithinExceptionHole( fh->neighbor(1)->vertex(1)->point(), fh->neighbor(1)->vertex(2)->point(), fh->neighbor(1)->vertex(0)->point(), distance(a, b)))
            {
                return true;
            }
        }
        else
        {
            if (fh->neighbor(1)->vertex(1)->point() != a && fh->neighbor(1)->vertex(1)->point() != b && fh->neighbor(1)->vertex(1)->point() != c)
            {
                if (liesWithinExceptionHole( fh->neighbor(1)->vertex(2)->point(), fh->neighbor(1)->vertex(0)->point(), fh->neighbor(1)->vertex(1)->point(), distance(a, b)))
                {
                    return true;
                }
            }
            else
            {
                if (fh->neighbor(1)->vertex(2)->point() != a && fh->neighbor(1)->vertex(2)->point() != b && fh->neighbor(1)->vertex(2)->point() != c)
                {
                    if (liesWithinExceptionHole( fh->neighbor(1)->vertex(0)->point(), fh->neighbor(1)->vertex(1)->point(), fh->neighbor(1)->vertex(2)->point(), distance(a, b)))
                    {
                        return true;
                    }
                }
            }
        }
    }

    if (isFinite( fh->neighbor(2)))
    {
        if (fh->neighbor(2)->vertex(0)->point() != a && fh->neighbor(2)->vertex(0)->point() != b && fh->neighbor(2)->vertex(0)->point() != c)
        {
            if (liesWithinExceptionHole( fh->neighbor(2)->vertex(1)->point(), fh->neighbor(2)->vertex(2)->point(), fh->neighbor(2)->vertex(0)->point(), distance(a, b)))
            {
                return true;
            }
        }
        else
        {
            if (fh->neighbor(2)->vertex(1)->point() != a && fh->neighbor(2)->vertex(1)->point() != b && fh->neighbor(2)->vertex(1)->point() != c)
            {
                if (liesWithinExceptionHole( fh->neighbor(2)->vertex(2)->point(), fh->neighbor(2)->vertex(0)->point(), fh->neighbor(2)->vertex(1)->point(), distance(a, b)))
                {
                    return true;
                }
            }
            else
            {
                if (fh->neighbor(2)->vertex(2)->point() != a && fh->neighbor(2)->vertex(2)->point() != b && fh->neighbor(2)->vertex(2)->point() != c)
                {
                    if (liesWithinExceptionHole( fh->neighbor(2)->vertex(0)->point(), fh->neighbor(2)->vertex(1)->point(), fh->neighbor(2)->vertex(2)->point(), distance(a, b)))
                    {
                        return true;
                    }
                }
            }
        }
    }

    return false;
}


void FindHole::outerBoundary( Points& exists, Edges& shape, PriorityQueue& pqOuterBoundary)
{
    for (Delaunay::Finite_faces_iterator ffi = this->m_delaunayTriangulation.finite_faces_begin(); ffi != this->m_delaunayTriangulation.finite_faces_end(); ffi++)
    {
        if (isInfinite( ffi->neighbor(0)) || isInfinite( ffi->neighbor(1)) || isInfinite( ffi->neighbor(2)))
        {
            if (isInfinite( ffi->neighbor(0)))
            {

                ffi->neighbor(0) = this->m_delaunayTriangulation.infinite_face(); //set the 0th neighbor as infinite face
                exists.push_back(ffi->vertex(1)->point());
                exists.push_back(ffi->vertex(2)->point());
                pqOuterBoundary.pqInsert(distance(ffi->vertex(1)->point(), ffi->vertex(2)->point()), ffi);
            }
            else
            {
                if (isInfinite( ffi->neighbor(1)))
                {
                    ffi->neighbor(1) = this->m_delaunayTriangulation.infinite_face();
                    exists.push_back(ffi->vertex(0)->point());
                    exists.push_back(ffi->vertex(2)->point());
                    pqOuterBoundary.pqInsert(distance(ffi->vertex(0)->point(), ffi->vertex(2)->point()), ffi);
                }
                else
                {
                    if (isInfinite( ffi->neighbor(2)))
                    {
                        ffi->neighbor(2) = this->m_delaunayTriangulation.infinite_face();
                        exists.push_back(ffi->vertex(1)->point());
                        exists.push_back(ffi->vertex(0)->point());
                        pqOuterBoundary.pqInsert(distance(ffi->vertex(0)->point(), ffi->vertex(1)->point()), ffi);
                    }
                }
            }
        }
    }

    while (!pqOuterBoundary.isEmpty())
    {
        allFaceItr currentFace = pqOuterBoundary.pqDelete();

        if (isInfinite( currentFace->neighbor(0)) && isFinite( currentFace->neighbor(1)) && isFinite( currentFace->neighbor(2)))
        {
            if (outerBoundaryCondition( currentFace->vertex(1)->point(), currentFace->vertex(2)->point(), currentFace->vertex(0)->point(), currentFace))/*If any of the circle is non-empty*/
            {
                if (!isInShape( exists, currentFace->vertex(0)->point()))
                {
                    exists.push_back(currentFace->vertex(0)->point());
                    updateNeighbors(pqOuterBoundary, currentFace, currentFace->neighbor(1), currentFace->neighbor(2), distance(currentFace->vertex(0)->point(), currentFace->vertex(2)->point()), distance(currentFace->vertex(1)->point(), currentFace->vertex(0)->point()));
                }
                else
                {
                    insertToShape( shape, currentFace->vertex(1)->point(), currentFace->vertex(2)->point());
                }
            }
            else
            {
                insertToShape( shape, currentFace->vertex(1)->point(), currentFace->vertex(2)->point());
            }
        }

        else
        {
            if (isInfinite( currentFace->neighbor(1)) && isFinite( currentFace->neighbor(0)) && isFinite( currentFace->neighbor(2)))
            {
                if (outerBoundaryCondition( currentFace->vertex(0)->point(), currentFace->vertex(2)->point(), currentFace->vertex(1)->point(), currentFace))
                {
                    if (!isInShape( exists, currentFace->vertex(1)->point()))
                    {
                        exists.push_back(currentFace->vertex(1)->point());
                        updateNeighbors( pqOuterBoundary, currentFace, currentFace->neighbor(2), currentFace->neighbor(0), distance(currentFace->vertex(0)->point(), currentFace->vertex(1)->point()), distance(currentFace->vertex(1)->point(), currentFace->vertex(2)->point()));
                    }
                    else
                    {
                        insertToShape( shape, currentFace->vertex(0)->point(), currentFace->vertex(2)->point());
                    }
                }
                else
                {
                    insertToShape( shape, currentFace->vertex(0)->point(), currentFace->vertex(2)->point());
                }
            }
            else
            {
                if (isInfinite( currentFace->neighbor(2)) && isFinite( currentFace->neighbor(1)) && isFinite( currentFace->neighbor(0)))
                {
                    if (outerBoundaryCondition( currentFace->vertex(1)->point(), currentFace->vertex(0)->point(), currentFace->vertex(2)->point(), currentFace))
                    {
                        if (!isInShape( exists, currentFace->vertex(2)->point()))
                        {
                            exists.push_back(currentFace->vertex(2)->point());
                            updateNeighbors( pqOuterBoundary, currentFace, currentFace->neighbor(1), currentFace->neighbor(0), distance(currentFace->vertex(0)->point(), currentFace->vertex(2)->point()), distance(currentFace->vertex(1)->point(), currentFace->vertex(2)->point()));
                        }
                        else
                        {
                            insertToShape( shape, currentFace->vertex(0)->point(), currentFace->vertex(1)->point());
                        }
                    }
                    else
                    {
                        insertToShape( shape, currentFace->vertex(0)->point(), currentFace->vertex(1)->point());
                    }
                }
                else
                {
                    if (isInfinite( currentFace->neighbor(2)) && isInfinite( currentFace->neighbor(1)) && isFinite( currentFace->neighbor(0)))
                    {
                        insertToShape( shape, currentFace->vertex(2)->point(), currentFace->vertex(0)->point());
                        insertToShape( shape, currentFace->vertex(0)->point(), currentFace->vertex(1)->point());
                    }
                    else
                    {
                        if (isInfinite( currentFace->neighbor(2)) && isFinite( currentFace->neighbor(1)) && isInfinite( currentFace->neighbor(0)))
                        {
                            insertToShape( shape, currentFace->vertex(2)->point(), currentFace->vertex(1)->point());
                            insertToShape( shape, currentFace->vertex(1)->point(), currentFace->vertex(0)->point());
                        }
                        else
                        {
                            if (isFinite( currentFace->neighbor(2)) && isInfinite( currentFace->neighbor(1)) && isInfinite( currentFace->neighbor(0)))
                            {
                                insertToShape( shape, currentFace->vertex(0)->point(), currentFace->vertex(2)->point());
                                insertToShape( shape, currentFace->vertex(2)->point(), currentFace->vertex(1)->point());
                            }
                        }
                    }
                }
            }
        }
    }
    std::string str = "";
    std::cout << "Boundary Edge:" << std::endl;
    for (const auto& ed : shape) {
        //std::cout << " {" << ed.source << ", " << ed.target << "}, ";

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
        //std::getline(std::cin, str);
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


void FindHole::holeDetection( Points exists, Edges& shape, Edges& holeEdges, PriorityQueue& pqInnerBoundary, std::vector< std::vector<Point>>& holes, int& numberOfHoles)
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
            if (holeEdges.size() == 0)
            {
                exit(0);
            }
            else
            {
                break;
            }
        }

        //checking any of the vertices of HighestAreaTriangle are in the outerbdry, if so do not put that it to Q
        for (int i = 0; i < holeEdges.size(); i++)
        {
            if ((highestAreaTriangleHandle->vertex(0)->point().x() == holeEdges[i].source.x() && highestAreaTriangleHandle->vertex(0)->point().y() == holeEdges[i].source.y()) ||
                (highestAreaTriangleHandle->vertex(0)->point().x() == holeEdges[i].target.x() && highestAreaTriangleHandle->vertex(0)->point().y() == holeEdges[i].target.y()) ||
                (highestAreaTriangleHandle->vertex(1)->point().x() == holeEdges[i].source.x() && highestAreaTriangleHandle->vertex(1)->point().y() == holeEdges[i].source.y()) ||
                (highestAreaTriangleHandle->vertex(1)->point().x() == holeEdges[i].target.x() && highestAreaTriangleHandle->vertex(1)->point().y() == holeEdges[i].target.y()) ||
                (highestAreaTriangleHandle->vertex(2)->point().x() == holeEdges[i].source.x() && highestAreaTriangleHandle->vertex(2)->point().y() == holeEdges[i].source.y()) ||
                (highestAreaTriangleHandle->vertex(2)->point().x() == holeEdges[i].target.x() && highestAreaTriangleHandle->vertex(2)->point().y() == holeEdges[i].target.y()))
            {
                inHole = true;
                break;
            }
        }

        //checking any of the vertices of Highest Area Triangle are in the outerbdry, if so do not put that it to Q
        for (int i = 0; i < shape.size(); i++)
        {
            if ((highestAreaTriangleHandle->vertex(0)->point().x() == shape[i].source.x() && highestAreaTriangleHandle->vertex(0)->point().y() == shape[i].source.y()) ||
                (highestAreaTriangleHandle->vertex(0)->point().x() == shape[i].target.x() && highestAreaTriangleHandle->vertex(0)->point().y() == shape[i].target.y()) ||
                (highestAreaTriangleHandle->vertex(1)->point().x() == shape[i].source.x() && highestAreaTriangleHandle->vertex(1)->point().y() == shape[i].source.y()) ||
                (highestAreaTriangleHandle->vertex(1)->point().x() == shape[i].target.x() && highestAreaTriangleHandle->vertex(1)->point().y() == shape[i].target.y()) ||
                (highestAreaTriangleHandle->vertex(2)->point().x() == shape[i].source.x() && highestAreaTriangleHandle->vertex(2)->point().y() == shape[i].source.y()) ||
                (highestAreaTriangleHandle->vertex(2)->point().x() == shape[i].target.x() && highestAreaTriangleHandle->vertex(2)->point().y() == shape[i].target.y()))
            {
                inOuterBdry = true;
                break;
            }
        }

        if (inOuterBdry)
        {
            std::cout << "highest area triangle is outer bdry triangle and it is an invalid input" << std::endl;
            if (holeEdges.size() == 0)
            {
                exit(0);
            }
            else
            {
                break;
            }
        }
        else if ((!inOuterBdry) && (!inHole))//if the hole is not along the outer bdry, initialize the queue1 with the highest area traingle
        {
            /* Initialize a Priority Queue(PQ) with the adjacent triangles of
             △i jk in the descending order of area.
             inserting all neigbhours*/

            faceHandle0 = highestAreaTriangleHandle->neighbor(0);
            pqInnerBoundary.pqInsert(area(faceHandle0->vertex(0)->point(), faceHandle0->vertex(1)->point(), faceHandle0->vertex(2)->point()), faceHandle0);

            faceHandle1 = highestAreaTriangleHandle->neighbor(1);
            pqInnerBoundary.pqInsert(area(faceHandle1->vertex(0)->point(), faceHandle1->vertex(1)->point(), faceHandle1->vertex(2)->point()), faceHandle1);

            faceHandle2 = highestAreaTriangleHandle->neighbor(2);
            pqInnerBoundary.pqInsert(area(faceHandle2->vertex(0)->point(), faceHandle2->vertex(1)->point(), faceHandle2->vertex(2)->point()), faceHandle2);

            zerothPointFaceHandle2 = firstPointFaceHandle2 = secondPointFaceHandle2 = false;

            /*the next three if loops are to find out the vertex number of the third vertex w.r.t facehandle2 when v0,v1 is processed. this is to use it in lies inside condition checking for diacircle
            to get the unvisited vertex of the neighboring triangle of 2nd vertex for cheking lies inside condition*/

            if (((faceHandle2->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x()) && (faceHandle2->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y())) ||
                ((faceHandle2->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x()) && (faceHandle2->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y())))
            {
                zerothPointFaceHandle2 = false;
            }
            else
            {
                zerothPointFaceHandle2 = 1;
            }

            if ((faceHandle2->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle2->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y()) ||
                (faceHandle2->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle2->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y()))
            {
                firstPointFaceHandle2 = false;
            }
            else
            {
                firstPointFaceHandle2 = true;
            }

            if ((faceHandle2->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle2->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y()) ||
                (faceHandle2->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle2->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y()))
            {
                secondPointFaceHandle2 = false;
            }
            else
            {
                secondPointFaceHandle2 = true;
            }

            if (zerothPointFaceHandle2)
            {
                faceHandle2->set_neighbor(0, this->m_delaunayTriangulation.infinite_face());
            }
            else if (firstPointFaceHandle2)
            {
                faceHandle2->set_neighbor(1, this->m_delaunayTriangulation.infinite_face());
            }


            else if (secondPointFaceHandle2)
            {
                faceHandle2->set_neighbor(2, this->m_delaunayTriangulation.infinite_face());
            }

            zerothPointFaceHandle0 = firstPointFaceHandle0 = secondPointFaceHandle0 = false;
            if (((faceHandle0->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x()) && (faceHandle0->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y())) ||
                ((faceHandle0->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x()) && (faceHandle0->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y())))
            {
                zerothPointFaceHandle0 = false;
            }
            else
            {
                zerothPointFaceHandle0 = true;
            }

            if ((faceHandle0->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle0->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y()) ||
                (faceHandle0->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle0->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y()))
            {
                firstPointFaceHandle0 = false;
            }
            else
            {
                firstPointFaceHandle0 = true;
            }

            if ((faceHandle0->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle0->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y()) ||
                (faceHandle0->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle0->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y()))
            {
                secondPointFaceHandle0 = false;
            }
            else
            {
                secondPointFaceHandle0 = true;
            }

            if (zerothPointFaceHandle0)
            {
                faceHandle0->set_neighbor(0, this->m_delaunayTriangulation.infinite_face());
            }
            else if (firstPointFaceHandle0)
            {
                faceHandle0->set_neighbor(1, this->m_delaunayTriangulation.infinite_face());
            }
            else if (secondPointFaceHandle0)
            {
                faceHandle0->set_neighbor(2, this->m_delaunayTriangulation.infinite_face());
            }

            zerothPointFaceHandle1 = firstPointFaceHandle1 = secondPointFaceHandle1 = false;

            if (((faceHandle1->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x()) && (faceHandle1->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y())) ||
                ((faceHandle1->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x()) && (faceHandle1->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y())))
            {
                zerothPointFaceHandle1 = false;
            }
            else
            {
                zerothPointFaceHandle1 = true;
            }

            if ((faceHandle1->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle1->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y()) ||
                (faceHandle1->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle1->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y()))
            {
                firstPointFaceHandle1 = false;
            }
            else
            {
                firstPointFaceHandle1 = true;
            }

            if ((faceHandle1->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle1->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y()) ||
                (faceHandle1->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle1->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y()))
            {
                secondPointFaceHandle1 = false;
            }
            else
            {
                secondPointFaceHandle1 = true;
            }


            if (zerothPointFaceHandle1)
            {
                faceHandle1->set_neighbor(0, this->m_delaunayTriangulation.infinite_face());
            }
            else if (firstPointFaceHandle1)
            {
                faceHandle1->set_neighbor(1, this->m_delaunayTriangulation.infinite_face());
            }
            else if (secondPointFaceHandle1)
            {
                faceHandle1->set_neighbor(2, this->m_delaunayTriangulation.infinite_face());
            }

            exists.push_back(highestAreaTriangleHandle->vertex(0)->point());
            exists.push_back(highestAreaTriangleHandle->vertex(1)->point());
            exists.push_back(highestAreaTriangleHandle->vertex(2)->point());
            this->m_delaunayTriangulation.delete_face(highestAreaTriangleHandle);
        }

        while (!pqInnerBoundary.isEmpty())
        {
            highestAreaTriangleHandle = pqInnerBoundary.pqDelete();

            if (isInfinite(highestAreaTriangleHandle->neighbor(0)) && isFinite( highestAreaTriangleHandle->neighbor(1)) && isFinite( highestAreaTriangleHandle->neighbor(2)))
            {
                if (holeDetectionCondition( highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle))
                {
                    if (!inBoundary( exists, highestAreaTriangleHandle->vertex(0)->point()))
                    {
                        exists.push_back(highestAreaTriangleHandle->vertex(0)->point());
                        updateHole( pqInnerBoundary, highestAreaTriangleHandle, highestAreaTriangleHandle->neighbor(1), highestAreaTriangleHandle->neighbor(2));
                    }
                    else
                    {
                        insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(2)->point());
                        hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                        //  hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                    }
                }
                else
                {
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(2)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                    // hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                }
            }
            else if (isInfinite( highestAreaTriangleHandle->neighbor(1)) && isFinite( highestAreaTriangleHandle->neighbor(0)) && isFinite( highestAreaTriangleHandle->neighbor(2)))
            {
                if (holeDetectionCondition( highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle))
                {
                    if (!inBoundary( exists, highestAreaTriangleHandle->vertex(1)->point()))
                    {
                        exists.push_back(highestAreaTriangleHandle->vertex(1)->point());
                        updateHole( pqInnerBoundary, highestAreaTriangleHandle, highestAreaTriangleHandle->neighbor(2), highestAreaTriangleHandle->neighbor(0));
                    }
                    else
                    {
                        insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point());
                        hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                        // hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                    }
                }
                else
                {
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    //hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                }
            }
            else if (isInfinite( highestAreaTriangleHandle->neighbor(2)) && isFinite( highestAreaTriangleHandle->neighbor(1)) && isFinite( highestAreaTriangleHandle->neighbor(0)))
            {
                if (holeDetectionCondition( highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle))
                {
                    if (!inBoundary( exists, highestAreaTriangleHandle->vertex(2)->point()))
                    {
                        exists.push_back(highestAreaTriangleHandle->vertex(2)->point());
                        updateHole( pqInnerBoundary, highestAreaTriangleHandle, highestAreaTriangleHandle->neighbor(1), highestAreaTriangleHandle->neighbor(0));
                    }
                    else
                    {
                        insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(1)->point());
                        hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                        //hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                    }
                }
                else
                {
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(1)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    // hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                }
            }
            else
            {
                if (isInfinite( highestAreaTriangleHandle->neighbor(2)) && isInfinite( highestAreaTriangleHandle->neighbor(1)) && isFinite( highestAreaTriangleHandle->neighbor(0)))
                {
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(0)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                    //hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(1)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                }
                else if (isInfinite( highestAreaTriangleHandle->neighbor(2)) && isFinite( highestAreaTriangleHandle->neighbor(1)) && isInfinite( highestAreaTriangleHandle->neighbor(0)))
                {
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(1)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(0)->point());
                    //hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                }
                else if (isFinite( highestAreaTriangleHandle->neighbor(2)) && isInfinite( highestAreaTriangleHandle->neighbor(1)) && isInfinite( highestAreaTriangleHandle->neighbor(0)))
                {
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    // hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(1)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                }
            }
        }
        if (hole.size() > 4) {
            holes.push_back(hole);
            ++numberOfHoles;
        }


    } while ((holeEdges.size() - 1 - tempVar) > 4);


    std::cout << "holedetection :No of Holes:" << numberOfHoles << "," << holeEdges.size() << std::endl;
}


void FindHole::singleHoleDetection( Edges& shape, Edges& holeEdges, std::vector< std::vector<Point>>& holes, int& numberOfHoles)
{
    Points  exists;
    // std::vector<std::vector<Point>> holes
    PriorityQueue pqOuterBoundary;
    PriorityQueue pqInnerBoundary;
    //Edges shape;
    //Edges holeEdges;
    //int numberOfHoles = 0;




    Delaunay::Finite_faces_iterator ffi = this->m_delaunayTriangulation.finite_faces_begin();

    Delaunay::Face_handle highestAreaTriangleHandle;//the face handle for the highest area triangle
    Delaunay::Face_handle faceHandle0, faceHandle1, faceHandle2;//facehandles for neighbouring triangles
    bool inHole;
    int tempVar = 0;//variable for for multiple hole detection
    bool zerothPointFaceHandle0, firstPointFaceHandle0, secondPointFaceHandle0;
    bool zerothPointFaceHandle1, firstPointFaceHandle1, secondPointFaceHandle1;
    bool zerothPointFaceHandle2, firstPointFaceHandle2, secondPointFaceHandle2;
    bool inOuterBdry = false;

    //do while loop for multiple hole detection

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
        if (holeEdges.size() == 0)
        {
            std::cout << "line 985 exit" << std::endl;
            exit(0);
        }
        else
        {
            //break;
        }
    }

    //checking any of the vertices of HAT are in the outerbdry, if so do not put that it to Q
    for (int i = 0; i < holeEdges.size(); i++)
    {
        if ((highestAreaTriangleHandle->vertex(0)->point().x() == holeEdges[i].source.x() && highestAreaTriangleHandle->vertex(0)->point().y() == holeEdges[i].source.y()) ||
            (highestAreaTriangleHandle->vertex(0)->point().x() == holeEdges[i].target.x() && highestAreaTriangleHandle->vertex(0)->point().y() == holeEdges[i].target.y()) ||
            (highestAreaTriangleHandle->vertex(1)->point().x() == holeEdges[i].source.x() && highestAreaTriangleHandle->vertex(1)->point().y() == holeEdges[i].source.y()) ||
            (highestAreaTriangleHandle->vertex(1)->point().x() == holeEdges[i].target.x() && highestAreaTriangleHandle->vertex(1)->point().y() == holeEdges[i].target.y()) ||
            (highestAreaTriangleHandle->vertex(2)->point().x() == holeEdges[i].source.x() && highestAreaTriangleHandle->vertex(2)->point().y() == holeEdges[i].source.y()) ||
            (highestAreaTriangleHandle->vertex(2)->point().x() == holeEdges[i].target.x() && highestAreaTriangleHandle->vertex(2)->point().y() == holeEdges[i].target.y()))
        {
            inHole = true;
            std::cout << "line 1005 break" << std::endl;
            break;
        }
    }

    //checking any of the vertices of Highest Area Triangle are in the outerbdry, if so do not put that it to Q
    for (int i = 0; i < shape.size(); i++)
    {
        if ((highestAreaTriangleHandle->vertex(0)->point().x() == shape[i].source.x() && highestAreaTriangleHandle->vertex(0)->point().y() == shape[i].source.y()) ||
            (highestAreaTriangleHandle->vertex(0)->point().x() == shape[i].target.x() && highestAreaTriangleHandle->vertex(0)->point().y() == shape[i].target.y()) ||
            (highestAreaTriangleHandle->vertex(1)->point().x() == shape[i].source.x() && highestAreaTriangleHandle->vertex(1)->point().y() == shape[i].source.y()) ||
            (highestAreaTriangleHandle->vertex(1)->point().x() == shape[i].target.x() && highestAreaTriangleHandle->vertex(1)->point().y() == shape[i].target.y()) ||
            (highestAreaTriangleHandle->vertex(2)->point().x() == shape[i].source.x() && highestAreaTriangleHandle->vertex(2)->point().y() == shape[i].source.y()) ||
            (highestAreaTriangleHandle->vertex(2)->point().x() == shape[i].target.x() && highestAreaTriangleHandle->vertex(2)->point().y() == shape[i].target.y()))
        {
            inOuterBdry = true;
            std::cout << "line 2021 break" << std::endl;
            break;
        }
    }

    if (inOuterBdry)
    {
        std::cout << "highest area triangle is outer bdry triangle and it is an invalid input" << std::endl;
        if (holeEdges.size() == 0)
        {
            std::cout << "line 1031 exit" << std::endl;
            exit(0);
        }
        else
        {
            //break;
        }
    }
    //inOuterBdry = false; inHole = false;
    if ((!inOuterBdry) && (!inHole))//if the hole is not along the outer bdry, initialize the queue1 with the highest area traingle
    {
        faceHandle0 = highestAreaTriangleHandle->neighbor(0);
        pqInnerBoundary.pqInsert(area(faceHandle0->vertex(0)->point(), faceHandle0->vertex(1)->point(), faceHandle0->vertex(2)->point()), faceHandle0);

        faceHandle1 = highestAreaTriangleHandle->neighbor(1);
        pqInnerBoundary.pqInsert(area(faceHandle1->vertex(0)->point(), faceHandle1->vertex(1)->point(), faceHandle1->vertex(2)->point()), faceHandle1);

        faceHandle2 = highestAreaTriangleHandle->neighbor(2);
        pqInnerBoundary.pqInsert(area(faceHandle2->vertex(0)->point(), faceHandle2->vertex(1)->point(), faceHandle2->vertex(2)->point()), faceHandle2);

        zerothPointFaceHandle2 = firstPointFaceHandle2 = secondPointFaceHandle2 = false;
        //the next three if loops are to find out the vertex number of the third vertex w.r.t facehandle2 when v0,v1 is processed. this is to use it in lies inside condition checking for diacircle
        //to get the unvisited vertex of the neighboring triangle of 2nd vertex for cheking lies inside condition

        if (((faceHandle2->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x()) && (faceHandle2->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y())) ||
            ((faceHandle2->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x()) && (faceHandle2->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y())))
        {
            zerothPointFaceHandle2 = false;
        }
        else
        {
            zerothPointFaceHandle2 = 1;
        }

        if ((faceHandle2->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle2->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y()) ||
            (faceHandle2->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle2->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y()))
        {
            firstPointFaceHandle2 = false;
        }
        else
        {
            firstPointFaceHandle2 = true;
        }

        if ((faceHandle2->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle2->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y()) ||
            (faceHandle2->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle2->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y()))
        {
            secondPointFaceHandle2 = false;
        }
        else
        {
            secondPointFaceHandle2 = true;
        }

        if (zerothPointFaceHandle2)
        {
            faceHandle2->set_neighbor(0, this->m_delaunayTriangulation.infinite_face());
        }
        else if (firstPointFaceHandle2)
        {
            faceHandle2->set_neighbor(1, this->m_delaunayTriangulation.infinite_face());
        }


        else if (secondPointFaceHandle2)
        {
            faceHandle2->set_neighbor(2, this->m_delaunayTriangulation.infinite_face());
        }

        zerothPointFaceHandle0 = firstPointFaceHandle0 = secondPointFaceHandle0 = false;
        if (((faceHandle0->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x()) && (faceHandle0->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y())) ||
            ((faceHandle0->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x()) && (faceHandle0->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y())))
        {
            zerothPointFaceHandle0 = false;
        }
        else
        {
            zerothPointFaceHandle0 = true;
        }

        if ((faceHandle0->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle0->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y()) ||
            (faceHandle0->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle0->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y()))
        {
            firstPointFaceHandle0 = false;
        }
        else
        {
            firstPointFaceHandle0 = true;
        }

        if ((faceHandle0->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle0->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(1)->point().y()) ||
            (faceHandle0->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle0->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y()))
        {
            secondPointFaceHandle0 = false;
        }
        else
        {
            secondPointFaceHandle0 = true;
        }

        if (zerothPointFaceHandle0)
        {
            faceHandle0->set_neighbor(0, this->m_delaunayTriangulation.infinite_face());
        }
        else if (firstPointFaceHandle0)
        {
            faceHandle0->set_neighbor(1, this->m_delaunayTriangulation.infinite_face());
        }
        else if (secondPointFaceHandle0)
        {
            faceHandle0->set_neighbor(2, this->m_delaunayTriangulation.infinite_face());
        }

        zerothPointFaceHandle1 = firstPointFaceHandle1 = secondPointFaceHandle1 = false;

        if (((faceHandle1->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x()) && (faceHandle1->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y())) ||
            ((faceHandle1->vertex(0)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x()) && (faceHandle1->vertex(0)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y())))
        {
            zerothPointFaceHandle1 = false;
        }
        else
        {
            zerothPointFaceHandle1 = true;
        }

        if ((faceHandle1->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle1->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y()) ||
            (faceHandle1->vertex(1)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle1->vertex(1)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y()))
        {
            firstPointFaceHandle1 = false;
        }
        else
        {
            firstPointFaceHandle1 = true;
        }

        if ((faceHandle1->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle1->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(2)->point().y()) ||
            (faceHandle1->vertex(2)->point().x() == highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle1->vertex(2)->point().y() == highestAreaTriangleHandle->vertex(0)->point().y()))
        {
            secondPointFaceHandle1 = false;
        }
        else
        {
            secondPointFaceHandle1 = true;
        }


        if (zerothPointFaceHandle1)
        {
            faceHandle1->set_neighbor(0, this->m_delaunayTriangulation.infinite_face());
        }
        else if (firstPointFaceHandle1)
        {
            faceHandle1->set_neighbor(1, this->m_delaunayTriangulation.infinite_face());
        }
        else if (secondPointFaceHandle1)
        {
            faceHandle1->set_neighbor(2, this->m_delaunayTriangulation.infinite_face());
        }

        exists.push_back(highestAreaTriangleHandle->vertex(0)->point());
        exists.push_back(highestAreaTriangleHandle->vertex(1)->point());
        exists.push_back(highestAreaTriangleHandle->vertex(2)->point());
        this->m_delaunayTriangulation.delete_face(highestAreaTriangleHandle);
    }

    while (!pqInnerBoundary.isEmpty())
    {
        highestAreaTriangleHandle = pqInnerBoundary.pqDelete();

        if (isInfinite( highestAreaTriangleHandle->neighbor(0)) && isFinite( highestAreaTriangleHandle->neighbor(1)) && isFinite( highestAreaTriangleHandle->neighbor(2)))
        {
            if (holeDetectionCondition( highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle))
            {
                if (!inBoundary( exists, highestAreaTriangleHandle->vertex(0)->point()))
                {
                    exists.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    updateHole( pqInnerBoundary, highestAreaTriangleHandle, highestAreaTriangleHandle->neighbor(1), highestAreaTriangleHandle->neighbor(2));
                }
                else
                {
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(2)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                    //  hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                }
            }
            else
            {
                insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(2)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                // hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
            }
        }
        else if (isInfinite( highestAreaTriangleHandle->neighbor(1)) && isFinite( highestAreaTriangleHandle->neighbor(0)) && isFinite( highestAreaTriangleHandle->neighbor(2)))
        {
            if (holeDetectionCondition( highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle))
            {
                if (!inBoundary( exists, highestAreaTriangleHandle->vertex(1)->point()))
                {
                    exists.push_back(highestAreaTriangleHandle->vertex(1)->point());
                    updateHole( pqInnerBoundary, highestAreaTriangleHandle, highestAreaTriangleHandle->neighbor(2), highestAreaTriangleHandle->neighbor(0));
                }
                else
                {
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    // hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                }
            }
            else
            {
                insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                //hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
            }
        }
        else if (isInfinite( highestAreaTriangleHandle->neighbor(2)) && isFinite( highestAreaTriangleHandle->neighbor(1)) && isFinite( highestAreaTriangleHandle->neighbor(0)))
        {
            if (holeDetectionCondition( highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle))
            {
                if (!inBoundary( exists, highestAreaTriangleHandle->vertex(2)->point()))
                {
                    exists.push_back(highestAreaTriangleHandle->vertex(2)->point());
                    updateHole( pqInnerBoundary, highestAreaTriangleHandle, highestAreaTriangleHandle->neighbor(1), highestAreaTriangleHandle->neighbor(0));
                }
                else
                {
                    insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(1)->point());
                    hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                    //hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                }
            }
            else
            {
                insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(1)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                // hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
            }
        }
        else
        {
            if (isInfinite( highestAreaTriangleHandle->neighbor(2)) && isInfinite( highestAreaTriangleHandle->neighbor(1)) && isFinite( highestAreaTriangleHandle->neighbor(0)))
            {
                insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(0)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                //hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(1)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
            }
            else if (isInfinite( highestAreaTriangleHandle->neighbor(2)) && isFinite( highestAreaTriangleHandle->neighbor(1)) && isInfinite( highestAreaTriangleHandle->neighbor(0)))
            {
                insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(1)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(0)->point());
                //hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
            }
            else if (isFinite( highestAreaTriangleHandle->neighbor(2)) && isInfinite( highestAreaTriangleHandle->neighbor(1)) && isInfinite( highestAreaTriangleHandle->neighbor(0)))
            {
                insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(0)->point());
                // hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                insertToHoleBoundary( holeEdges, highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(1)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(2)->point());
                hole.push_back(highestAreaTriangleHandle->vertex(1)->point());
            }
        }

    }

    holes.push_back(hole);
    numberOfHoles++;

    std::cout << "No of Holes:" << numberOfHoles << std::endl;
}


int FindHole::getholes(std::vector<std::vector<Point>>& holes, Edges& holeEdges) {

    // std::vector<std::vector<Point>> holes
    PriorityQueue pqOuterBoundary;
    PriorityQueue pqInnerBoundary;
    Edges shape;
  //  Edges holeEdges;


    double minX = std::numeric_limits<int>::max(), maxX = std::numeric_limits<int>::min();
    double minY = std::numeric_limits<int>::max(), maxY = std::numeric_limits<int>::min();
    Points inputPoints, exists;

    int numberOfHoles = 0;
    //saveDelaunayTriangulationTxtFile(dt, points);

    for (Delaunay::Vertex_iterator vi = this->m_delaunayTriangulation.vertices_begin(); vi != this->m_delaunayTriangulation.vertices_end(); vi++)
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

    //diagonalDistance = distance(Point(minX, minY, 1), Point(maxX, maxY, 1));

   // outerBoundary(dt, exists, shape, pqOuterBoundary);
  //  holeDetection(dt, exists, shape, holeEdges, pqInnerBoundary, holes, numberOfHoles);
    
    FindHole::singleHoleDetection( shape, holeEdges, holes, numberOfHoles);
    //  singleHoleDetection(dt, shape, holeEdges, holes);
    //  singleHoleDetection(dt, shape, holeEdges, holes);
    //  singleHoleDetection(dt, shape, holeEdges, holes);
    //  singleHoleDetection(dt, shape, holeEdges, holes);
    //  singleHoleDetection(dt, shape, holeEdges, holes);
    //  singleHoleDetection(dt, shape, holeEdges, holes);
    //saveFilteredDelaunayTriangulationTxtFile(dt, points);

    return numberOfHoles;

}





//void chechFreq(std::vector<std::pair<double, int>>& freqOfMinDist, double& m) {
//
//    for (int j = 0; j < freqOfMinDist.size(); ++j) {
//        auto f = freqOfMinDist.at(j).first;
//        auto d = f - m;
//        //std::cout << freqOfMinDist.at(j).first<<": "<<m << std::endl;
//
//        if (d <= 0.00001 && d > -0.00001) {
//            freqOfMinDist.at(j).first = m;
//            freqOfMinDist.at(j).second = ++freqOfMinDist.at(j).second;
//            /// break;
//             //std::pair<double, int> tem = { m,1 };
//             //freqOfMinDist.push_back(tem);
//             // std::cout << freqOfMinDist.at(j).first << ": " << freqOfMinDist.at(j).second << std::endl;
//        }
//        else {
//            std::pair<double, int> tem = { m,1 };
//            freqOfMinDist.push_back(tem);
//            //break;
//           // freqOfMinDist.at(j).second = freqOfMinDist.at(j).second + 1;
//            // std::cout << freqOfMinDist.at(j).first << ": " << freqOfMinDist.at(j).second << std::endl;
//
//        }
//    }
//
//}
//
//
//void minDist(const std::vector<Point>& points) {
//
//
//
//
//
//
//    Delaunay dt;
//    dt.insert(points.begin(), points.end());
//    double minX = std::numeric_limits<int>::max(), maxX = std::numeric_limits<int>::min();
//    double minY = std::numeric_limits<int>::max(), maxY = std::numeric_limits<int>::min();
//    std::vector<std::pair<double, int>> freqOfMinDist;
//    int cc = 0;
//    for (Delaunay::Vertex_iterator vi = dt.vertices_begin(); vi != dt.vertices_end(); vi++)
//    {
//        if (vi->point().x() > maxX)
//        {
//            maxX = vi->point().x();
//        }
//        if (vi->point().y() > maxY)
//        {
//            maxY = vi->point().y();
//        }
//        if (vi->point().x() < minX)
//        {
//            minX = vi->point().x();
//        }
//        if (vi->point().y() < minY)
//        {
//            minY = vi->point().y();
//        }
//
//        // inputPoints.push_back(vi->point());
//         //std::string temp = std::to_string(vi->point().x()) + " " + std::to_string(vi->point().y()) + "\n";
//        // str = str + temp;
//    }
//
//
//
//
//
//
//    auto ddy = maxY - minY;
//    ddy = ddy * ddy;
//    auto ddx = maxX - minX;
//    ddx = ddx * ddx;
//    auto dimtr = std::sqrt(ddx + ddy);
//    std::vector<double> minDis;
//    double pixelSize = dimtr;
//    for (int i = 0; i < points.size(); ++i) {
//        double m = dimtr;
//        for (int j = 0; j < points.size(); ++j) {
//            auto d = distance(points.at(i), points.at(j));
//            if (!d == 0 && d < m) {
//                m = d;
//                if (m < pixelSize) {
//                    pixelSize = m;
//                }
//            }
//        }
//        minDis.push_back(m);
//
//        //if (i == 0) {
//        //    std::pair<double, int> tem = {m,1};
//        //    freqOfMinDist.push_back(tem);
//        //   // std::cout << freqOfMinDist.at(j).first << ": " << freqOfMinDist.at(j).second << std::endl;
//        //    
//        //}
//        //else{
//        //    for (int j = 0; j < freqOfMinDist.size(); ++j) {
//        //        double t = freqOfMinDist.at(j).first - m;
//        //        if (std::abs(t)<=0.1) {
//        //            //freqOfMinDist.at(j).first = m;
//        //            //freqOfMinDist.at(j).second =  1;
//        //            std::pair<double, int> tem = { m,freqOfMinDist.at(j).second + 1 };
//        //            freqOfMinDist.push_back(tem);
//        //           // std::cout << freqOfMinDist.at(j).first << ": " << freqOfMinDist.at(j).second << std::endl;
//        //            break;
//        //        }
//        //        else {
//        //            std::pair<double, int> tem = { m,1 };
//        //            freqOfMinDist.push_back(tem);
//        //           // std::cout << freqOfMinDist.at(j).first << ": " << freqOfMinDist.at(j).second << std::endl;
//        //            break;
//
//        //        }
//        //}
//        //   // chechFreq(freqOfMinDist, m);
//        //
//        //
//        //}
//       // std::cout <<cc<<"min dist of Point:" << points.at(i) <<" " << m << std::endl;
//        ++cc;
//    }
//    int cf = 0;
//
//    std::vector<std::pair<double, int>> freqOfMinDist1;
//    std::string str = "";
//    std::vector<int>fd;
//
//    for (auto pi : minDis) {
//        if (pi < pixelSize) {
//            pixelSize = pi;
//        }
//        str = str + std::to_string(pi) + "\n";
//        auto digit = pi;
//        // Specify the digit as a input to find the
//        // occurrence of digit
//        //std::pair<double, int> tem = { m,freqOfMinDist.at(j).second + 1 };
//        //  auto itr= partition(minDis.begin(), minDis.end(), [&digit](double x) { return x == digit; });
//        // int count = distance(minDis.begin(), itr);
//        // std::cout <<cf<<" " << distance(minDis.begin(), itr) << std::endl;
//        // if(count>0)++cf;
//        // 
//        // 
//        // 
//        std::cout << minDis.size() << ":" << pi << std::endl;
//        if (cf <= 2) {
//            std::pair<double, int> tem = { pi,1 };
//            freqOfMinDist.push_back(tem);
//        }
//        else {
//            chechFreq(freqOfMinDist, pi);
//
//        }
//
//
//
//
//        //for (int j = 0; j < freqOfMinDist1.size(); ++j) {
//        //    double t = freqOfMinDist1.at(j).first - pi;
//        //    if (std::abs(t) < 5) {
//        //        //freqOfMinDist.at(j).first = m;
//        //        //freqOfMinDist.at(j).second =  1;
//        //        std::pair<double, int> tem = { pi,freqOfMinDist.at(j).second + 1 };
//        //        freqOfMinDist1.push_back(tem);
//        //        // std::cout << freqOfMinDist.at(j).first << ": " << freqOfMinDist.at(j).second << std::endl;
//        //        break;
//        //    }
//        //    else {
//        //        std::pair<double, int> tem = { pi,1 };
//        //        freqOfMinDist1.push_back(tem);
//        //        // std::cout << freqOfMinDist.at(j).first << ": " << freqOfMinDist.at(j).second << std::endl;
//        //        break;
//
//        //    }
//        //}
//
//
//
//
//
//
//
//
//        ++cf;
//    }
//
//
//    //for (int i = 0; i < minDis.size(); ++i) {
//    //    int f = 1;
//    //    double t1 = minDis.at(i);
//    //    for (int j = 0; j < minDis.size(); ++j) {
//    //        double t2 = minDis.at(j);
//    //        auto td =std::abs( t1 - t2);
//    //        if ((i != j) && td<0.001) {
//    //            ++f;
//    //        }
//
//
//
//
//    //    }
//
//    //    fd.push_back(f);
//    //}
//
//
//
//
//
//   // str = str + "\n";
//
//
//
//
//
//
//
//    std::sort(freqOfMinDist.begin(), freqOfMinDist.end(), [](const auto& lhs, const auto& rhs) {
//        return lhs.second < rhs.second; });
//
//    std::sort(freqOfMinDist.begin(), freqOfMinDist.end(), [](const auto& lhs, const auto& rhs) {
//        return lhs.first > rhs.first; });
//
//    //for (const auto& pair : freqOfMinDist) {
//    //    std::cout << pair.first << ": " << pair.second << std::endl;
//    //    //str = str + std::to_string(pair.first)+" "+ std::to_string(pair.second) + "\n";
//
//    //}
//    ////str = str + "\n";
//    //std::cout << "freq size:" << freqOfMinDist.size() << std::endl;
//    //std::cout << "pixelSize radius:" << pixelSize << std::endl;
//
//    //auto cp = CGAL::centroid(points.begin(), points.end(), CGAL::Dimension_tag<0>());
//
//    //std::cout << "Centriod:" << cp << std::endl;
//
//    //std::cout << "horizontal length:" << maxX - minX << std::endl;
//    //std::cout << "vertical length:" << maxY - minY << std::endl;
//    //std::cout << "horizontal length pix no:" << (maxX - minX) / pixelSize << std::endl;
//    //std::cout << "vertical length pix no:" << (maxY - minY) / pixelSize << std::endl;
//
//
//
//    try {
//        std::string tfilename = "E:\\Avinash C drive\\Downloads\\gp610-20230826-win64-mingw\\gnuplot\\bin\\freqDist.txt";
//        std::ofstream tfile(tfilename);
//
//        if (!tfile.is_open())
//            throw "Error to Open File !!!" + tfilename;
//        /*cout << "Error to Open File !!!" << endl;
//        exit(1);*/
//        //std::vector<K::Point_2> tconvexhull;
//        // CGAL::convex_hull_2(hole.begin(), hole.end(), std::back_inserter(tconvexhull));
//         //std::getline(std::cin, str);
//        if (str.empty())
//            throw "Write the file atleast one char !!!";
//        tfile.write(str.c_str(), str.length());
//        tfile.close();
//
//        std::cout << "\n hole file saved to " << tfilename << std::endl;
//
//
//    }
//    catch (...) {
//        std::cout << "\n ERROR!!! 1208" << std::endl;
//    }
//
//
//}
//
//
//void mostProbableHoles(std::vector<std::vector<Point>>& holes) {
//    int nh = 1;
//    Edges holeEdges;
//    for (auto hole : holes) {
//        std::cout << " ............................................................................... " << std::endl;
//        std::cout << "Hole " << nh << ":" << std::endl;
//        std::string str = "";
//
//        std::vector<Point> tp;
//        std::vector< std::vector<Point>> tpholes;
//        for (auto ptInHole : hole) {
//            //std::cout << ptInHole << std::endl;
//            //std::string temp = std::to_string(ptInHole.x()) + " " + std::to_string(ptInHole.y()) + "\n";
//            //str = str + temp;
//            tp.push_back(ptInHole);
//        }
//        //str = str + "\n";
//        //std::cout <<str<< std::endl;
//        Delaunay tdt;
//        tdt.insert(tp.begin(), tp.end());
//        double minX = std::numeric_limits<int>::max(), maxX = std::numeric_limits<int>::min();
//        double minY = std::numeric_limits<int>::max(), maxY = std::numeric_limits<int>::min();
//        Points inputPoints, exists;
//        // std::vector<std::vector<Point>> holes
//        PriorityQueue pqOuterBoundary;
//        PriorityQueue pqInnerBoundary;
//        Edges shape;
//        Edges holeEdges;
//        int numberOfHoles = 0;
//
//        for (Delaunay::Vertex_iterator vi = tdt.vertices_begin(); vi != tdt.vertices_end(); vi++)
//        {
//            if (vi->point().x() > maxX)
//            {
//                maxX = vi->point().x();
//            }
//            if (vi->point().y() > maxY)
//            {
//                maxY = vi->point().y();
//            }
//            if (vi->point().x() < minX)
//            {
//                minX = vi->point().x();
//            }
//            if (vi->point().y() < minY)
//            {
//                minY = vi->point().y();
//            }
//
//            inputPoints.push_back(vi->point());
//            //std::string temp = std::to_string(vi->point().x()) + " " + std::to_string(vi->point().y()) + "\n";
//            //str = str + temp;
//        }
//        //str = str + "\n";
//        //std::cout << str << std::endl;
//
//       // outerBoundary(tdt, exists, shape, pqOuterBoundary);
//        // holeDetection(tdt, exists, shape, holeEdges, pqInnerBoundary, tpholes, numberOfHoles);
//
//
//
//        for (Delaunay::Vertex_iterator vi = tdt.vertices_begin(); vi != tdt.vertices_end(); vi++)
//        {
//            if (vi->point().x() > maxX)
//            {
//                maxX = vi->point().x();
//            }
//            if (vi->point().y() > maxY)
//            {
//                maxY = vi->point().y();
//            }
//            if (vi->point().x() < minX)
//            {
//                minX = vi->point().x();
//            }
//            if (vi->point().y() < minY)
//            {
//                minY = vi->point().y();
//            }
//
//            // inputPoints.push_back(vi->point());
//             //std::string temp = std::to_string(vi->point().x()) + " " + std::to_string(vi->point().y()) + "\n";
//            // str = str + temp;
//        }
//        //str = str + "\n";
//        //std::cout << str << std::endl;
//        //holeDetection(tdt, exists, shape, holeEdges, pqInnerBoundary, tpholes, numberOfHoles);
//
//
//
//
//
//
//        str = "";
//        int trnls = 0;
//        for (Delaunay::Finite_faces_iterator fit = tdt.finite_faces_begin(); fit != tdt.finite_faces_end(); ++fit) {
//            if (tdt.is_infinite(fit)) continue;  // Skip infinite faces
//            Delaunay::Face_handle face = fit;
//            // std::cout << "Triangle: ";
//            ++trnls;
//            Delaunay::Point p1 = fit->vertex(0)->point();
//            Delaunay::Point p2 = fit->vertex(1)->point();
//            Delaunay::Point p3 = fit->vertex(2)->point();
//            //std::cout << "(" << p1.x() << ", " << p1.y() << ") ";
//            //std::cout << "(" << p2.x() << ", " << p2.y() << ") ";
//            //std::cout << "(" << p3.x() << ", " << p3.y() << ") ";
//            std::string temp = std::to_string(p1.x()) + " " + std::to_string(p1.y()) + "\n" +
//                std::to_string(p2.x()) + " " + std::to_string(p2.y()) + "\n" +
//                std::to_string(p3.x()) + " " + std::to_string(p3.y()) + "\n" +
//                std::to_string(p1.x()) + " " + std::to_string(p1.y()) + "\n";
//            //str = str + temp;
//            //str = str + "\n\n";
//            //double area = 0.5 * std::abs((p2 - p1) * (p3 - p1));
//           // std::cout << " area:" << areaOfTriangle(p1, p2, p3);
//            /*for (int i = 0; i < 3; ++i) {
//                Point vertex = face->vertex(i)->point();
//                std::cout << "(" << vertex.x() << ", " << vertex.y() << ") ";
//
//                std::string temp= to_string(vertex.x())+" "+to_string(vertex.y())+"\n";
//                str = str + temp;
//
//
//            }*/
//
//            // str = str + "\n\n";
//           // std::cout << std::endl;
//
//
//        }
//
//
//        // outer boundary
//        for (const auto& ed : shape) {
//            //std::cout << " {" << ed.source << ", " << ed.target << "}, ";
//
//            std::string temp = std::to_string(ed.source.x()) + " " + std::to_string(ed.source.y()) + "\n" +
//                std::to_string(ed.target.x()) + " " + std::to_string(ed.target.y()) + "\n";
//            // str = str + temp;
//            // str = str + "\n\n";
//
//        }
//        //str = str + std::to_string(shape.at(0).source.x()) + " " + std::to_string(shape.at(0).source.y()) + "\n" +
//        //    std::to_string(shape.at(0).target.x()) + " " + std::to_string(shape.at(0).target.y()) + "\n";
//       // str = str + "\n\n";
//
//        std::vector<K::Point_2> tconvexhull;
//        CGAL::convex_hull_2(hole.begin(), hole.end(), std::back_inserter(tconvexhull));
//        for (const auto& con : tconvexhull) {
//
//            std::string temp = std::to_string(con.x()) + " " + std::to_string(con.y()) + "\n";
//            // str = str + temp;
//        }
//        // str = str + "\n\n";
//
//
//
//        // std::list<Point> hole;
//
//        Alpha_shape_2 A(hole.begin(), hole.end(), FT(1000), Alpha_shape_2::GENERAL);
//        std::vector<Segment> segments;
//        //alpha_edges(A, segments);
//        ////alpha_edges(A, std::back_inserter(segments));
//        //std::cout << "Alpha Shape computed" << std::endl;
//        //std::cout << segments.size() << " alpha shape edges" << std::endl;
//        //std::cout << "Optimal alpha: " << *A.find_optimal_alpha(1) << std::endl;
//
//        //str = "";
//        //for (auto s : segments) {
//        //    std::string temp = std::to_string(s.source().x()) + " " + std::to_string(s.source().y()) + "\n" +
//        //        std::to_string(s.end().x()) + " " + std::to_string(s.end().y()) + "\n";
//        //    str = str + temp;
//        //}
//        //str = str + "\n\n";
//
//        for (auto f = A.finite_faces_begin(); f != A.finite_faces_end(); ++f) {
//            //std::cout<< *f  << std::endl;
//                        // Delaunay::Face_handle face = f;
//             // std::cout << "Triangle: ";
//            // ++trnls;
//            Delaunay::Point p1 = f->vertex(0)->point();
//            Delaunay::Point p2 = f->vertex(1)->point();
//            Delaunay::Point p3 = f->vertex(2)->point();
//            //std::cout << "(" << p1.x() << ", " << p1.y() << ") ";
//            //std::cout << "(" << p2.x() << ", " << p2.y() << ") ";
//            //std::cout << "(" << p3.x() << ", " << p3.y() << ") ";
//            std::string temp = std::to_string(p1.x()) + " " + std::to_string(p1.y()) + "\n" +
//                std::to_string(p2.x()) + " " + std::to_string(p2.y()) + "\n" +
//                std::to_string(p3.x()) + " " + std::to_string(p3.y()) + "\n" +
//                std::to_string(p1.x()) + " " + std::to_string(p1.y()) + "\n";
//            // str = str + temp;
//             //str = str + "\n\n";
//            //std::string temp = std::to_string(f->point().x()) + " " + std::to_string(f->point().y()) + "\n";
//            //str = str + temp;
//
//        }
//        // str = str + "\n\n";
//        std::vector<Point> newInput;
//        for (auto vh = A.all_vertices_begin(); vh != A.all_vertices_end(); ++vh) {
//            newInput.push_back(vh->point());
//        }
//
//        Delaunay dtN;
//        dtN.insert(newInput.begin(), newInput.end());
//        outerBoundary(dtN, exists, shape, pqOuterBoundary);
//        //holeDetection(dtN, exists, shape, holeEdges, pqInnerBoundary, tpholes, numberOfHoles);
//
//
//        //auto cp = CGAL::centroid(hole.begin(),hole.end(), CGAL::Dimension_tag<0>());
//        //CGAL::centroid()
//
//        //str = str + std::to_string(cp.x()) + " " + std::to_string(cp.y()) + "\n";
//       // str = str + std::to_string(cp.x()) + " " + std::to_string(cp.y()) + "\n";
//       // str = str + "\n\n";
//
//
//        //holedetection for one hole only
//        if (nh < 2)
//        {
//            singleHoleDetection(dtN, shape, holeEdges, tpholes);
//
//            std::string tstr = "";
//            for (auto a : tpholes) {
//                for (auto b : a) {
//                    std::string temp = std::to_string(b.x()) + " " + std::to_string(b.y()) + "\n";
//                    tstr = tstr + temp;
//                }
//                tstr = tstr + "\n";
//            }
//            try {
//                std::string tfilename = "E:\\Avinash C drive\\Downloads\\gp610-20230826-win64-mingw\\gnuplot\\bin\\sampleHoleIn" + std::to_string(nh) + ".txt";
//                std::ofstream tfile(tfilename);
//
//                if (!tfile.is_open())
//                    throw "Error to Open File !!!" + tfilename;
//                /*cout << "Error to Open File !!!" << endl;
//                exit(1);*/
//                //std::vector<K::Point_2> tconvexhull;
//                // CGAL::convex_hull_2(hole.begin(), hole.end(), std::back_inserter(tconvexhull));
//                 //std::getline(std::cin, str);
//                if (tstr.empty())
//                    throw "Write the file atleast one char !!!";
//                tfile.write(tstr.c_str(), tstr.length());
//                tfile.close();
//
//                std::cout << "\n hole file saved to " << tfilename << std::endl;
//
//
//            }
//            catch (...) {
//                std::cout << "\n ERROR!!! 1208" << std::endl;
//            }
//
//        }
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//        for (auto f = dtN.finite_faces_begin(); f != dtN.finite_faces_end(); ++f) {
//            //std::cout<< *f  << std::endl;
//                        // Delaunay::Face_handle face = f;
//             // std::cout << "Triangle: ";
//            // ++trnls;
//            Delaunay::Point p1 = f->vertex(0)->point();
//            Delaunay::Point p2 = f->vertex(1)->point();
//            Delaunay::Point p3 = f->vertex(2)->point();
//            //std::cout << "(" << p1.x() << ", " << p1.y() << ") ";
//            //std::cout << "(" << p2.x() << ", " << p2.y() << ") ";
//            //std::cout << "(" << p3.x() << ", " << p3.y() << ") ";
//            std::string temp = std::to_string(p1.x()) + " " + std::to_string(p1.y()) + "\n" +
//                std::to_string(p2.x()) + " " + std::to_string(p2.y()) + "\n" +
//                std::to_string(p3.x()) + " " + std::to_string(p3.y()) + "\n" +
//                std::to_string(p1.x()) + " " + std::to_string(p1.y()) + "\n";
//            str = str + temp;
//            str = str + "\n\n";
//            //std::string temp = std::to_string(f->point().x()) + " " + std::to_string(f->point().y()) + "\n";
//            //str = str + temp;
//
//        }
//        str = str + "\n\n";
//
//
//
//        try {
//            std::string tfilename = "E:\\Avinash C drive\\Downloads\\gp610-20230826-win64-mingw\\gnuplot\\bin\\sampleHole" + std::to_string(nh) + ".txt";
//            std::ofstream tfile(tfilename);
//
//            if (!tfile.is_open())
//                throw "Error to Open File !!!" + tfilename;
//            /*cout << "Error to Open File !!!" << endl;
//            exit(1);*/
//            //std::vector<K::Point_2> tconvexhull;
//            // CGAL::convex_hull_2(hole.begin(), hole.end(), std::back_inserter(tconvexhull));
//             //std::getline(std::cin, str);
//            if (str.empty())
//                throw "Write the file atleast one char !!!";
//            tfile.write(str.c_str(), str.length());
//            tfile.close();
//
//            std::cout << "\n hole file saved to " << tfilename << std::endl;
//
//
//        }
//        catch (...) {
//            std::cout << "\n ERROR!!! 1208" << std::endl;
//        }
//        ++nh;
//        std::cout << " ............................................................................... " << std::endl;
//
//    }
//
//}
//
//
//void binarayMatrixOfEstimatedGridPixel(const std::vector<Point>& points, double minX, double minY, double maxX, double maxY, double m, bool** binMat, int row, int col) {
//    //int** binMat{ nullptr };
//    //double**** mapingToPoint;
//    //binMat = new int* [row];
//
//
//    for (auto p : points) {
//        for (int i = 0; i < row; ++i) {
//
//
//            for (int j = 0; j < col && i != row - 1 && j != col - 1; ++j) {
//                auto x = (i * m) + minX;
//                auto y = maxY - (j * m);
//                if (((x <= p.x()) && (x + m >= p.x())) && ((y - m <= p.y()) && (y >= p.y()))) {
//                    binMat[i][j] = 1;
//                    // pt.push_back(p);
//                }
//            }
//
//
//        }
//    }
//
//
//
//
//
//
//
//
//
//
//    //mapingToPoint = new double*** [row];
//    //for (int i = 0; i < row; ++i) {
//    //    mapingToPoint[i] = new double** [col];
//    //    for (int j = 0; j < col; j++) { binMat[i][j] = 0; }
//
//
//    //    for (int j = 0; j < col; ++j) {
//
//    //    }
//    //}
//}
//
//
//void nearestNeighbour(const std::vector<Point>& points) {
//
//    Delaunay dt;
//    std::vector<double> minDistances;
//    dt.insert(points.begin(), points.end());
//    double minX = std::numeric_limits<int>::max(), maxX = std::numeric_limits<int>::min();
//    double minY = std::numeric_limits<int>::max(), maxY = std::numeric_limits<int>::min();
//    for (Delaunay::Vertex_iterator vi = dt.vertices_begin(); vi != dt.vertices_end(); vi++)
//    {
//        if (vi->point().x() > maxX)
//        {
//            maxX = vi->point().x();
//        }
//        if (vi->point().y() > maxY)
//        {
//            maxY = vi->point().y();
//        }
//        if (vi->point().x() < minX)
//        {
//            minX = vi->point().x();
//        }
//        if (vi->point().y() < minY)
//        {
//            minY = vi->point().y();
//        }
//    }
//
//    std::string str = "";
//    for (Delaunay::Vertex_iterator vi = dt.vertices_begin(); vi != dt.vertices_end(); vi++)
//    {
//        //Delaunay::Vertex_iterator nbv=Delaunay::nearest_vertex(vi,nbv);
//        auto minDist = distance(Point(minX, minY), Point(maxX, maxY));
//        for (Delaunay::Vertex_iterator vj = dt.vertices_begin(); vj != dt.vertices_end(); vj++)
//        {
//            if (vi != vj) {
//                auto d = distance(vi->point(), vj->point());
//
//                //if()
//                auto t = minDist - d;
//                if (d < minDist) {
//                    minDist = d;
//                }
//
//            }
//        }
//        //std::cout << vi->point()<<":" << minDist << std::endl;
//        minDistances.push_back(minDist);
//        str = str + std::to_string(minDist) + "\n";
//
//    }
//
//
//
//    try {
//        std::string tfilename = "E:\\Avinash C drive\\Downloads\\gp610-20230826-win64-mingw\\gnuplot\\bin\\freqDist.txt";
//        std::ofstream tfile(tfilename);
//
//        if (!tfile.is_open())
//            throw "Error to Open File !!!" + tfilename;
//        /*cout << "Error to Open File !!!" << endl;
//        exit(1);*/
//        //std::vector<K::Point_2> tconvexhull;
//        // CGAL::convex_hull_2(hole.begin(), hole.end(), std::back_inserter(tconvexhull));
//         //std::getline(std::cin, str);
//        if (str.empty())
//            throw "Write the file atleast one char !!!";
//        tfile.write(str.c_str(), str.length());
//        tfile.close();
//
//        std::cout << "\n hole file saved to " << tfilename << std::endl;
//
//
//    }
//    catch (...) {
//        std::cout << "\n ERROR!!! 1208" << std::endl;
//    }
//
//
//
//    auto m = distance(Point(minX, minY), Point(maxX, maxY));
//    auto min = std::min(minDistances.begin(), minDistances.end());
//    for (auto it = minDistances.begin(); it != minDistances.end(); ++it) {
//        auto d = double(*it) - m;
//        if (*it < m) {
//            m = *it;
//
//        }
//
//    }
//
//    std::cout << *min << " :" << m << std::endl;
//    // m = *min;
//     //2d binary matrix
//
//    auto r = CGAL::abs(maxY - minY) / m;
//    // ceil(2.3);
//    // std::ceil(2.3);
//    const int row = std::ceil(CGAL::abs(maxY - minY) / m);
//    const int col = std::ceil(CGAL::abs(maxX - minX) / m);
//    //Matrix M(row, col);
//
//    std::cout << "row: " << CGAL::abs(maxY - minY) / m << " col:" << CGAL::abs(maxX - minX) / m << std::endl;
//    std::cout << "row: " << std::ceil(CGAL::abs(maxY - minY) / m) << " col:" << std::ceil(CGAL::abs(maxX - minX) / m) << std::endl;
//
//
//    // std::vector<std::pair<std::pair<int, int>, std::vector<int>>> matToPoints;
//     //populate the binary matrix
//
//    bool** binMat{ nullptr };
//    binMat = new bool* [row];
//    for (int i = 0; i < row; ++i) {
//        binMat[i] = new bool[col];
//    }
//
//    for (int i = 0; i < row; ++i) {
//        for (int j = 0; j < col; ++j) {
//            binMat[i][j] = 0;
//        }
//    }
//
//
//    binarayMatrixOfEstimatedGridPixel(points, minX, minY, maxX, maxY, m, binMat, row, col);
//
//    for (int i = 0; i < row; ++i) {
//        for (int j = 0; j < col; ++j) {
//
//            if (binMat[i][j]) { std::cout << "+"; }
//            else {
//                std::cout << " ";
//            }
//        }
//        std::cout << "\n";
//    }
//}
//
//
//int indxPoint(const std::vector<Point>& points, Point p) {
//    int inx = 0;
//    int cnt = 0;
//    for (auto a : points) {
//        auto d = distance(p, a);
//        if (d == 0) {
//            inx = cnt;
//            break;
//        }
//        ++cnt;
//    }
//    return inx;
//}
