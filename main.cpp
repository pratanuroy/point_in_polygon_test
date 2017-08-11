//=========================================================================================================
//*********************************************************************************************************
//
//  Purpose:
//
//    MAIN is the main program for polygon_in_point_test
//
//  Discussion:
//
//    A rectangular two-dimensional domain is populated by a Cartesian mesh with identical spacing in the X
//    and Y directions. The left bottom corner of the domain is located at (0.0, 0.0). An arbitrary shaped
//    polygon described as a sequence of vertices ordered in anti-clockwise direction is overlaid on this
//    Cartesian mesh. The polygon is completely contained in the rectangular domain. This program find the
//    location of the points with respect to the polygon and mark the cells of the Cartesian mesh with a
//    location attribute described in the following manner.
//		- Completely inside the polygon: I
//		- Intersecting the polygon segments: X
//		- Completely outside the polygon: O
//
//  Modified:
//
//    August 5, 2017
//
//  Author:
//
//    Pratanu Roy
//
//*********************************************************************************************************
//=========================================================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include "point.h"

using std::vector;
using std::cout;
using std::cin;
using std::endl;
using std::min;
using std::max;

#define DEBUG 0
#define ONSEGMENT 999

#define TOL 1e-8
#define INTERSECTION 1
#define NO_INTERSECTION 0
#define INSIDE 2
#define OUTSIDE 3


// getCoords function calculates the nodal and cell center values of the coordinates

void getCoords(int nx, int ny, double spacing,
               vector<double> &x_coord, vector<double> &y_coord,
               vector<double> &x_cellcoord, vector<double> &y_cellcoord)
{
    
    x_cellcoord.resize(nx+2);
    y_cellcoord.resize(ny+2);
    
    x_coord[0] = 0.0;
    y_coord[0] = 0.0;
    
    for (int i=1; i<nx+1; ++i)
    {
        x_coord[i] = x_coord[i-1] + spacing;
    }
    
    for (int j=1; j<ny+1; ++j)
    {
        y_coord[j] = y_coord[j-1] + spacing;
    }
    
    x_cellcoord[0] = 0.0;
    y_cellcoord[0] = 0.0;
    
    x_cellcoord[1] = x_cellcoord[0] + spacing/2.0;
    y_cellcoord[1] = y_cellcoord[0] + spacing/2.0;
    
    for (int i=2; i<nx+1; i++)
    {
        x_cellcoord[i] = x_cellcoord[i-1] + spacing;
    }
    
    for (int j=2; j<ny+1; j++)
    {
        y_cellcoord[j] = y_cellcoord[j-1] + spacing;
    }
    
    x_cellcoord[nx+1] = x_cellcoord[nx] + spacing/2.0;
    y_cellcoord[ny+1] = y_cellcoord[ny] + spacing/2.0;
    
}

// isOutsideBoundingBox determines if a point is outside the bounding box of polygon
bool isOutsideBoundingBox(const Point &point, const Point &minPoint, const Point &maxPoint)
{
    double xvalue = point.X;
    double yvalue = point.Y;
    
    if (DEBUG == 2) cout << xvalue << " " << yvalue << " " << minPoint.X<<" "<<minPoint.Y<< " "<< maxPoint.X<<" "<<maxPoint.Y<<" ";
    
    if(xvalue < minPoint.X || yvalue < minPoint.Y || xvalue > maxPoint.X || yvalue > maxPoint.Y)
    {
        return true;
    }
    return false;
}

// Compare the x or y values of two points and returns the minimum value

bool _compare_min_x(Point const &p1, Point const &p2) { return p1.X < p2.X; }
bool _compare_min_y(Point const &p1, Point const &p2) { return p1.Y < p2.Y; }


// getBoundingBox gives the minimum and maximum x and y-values for coordinates to determine the bounding box of polygon

void getBoundingBox(vector<Point> &points, Point &minPoint, Point &maxPoint)
{
    if (DEBUG) cout << "size of Polygon : " << points.size() << endl;
    
    if(points.size()>1) {
        
        double min_x = (*std::min_element(points.begin(), points.end(), &_compare_min_x)).X;
        double min_y = (*std::min_element(points.begin(), points.end(), &_compare_min_y)).Y;
        
        double max_x = (*std::max_element(points.begin(), points.end(), &_compare_min_x)).X;
        double max_y = (*std::max_element(points.begin(), points.end(), &_compare_min_y)).Y;
        
        minPoint.setvalues(min_x,min_y);
        maxPoint.setvalues(max_x,max_y);
    }
    
}

// orientation function determines the orientation of an ordered triplet (p1, p2, p3).
// The function returns following values
// 0 if p1, p2 and p3  are Colinear
// 1 if Clockwise
// -1 if Counterclockwise
// Theory:
// Slope of segment p1,p2: s12= (p2.Y-p1.Y)/(p2.X - p1.X)
// Slope of segment p2,p3: s23 = (p3.Y - p2.Y)/(p3.X - p2.X)
// Counterclockwise => s12 < s23
// Clockwise => s12  > s23
// Colinear => s12 = s23

int orientation(const Point& p1, const Point& p2, const Point& p3) {
    int val = (p2.Y - p1.Y) * (p3.X - p2.X) - (p3.Y - p2.Y) * (p2.X - p1.X);
    if (val == 0)
        return 0;
    else
        return (val < 0) ? -1 : 1;
}

// Returns true if q lies on line p1p2
bool onSegment(const Point& p1, const Point& p2, const Point& q) {
    if (min(p1.X, p2.X) <= q.X+TOL && q.X-TOL <= max(p1.X, p2.X)
        && min(p1.Y, p2.Y) <= q.Y+TOL && q.Y-TOL <= max(p1.Y, p2.Y))
        return true;
    else
        return false;
}

// The intersectionTest function returns true if line p1p2 intersects with line p3p4

bool intersectionTest(const Point& p1, const Point& p2,
                      const Point& p3, const Point& p4) {
    
    // Four orientations required for the general and special cases
    int o1 = orientation(p1, p2, p3);
    int o2 = orientation(p1, p2, p4);
    int o3 = orientation(p3, p4, p1);
    int o4 = orientation(p3, p4, p2);
    
    // General case
    if (o1 != o2 && o3 != o4)
        return true;
    
    // Special cases
    if (o1 == 0 && onSegment(p1, p2, p3))
        return true;
    if (o2 == 0 && onSegment(p1, p2, p4))
        return true;
    if (o3 == 0 && onSegment(p3, p4, p1))
        return true;
    if (o4 == 0 && onSegment(p3, p4, p2))
        return true;
    
    return false;
}

// Liang-Barsky function for line clipping algorithm
// This function inputs 4 edges of a rectangular box, and 2 points of a line (outputs a boolean value to say whether the line is clipped or not).
// https://en.wikipedia.org/wiki/Liang%E2%80%93Barsky_algorithm

bool LiangBarsky(double edgeLeft, double edgeRight, double edgeBottom, double edgeTop,   // Define the x/y clipping values for the border.
                 double x0src, double y0src, double x1src, double y1src)                 // Define the start and end points of the line.
{
    
    double t0 = 0.0;    double t1 = 1.0;
    double xdelta = x1src-x0src;
    double ydelta = y1src-y0src;
    double p,q,r;
    
    for(int edge=0; edge<4; edge++) {   // Traverse through left, right, bottom, top edges.
        if (edge==0) {  p = -xdelta;    q = -(edgeLeft-x0src);  }
        if (edge==1) {  p = xdelta;     q =  (edgeRight-x0src); }
        if (edge==2) {  p = -ydelta;    q = -(edgeBottom-y0src);}
        if (edge==3) {  p = ydelta;     q =  (edgeTop-y0src);   }
        r = q/p;
        if(p==0.0 && q<0.0) return false;   // line is outside
        
        if(p<0.0) {
            if(r>t1) return false;         // line is outside
            else if(r>t0) t0=r;            // Line is clipped!
        } else if(p>0.0) {
            if(r<t0) return false;      // line is outside
            else if(r<t1) t1=r;         // Line is clipped!
        }
    }
    
    return true; // Line intersects the box
}

// This function returns ONSEGMENT if the point is on the polygon segment, 1 if the point is
// completely inside the polygon, and 0 if the point is completely outside the polygon
// The concepts of orientation and intersection can be found in:
// http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf

int checkLineIntersection(Point& p_lb, Point& p_lu, Point& p_rb, Point& p_ru, vector<Point>& polygon)
{
    bool checkIntersection = 0;
    int cellFlag = 0;
    int i = 0;
    int j = i+1;
    
    do{
        double edgeLeft = p_lb.X;
        double edgeRight = p_rb.X;
        double edgeBottom = p_lb.Y;
        double edgeTop = p_lu.Y;
        
        double x0src = polygon[i].X;
        double y0src = polygon[i].Y;
        double x1src = polygon[j].X;
        double y1src = polygon[j].Y;
        
        checkIntersection = LiangBarsky (edgeLeft, edgeRight, edgeBottom, edgeTop,x0src, y0src, x1src, y1src);
        
        if(checkIntersection) {cellFlag = INTERSECTION; return cellFlag;}
        
        i = ((i + 1) % polygon.size());
        j = ((j + 1) % polygon.size());
    }while(i!=0);
    
    // Calculate cell center and check if it lies inside or outside the polygon
    Point p_cent;
    p_cent.X = (p_lu.X+p_ru.X+p_lb.X+p_rb.X)/4.0;
    p_cent.Y = (p_lu.Y+p_ru.Y+p_lb.Y+p_rb.Y)/4.0;
    double inf = 1000000.0;
    Point PtoInfty(inf , p_cent.Y);
    int intersectionsCount = 0;
    i = 0; j = i + 1;
    
    if(!checkIntersection){
        do{
            if (intersectionTest(p_cent, PtoInfty, polygon[i], polygon[j]) == true) {
                ++intersectionsCount;
                if (orientation(polygon[i], polygon[j], p_cent) == 0) { // Colinear
                    if (onSegment(polygon[i], polygon[j], p_cent) == true){
                        
                        cellFlag = ONSEGMENT;
                     
                    }
                    else {
                        // Exception case when point is colinear but not on segment
                        // The colinear segment is worth nothing if they have the same
                        // vertical direction
                        
                        int k = (((i - 1) >= 0) ? // Negative wrap-around
                                 (i - 1) % static_cast<int>(polygon.size()):static_cast<int>(polygon.size()) + (i - 1));
                        
                        int w = ((j + 1) % polygon.size());
                        
                        if ((polygon[k].Y <= polygon[i].Y && polygon[w].Y <= polygon[j].Y)
                            || (polygon[k].Y >= polygon[i].Y && polygon[w].Y >= polygon[j].Y))
                            --intersectionsCount;
                    }
                }
            }
            i = ((i + 1) % polygon.size());
            j = ((j + 1) % polygon.size());
            
        } while (i != 0);
        if(intersectionsCount % 2 != 0) cellFlag = INSIDE;
        else cellFlag = OUTSIDE;
    }
    
    return cellFlag;
}


void writeOutput(vector<Point> &polygon, int nx, int ny, double s, vector<double> &xcoord,vector<double> &ycoord,
                 vector<double> &xcellcoord, vector<double> &ycellcoord, Point &minPoint, Point  &maxPoint)
{
    int numVertices = polygon.size();
    
    if (DEBUG >=1)
    {
        std::cout << s <<" "<< nx <<" "<< ny <<" "<< std::endl;
        std::cout << numVertices << std::endl;
        for (int i=0; i<numVertices; ++i)
        {
            polygon[i].write();
        }
    }
    
    if(DEBUG>=1)
    {
        cout<<"Bounding box : "<<endl;
        minPoint.write();
        maxPoint.write();
    }
    
    bool isOutsideBB=false;
    int isInsidePolygon=0;
    
    Point p_lb, p_rb, p_lu, p_ru;
    
    // Loop through the cell coordinates in x and y directions and determine the
    // location of each cell coordinate with respect to the polygon.
    // The cells are marked in following manner:
    // Completely inside the polygon: I
    // Intersecting the polygon segments: X
    // Completely outside the polygon: O
    // For each case there are ny lines in the output and each line contains nx entries,
    // where each entry represents the location attribute of a cell.
    
    
    for(int j=ny; j>0; j--)
    {
        for(int i=1; i<nx+1; i++)
        {
            // Get four vertices of cell (i,j)
            
            p_lb.setvalues(xcoord[i-1],ycoord[j-1]);
            p_rb.setvalues(xcoord[i],ycoord[j-1]);
            p_lu.setvalues(xcoord[i-1],ycoord[j]);
            p_ru.setvalues(xcoord[i],ycoord[j]);
            
            // First test if the cell is completely outside the bounding box.
            //If not, then check if the cell is inside or outside or intersected by the polygon.
            
            isOutsideBB = isOutsideBoundingBox(p_lb, minPoint, maxPoint) &&
                          isOutsideBoundingBox(p_rb, minPoint, maxPoint) &&
                          isOutsideBoundingBox(p_lu, minPoint, maxPoint) &&
                          isOutsideBoundingBox(p_ru, minPoint, maxPoint);
            
            if(isOutsideBB){cout <<"O";}
            
            else
            {
                int cellFlag = checkLineIntersection(p_lb, p_lu, p_rb, p_ru, polygon);
                
                if(cellFlag == INTERSECTION) {cout<< "X";}
                else if(cellFlag == INSIDE) {cout<< "I";}
                else {cout<< "O";}
            }
        }
        cout<<endl;
        
    }
    
}

int main()
{
    // Number of tests
    int numTest = 0;
    
    // Number of cells in x and y directions
    int nx=20;
    int ny=20;
    
    // Number of vertices in polygon
    int nv;
    
    // Spacing of Cartesian mesh along X and Y directions
    double s=1.0;
    
    // Vector arrays for node coordinates and cell coordinates
    vector<double> xcoord;
    vector<double> ycoord;
    
    vector<double> xcellcoord;
    vector<double> ycellcoord;
    
    // Points for constructing bounding box with minimum and maximum coordinates of polygon
    Point minPoint;
    Point maxPoint;
    
    // Reading input file
    
    std::ifstream infile("input.txt");
    
    // check if file opened successfully
    if (!infile) {
        std::cerr << "Failure!! can't open file" << endl;
        cin.get();
        return EXIT_FAILURE;
    }
    
    // the container in which we will store all the points for a polygon
    vector<Point> polygonPoints;
    
    // a temporary point to hold the coords when we read them from the input file
    Point tmp;
    
    if (DEBUG) cout << "Reading input file..." << endl;
    
    // Read number of tests
    infile >> numTest;
    
    if (DEBUG) cout << numTest << endl;
    
    // Read inputs for specified number of tests
    
    for(int i = 0; i<numTest;++i)
    {
        // Read spacing and number of cells of the rectangular domain
        infile >> s >> nx >> ny;
        
        // Read the number of vertices of the polygon
        infile >> nv;
        
        // Read Nv lines of location of vertices of polygon and fill up the polygonPoints container accordingly
        for(int ivert = 0; ivert<nv; ++ivert)
            
        {
            infile >> tmp.X >> tmp.Y;
            
            polygonPoints.push_back(tmp);
        }
        
        // Resize node and cell coordinates with rectangular domain size information
        
        xcoord.resize(nx+1);
        ycoord.resize(ny+1);
        
        xcellcoord.resize(nx+2);
        ycellcoord.resize(ny+2);
        
        // Get node and cell coordinates of rectangular domain.
        
        getCoords(nx,ny,s,xcoord,ycoord,xcellcoord,ycellcoord);
        
        // Get minimum and maximum points of polygon
        getBoundingBox(polygonPoints, minPoint, maxPoint);
        
        // Determine the location of cell coordinates and write output
        writeOutput(polygonPoints, nx, ny, s, xcoord, ycoord, xcellcoord, ycellcoord, minPoint, maxPoint);
        
        // Set the size of the polygon to zero
        polygonPoints.resize(0);
    }
    
    return 0;
}





