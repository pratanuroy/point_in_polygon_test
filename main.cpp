//=========================================================================================================
//*********************************************************************************************************
//
//  Purpose:
//
//    MAIN is the main program for polygon_in_point_test
//
//  Discussion:
//
// 	A rectangular two-dimensional domain is populated by a Cartesian mesh with identical spacing in the X
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
//    July 31, 2017
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
bool isOutsideBoundingBox(double xvalue, double yvalue, const Point &minPoint, const Point &maxPoint)
{
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
// Slope of segment p1,p2: s23 = (p3.Y - p2.Y)/(p3.X - p2.X)
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

// This function returns ONSEGMENT if the point is on the polygon segment, 1 if the point is
// completely inside the polygon, and 0 if the point is completely outside the polygon
// The concepts of orientation and intersection can be found in:
// http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf

int checkPointInPolygon(const Point& p, const vector<Point>& polygon) {

	if (polygon.size() < 3)
		return 0; // As the minimum number of vertiecs should be 3 to construct a polygon

	// Constructing a horizontal line of inifinite length
	double inf = 10000.0;
	Point PtoInfty(inf , p.Y);

	int intersectionsCount = 0;
	int i = 0, j = i + 1;

	do {

		if (intersectionTest(p, PtoInfty, polygon[i], polygon[j]) == true) {

			++intersectionsCount;

			if (orientation(polygon[i], polygon[j], p) == 0) { // Colinear
				if (onSegment(polygon[i], polygon[j], p) == true)
					return ONSEGMENT;
				else {
					// Exception case when point is colinear but not on segment
					// for example:
					//					   k      w
					//               /        \
					//           *  ************
					// The colinear segment is worth nothing if k and w have the same
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

	// Return 1 if intersectionsCount is odd, otherwise return 0
	if (intersectionsCount % 2 != 0) return 1;
	else return 0;

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
			Point p;
			p.setvalues(xcellcoord[i],ycellcoord[j]);

			// First test if the coordinate is outside the bounding box. 
			//If not, then check if the point is inside or on the polygon.

			isOutsideBB = isOutsideBoundingBox(xcellcoord[i], ycellcoord[j], minPoint, maxPoint);

			if (!isOutsideBB) isInsidePolygon = checkPointInPolygon( p, polygon);

			if(isOutsideBB){cout <<"O";}
			else if (isInsidePolygon == ONSEGMENT){cout << "X";}
			else if (isInsidePolygon){cout <<"I";}
			else {cout << "O";}

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





