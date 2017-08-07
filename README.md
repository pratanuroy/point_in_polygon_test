# point_in_polygon_test
A C++ program to test if a cell in 2D cartesian domain resides inside a polygon, outside a polygon or on the polygon. 

Description:

A rectangular two-dimensional domain is populated by a Cartesian mesh with identical spacing in the X
and Y directions. The left bottom corner of the domain is located at (0.0, 0.0). An arbitrary shaped
polygon described as a sequence of vertices ordered in anti-clockwise direction is overlaid on this
Cartesian mesh. The polygon is completely contained in the rectangular domain. This program finds the
location of the points with respect to the polygon and mark the cells of the Cartesian mesh with location attribute described in the following manner.
- Completely inside the polygon: I
- Intersecting the polygon segments: X
- Completely outside the polygon: O

Algorithm:

1. First find out the bounding box using the maximum and minimum (x,y) values of the Polygon coordinates.
2. Traverse through all the 2D Cartesian cells within the bounding box. For each cell do the following:
a. Use a line clipping algorithm to find out if any edge of the Polygon is intersecting the cell. I used Liang-Barsky line clipping algorithm in this case.
b. If the cell is not intersecting any of the edges, check whether the cell center lies inside the polygon or outside the polygon based on the "Point in Polygon" algorithm. If the cell center lies inside the polygon that means the cell is completely inside the polygon. Otherwise, the cell is completely outside. Flag the cells accordingly.
3. Write output using the cell flags.
