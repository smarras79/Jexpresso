// Mesh size - adjusted for the new scale
lc = 0.1; 

// Domain boundaries
xmin =  0.0;
xmax =  3.0;
ymin = -1.0;
ymax =  1.0;

Point(1) = {xmin, ymin, 0, lc};
Point(2) = {xmax, ymin, 0, lc};
Point(3) = {xmax, ymax, 0, lc};
Point(4) = {xmin, ymax, 0, lc};

// Circle parameters (scaled to fit inside the 2x2 box)
radius = 0.2;
xc = 1.0;
yc = 0.0;

Point(5) = {xc, yc, 0, lc};      // Center
Point(6) = {xc + radius, yc, 0, lc};
Point(7) = {xc, yc + radius, 0, lc};
Point(8) = {xc - radius, yc, 0, lc};
Point(9) = {xc, yc - radius, 0, lc};

// Square Boundary
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Circle Boundary
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// Loops
// Outer loop should be counter-clockwise; inner loop (hole) clockwise
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1, 2};

// Recombine to create unstructured quads
Recombine Surface {1};

// Physical Groups
Physical Curve("bottom", 1) = {1};
Physical Curve("outflow", 2) = {2};
Physical Curve("top", 3) = {3};	
Physical Curve("inflow", 4) = {4};
Physical Curve("circle_boundary", 5) = {5, 6, 7, 8};
Physical Surface("domain") = {1};