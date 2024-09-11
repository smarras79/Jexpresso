lc = 400;

xmin = -3000.0;
xmax =  3000.0;
ymin =     0.0;
ymax = 10000.0;

Point(1) = {xmin, ymin, 0, lc};
Point(2) = {xmax, ymin, 0, lc};
Point(3) = {xmax, ymax, 0, lc};
Point(4) = {xmin, ymax, 0, lc};

radius = 1000.0;

xc =  0.0;
yc =  5000.0;

x6 = xc + radius;
y6 = yc;

x7 =  xc;
y7 =  yc + radius;

x8 =  xc - radius;
y8 =  yc;

x9 =  xc;
y9 =  yc - radius;

Point(5) = {xc, yc, 0, lc};
Point(6) = {x6, y6, 0, lc};
Point(7) = {x7, y7, 0, lc};
Point(8) = {x8, y8, 0, lc};
Point(9) = {x9, y9, 0, lc};


Line(1) = {1,2};
Line(2) = {3,2};
Line(3) = {3,4};
Line(4) = {4,1};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Line Loop(1) = {4,1,-2,3};
Line Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1,2} ;
Recombine Surface {1};

Physical Point("free_slip", 1) = {1, 2, 3, 4, 6, 7, 8, 9};
Physical Curve("free_slip", 2) = {1, 2, 3, 4};
//Physical Point("free_slip", 1) = {1, 2, 3, 4, 6, 7, 8, 9};
Physical Curve("circle", 3) = {5, 6, 7, 8};
Physical Surface("domain") = {1};