lc = 500;

xmin = -5000.0;
xmax =  5000.0;
ymin =     0.0;
ymax = 10000.0;

Point(1) = {xmin, ymin, 0, lc};
Point(2) = {xmax, ymin, 0, lc};
Point(3) = {xmax, ymax, 0, lc};
Point(4) = {xmin, ymax, 0, lc};

x5 = -1000;
y5 =  5000;

x6 = 1000;
y6 = 5000;

x7 = 1000;
y7 = 5500;

x8 = -1000;
y8 =  5500;

Point(5) = {x5, y5, 0, lc};
Point(6) = {x6, y6, 0, lc};
Point(7) = {x7, y7, 0, lc};
Point(8) = {x8, y8, 0, lc};


Line(1) = {1,2};
Line(2) = {3,2};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line Loop(1) = {4,1,-2,3} ;
Line Loop(2) = {5,6,7,8} ;

Plane Surface(1) = {1,2} ;
Recombine Surface {1};

Physical Point("free_slip", 1) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Curve("free_slip", 2) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Surface("domain") = {1};