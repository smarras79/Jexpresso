lc = 250;

xmin = 0.0;
xmax = 10000.0;
ymin = 0.0;
ymax = 6400.0;

Point(1) = {xmin, ymin, 0, lc};
Point(2) = {xmax, ymin, 0, lc};
Point(3) = {xmax, ymax, 0, lc};
Point(4) = {xmin, ymax, 0, lc};

x5 = 1000;
y5 = 1000;
x6 = 2500;
y6 = 1000;
x7 = 1250;
y7 = 2000;

Point(5) = {x5, y5, 0, lc};
Point(6) = {x6, y6, 0, lc};
Point(7) = {x7, y7, 0, lc};

Line(1) = {1,2};
Line(2) = {3,2};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,5};

Line Loop(1) = {4,1,-2,3} ;
Line Loop(2) = {5,6,7} ;

Plane Surface(1) = {1,2} ;
Recombine Surface {1};

Physical Point("zero_all",   1) = {1, 2, 3, 4, 5, 6, 7};
Physical Curve("zero_all",   2) = {1, 2, 3, 4, 5, 6, 7};
Physical Surface("domain") = {1};