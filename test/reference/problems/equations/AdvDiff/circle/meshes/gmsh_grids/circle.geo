rx = 2.0;
ry = 2.0;
xmin = -1.0;
xmax =  1.0;
ymin = -1.0;
ymax =  1.0;

lc = 0.2;
Point(100) = {0, 0, 0, lc};
Point(1) = {xmin, ymin, 0, lc};
Point(2) = {xmax, ymin, 0, lc};
Point(3) = {xmax, ymax, 0, lc};
Point(4) = {xmin, ymax, 0, lc};
Point(5) = {-rx, -ry, 0, lc};
Point(6) = { rx, -ry, 0, lc};
Point(7) = { rx,  ry, 0, lc};
Point(8) = {-rx,  ry, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {5, 100, 6};
Circle(6) = {6, 100, 7};
Circle(7) = {7, 100, 8};
Circle(8) = {8, 100, 5};
Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Transfinite Line {1,2,3,4, 5,6,7,8,9,10,11,12} = 10 Using Progression 1.0;

//+
Curve Loop(13) = {4, 1, 2, 3};
Plane Surface(14) = {13};

//+
Curve Loop(15) = {-1, 9, 5, -10};
Plane Surface(16) = {15};

//+
Curve Loop(17) = {-2, 10, 6, -11};
Plane Surface(18) = {17};

//+
Curve Loop(19) = {-3, 11, 7, -12};
Plane Surface(20) = {19};

//+
Curve Loop(21) = {-4, 12, 8, -9};
Plane Surface(22) = {21};

Transfinite Surface {14,16,18,20,22};
Recombine Surface {14,16,18,20,22};
Coherence;

Physical Point("boundary",   1) = {5, 6, 7, 8};
Physical Curve("nothing",  2) = {7, 5};
Physical Curve("nothing",  3) = {8, 6};
Physical Surface("domain") = {14,16,18,20,22};