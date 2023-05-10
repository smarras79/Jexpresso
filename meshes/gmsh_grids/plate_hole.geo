lc = 0.25;

Point(1) = {0, 0, 0, lc};
Point(2) = {10, 0,  0, lc};
Point(3) = {10, 4, 0, lc};
Point(4) = {0,  4, 0, lc};

Point(5) = {1,  1, 0, lc};
Point(6) = {3,  1, 0, lc};
Point(7) = {2,  2, 0, lc};

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