rx = 2.0;
ry = 2.0;
xmin = -1.0;
xmax =  1.0;
ymin = -1.0;
ymax =  1.0;

lc = 0.2;
//Point(100) = {0, 0, 0, lc};
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
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

//+
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
//Transfinite Line {1} = 10 Using Progression 1.;
//Transfinite Line {3} = 10 Using Progression 1.;
//Transfinite Line {2} = 10;
//Transfinite Line {4} = 10;
//Transfinite Surface {1} = {1, 2, 3, 4};

//+
Curve Loop(2) = {-1, 9, 5, -10};
Plane Surface(2) = {2};
//Transfinite Line {1} = 10 Using Progression 1.;
//Transfinite Line {5} = 10 Using Progression 1.;
//Transfinite Line {9} = 10;
//Transfinite Line {10} = 10;
//Transfinite Surface {2} = {1, 5, 9, 10};


//+
Curve Loop(3) = {-2, 10, 6, -11};
Plane Surface(3) = {3};
//Transfinite Line {5} = 10 Using Progression 1.;
//Transfinite Line {7} = 10 Using Progression 1.;
//Transfinite Line {2} = 10;
//Transfinite Line {8} = 10;
//Transfinite Surface {3} = {5, 7, 2, 8};



//+
Curve Loop(4) = {-3, 11, 7, -12};
Plane Surface(4) = {4};
Transfinite Line {11} = 10Using Progression 1.;
//Transfinite Line {7} = 10Using Progression 1.;
//Transfinite Line {3} = 10;
//Transfinite Line {12} = 10;
//Transfinite Surface {4} = {3, 7, 12, 11};

//+
Curve Loop(5) = {-4, 12, 8, -9};
Plane Surface(5) = {5};
//Transfinite Line {8} = 10 Using Progression 1.;
//Transfinite Line {4} = 10;
//Transfinite Line {10} = 10;
//Transfinite Surface {3} = {3,6,4,2};
//Transfinite Surface {5} = {4, 12, 8, 9};


Transfinite Line {1} = 10 Using Progression 1.;
Transfinite Line {2} = 10;
Transfinite Line {3} = 10;
Transfinite Line {4} = 10;

Transfinite Line {5} = 10 Using Progression 1.;
Transfinite Line {6} = 10;
Transfinite Line {7} = 10;
Transfinite Line {8} = 10;

Transfinite Line {9} = 10 Using Progression 1.;
Transfinite Line {10} = 10;
Transfinite Line {11} = 10;
Transfinite Line {12} = 10;

Transfinite Surface (13) = {1};
Transfinite Surface (14) = {2};
Transfinite Surface (15) = {3};
Transfinite Surface (16) = {4};
Transfinite Surface (17) = {5};

Recombine Surface {13, 14, 15, 16, 17};
Coherence;

//Physical Point("nothing",  1) = {2, 7, 3, 8};
//Physical Curve("nothing", 2) = {8, 9, 10, 11};
//Physical Surface("domain") = {1};//+
Show "*";
