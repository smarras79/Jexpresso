rx = 2.0;
ry = 2.0;
xmin = -1.0;
xmax =  1.0;
ymin = -1.0;
ymax =  1.0;

lc = 0.2;
Point(100) = {0, 0, 0, lc};
Point(1) = {xmin, 0, 0, lc};
Point(2) = {xmax, 0, 0, lc};
Point(3) = {xmax, ymax, 0, lc};
Point(4) = {xmin, ymax, 0, lc};

Point(5) = {-Sqrt(8), 0, 0, lc};
Point(6) = { Sqrt(8), 0, 0, lc};

Point(7) = {-rx, ry, 0, lc};
Point(8) = {rx, ry, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 1};
Line(6) = {2, 6};
Line(7) = {4, 7};
Line(8) = {3, 8};

Circle(9) = {5, 100, 7};
Circle(10) = {7, 100, 8};
Circle(11) = {8, 100, 6};

//+
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
//Transfinite Line {3} = 10Using Progression 1.;
//Transfinite Line {4} = 10Using Progression 1.;
//Transfinite Line {1} = 10;
//Transfinite Line {2} = 10;
//Transfinite Surface {1} = {1, 4, 6, 5};
//Recombine Surface {1};

//+
Curve Loop(2) = {2, 6, 11, 8};
Plane Surface(2) = {2};

//Transfinite Line {7} = 10Using Progression 1.;
//Transfinite Line {9} = 10Using Progression 1.;
//Transfinite Line {6} = 10;
//Transfinite Line {1} = 10;
//Transfinite Surface {2} = {5, 6, 7, 3};
//Recombine Surface {2};

//+
Curve Loop(3) = {3, 8, 10, 7};
Plane Surface(3) = {3};
//Transfinite Line {5} = 10Using Progression 1.;
//Transfinite Line {7} = 10Using Progression 1.;
//Transfinite Line {2} = 10;
//Transfinite Line {8} = 10;
//Transfinite Surface {3} = {3,6,4,2};
//Transfinite Surface {3} = {4, 2, 7, 6};
//Recombine Surface {3};


//+
Curve Loop(4) = {4, 7, 9, 5};
Plane Surface(4) = {4};
//Transfinite Line {5} = 10Using Progression 1.;
//Transfinite Line {88} = 10Using Progression 1.;
//Transfinite Line {3} = 10;
//Transfinite Line {11} = 10;
//Transfinite Surface {3} = {3,6,4,2};
//Transfinite Surface {4} = {2, 4, 1, 8};
//Recombine Surface {4};
/*
//-------------------------------------------------------------------------------
//Boundary tagging
//-------------------------------------------------------------------------------
// At this level, Gmsh knows everything to display the rectangular surface 1 and
// to mesh it. An optional step is needed if we want to group elementary
// geometrical entities into more meaningful groups, e.g. to define some
// mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
// material ("steel", "carbon") properties.
//
// Such groups are called "Physical Groups" in Gmsh. By default, if physical
// groups are defined, Gmsh will export in output files only mesh elements that
// belong to at least one physical group. (To force Gmsh to save all elements,
// whether they belong to physical groups or not, set `Mesh.SaveAll=1;', or
// specify `-save_all' on the command line.) Physical groups are also identified
// by tags, i.e. strictly positive integers, that should be unique per dimension
// (0D, 1D, 2D or 3D). Physical groups can also be given names.
//
// Here we define a physical curve that groups the left, bottom and right curves
// in a single group (with prescribed tag 5); and a physical surface with name
// "My surface" (with an automatic tag) containing the geometrical surface 1:
//
Physical Point("boundary",  1) = {2, 7, 3, 8};
Physical Curve("free_slip", 2) = {8, 9, 10, 11};
Physical Surface("domain") = {1};
*/