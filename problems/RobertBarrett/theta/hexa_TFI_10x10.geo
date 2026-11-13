nelemx = 10;
nelemy = 10;
nelemz = 1;

xmin = -10000;
xmax = 10000;
ymin = 0;
ymax = 20000;
ymid = (ymax+ymin)/2;
gridsize = (xmax-xmin) / nelemx;

Point(1) = {xmin, ymin, 0};
Point(2) = {0, ymin, 0};
Point(3) = {xmax, ymin, 0};
Point(4) = {xmax, ymid, 0};
Point(5) = {xmax, ymax, 0};
Point(6) = {0, ymax, 0};
Point(7) = {xmin, ymax, 0};
Point(8) = {xmin, ymid, 0};
Point(9) = {0, ymid, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {2, 9};
Line(10) = {9, 6};
Line(11) = {4, 9};
Line(12) = {9, 8};

npx = nelemx + 1;
npy = nelemy + 1;

//Horizontal sides
Transfinite Line {1, 12} = npx Using Progression 1;
Transfinite Line {12, 6} = npx Using Progression 1;
Transfinite Line {2, 11} = npx Using Progression 1;
Transfinite Line {11, 5} = npx Using Progression 1;

//Vertical sides
Transfinite Line {8, 9} = npy Using Progression 1.0;
Transfinite Line {7, 10} = npy Using Progression 1.0;
Transfinite Line {9, 3} = npy Using Progression 1.0;
Transfinite Line {10, 4} = npy Using Progression 1.0;

Line Loop(11) = {8, 1, 9, 12};
Plane Surface(12) = {11};

Transfinite Surface {12};
Recombine Surface {12};

Line Loop(21) = {7, -12, 10, 6};
Plane Surface(22) = {21};

Transfinite Surface {22};
Recombine Surface {22};

Line Loop(31) = {-9, 2, 3, 11};
Plane Surface(32) = {31};

Transfinite Surface {32};
Recombine Surface {32};

Line Loop(41) = {-10, -11, 4, 5};
Plane Surface(42) = {41};

Transfinite Surface {42};
Recombine Surface {42};

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
Physical Point("boundary",   1) = {1, 2, 3, 4, 5, 6, 7, 8};
//Physical Curve("solid",  2) = {3};
Physical Curve("solid",  3) = {3, 4, 5, 6, 7, 8};
Physical Curve("bottom",  4) = {1, 2};
Physical Surface("domain") = {12, 22, 32, 42};
