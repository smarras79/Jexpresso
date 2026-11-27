// Mesh for problems/AdvDiff/2d_robin: the rectangle [0, 2] x [0, 1], periodic in x, with the two
// horizontal walls tagged "bottom" and "top" (the Robin boundaries).
//
//   gmsh robin_2d.geo -2 -format msh4 -o robin_2d_10x5.msh
//   gmsh robin_2d.geo -2 -format msh4 -setnumber nelemx 20 -setnumber nelemy 10 -o robin_2d_20x10.msh
//
// nelemx/nelemy set the resolution (keep them in a 2:1 ratio for square elements).

DefineConstant[ nelemx = 10, nelemy = 5 ];

xmin = 0.0;
xmax = 2.0;
ymin = 0.0;
ymax = 1.0;

Point(1) = {xmin, ymin, 0};
Point(2) = {xmax, ymin, 0};
Point(3) = {xmax, ymax, 0};
Point(4) = {xmin, ymax, 0};

Line(1) = {1, 2};   // bottom
Line(2) = {2, 3};   // right
Line(3) = {3, 4};   // top
Line(4) = {4, 1};   // left

Transfinite Line {1, 3} = nelemx + 1;
Transfinite Line {2, 4} = nelemy + 1;

Line Loop(11) = {1, 2, 3, 4};
Plane Surface(12) = {11};
Transfinite Surface {12};
Recombine Surface {12};

// Tags read by Jexpresso. "periodicx" pairs the left and right sides; "bottom" and "top" are
// ordinary tags, so user_bc.jl decides what happens there (here: a Robin flux).
Physical Point("boundary", 1)  = {1, 2, 3, 4};
Physical Curve("periodicx", 2) = {2, 4};
Physical Curve("bottom", 3)    = {1};
Physical Curve("top", 4)       = {3};
Physical Surface("domain", 5)  = {12};
