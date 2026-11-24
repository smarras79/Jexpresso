// ============================================================
// SoliWaveIsland: 2D non-linear shallow water -- closed basin
//
// Domain: [0, 25] x [-15, 15] m  (25 m by 30 m)
// 25 x 30 structured quadrilaterals -> Δx = Δy = 1 m at N=4.
//
// All four boundaries are reflective (free-slip) walls.
//
// Generate with:
//   gmsh -2 SoliWaveIsland.geo -o ../../../meshes/gmsh_grids/SoliWaveIsland.msh
// ============================================================

xmin =   0.0;
xmax =  25.0;
ymin = -15.0;
ymax =  15.0;
nx   =  25;
ny   =  30;

Point(1) = {xmin, ymin, 0};
Point(2) = {xmax, ymin, 0};
Point(3) = {xmax, ymax, 0};
Point(4) = {xmin, ymax, 0};

Line(1) = {1, 2};  // bottom
Line(2) = {2, 3};  // right
Line(3) = {3, 4};  // top
Line(4) = {4, 1};  // left

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve{1, 3} = nx + 1;
Transfinite Curve{2, 4} = ny + 1;
Transfinite Surface{1};
Recombine Surface{1};

Physical Surface("domain") = {1};
Physical Curve("bottom") = {1};
Physical Curve("right")  = {2};
Physical Curve("top")    = {3};
Physical Curve("left")   = {4};

Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;
