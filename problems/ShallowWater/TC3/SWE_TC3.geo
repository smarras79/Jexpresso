// ============================================================
// TC3: Planetary Rossby Wave -- Non-periodic domain
// Domain: [0, 1e6] x [0, 1e6] m  (1000 km x 1000 km)
// 16x16 structured quad elements
//
// Generate with:  gmsh -2 SWE_TC3.geo
// ============================================================

Lx = 1.0e6;
Ly = 1.0e6;
nx = 16;
ny = 16;

Point(1) = {0,  0,  0};
Point(2) = {Lx, 0,  0};
Point(3) = {Lx, Ly, 0};
Point(4) = {0,  Ly, 0};

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
