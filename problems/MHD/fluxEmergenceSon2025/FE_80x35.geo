// ============================================================
// Two-dimensional flux-emergence domain of
//
//   D. Son, Y. Jang, T. Magara, "A Comparative Analysis of
//   High-resolution Shock-capturing Schemes for Two-dimensional
//   Magnetohydrodynamic Simulation of Flux Emergence in the
//   Solar Atmosphere", ApJS 277:46 (2025), Sec. 2.1.
//
//   [0, Xmax] x [0, Zmax] = [0, 80 H0] x [0, 35 H0]
//
// (Jexpresso's y is the paper's vertical z.)
//
// nx x ny structured (TFI) quad elements of size (Xmax/nx, Zmax/ny).
// With nx = 80, ny = 35 every element is 1 H0 x 1 H0 and the case's
// :nop => 4 gives 320 x 140 unique LGL points, i.e. a mean nodal
// spacing of 0.25 H0 (smallest LGL gap 0.17 H0). This is the
// coarsest resolution that still puts ~3 nodes across the
// w = 0.5-0.6 H0 tanh transitions of the flux sheet and of the
// chromosphere-corona interface; the paper's coarsest grid is
// 300^2 cells (dx = 0.27 H0, dz = 0.12 H0). Refine with ny = 70
// (dz_elem = 0.5 H0) to match the paper's vertical spacing.
//
// Boundary tags:
//   "periodicx" left/right  -> Jexpresso periodicity tags (the mesh
//                              reader builds the point pairing itself)
//   "bottom"                -> symmetric (free-slip) wall, see user_bc.jl
//   "top"                   -> free-slip + absorbing layer (user_source.jl)
//
// Generate with:
//   gmsh -2 FE_80x35.geo -o FE_80x35.msh
// ============================================================

Xmax = 80.0;
Zmax = 35.0;
nx   = 80;
ny   = 35;

Point(1) = {0,    0,    0};
Point(2) = {Xmax, 0,    0};
Point(3) = {Xmax, Zmax, 0};
Point(4) = {0,    Zmax, 0};

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

Physical Surface("domain")  = {1};
Physical Curve("periodicx") = {2, 4};  // left <-> right
Physical Curve("bottom")    = {1};
Physical Curve("top")       = {3};

Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;
