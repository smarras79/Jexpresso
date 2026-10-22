// ============================================================
// Domain of the Mach 15 shock / round cloud test of Dai & Woodward (1998, JCP 142:331, p. 365, Fig. 23):
//
//   [0, 1.5] x [0, 1], open on all sides (no periodic direction)
//
// nx x ny square quad elements. DW use 384 x 256 cells; at :nop => 4 the
// 96 x 64 mesh carries 385 x 257 nodes, the same spacing.
//
//   ref = 1  ->  48 x 32 elements = SC15_48x32.msh   (193 x 129 nodes at nop 4)
//   ref = 2  ->  96 x 64 elements = SC15_96x64.msh   (385 x 257 nodes, DW's resolution)
//
// The surface must be named "domain" (see problems/MHD/astroJetWuShu2018/AJ.geo).
//
// Boundary tags (the strings that reach user_bc_dirichlet!):
//   "right"   x = 1.5   super-fast inflow of the pre-shock state (cloud ahead of the shock)
//   "left"    x = 0     open
//   "bottom"  y = 0     open
//   "top"     y = 1     open
//
// Generate with:
//   gmsh -2 SC15.geo -setnumber ref 1 -o SC15_48x32.msh
//   gmsh -2 SC15.geo -setnumber ref 2 -o SC15_96x64.msh
// ============================================================

Lx = 1.5;
Ly = 1.0;

DefineConstant[ ref = 2 ];
If (ref == 1)
  nx = 48;  ny = 32;
EndIf
If (ref == 2)
  nx = 96;  ny = 64;
EndIf

Point(1) = {0,  0,  0};
Point(2) = {Lx, 0,  0};
Point(3) = {Lx, Ly, 0};
Point(4) = {0,  Ly, 0};

// Opposite edges run the same way, as in the periodic meshes of the other MHD cases.
Line(1) = {1, 2};  // bottom
Line(2) = {2, 3};  // right
Line(3) = {4, 3};  // top
Line(4) = {1, 4};  // left

Curve Loop(1) = {1, 2, -3, -4};
Plane Surface(1) = {1};

Transfinite Curve{1, 3} = nx + 1;
Transfinite Curve{2, 4} = ny + 1;
Transfinite Surface{1};
Recombine Surface{1};

Physical Surface("domain") = {1};
Physical Curve("bottom")   = {1};
Physical Curve("right")    = {2};
Physical Curve("top")      = {3};
Physical Curve("left")     = {4};

Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;
