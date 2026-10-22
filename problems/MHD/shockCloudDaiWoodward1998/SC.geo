// ============================================================
// Domain of the shock / "MHD" clouds interaction of Dai & Woodward (1998, JCP 142:331, Fig. 20):
//
//   [0, 2] x [0, 1.5], open on all sides (no periodic direction)
//
// nx x ny square quad elements. DW use 512 x 384 cells; at :nop => 4 the
// 128 x 96 mesh carries 513 x 385 nodes, the same spacing.
//
//   ref = 1  ->   64 x 48 elements = SC_64x48.msh    (257 x 193 nodes at nop 4)
//   ref = 2  ->  128 x 96 elements = SC_128x96.msh   (513 x 385 nodes, DW's resolution)
//
// The surface must be named "domain" (see problems/MHD/astroJetWuShu2018/AJ.geo:
// Geom.jl filters that name out of the boundary tags).
//
// Boundary tags (the strings that reach user_bc_dirichlet!):
//   "left"    x = 0     super-fast inflow of the post-shock state
//   "right"   x = 2     open
//   "bottom"  y = 0     open
//   "top"     y = 1.5   open
//
// Generate with:
//   gmsh -2 SC.geo -setnumber ref 1 -o SC_64x48.msh
//   gmsh -2 SC.geo -setnumber ref 2 -o SC_128x96.msh
// ============================================================

Lx = 2.0;
Ly = 1.5;

DefineConstant[ ref = 2 ];
If (ref == 1)
  nx = 64;  ny = 48;
EndIf
If (ref == 2)
  nx = 128; ny = 96;
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
