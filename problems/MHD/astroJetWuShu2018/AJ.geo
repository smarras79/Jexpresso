// ============================================================
// Two-dimensional domain of the magnetized astrophysical jet of
//
//   K. Wu, C.-W. Shu, "A provably positive discontinuous Galerkin
//   method for multidimensional ideal magnetohydrodynamics",
//   SIAM J. Sci. Comput. 40(5) (2018) B1302-B1329, Example 5.6,
//
//   [Xmin, Xmax] x [Ymin, Ymax] = [-0.5, 0.5] x [0, 1.5]
//
// nx x ny structured (TFI) quad elements of size (1/nx, 1.5/ny).
//
// THE ONE CONSTRAINT ON nx. The nozzle is {y = 0, |x| <= 0.05} and the
// inflow is imposed at NODES (user_bc.jl), so the lip |x| = 0.05 must be
// an element boundary: 0.05·nx must be an integer, i.e. nx must be a
// multiple of 20. nx = 40 (h = 0.025, 2 elements per nozzle half-width)
// and nx = 100 (h = 0.01, 5 per half-width) are the two that ship with
// the case; anything else and the beam edge falls inside an element,
// where the top-hat datum cannot be represented at all.
//
// Keep the elements square (ny = 1.5·nx) so that the DynSGS element
// length scale Δ = min(edge)/(N+1) is isotropic; the reference
// computations are all on square cells.
//
//   ref = 1  ->  40 x  60 elements = AJ_40x60.msh    (160 x 240 LGL points at nop 4)
//   ref = 2  -> 100 x 150 elements = AJ_100x150.msh  (400 x 600 LGL points)
//
// ref = 2 is the papers' own resolution: Wu & Shu compute the right half
// [0, 0.5] x [0, 1.5] on 200 x 600 cells, i.e. Δ = 2.5e-3, and 100 x 150
// elements at nop = 4 give the same nodal spacing over the full width.
//
// Boundary tags (the strings that reach user_bc_dirichlet!):
//   "bottom"  y = 0      the nozzle for |x| <= 0.05, outflow outside it
//   "right"   x = +0.5   outflow
//   "top"     y = 1.5    outflow
//   "left"    x = -0.5   outflow
// There is no periodic direction in this problem.
//
// Generate with:
//   gmsh -2 AJ.geo -o AJ_40x60.msh
//
// The two .msh files that ship with the case were written directly, with
// the same entity/node/element layout gmsh produces for this .geo, so the
// case runs without a gmsh installation:
//
//   tools/astro_jet_mesh.py            // regenerates both, byte-identical
//   tools/astro_jet_mesh.py 60 90      // any nx (multiple of 20) and ny
//
// This file is the reference definition and the gmsh path to the same
// meshes.
// ============================================================

Xmin = -0.5;
Xmax =  0.5;
Ymin =  0.0;
Ymax =  1.5;

ref = 1;                        // 1 -> 40x60, 2 -> 100x150
If (ref == 1)
  nx = 40;  ny = 60;
EndIf
If (ref == 2)
  nx = 100; ny = 150;
EndIf

Point(1) = {Xmin, Ymin, 0};
Point(2) = {Xmax, Ymin, 0};
Point(3) = {Xmax, Ymax, 0};
Point(4) = {Xmin, Ymax, 0};

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
Physical Curve("bottom")   = {1};
Physical Curve("right")    = {2};
Physical Curve("top")      = {3};
Physical Curve("left")     = {4};

Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;
