// ============================================================
// Doubly periodic UNIT SQUARE [0,1]^2 quad grid for the
// Orszag-Tang vortex (Bormanis, Leon & Scheinker,
// Phys. Plasmas 31, 012101 (2024), Sec. I A / III).
//
// 32x32 structured (TFI) quad elements. Combined with the
// case's :nop => 4 this gives 32*4 = 128 unique points per
// direction, i.e. exactly the 128 x 128 resolution of the
// reference simulation.
//
// The physical-curve names "periodicx"/"periodicy" are
// Jexpresso's periodicity tags: the mesh reader builds the
// periodic point pairing itself, so no gmsh-level
// "Periodic Curve" constraints are needed.
//
// Generate with:
//   gmsh -2 OT_32x32_periodic.geo -o OT_32x32_periodic.msh
// ============================================================

L  = 1.0;
nx = 32;
ny = 32;

Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, L, 0};
Point(4) = {0, L, 0};

// OPPOSITE EDGES RUN THE SAME WAY (top like bottom, left to right; left like
// right, bottom to top). gmsh places the nodes of an edge with a small smooth
// error (~1e-12 of the edge length) that follows the edge's direction; on
// edges running opposite ways the errors are mirrored, the nodes Jexpresso
// pairs across a periodic seam no longer coincide, and the flux through the
// two sides of the seam does not cancel: mass and energy then leak at ~1e-13
// wherever the flow varies along the seam. Run the same way, the paired nodes
// coincide to the last bit (see problems/MHD/alfvenWaveBabbarRanocha2026).
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

Physical Surface("domain")   = {1};
Physical Curve("periodicx")  = {2, 4};  // left <-> right
Physical Curve("periodicy")  = {1, 3};  // bottom <-> top

Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;
