// ============================================================
// Doubly periodic [-1,1]^2 quad grid for the magnetized
// Kelvin-Helmholtz instability (Chan et al. 2022, Sec. 3.2.1).
//
// 32x32 structured (TFI) quad elements. The physical-curve names
// "periodicx"/"periodicy" are Jexpresso's periodicity tags: the
// mesh reader builds the periodic point pairing itself, so no
// gmsh-level "Periodic Curve" constraints are needed.
//
// Generate with:  gmsh -2 KHI_32x32_periodic.geo -o hexa_TFI_32x32_unitsquare.msh
// ============================================================

L  = 1.0;
nx = 32;
ny = 32;

Point(1) = {-L, -L, 0};
Point(2) = { L, -L, 0};
Point(3) = { L,  L, 0};
Point(4) = {-L,  L, 0};

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
