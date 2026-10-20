// ============================================================
// Doubly periodic square [-L/2, L/2]^2 quad grid for the smooth
// MHD vortex (Dao & Nazarov, J. Sci. Comput. 92:77 (2022), §5.1.1;
// Balsara, ApJS 151:149 (2004)). L = 10 is the box the vortex is
// written for.
//
// The physical-curve names "periodicx"/"periodicy" are Jexpresso's
// periodicity tags: the mesh reader builds the periodic point
// pairing itself, so no gmsh-level "Periodic Curve" constraints
// are needed.
//
// nx is a command-line constant, so one .geo makes every
// resolution of a convergence sweep:
//
//   gmsh -2 -setnumber nx 16 vortex_periodic.geo -o vortex_16x16.msh
//
// or, for the whole set, tools/smooth_vortex_mesh.sh
// ============================================================

DefineConstant[ nx = {16, Name "nx"} ];
DefineConstant[ ny = {nx, Name "ny"} ];
DefineConstant[ L  = {10.0, Name "L"} ];

h = L/2;

Point(1) = {-h, -h, 0};
Point(2) = { h, -h, 0};
Point(3) = { h,  h, 0};
Point(4) = {-h,  h, 0};

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
