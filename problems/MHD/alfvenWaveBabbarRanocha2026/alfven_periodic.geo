// ============================================================
// Doubly periodic box [0, 1/cos(alpha)] x [0, 1/sin(alpha)] quad grid for
// the smooth Alfven wave of Babbar & Ranocha, J. Comput. Phys. 563 (2026)
// 115060, Sec. 5.5.1 (Toth, J. Comput. Phys. 161 (2000) 605). With the
// propagation angle alpha = pi/4 of the paper the box is [0, sqrt(2)]^2,
// one wavelength of the wave along each side.
//
// The physical-curve names "periodicx"/"periodicy" are Jexpresso's
// periodicity tags: the mesh reader builds the periodic point
// pairing itself, so no gmsh-level "Periodic Curve" constraints
// are needed.
//
// nx is a command-line constant, so one .geo makes every
// resolution of a convergence sweep:
//
//   gmsh -2 -setnumber nx 16 alfven_periodic.geo -o alfven_16x16.msh
//
// or, for the whole set, tools/alfven_wave_mesh.sh
// ============================================================

DefineConstant[ nx    = {16, Name "nx"} ];
DefineConstant[ ny    = {nx, Name "ny"} ];
DefineConstant[ alpha = {Pi/4, Name "alpha"} ];

Lx = 1/Cos(alpha);
Ly = 1/Sin(alpha);

Point(1) = { 0,  0, 0};
Point(2) = {Lx,  0, 0};
Point(3) = {Lx, Ly, 0};
Point(4) = { 0, Ly, 0};

// OPPOSITE EDGES RUN THE SAME WAY (top like bottom, left to right; left like
// right, bottom to top). gmsh places the nodes of an edge with a smooth error
// of up to 2e-12 that follows the edge's direction; on edges running opposite
// ways the errors are mirrored, and the nodes Jexpresso pairs across a seam
// were up to 4e-12 apart. The faces on the two sides of the seam then differ
// in length, the flux through them does not cancel, and mass and energy
// drifted by 1e-13 with the period of the wave. Run the same way, the paired
// nodes coincide to the last bit and both are conserved to round-off.
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
