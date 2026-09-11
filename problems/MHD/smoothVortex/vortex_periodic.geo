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

Physical Surface("domain")   = {1};
Physical Curve("periodicx")  = {2, 4};  // left <-> right
Physical Curve("periodicy")  = {1, 3};  // bottom <-> top

Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;
