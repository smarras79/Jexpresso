// ===========================================================================
//  rtb2d_dsgs_10x10.geo -- the 2D half of the DynSGS ladder.
//
//  10000 x 10000 m, 10 x 10 quad elements. At :nop => 4 that is 41 x 41 =
//  1681 gridpoints and it runs in seconds.
//
//  WHY THIS MESH EXISTS WHEN CompEuler/theta_dsgs ALREADY RUNS THE SAME CASE.
//  theta_dsgs points at ./meshes/gmsh_grids/hexa_TFI_10x10.msh, and `meshes`
//  is a SYMLINK to the separate smarras79/JexpressoMeshes repository -- it is
//  in .gitignore and does not exist in a fresh clone or on a CI runner. So the
//  one 2D DynSGS case in the repository cannot be run by someone who has only
//  this repository, which is exactly the wrong property for the case you reach
//  for when the 3D one misbehaves. This one is committed next to its deck.
//
//  IT IS THE SAME GEOMETRY AS rtb3d_dsgs, ON PURPOSE. Same extents, same
//  element count in the x-z plane, same nop, same bubble (initialize.jl:
//  centre (0, 2500) m, radius 2000 m, 2 K, linear taper), same dt, same tend.
//  A 3D slab of a 2D problem IS the 2D problem, so the pair is a ladder: if
//  2D is right and 3D is not, the difference is in the 3D path and nowhere
//  else. The one thing that is NOT the same is the closure's theta
//  coefficient -- the 2D DSGS path applies Nazarov & Hoffman's Pr/(gamma-1)*mu
//  and the 3D path applies mu/Pr_t -- and rtb2d_dsgs/user_inputs.jl says how
//  to line the two up when that matters.
//
//  Generate with:
//    gmsh -2 problems/CompEuler/rtb2d_dsgs/rtb2d_dsgs_10x10.geo \
//            -o problems/CompEuler/rtb2d_dsgs/rtb2d_dsgs_10x10.msh
//
//  The boundary names are free_slip on all four sides. user_bc_dirichlet! in
//  this case ignores the tag and removes the normal momentum component
//  whatever it is called, so the names are documentation; what matters is that
//  none of them starts with "periodic", which is the prefix BCs.jl reserves.
// ===========================================================================

xmin = -5000.0;
xmax =  5000.0;
ymin =     0.0;
ymax = 10000.0;

nx = 10;
ny = 10;

Point(1) = {xmin, ymin, 0};
Point(2) = {xmax, ymin, 0};
Point(3) = {xmax, ymax, 0};
Point(4) = {xmin, ymax, 0};

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
