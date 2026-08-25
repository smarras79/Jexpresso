// ===========================================================================
//  rtb2d_20x1x80.geo -- the 2D analogue of rtb3d_20x20x80.geo, as a SLAB.
//
//  WHY A SLAB AND NOT AN nsd == 2 MESH
//  -----------------------------------
//  HEVI and IMEX3D are three-dimensional by construction, not by accident, and
//  refuse a 2D mesh outright (params_setup.jl: "HEVI is 3D only" /
//  "IMEX3D is 3D only"). The refusal is protective rather than lazy:
//
//    * build_column_topology reads coords[3,:] and builds the horizontal
//      column catalogue as the cross product of distinct x and y. A 2D mesh
//      allocates coords as (2, npoin), so that is a BoundsError -- and in 2D
//      the vertical is y, not z, so even a padded array would identify the
//      wrong axis.
//    * every element kernel is a triple i,j,k loop over ngl^3 nodes of
//      connijk[iel,i,j,k], which a 2D mesh shapes (nelem, ngl, ngl, 1).
//    * the operator is contracted with dzeta/d{x,y,z}. A 2D mesh ALLOCATES
//      those and never fills them, so without the guard the acoustic operator
//      would assemble to zero -- a silent wrong answer, not a crash. That is
//      the one worth knowing about.
//
//  So the 2D analogue is the idiom the repository already uses for rtb_hevi
//  and rtb_imex: a three-dimensional mesh ONE ELEMENT THICK in y, with
//  free-slip walls front and back. The bubble in initialize.jl is a CYLINDER
//  rather than a sphere -- no y dependence -- and free slip sets v = 0 on both
//  y faces, so v stays identically zero and the solution is exactly the 2D one.
//  Everything downstream (the column topology, the Schur reduction, the
//  preconditioner, the profile) runs the same code as the 3D case.
//
//  WHAT IS HELD FIXED FROM THE 3D CASE, so the two are comparable
//  --------------------------------------------------------------
//  Same 10 km height, same 10 km width, same 20 x 80 element count in the x-z
//  plane, therefore the SAME mesh spacing and the SAME acoustic anisotropy:
//
//       h_x = 500 m,  h_z = 125 m   ->  4:1
//
//  which is the ratio the Schur reduction is most worth having at (the table
//  in rtb3d_20x20x80.geo). Same Delta t, same CFL_h ~ 1.05, same bubble. The
//  only thing that changes is that 20 elements of y collapse to 1.
//
//  y is one element 500 m wide -- ISOTROPIC WITH h_x, deliberately. A very
//  thin slab (say 50 m) would make h_x/h_y = 10 and hand the horizontal
//  operator a stiffness the 3D case does not have, which is exactly the thing
//  this benchmark is supposed to hold fixed.
//
//  SIZE. 81 x 5 x 321 = 130,005 gridpoints, 1/16 of the 3D case's 2,106,081.
//  That is a desktop run: minutes on a handful of cores, not a queue slot.
//
//  Regenerate with either of
//
//      gmsh -3 problems/CompEuler/rtb2d_schur/rtb2d_20x1x80.geo
//      julia --project=. problems/CompEuler/rtb2d_schur/generate_mesh.jl
//
//  The second needs no gmsh on PATH.
// ===========================================================================

nelemx = 20;
nelemy =  1;
nelemz = 80;

xmin = -5000.0;
xmax =  5000.0;
ymin =  -250.0;
ymax =   250.0;
zmin =     0.0;
zmax = 10000.0;

gridsize = (xmax-xmin) / nelemx;

Point(1) = {xmin, ymin, zmin, gridsize};
Point(2) = {xmax, ymin, zmin, gridsize};
Point(3) = {xmax, ymin, zmax, gridsize};
Point(4) = {xmin, ymin, zmax, gridsize};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

npx = nelemx + 1;
npz = nelemz + 1;

// Horizontal sides
Transfinite Line {1, 3} = npx;
// Vertical sides
Transfinite Line {4, -2} = npz Using Progression 1.0;

Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {11};

Transfinite Surface {12};
Recombine Surface {12};

surfaceVector = Extrude {0,(ymax-ymin),0} {
  Surface{12};
  Layers{nelemy};
  Recombine;
};

// The surface IDs below are gmsh's EXTRUSION numbering, not anything chosen
// here. They are copied from rtb3d_20x20x80.geo, which has the same
// surface-in-x-z-extruded-in-y structure, and they are what the free-slip
// boundary tags hang off.
    Physical Surface("front")  = {12};
    Physical Surface("back")   = {34};
    Physical Volume("internal") = {1};
    Physical Surface("top")    = {33};
    Physical Surface("bottom") = {25};
    Physical Surface("left")   = {21};
    Physical Surface("right")  = {29};
