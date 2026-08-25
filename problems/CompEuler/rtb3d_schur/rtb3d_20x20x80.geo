// ===========================================================================
//  rtb3d_20x20x80.geo -- a 10 km cube for the Schur A/B.
//
//  THE ELEMENT ASPECT RATIO IS THE POINT, not an accident of resolution.
//  The column preconditioner is exact for the vertical acoustic operator and
//  does nothing for the horizontal one, so what it leaves for the Krylov
//  iteration -- and therefore how much the Schur reduction can save -- scales
//  with the grid's ACOUSTIC ANISOTROPY h_x/h_z. Measured on the mock at
//  production stiffness (test/hevi/test_schur_precond.jl):
//
//       h_x/h_z    H iterations   full   ratio
//         1:1           11          23   1.66x
//         4:1            3           8   2.67x    <-- this mesh
//        16:1            2           3   1.50x
//
//  20 x 20 x 80 elements over 10 km cubed gives
//
//       h_x = h_y = 500 m,  h_z = 125 m   ->  4:1
//
//  and at nop = 4 the smallest LGL gaps are 86.4 m horizontally and 21.6 m
//  vertically, which is the ratio that actually matters (the gaps scale with
//  the element, so the 4:1 carries through).
//
//  SIZE. 81 x 81 x 321 = 2,106,081 gridpoints. tools/pick_nranks.jl recommends
//  25 ranks (a 5 x 5 rank grid over the 20 x 20 element columns, 16 columns and
//  84k points each); 16 ranks is marginally more efficient per rank. Do not use
//  a rank count that tool does not list -- with :lxy_partition the usable
//  parallelism is bounded by nelemx*nelemy and the rest leave ranks empty.
//
//  Regenerate with either of
//
//      gmsh -3 problems/CompEuler/rtb3d_schur/rtb3d_20x20x80.geo
//      julia --project=. problems/CompEuler/rtb3d_schur/generate_mesh.jl
//
//  The second needs no gmsh on PATH: it drives the SDK that GridapGmsh already
//  ships, which is the same mesher.
//
//  Structure copied from CompEuler/rtb_imex/rtb_10x1x10.geo -- a transfinite
//  surface in x-z, extruded in y -- because the surface IDs the Physical
//  Surface lines below name are gmsh's extrusion numbering, and they are what
//  the free-slip boundary tags hang off.
// ===========================================================================

nelemx = 20;
nelemy = 20;
nelemz = 80;

xmin = -5000.0;
xmax =  5000.0;
ymin = -5000.0;
ymax =  5000.0;
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

    Physical Surface("front")  = {12};
    Physical Surface("back")   = {34};
    Physical Volume("internal") = {1};
    Physical Surface("top")    = {33};
    Physical Surface("bottom") = {25};
    Physical Surface("left")   = {21};
    Physical Surface("right")  = {29};
