// ===========================================================================
//  rtb3d_dsgs_10x1x10.geo -- the 2D rising thermal bubble as a 3D slab,
//                            one element thick, with ISOTROPIC cells.
//
//  10000 x 1000 x 10000 m, 10 x 1 x 10 elements  ->  dx = dy = dz = 1000 m.
//  At nop = 4 that is 41 x 5 x 41 = 8405 gridpoints and it runs on one core.
//
//  WHY THE y EXTENT IS 1000 m -- IT IS ONE ELEMENT WIDE, NOT A SLAB
//  ----------------------------------------------------------------
//  A dummy direction that is many elements WIDE but one element DEEP is
//  harmless for an inviscid or constant-viscosity run and is not harmless for
//  an LES closure. The LES filter width is chosen per element by
//  compute_element_size_driver from :les_filter_width, whose default (:max)
//  is max(dx, dy, dz), and the eddy viscosity goes as Delta SQUARED. A 5000 m
//  y extent on this mesh would hand the whole domain Delta = 1250 m instead
//  of 250 m -- 25x on nu -- from a direction the solution does not vary in.
//  `mesh.jl` warns about exactly this, and :geometric would only soften it
//  (cbrt(1000*5000*1000) = 1710 m), not fix it.
//
//  Making the dummy direction exactly one real element wide removes the
//  question instead of answering it: :max, :geometric and :min all give
//  1000 m here, so no filter-width setting can be wrong and the deck does not
//  have to carry a key whose only purpose is to work around the mesh.
//
//  NOTE ON THE FILE THIS IS DERIVED FROM. rtb_hevi/rtb_10x1x10.geo now reads
//  ymax = 5000, nelemz = 50; the rtb_10x1x10.msh committed beside it is
//  10000 x 1000 x 10000 with 10 x 1 x 10 -- the two have drifted, and it is
//  the .msh that matches this file. Regenerating that .geo would produce a
//  DIFFERENT mesh from the one rtb_hevi has been running, which is why this
//  case keeps its own .msh and its own .geo in step rather than pointing at
//  that one. (Verified: the mesh generated from this .geo is identical to
//  rtb_10x1x10.msh except in the last digit of gmsh's transfinite rounding.)
//
//  ONE ELEMENT IN y, AND WHAT IT BUYS. The bubble is a cylinder along y
//  (initialize.jl measures r in the x-z plane only), so the exact solution is
//  y-invariant and the third dimension is there to make this a 3D case, not
//  to resolve anything. |v| staying at round-off is then a free correctness
//  check on the whole 3D path -- the stress tensor, the metrics, the DSS.
//
//  Generate with either of
//
//      gmsh -3 problems/CompEuler/rtb3d_dsgs/rtb3d_dsgs_10x1x10.geo
//      julia --project=. problems/CompEuler/rtb3d_dsgs/generate_mesh.jl
//
//  The second needs no gmsh on PATH: it drives the SDK GridapGmsh ships.
//
//  Geometry and extrude order follow rtb_hevi/rtb_10x1x10.geo, because the
//  surface IDs the Physical Surface lines name are gmsh's extrusion numbering
//  and they are what the free-slip boundary tags hang off. Only the extents
//  and the element counts differ.
// ===========================================================================

nelemx = 10;
nelemy =  1;
nelemz = 10;

xmin = -5000.0;
xmax =  5000.0;
ymin =     0.0;
ymax =  1000.0;
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
