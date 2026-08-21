// ---------------------------------------------------------------------------
// COARSE mesh for the p4est path (:linitial_refine => true).
//
// 16 x 16 x 15 elements over 10240 x 10240 x 5000 m.
// With :init_refine_lvl => 1 p4est refines this to 32 x 32 x 30, which is the
// grid the simulation actually runs on.
//
// The point of going through p4est is that ONLY this coarse mesh is read
// serially and broadcast to every rank; the refinement happens in parallel.
// 3,840 elements instead of 30,720 -- and for a 128^2 target, 3,840 instead
// of 1,966,080, which is the size that OOM-killed nodes.
//
// Generate with:   gmsh -3 LESICP_coarse_16x16x15.geo
//
// Geometry, extrude order and Physical Surface ids are identical to
// LESICP2-coarse/LESICP.geo -- only the element counts and extents differ --
// so "MOST", "top_wall", "periodicx" and "periodicy" keep their meaning.
// ---------------------------------------------------------------------------
nelemx = 16;
nelemy = 16;
nelemz = 15;

xmin =      0.0;
xmax =  10240.0;
ymin =      0.0;
ymax =  10240.0;
zmin =      0.0;
zmax =   5000.0;

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
npy = nelemy + 1;
npz = nelemz + 1;

//Horizontal sides
Transfinite Line {1, 3} = npx;
//Vertical sides
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

    Physical Surface("periodicy") = {12,34};
    Physical Volume("internal") = {1};
    Physical Surface("top_wall") = {33};
    Physical Surface("MOST") = {25};
    Physical Surface("periodicx") = {21,29};
