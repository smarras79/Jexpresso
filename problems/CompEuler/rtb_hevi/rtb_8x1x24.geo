// ---------------------------------------------------------------------------
// Rising thermal bubble, deliberately anisotropic.
//
// 1000 x 1000 x 1000 m, 8 x 1 x 24 elements. Δx = 125 m against Δz = 41.7 m,
// so the vertical spacing is three times finer than the horizontal -- which is
// the whole reason to run this case with HEVI. On an isotropic bubble mesh a
// vertically-implicit split buys nothing, because there is no vertical
// acoustic term to remove that the horizontal one does not match.
//
// One element in y: the bubble is a cylinder along y, so the solution is
// y-invariant and the third dimension is there to make it a 3D case (HEVI is
// 3D only), not to resolve anything. |v| staying at round-off is a useful
// check that nothing has gone sideways.
//
// Generate with:   gmsh -3 rtb_8x1x24.geo
//
// Geometry and extrude order follow LESICP_coarse_16x16x15.geo; only the
// counts, extents and Physical Surface names differ. The names matter: BCs.jl
// sends every non-"periodic*" face to user_bc_dirichlet! with its name as the
// tag, and user_bc_neumann! in this case keys on "bottom" and "top".
// ---------------------------------------------------------------------------
nelemx = 8;
nelemy = 1;
nelemz = 24;

xmin =    0.0;
xmax = 1000.0;
ymin =    0.0;
ymax = 1000.0;
zmin =    0.0;
zmax = 1000.0;

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
