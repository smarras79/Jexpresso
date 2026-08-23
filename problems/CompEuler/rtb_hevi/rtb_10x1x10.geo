// ---------------------------------------------------------------------------
// Rising thermal bubble on the SAME domain and the SAME resolution as the 2D
// CompEuler/theta case, extruded one element in y.
//
// 10000 x 1000 x 10000 m, 10 x 1 x 10 elements. theta's hexa_TFI_10x10.msh is
// x in [-5000, 5000] by y in [0, 10000] with 10 x 10 elements, so Δx = Δz =
// 1000 m here matches it exactly and (with the matching bubble and viscosity
// in initialize.jl / user_inputs.jl) this case is theta as a 3D slab. That is
// the point: run it explicitly and it must reproduce theta, which makes any
// difference under HEVI attributable to the time integrator alone.
//
// One element in y: the bubble is a cylinder along y, so the solution is
// y-invariant and the third dimension is there only to make it a 3D case
// (HEVI is 3D only). |v| staying at round-off is a useful check that nothing
// has gone sideways. All six faces are free-slip walls (user_bc.jl), so the
// spanwise faces do not disturb that.
//
// NOTE this mesh is ISOTROPIC (Δx = Δz), which is the one geometry where a
// vertically-implicit split buys nothing -- there is no vertical acoustic
// term to remove that the horizontal one does not match. It exists to make
// HEVI COMPARABLE to a known 2D answer, not to show it off. For the
// performance claim use the deliberately anisotropic rtb_8x1x24.geo, kept
// alongside this file.
//
// Generate with:   gmsh -3 rtb_10x1x10.geo
//
// Geometry and extrude order follow rtb_8x1x24.geo; only the counts and
// extents differ. The Physical Surface names matter: BCs.jl sends every
// non-"periodic*" face to user_bc_dirichlet! with its name as the tag, and
// user_bc_neumann! in this case keys on "bottom" and "top".
// ---------------------------------------------------------------------------
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
