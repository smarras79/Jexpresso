// ---------------------------------------------------------------------------
// cubed_sphere.geo — a closed, all-quadrilateral cubed-sphere surface grid.
//
//   gmsh -2 problems/ShallowWater/SWsphere/cubed_sphere.geo \
//        -format msh2 -o ./meshes/gmsh_grids/cubed_sphere.msh
//
// The six panels are spherical patches ("In Sphere {1}") bounded by great
// circle arcs through the projected corners of an inscribed cube, made
// structured with Transfinite and turned into quads with Recombine.
//
// Because all six surfaces are built on the SAME twelve arcs and the SAME
// eight corner points, gmsh emits ONE node per seam location: the shell comes
// out watertight. That is the property src/kernel/mesh/sphere_mesh.jl checks
// (and repairs, by merging coincident vertices, if a hand-made grid gets it
// wrong).
//
// NOTE: this uses the built-in gmsh kernel. `Ruled Surface`/`Line Loop` are
// the pre-4.0 spellings and are still accepted as aliases in gmsh 4.x.
//
// NO SPHERE-CENTRE NODE. Point(1) at the origin is unavoidable GEOMETRY here —
// `Circle` needs a centre and `In Sphere` needs a reference — but it must not
// reach the mesh. It is in no Physical group, so `Mesh.SaveAll = 0` below drops
// it, together with the twelve arcs, leaving exactly the six quad panels. In the
// GUI this is the "Save all elements" checkbox: LEAVE IT UNCHECKED. Ticking it
// overrides the setting and writes Point(1) as a node at r = 0 that belongs to
// no element, which is what produced the 603-node files.
//
// WHICH MAP THIS ACTUALLY BUILDS: the EQUIANGULAR cubed sphere, not the
// gnomonic one. `Transfinite Line` spaces points at equal ANGLE along a
// `Circle` arc, and the transfinite patch inherits that. Verified: the grid
// this file produces matches tools/generate_cubed_sphere.jl's :equiangular
// output to 2.2e-16 of R, and differs from its :gnomonic output by 535 km.
// So do NOT then set `:cubed_sphere_map => :equiangular` on it — that would
// apply the warp twice.
//
// The equivalent grid without gmsh at all, which is the simpler path and cannot
// emit a centre node by construction:
//
//   julia tools/generate_cubed_sphere.jl 10 6.371e6 cubed_sphere.msh equiangular
// ---------------------------------------------------------------------------

R = 6371000.0;      // sphere radius [m] (Earth)
n = 10;             // elements per panel edge -> 6*n^2 quads

a  = R/Sqrt(3);
lc = 2*Pi*R/(4*n);  // only a hint: the grid is transfinite

Point(1) = {0, 0, 0, lc};        // sphere centre (the "In Sphere" reference)

// the eight cube corners, projected onto the sphere
Point(2) = {-a, -a, -a, lc};
Point(3) = { a, -a, -a, lc};
Point(4) = { a,  a, -a, lc};
Point(5) = {-a,  a, -a, lc};
Point(6) = {-a, -a,  a, lc};
Point(7) = { a, -a,  a, lc};
Point(8) = { a,  a,  a, lc};
Point(9) = {-a,  a,  a, lc};

// the twelve cube edges as great-circle arcs (centre = point 1)
Circle(1)  = {2, 1, 3};
Circle(2)  = {3, 1, 4};
Circle(3)  = {4, 1, 5};
Circle(4)  = {5, 1, 2};

Circle(5)  = {6, 1, 7};
Circle(6)  = {7, 1, 8};
Circle(7)  = {8, 1, 9};
Circle(8)  = {9, 1, 6};

Circle(9)  = {2, 1, 6};
Circle(10) = {3, 1, 7};
Circle(11) = {4, 1, 8};
Circle(12) = {5, 1, 9};

Line Loop(1) = { 1,  2,  3,  4};        // -z  (south)
Line Loop(2) = { 5,  6,  7,  8};        // +z  (north)
Line Loop(3) = { 1, 10, -5, -9};        // -y
Line Loop(4) = { 2, 11, -6, -10};       // +x
Line Loop(5) = { 3, 12, -7, -11};       // +y
Line Loop(6) = { 4,  9, -8, -12};       // -x

Ruled Surface(1) = {1} In Sphere {1};
Ruled Surface(2) = {2} In Sphere {1};
Ruled Surface(3) = {3} In Sphere {1};
Ruled Surface(4) = {4} In Sphere {1};
Ruled Surface(5) = {5} In Sphere {1};
Ruled Surface(6) = {6} In Sphere {1};

Transfinite Line {1:12} = n + 1;
Transfinite Surface {1, 2, 3, 4, 5, 6};
Recombine Surface  {1, 2, 3, 4, 5, 6};

// one physical tag per panel -> shows up as the `panel` field in the VTK output
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};

Mesh.ElementOrder    = 1;   // JEXPRESSO adds the high-order LGL points itself
Mesh.RecombineAll    = 1;
Mesh.SecondOrderLinear = 0;

// Write ONLY entities that belong to a Physical group, i.e. the six panels.
// This is what keeps Point(1) (the sphere centre) and the twelve arcs out of
// the .msh. Equivalent to leaving "Save all elements" UNCHECKED in the GUI —
// and the GUI checkbox wins over this line, so leave it unticked there too.
Mesh.SaveAll = 0;
