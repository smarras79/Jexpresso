// ============================================================
// city2d: 2D flow around a single building
//
// Domain:   [0, 1000] x [0, 1000] m
// Building: 100 x 100 m, positioned with lower-left at (x,y) = (350, 0)
//
// Physical groups:
//   "bottom"    : full bottom of the domain including the building profile
//   "top"       : top wall
//   "periodicx" : both vertical walls (left and right)
//   "domain"    : interior surface
//
// Mesh: unstructured, QUADS only.
//
// Generate with:
//   gmsh -2 city2d.geo -o ../../../meshes/gmsh_grids/city2d.msh
// ============================================================

xmin = 0.0;
xmax = 1000.0;
ymin = 0.0;
ymax = 1000.0;

// Building extent
bx0 = 350.0;        // left foot
bx1 = 450.0;        // right foot  (length = 100)
by0 = 0.0;
by1 = 100.0;        // top         (height = 100)

// Characteristic mesh sizes
lc_far = 40.0;      // far field
lc_bld = 10.0;      // near the building

Point(1) = {xmin, ymin, 0, lc_far};
Point(2) = {bx0,  ymin, 0, lc_bld};
Point(3) = {bx0,  by1,  0, lc_bld};
Point(4) = {bx1,  by1,  0, lc_bld};
Point(5) = {bx1,  ymin, 0, lc_bld};
Point(6) = {xmax, ymin, 0, lc_far};
Point(7) = {xmax, ymax, 0, lc_far};
Point(8) = {xmin, ymax, 0, lc_far};

Line(1) = {1, 2};   // bottom, left strip
Line(2) = {2, 3};   // building, left wall
Line(3) = {3, 4};   // building, top
Line(4) = {4, 5};   // building, right wall
Line(5) = {5, 6};   // bottom, right strip
Line(6) = {6, 7};   // right vertical wall (periodic)
Line(7) = {7, 8};   // top
Line(8) = {8, 1};   // left vertical wall  (periodic)

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// Force a periodic node distribution between the two vertical walls
Periodic Curve{6} = {-8} Translate{xmax - xmin, 0, 0};

// Physical groups (these names become the boundary tags in Jexpresso)
Physical Surface("domain")    = {1};
Physical Curve("bottom")      = {1, 2, 3, 4, 5};
Physical Curve("top")         = {7};
Physical Curve("periodicx")   = {6, 8};

// ---- Meshing options: unstructured, all-quad ----
Recombine Surface{1};

Mesh.Algorithm              = 8;   // Frontal-Delaunay for Quads
Mesh.RecombinationAlgorithm = 3;   // Blossom full-quad
Mesh.RecombineAll           = 1;
Mesh.SubdivisionAlgorithm   = 1;   // Guarantees an all-quad mesh
Mesh.ElementOrder           = 1;
