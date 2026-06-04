// ============================================================
// city2d_unstr: 2D flow around a building -- UNSTRUCTURED quad mesh
//               with EXPLICIT lateral periodicity preserved.
//
// Domain:   [0, 1000] x [0, 300] m
// Building: 100 x 100 m, lower-left at (x,y) = (350, 0)
//
// Single Plane Surface with the building cut out of the bottom.
// All-quad output via Algorithm 8 + Blossom + SubdivisionAlgorithm = 1.
//
// LATERAL PERIODICITY (two safeguards, both active):
//
//   (a) Transfinite Line with MATCHED counts on the two periodic walls
//       (ny_low on the lower segment, ny_high on the upper).  Equal
//       segment counts on equal-length walls force identical node
//       y-coordinates on the left and right — i.e. geometric
//       periodicity of the boundary discretisation.
//
//   (b) Periodic Line{...} Translate{...} -- gmsh's explicit
//       master/slave directive.  Every node on lines 6 and 7 is
//       slaved to its translate on lines 12 and 11.  After this
//       directive gmsh writes the right-wall nodes as exact
//       (x + xmax, y) copies of the left-wall nodes.
//
// Physical groups:
//   "bottom"    : full bottom of the domain including the building
//   "top"       : top wall
//   "periodicx" : both vertical walls (left and right)
//   "domain"    : interior surface
//
// Generate with:
//   gmsh -2 city2d_unstr.geo -o ../../../meshes/gmsh_grids/city2d_unstr.msh
// ============================================================

xmin = 0.0;
xmax = 1000.0;
ymin = 0.0;
ymax = 300.0;

// Building extent
bx0 = 350.0;
bx1 = 450.0;
by1 = 100.0;

// Element-count knobs on each perimeter segment.
// All counts are even (Blossom-friendly).  Smallest element is ~10 m
// before subdivision, ~5 m after — right at the user's floor.
nx_A    = 18;   // bottom strip 0   ... 350    (Δx ≈ 19.4 m)
nx_C    = 6;    // building top   350 ... 450    (Δx ≈ 16.7 m)
nx_E    = 28;   // bottom strip 450 ... 1000   (Δx ≈ 19.6 m)
ny_low  = 10;   // y in [0,   100]              (Δy = 10 m)
ny_high = 24;   // y in [100, 300]              (Δy ≈ 8.3 m)

// -------- Points --------
Point(1)  = {xmin, ymin, 0};   // bottom-left
Point(2)  = {bx0,  ymin, 0};   // building left foot
Point(3)  = {bx0,  by1,  0};   // building top-left
Point(4)  = {bx1,  by1,  0};   // building top-right
Point(5)  = {bx1,  ymin, 0};   // building right foot
Point(6)  = {xmax, ymin, 0};   // bottom-right
Point(7)  = {xmax, by1,  0};   // right-wall split
Point(8)  = {xmax, ymax, 0};   // top-right
Point(9)  = {xmin, ymax, 0};   // top-left
Point(10) = {xmin, by1,  0};   // left-wall split

// -------- Lines --------
Line(1)  = {1,  2};    // bottom, A-strip
Line(2)  = {2,  3};    // building, left wall
Line(3)  = {3,  4};    // building, top
Line(4)  = {4,  5};    // building, right wall
Line(5)  = {5,  6};    // bottom, D-strip
Line(6)  = {6,  7};    // right wall, lower  (periodicx)
Line(7)  = {7,  8};    // right wall, upper  (periodicx)
Line(8)  = {8,  9};    // top
Line(9)  = {9, 10};    // left wall, upper   (periodicx)
Line(10) = {10, 1};    // left wall, lower   (periodicx)

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Plane Surface(1) = {1};
Recombine Surface{1};

// -------- (a) Transfinite Line: matched discretisation on opposite walls --------
Transfinite Line{1}  = nx_A    + 1;   // bottom A
Transfinite Line{2}  = ny_low  + 1;   // building left wall
Transfinite Line{3}  = nx_C    + 1;   // building top
Transfinite Line{4}  = ny_low  + 1;   // building right wall
Transfinite Line{5}  = nx_E    + 1;   // bottom D
Transfinite Line{6}  = ny_low  + 1;   // right wall, lower  (matches Line 10)
Transfinite Line{7}  = ny_high + 1;   // right wall, upper  (matches Line 9)
Transfinite Line{8}  = (nx_A + nx_C + nx_E) + 1;  // top
Transfinite Line{9}  = ny_high + 1;   // left  wall, upper
Transfinite Line{10} = ny_low  + 1;   // left  wall, lower

// -------- (b) gmsh-level periodicity directive --------
// Slave the right-wall lines to the (reversed) left-wall lines so the
// right-wall boundary nodes are exact (+xmax, 0) translates of the
// left-wall nodes.  Reversal is needed because Line 6 goes upward
// (P6->P7) while Line 10 also goes upward (P10->P1)... no — Line 10
// goes from P10 (top of lower-left split) down to P1, so the matching
// direction is opposite.  Likewise Line 7 (upward) matches the
// reversal of Line 9 (downward).
Periodic Line{6} = {-10} Translate{xmax - xmin, 0, 0};
Periodic Line{7} = {-9}  Translate{xmax - xmin, 0, 0};

// -------- Physical groups --------
Physical Surface("domain")    = {1};
Physical Curve("bottom")      = {1, 2, 3, 4, 5};
Physical Curve("top")         = {8};
Physical Curve("periodicx")   = {6, 7, 9, 10};

// -------- Meshing knobs --------
Mesh.ElementOrder           = 1;
Mesh.Algorithm              = 8;   // Frontal-Delaunay for Quads (quad-friendly triangulation)
Mesh.RecombinationAlgorithm = 1;   // Blossom (compatible with Transfinite Line and Periodic Line)
Mesh.SubdivisionAlgorithm   = 1;   // Post-split any residual triangle into 3 quads
                                   // → 100% quads in the .msh.  Conforming because
                                   // edge midpoints are shared between neighbours.
                                   // Safe here because explicit Periodic Line is set
                                   // BEFORE subdivision, so the periodic remap covers
                                   // the new midpoint nodes too.
