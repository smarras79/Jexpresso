// ============================================================
// city2d: 2D flow around a single building
//
// Domain:   [0, 1000] x [0, 1000] m
// Building: 100 x 100 m, lower-left at (x,y) = (350, 0)
//
// Single unstructured surface, all-quad via Recombine, with
// Transfinite Line on the perimeter to (a) control element size,
// and (b) force matching node counts on the two periodic vertical
// walls — required by Jexpresso's "periodicx" pairing.
//
// Physical groups:
//   "bottom"    : full bottom of the domain including the building profile
//   "top"       : top wall
//   "periodicx" : both vertical walls (left and right)
//   "domain"    : interior surface
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
by1 = 100.0;        // top         (height = 100)

// Element-count "knobs" on each perimeter segment.
// Left/right walls are split at y = by1 so that the periodic counts
// (lines 6 ↔ 10 lower, lines 7 ↔ 9 upper) match by construction.
// Coarse mesh: smallest element ~17 m near the building, ~37 m in the
// far field (well above the requested 5 m floor).
//
// IMPORTANT: with Mesh.RecombinationAlgorithm = 3 (Blossom Full-Quad)
// each Transfinite Line MUST have an even number of segments,
// otherwise gmsh emits "1D mesh cannot be divided by 2" and the
// boundary loop fails to close.  Every count below is even.
nx_A    = 18;   // bottom strip 0   ... 350   (Δx ≈ 19.4 m)
nx_C    = 6;    // building top   350 ... 450   (Δx ≈ 16.7 m)
nx_E    = 28;   // bottom strip 450 ... 1000  (Δx ≈ 19.6 m)
ny_low  = 6;    // y in [0,   100]              (Δy ≈ 16.7 m)
ny_high = 24;   // y in [100, 1000]             (Δy = 37.5 m)

// Point size only used as a fallback away from transfinite lines
lc = (xmax - xmin) / (nx_A + nx_C + nx_E);

Point(1)  = {xmin, ymin, 0, lc};   // bottom-left
Point(2)  = {bx0,  ymin, 0, lc};   // building left foot
Point(3)  = {bx0,  by1,  0, lc};   // building top-left
Point(4)  = {bx1,  by1,  0, lc};   // building top-right
Point(5)  = {bx1,  ymin, 0, lc};   // building right foot
Point(6)  = {xmax, ymin, 0, lc};   // bottom-right
Point(7)  = {xmax, by1,  0, lc};   // right-wall split
Point(8)  = {xmax, ymax, 0, lc};   // top-right
Point(9)  = {xmin, ymax, 0, lc};   // top-left
Point(10) = {xmin, by1,  0, lc};   // left-wall split

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

// -------- Transfinite Line: boundary node distribution --------
// Bottom + building profile
Transfinite Line{1}  = nx_A    + 1;
Transfinite Line{2}  = ny_low  + 1;
Transfinite Line{3}  = nx_C    + 1;
Transfinite Line{4}  = ny_low  + 1;
Transfinite Line{5}  = nx_E    + 1;

// Periodic right wall (lower + upper)
Transfinite Line{6}  = ny_low  + 1;
Transfinite Line{7}  = ny_high + 1;

// Top
Transfinite Line{8}  = (nx_A + nx_C + nx_E) + 1;

// Periodic left wall (upper + lower) — match the right wall
Transfinite Line{9}  = ny_high + 1;
Transfinite Line{10} = ny_low  + 1;

// -------- Physical groups --------
Physical Surface("domain")    = {1};
Physical Curve("bottom")      = {1, 2, 3, 4, 5};
Physical Curve("top")         = {8};
Physical Curve("periodicx")   = {6, 7, 9, 10};

Mesh.ElementOrder           = 1;
Mesh.Algorithm              = 8;   // Frontal-Delaunay for Quads (quad-friendly triangulation)
Mesh.RecombinationAlgorithm = 1;   // Blossom (NOT Full-Quad).
// Full-Quad (= 3) is incompatible with Transfinite Line: gmsh issues
// "Full-quad recombination only compatible with transfinite meshes if
// those are performed first" and leaves stubborn boundary triangles in
// the .msh even when the summary reports "0 triangles", which trips
// GridapGmsh's "only one element type per dimension" check.
// Plain Blossom (= 1) pairs every triangle as long as the closed-loop
// boundary has an even total segment count (it does: 176).
//
// NB: do NOT set Mesh.SubdivisionAlgorithm here — that one produced a
// non-conforming mesh (the earlier "wedges from the origin" artefact).
