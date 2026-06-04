// ============================================================
// city2d: 2D flow around a single building
//
// Domain:   [0, 1000] x [0, 1000] m
// Building: 100 x 100 m, lower-left at (x,y) = (350, 0)
//
// Multi-block transfinite layout (5 conforming quad blocks):
//
//     P9 ------- P8 ----- P7 -------- P6   y = 1000
//     |    B     |   C    |     E     |
//     |          |        |           |
//     P10 ----- P11 ---- P12 -------- P5   y = 100
//     |    A     | bldg.  |     D     |
//     |          |(solid) |           |
//     P1 ------- P2 ----- P3 -------- P4   y = 0
//     x = 0    x=350    x=450        x = 1000
//
// Physical groups (matched to Jexpresso conventions):
//   "bottom"    : full bottom of the domain including the building profile
//   "top"       : top wall
//   "periodicx" : both vertical walls (left and right)
//   "domain"    : interior
//
// Mesh: all-quad, conforming (transfinite + Recombine).
// Node counts on the left and right periodic walls match by construction
// (nA_y = nD_y on the lower segment, nB_y = nE_y on the upper segment),
// which is required by Jexpresso's "periodicx" pairing.
//
// Generate with:
//   gmsh -2 city2d.geo -o ../../../meshes/gmsh_grids/city2d.msh
// ============================================================

xmin = 0.0;
xmax = 1000.0;
ymin = 0.0;
ymax = 1000.0;

// Building extent
bx0 = 350.0;
bx1 = 450.0;
by1 = 100.0;

// Number of element segments along each strip.
// (Δx, Δy ≈ 10 m near the building, ≈ 20 m far from it.)
nx_A = 35;   // x in [0,   350]
nx_C = 10;   // x in [350, 450]   (above the building)
nx_E = 55;   // x in [450, 1000]
ny_low  = 10;   // y in [0,   100]
ny_high = 45;   // y in [100, 1000]

// -------- Points --------
Point(1)  = {xmin, ymin, 0};
Point(2)  = {bx0,  ymin, 0};
Point(3)  = {bx1,  ymin, 0};
Point(4)  = {xmax, ymin, 0};
Point(5)  = {xmax, by1,  0};
Point(6)  = {xmax, ymax, 0};
Point(7)  = {bx1,  ymax, 0};
Point(8)  = {bx0,  ymax, 0};
Point(9)  = {xmin, ymax, 0};
Point(10) = {xmin, by1,  0};
Point(11) = {bx0,  by1,  0};
Point(12) = {bx1,  by1,  0};

// -------- Outer boundary lines --------
Line(1)  = {1,  2};   // bottom, A-strip
Line(2)  = {2, 11};   // building left wall  (bottom)
Line(3)  = {11,12};   // building top        (bottom)
Line(4)  = {12, 3};   // building right wall (bottom)
Line(5)  = {3,  4};   // bottom, D-strip
Line(6)  = {4,  5};   // right wall, lower  (periodicx)
Line(7)  = {5,  6};   // right wall, upper  (periodicx)
Line(8)  = {6,  7};   // top, E-strip
Line(9)  = {7,  8};   // top, C-strip
Line(10) = {8,  9};   // top, B-strip
Line(11) = {9, 10};   // left wall, upper   (periodicx)
Line(12) = {10, 1};   // left wall, lower   (periodicx)

// -------- Interior (block-interface) lines --------
Line(13) = {10, 11};  // horizontal A-B interface
Line(14) = {11,  8};  // vertical   AB-C interface (left of building, above)
Line(15) = {12,  7};  // vertical   C-DE interface (right of building, above)
Line(16) = {12,  5};  // horizontal D-E interface

// -------- Curve loops (CCW) --------
Curve Loop(1) = { 1,  2, -13,  12};         // A : lower-left
Curve Loop(2) = {13, 14,  10,  11};         // B : upper-left
Curve Loop(3) = { 3, 15,   9, -14};         // C : above the building
Curve Loop(4) = { 5,  6, -16,  4};          // D : lower-right
Curve Loop(5) = {16,  7,   8, -15};         // E : upper-right

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

// -------- Transfinite distributions --------
// Lower strips (y in [0, 100]):
Transfinite Curve{1}   = nx_A + 1;
Transfinite Curve{5}   = nx_E + 1;
Transfinite Curve{2}   = ny_low + 1;     // building left wall
Transfinite Curve{4}   = ny_low + 1;     // building right wall
Transfinite Curve{12}  = ny_low + 1;     // left  wall, lower
Transfinite Curve{6}   = ny_low + 1;     // right wall, lower  (matches Line 12)
Transfinite Curve{13}  = nx_A + 1;       // A-B interface
Transfinite Curve{16}  = nx_E + 1;       // D-E interface

// Upper strips (y in [100, 1000]):
Transfinite Curve{10}  = nx_A + 1;
Transfinite Curve{9}   = nx_C + 1;
Transfinite Curve{8}   = nx_E + 1;
Transfinite Curve{3}   = nx_C + 1;       // building top
Transfinite Curve{11}  = ny_high + 1;    // left  wall, upper
Transfinite Curve{7}   = ny_high + 1;    // right wall, upper  (matches Line 11)
Transfinite Curve{14}  = ny_high + 1;    // AB-C interface
Transfinite Curve{15}  = ny_high + 1;    // C-DE interface

Transfinite Surface{1};
Transfinite Surface{2};
Transfinite Surface{3};
Transfinite Surface{4};
Transfinite Surface{5};

Recombine Surface{1, 2, 3, 4, 5};

// -------- Physical groups --------
Physical Surface("domain")    = {1, 2, 3, 4, 5};
Physical Curve("bottom")      = {1, 2, 3, 4, 5};   // bottom + building profile
Physical Curve("top")         = {8, 9, 10};
Physical Curve("periodicx")   = {6, 7, 11, 12};

// -------- Output options --------
Mesh.ElementOrder     = 1;
Mesh.MshFileVersion   = 2.2;    // Compatible with Jexpresso's GMSH reader
Mesh.RecombineAll     = 1;
