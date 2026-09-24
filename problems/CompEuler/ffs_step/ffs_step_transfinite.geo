// ============================================================
// ffs_step_transfinite: structured multi-block transfinite quad
// mesh for the Mach-3 forward-facing-step wind tunnel
// (Woodward & Colella 1984; Section 5.1 of Nazarov & Hoffman,
//  Int. J. Numer. Meth. Fluids 2013; 71:339-357).
//
// Geometry (Section 5.1):
//   Tunnel:  length 3, height 1   ->  [0, 3] x [0, 1]
//   Step:    height 0.2, located 0.6 from the inflow
//            solid region = [0.6, 3] x [0, 0.2]  (EXCLUDED)
//   The flow domain is therefore the L-shape
//      ( [0,3]x[0,1] )  \  ( [0.6,3]x[0,0.2] )
//
// Flow conditions (Section 5.1):
//   Supersonic inflow  M = 3 : rho = 1.4, m = (4.2, 0), e = 8.8
//   All characteristics leave at outflow -> nothing imposed there
//   Slip (reflecting) condition on every solid wall
//
//   P7 -------------- P6 ------------------------- P5   y = 1.0  (top)
//   |        B        |              C              |
//   |                 |                             |
//   P8 -------------- P3 ------------------------- P4   y = 0.2  (step top)
//   |        A        | <--- solid step (excluded) --->
//   |                 |
//   P1 -------------- P2                                 y = 0.0  (floor)
//   x = 0           x = 0.6                       x = 3.0
//
// Three conforming transfinite quad blocks:
//   A : [0,   0.6] x [0,   0.2]   (inflow channel, below step level)
//   B : [0,   0.6] x [0.2, 1.0]   (inflow channel, above step level)
//   C : [0.6, 3.0] x [0.2, 1.0]   (channel above the step)
//
// Generate with:
//   gmsh -2 ffs_step_transfinite.geo -o ffs_step_transfinite.msh
// ============================================================

// -------- Domain / step extents --------
xmin   = 0.0;
xmax   = 3.0;
ymin   = 0.0;
ymax   = 1.0;
xstep  = 0.6;   // streamwise location of the step face
ystep  = 0.2;   // step height

// -------- Resolution --------
// ref = 1 : coarse, h = 0.025   (120 x 40 over the bounding box)
// ref = 2 : medium, h = 0.0125  (240 x 80, ~ the 80x240 mesh of the paper)
ref = 1;

nx_in   = 24 * ref;   // x in [0,   0.6]   (Dx = 0.025/ref)
nx_out  = 96 * ref;   // x in [0.6, 3.0]   (Dx = 0.025/ref)
ny_low  =  8 * ref;   // y in [0,   0.2]   (Dy = 0.025/ref)
ny_high = 32 * ref;   // y in [0.2, 1.0]   (Dy = 0.025/ref)

// -------- Points --------
Point(1) = {xmin,  ymin,  0};   // floor, inflow corner
Point(2) = {xstep, ymin,  0};   // floor, base of the step
Point(3) = {xstep, ystep, 0};   // CONVEX step corner (top of the step face)
Point(4) = {xmax,  ystep, 0};   // outflow, bottom (on the step top)
Point(5) = {xmax,  ymax,  0};   // outflow, top
Point(6) = {xstep, ymax,  0};   // top wall, above the step corner
Point(7) = {xmin,  ymax,  0};   // top wall, inflow corner
Point(8) = {xmin,  ystep, 0};   // inflow wall, at step level

// -------- Boundary lines --------
Line(1) = {1, 2};   // floor  (wall)          x in [0,   0.6]
Line(2) = {2, 3};   // step vertical face (wall)  y in [0, 0.2]
Line(3) = {3, 4};   // step top (wall)        x in [0.6, 3.0]
Line(4) = {4, 5};   // outflow                y in [0.2, 1.0]
Line(5) = {5, 6};   // top wall (wall)        x in [0.6, 3.0]
Line(6) = {6, 7};   // top wall (wall)        x in [0,   0.6]
Line(7) = {7, 8};   // inflow, upper          y in [0.2, 1.0]
Line(8) = {8, 1};   // inflow, lower          y in [0,   0.2]

// -------- Interior (block-interface) lines --------
Line(9)  = {8, 3};   // horizontal A | B interface  (y = 0.2, x in [0,0.6])
Line(10) = {3, 6};   // vertical   B | C interface  (x = 0.6, y in [0.2,1])

// -------- Curve loops (CCW) --------
Curve Loop(1) = { 1,  2, -9,  8};   // A : lower-left
Curve Loop(2) = { 9, 10,  6,  7};   // B : upper-left
Curve Loop(3) = { 3,  4,  5, -10};  // C : upper-right (above step)

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// -------- Transfinite line distributions --------
// streamwise (x) divisions
Transfinite Curve{1} = nx_in  + 1;   // floor, left
Transfinite Curve{9} = nx_in  + 1;   // A | B interface
Transfinite Curve{6} = nx_in  + 1;   // top, left
Transfinite Curve{3} = nx_out + 1;   // step top
Transfinite Curve{5} = nx_out + 1;   // top, right

// vertical (y) divisions
Transfinite Curve{2}  = ny_low  + 1;   // step face
Transfinite Curve{8}  = ny_low  + 1;   // inflow, lower
Transfinite Curve{4}  = ny_high + 1;   // outflow
Transfinite Curve{7}  = ny_high + 1;   // inflow, upper
Transfinite Curve{10} = ny_high + 1;   // B | C interface

// -------- Transfinite surfaces + recombine to quads --------
Transfinite Surface{1};
Transfinite Surface{2};
Transfinite Surface{3};
Recombine Surface{1, 2, 3};

// -------- Physical groups --------
// All solid walls use a slip / reflecting BC (paper Sec. 3.3 & 5.1).
Physical Surface("domain")  = {1, 2, 3};
Physical Curve("inflow")    = {7, 8};            // left wall, x = 0
Physical Curve("outflow")   = {4};               // right wall, x = 3
Physical Curve("wall")      = {1, 2, 3, 5, 6};   // floor, step face, step top, top wall

Mesh.ElementOrder = 1;
