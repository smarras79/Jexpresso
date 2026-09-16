// ============================================================
// ffs_step_round: the forward-facing-step wind tunnel of
// ffs_step / ffs_step_M7, with the CONVEX STEP CORNER FILLETED.
//
// Everything else is the geometry of
// problems/CompEuler/ffs_step_M7/ffs_step_transfinite.geo:
//   Tunnel:  [0, 3] x [0, 1]
//   Step:    solid region [0.6, 3] x [0, 0.2]  (EXCLUDED)
//
// WHY. At (0.6, 0.2) the fluid's interior angle is 270 degrees —
// a reentrant corner, i.e. a POINT SINGULARITY of the Euler
// solution. The expansion there is a Prandtl-Meyer fan centred on
// that single point, and what it sheds downstream is an entropy
// layer. Woodward & Colella (1984) section IV say so outright
// ("the corner of the step is a singular point of the flow") and
// reset the state in the cells beside it for exactly this reason,
// at Mach 3. At Mach 7 the fan is far stronger and expands into
// far lower density, and a residual-based viscosity cannot help:
// an expansion fan is a SMOOTH solution, so the sensor correctly
// returns almost nothing in it.
//
// Filleting removes the singularity from the GEOMETRY instead of
// trying to patch the solution. The 90-degree turn is spread over
// an arc of radius r = rfac*h, so the fan is no longer centred on
// a point and every boundary node has a well-defined normal.
//
//   ...........                        ...........
//             |                                  \
//             |  <- 270 deg point       r         )  <- turn spread
//   __________|     singularity         __________/     over an arc
//
// The ramp's sharp leading edge is the same kind of point, so
// whatever this run shows transfers to rampCaoEtAl2021.
//
//   P8 -------------- P7 ------------------------- P6   y = 1.0  (top)
//   |        B        |              C              |
//   P9 -------------- P4 ------------------------- P5   y = 0.2
//   |        A       /  <- solid step (excluded) ->
//   |               P3   (0.6, 0.2-r), tangent point
//   |               |
//   P1 ------------ P2                                  y = 0.0  (floor)
//   x = 0         x = 0.6                        x = 3.0
//
// The three transfinite blocks of the sharp-corner mesh survive:
// the block corner simply moves from (0.6, 0.2) to (0.6+r, 0.2),
// and block A's right-hand side becomes the step face plus the
// arc. Element count 4033 against the sharp mesh's 4032, and
// the edge lengths run 0.0189-0.0260 m against that mesh's
// uniform 0.025, a max/min of 1.38 — this is NOT a stretched
// grid, but Dx_min IS 1.32x smaller, so the printed CFLs will be
// that much higher at the same :Dt. (n_arc = 3 was chosen over 4
// for exactly this: 4 gives max/min 1.70.)
//
// Generate with:
//   gmsh -2 ffs_step_round.geo -o ffs_step_round.msh
// ============================================================

// -------- Domain / step extents --------
xmin   = 0.0;
xmax   = 3.0;
ymin   = 0.0;
ymax   = 1.0;
xstep  = 0.6;   // streamwise location of the step face
ystep  = 0.2;   // step height

// -------- Resolution --------
// ref = 1 : h = 0.025   (as the sharp-corner mesh at ref = 1)
// ref = 2 : h = 0.0125
ref  = 1;
h    = 0.025 / ref;

// -------- Fillet radius --------
// rfac = how many elements the 90-degree turn is spread over.
// 2 is the smallest radius that is still resolved; raise it to
// separate "the corner is singular" from "the corner is sharp".
rfac = 2;
r    = rfac * h;

// -------- Divisions --------
nx_in   = 25 * ref;   // floor and the A|B interface
nx_out  = 94 * ref;   // step top and the top wall, right of x = 0.6+r
n_face  =  6 * ref;   // step face, y in [0, 0.2-r]
n_arc   =  3 * ref;   // the fillet arc itself
ny_low  = n_face + n_arc;   // inflow, y in [0, 0.2]
ny_high = 32 * ref;   // y in [0.2, 1.0]

// -------- Points --------
Point(1)  = {xmin,       ymin,        0};   // floor, inflow corner
Point(2)  = {xstep,      ymin,        0};   // floor, base of the step
Point(3)  = {xstep,      ystep - r,   0};   // tangent point on the step face
Point(4)  = {xstep + r,  ystep,       0};   // tangent point on the step top
Point(5)  = {xmax,       ystep,       0};   // outflow, bottom
Point(6)  = {xmax,       ymax,        0};   // outflow, top
Point(7)  = {xstep + r,  ymax,        0};   // top wall, above the block corner
Point(8)  = {xmin,       ymax,        0};   // top wall, inflow corner
Point(9)  = {xmin,       ystep,       0};   // inflow wall, at step-top level
Point(10) = {xstep + r,  ystep - r,   0};   // ARC CENTRE (not on the boundary)

// -------- Boundary curves --------
Line(1)   = {1, 2};          // floor                    (wall)
Line(2)   = {2, 3};          // step vertical face       (wall)
Circle(3) = {3, 10, 4};      // THE FILLET               (wall)
Line(4)   = {4, 5};          // step top                 (wall)
Line(5)   = {5, 6};          // outflow
Line(6)   = {6, 7};          // top wall, right          (wall)
Line(7)   = {7, 8};          // top wall, left           (wall)
Line(8)   = {8, 9};          // inflow, upper
Line(9)   = {9, 1};          // inflow, lower

// -------- Interior (block-interface) curves --------
Line(10) = {9, 4};   // horizontal A | B interface  (y = 0.2, x in [0, 0.6+r])
Line(11) = {4, 7};   // vertical   B | C interface  (x = 0.6+r, y in [0.2, 1])

// -------- Curve loops (CCW) --------
Curve Loop(1) = { 1,  2,  3, -10,  9};   // A : lower-left, carries the fillet
Curve Loop(2) = {10, 11,  7,   8};       // B : upper-left
Curve Loop(3) = { 4,  5,  6, -11};       // C : upper-right (above the step)

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// -------- Transfinite curve distributions --------
// streamwise (x)
Transfinite Curve{1}  = nx_in  + 1;   // floor
Transfinite Curve{10} = nx_in  + 1;   // A | B interface
Transfinite Curve{7}  = nx_in  + 1;   // top wall, left
Transfinite Curve{4}  = nx_out + 1;   // step top
Transfinite Curve{6}  = nx_out + 1;   // top wall, right

// vertical (y).  Block A's right-hand side is curve 2 PLUS curve 3,
// so n_face + n_arc must equal ny_low on the opposite side (curve 9).
Transfinite Curve{2}  = n_face  + 1;   // step face
Transfinite Curve{3}  = n_arc   + 1;   // fillet arc
Transfinite Curve{9}  = ny_low  + 1;   // inflow, lower
Transfinite Curve{5}  = ny_high + 1;   // outflow
Transfinite Curve{8}  = ny_high + 1;   // inflow, upper
Transfinite Curve{11} = ny_high + 1;   // B | C interface

// -------- Transfinite surfaces + recombine to quads --------
// The corners are named explicitly because block A has five curves
// on four sides (the step face and the arc share one side).
Transfinite Surface{1} = {1, 2, 4, 9};
Transfinite Surface{2} = {9, 4, 7, 8};
Transfinite Surface{3} = {4, 5, 6, 7};
Recombine Surface{1, 2, 3};

// -------- Physical groups --------
// Same three names as the sharp-corner mesh, so user_bc.jl is unchanged
// apart from the corner special case, which the fillet makes unnecessary.
Physical Surface("domain")  = {1, 2, 3};
Physical Curve("inflow")    = {8, 9};                  // left wall, x = 0
Physical Curve("outflow")   = {5};                     // right wall, x = 3
Physical Curve("wall")      = {1, 2, 3, 4, 6, 7};      // floor, face, FILLET, step top, roof

Mesh.ElementOrder = 1;
