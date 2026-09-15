// ============================================================================
// ramp15.geo — 2D hypersonic compression-ramp grid, case G1 of
//
//   S. Cao, J. Hao, I. Klioutchnikov, H. Olivier, C.-Y. Wen,
//   "Unsteady effects in a hypersonic compression ramp flow with laminar
//    separation", J. Fluid Mech. 912, A3 (2021), doi:10.1017/jfm.2020.1093
//
// reduced to the x-y plane (the paper's DNS is 3D and periodic in z; the
// paper's own 2D run is the base flow of its global stability analysis,
// Section 3.2, and that is what this grid is for).
//
// GEOMETRY (Section 2.2; model of Roghelia et al., Exp. Fluids 58:139, 2017)
//
//   sharp leading edge at the origin, flat plate of L = 100 mm, then a ramp
//   deflected by 15 deg, also 100 mm long MEASURED ALONG THE RAMP SURFACE
//   (set ramp_len below to 0.1/Cos(alpha) if you read the paper's "both
//   100 mm" as a horizontal extent instead — that moves the outflow from
//   x/L = 1.966 to 2.0 and changes nothing else).
//
//   20 grid points are placed in the 1 mm ahead of the leading edge, as in
//   Section 2.3, so that the free stream is established before the plate
//   starts and the leading-edge singularity does not sit on the inflow.
//
//   P8 ------- P7 ----------------------- P6                 y = H
//   |    A     |             B            |    \ C
//   |          |                          |         \
//   |          |                          |              \  P5
//   P1 ------- P2 ======================= P3                 y = 0
//                                            \  ramp 15 deg
//   x = -1 mm  x = 0                x = 100 mm    \  P4
//
//   A : free-stream strip ahead of the leading edge (bottom = "symmetry")
//   B : above the flat plate                        (bottom = "wall")
//   C : above the ramp, a parallelogram — the upper boundary is the wall
//       contour shifted VERTICALLY by H, so every grid line of constant
//       streamwise index is vertical and block C is a uniform shear of a
//       rectangle (no metric distortion, no skew beyond the 15 deg).
//
//   H = 60 mm.  The separation shock leaves the plate near x/L = 0.59 at
//   roughly 20 deg, so it crosses x = x_end some 36 mm below the upper
//   boundary; the whole shock system stays inside the domain and the
//   free-stream Dirichlet condition on "top" is never asked to swallow a
//   discontinuity.
//
// RESOLUTION (case G1 of Section 2.2: 1080 x 240 points)
//
//   Jexpresso is a spectral-element code, so "points" are the LGL nodes of
//   the elements: at :nop => 4 a curve of N elements carries 4N+1 nodes.
//        x :   5 + 134 + 130 = 269 elements -> 1077 nodes  (paper 1080)
//        y :               60 elements      ->  241 nodes  (paper  240)
//
//   The wall-normal progression 1.0808 over H = 60 mm puts the first
//   element at 4.62e-5 m; the first LGL interval of a 4th-order element is
//   0.17267 of it, i.e. Dy_wall = 7.97e-6 m — the paper's 8e-6 m, hence its
//   y+ ~ 0.3 just upstream of separation.  That places ~64 nodes inside the
//   1.38 mm boundary layer at separation against the paper's 75; the LGL
//   clustering at every element face recovers most of the difference.
//
// This file is the RECORD of the grid, not the way it is built: no gmsh is
// needed, generate_mesh.py writes ramp15.msh directly from the same
// numbers (identical Progression formula, s_i = (p^i-1)/(p^N-1)).  Either
// route gives the same mesh:
//
//   gmsh -2 ramp15.geo -o ramp15.msh
//   python3 generate_mesh.py
// ============================================================================

// -------- Geometry --------
L         = 0.1;                 // flat plate length [m]
alpha     = 15.0*Pi/180.0;       // ramp deflection [rad]
ramp_len  = 0.1;                 // ramp length ALONG the ramp surface [m]
x_up      = 0.001;               // free-stream strip ahead of the leading edge
H         = 0.06;                // domain height above the wall contour [m]

x_end     = L + ramp_len*Cos(alpha);
y_end     =     ramp_len*Sin(alpha);

// -------- Resolution --------
nx_up     =   5;   // uniform,          [-1 mm, 0]
nx_plate  = 134;   // Progression 1.006, clustered at the leading edge
nx_ramp   = 130;   // uniform,          along the ramp
ny        =  60;   // Progression 1.0808, clustered at the wall

px        = 1.006;
py        = 1.0808;

// -------- Points --------
Point(1) = {-x_up,  0.0,       0};   // inflow / symmetry corner
Point(2) = { 0.0,   0.0,       0};   // SHARP LEADING EDGE
Point(3) = { L,     0.0,       0};   // compression corner
Point(4) = { x_end, y_end,     0};   // ramp end, outflow bottom
Point(5) = { x_end, y_end + H, 0};   // outflow top
Point(6) = { L,     H,         0};   // top, above the corner
Point(7) = { 0.0,   H,         0};   // top, above the leading edge
Point(8) = {-x_up,  H,         0};   // inflow top

// -------- Curves --------
// Curves that face each other are given the SAME orientation so that one
// Progression applies to both and the blocks stay conforming.
Line(1)  = {1, 2};   // symmetry (ahead of the leading edge)
Line(2)  = {2, 3};   // wall, flat plate
Line(3)  = {3, 4};   // wall, ramp
Line(4)  = {4, 5};   // outflow
Line(5)  = {6, 5};   // top, above the ramp
Line(6)  = {7, 6};   // top, above the plate
Line(7)  = {8, 7};   // top, above the strip
Line(8)  = {1, 8};   // inflow
Line(9)  = {2, 7};   // block A | B interface
Line(10) = {3, 6};   // block B | C interface

// -------- Surfaces --------
Curve Loop(1) = { 1,  9, -7, -8};   Plane Surface(1) = {1};   // A
Curve Loop(2) = { 2, 10, -6, -9};   Plane Surface(2) = {2};   // B
Curve Loop(3) = { 3,  4, -5, -10};  Plane Surface(3) = {3};   // C

// -------- Transfinite distributions --------
Transfinite Curve{1, 7}  = nx_up    + 1;
Transfinite Curve{2, 6}  = nx_plate + 1 Using Progression px;
Transfinite Curve{3, 5}  = nx_ramp  + 1;
Transfinite Curve{4, 8, 9, 10} = ny + 1 Using Progression py;

Transfinite Surface{1};
Transfinite Surface{2};
Transfinite Surface{3};
Recombine Surface{1, 2, 3};

// -------- Physical groups (these names reach user_bc_dirichlet!) --------
Physical Curve("inflow",   2) = {8};
Physical Curve("outflow",  3) = {4};
Physical Curve("wall",     4) = {2, 3};
Physical Curve("top",      5) = {5, 6, 7};
Physical Curve("symmetry", 6) = {1};
Physical Surface("domain", 1) = {1, 2, 3};

Mesh.ElementOrder = 1;
