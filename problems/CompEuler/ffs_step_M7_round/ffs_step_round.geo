// ============================================================
// ffs_step_round: the forward-facing-step wind tunnel of
// ffs_step / ffs_step_M7, with the CONVEX STEP CORNER FILLETED,
// meshed as UNSTRUCTURED QUADS and curved by :exact_geometry.
//
// Geometry, otherwise as ffs_step_M7:
//   Tunnel:  [0, 3] x [0, 1]
//   Step:    solid region [0.6, 3] x [0, 0.2]  (EXCLUDED)
//   Fillet:  arc of radius r = 0.05 m centred at (0.65, 0.15),
//            tangent to the step face at (0.6, 0.15) and to the
//            step top at (0.65, 0.2).
//
// WHY FILLET. At a sharp (0.6, 0.2) the fluid's interior angle is
// 270 degrees — a POINT SINGULARITY of the Euler solution. The
// expansion there is a Prandtl-Meyer fan centred on that one
// point, and what it sheds downstream is an entropy layer, a
// contact-type feature that convects and never self-heals. A
// residual-based viscosity cannot help: an expansion fan is a
// SMOOTH solution, so the sensor correctly returns almost nothing
// in it. Woodward & Colella (1984) section IV call this corner "a
// singular point of the flow" and patch the cells beside it, at
// Mach 3.
//
// WHY UNSTRUCTURED, AND WHY :exact_geometry.
//
// Two separate mistakes are easy to make here, and the first
// version of this mesh made both.
//
//   1. FORCING A TRANSFINITE BLOCK AROUND THE FILLET. Putting the
//      step face and the arc on one side of a structured block
//      and a straight line on the opposite side makes the
//      transfinite map shear the cells near the corner: measured
//      minSICN 0.235, against 1.000 for the sharp rectangular
//      mesh. A badly distorted element at the exact place the
//      case is failing is worse than the sharp corner it was
//      meant to cure. Unstructured quads (gmsh Algorithm 11,
//      quasi-structured) give minSICN 0.711 here and stay
//      near-Cartesian away from the fillet.
//
//   2. LEAVING THE ARC AS A POLYLINE. gmsh writes a LINEAR grid:
//      the Circle below comes back as a few straight segments
//      whose endpoints happen to sit on the circle. Filling those
//      elements with LGL nodes puts every high-order node on the
//      CHORD, so however large :nop is, the wall the solver sees
//      is a polygon — and a slip wall spuriously generates
//      vorticity at every polygon corner. Rounding the corner and
//      then discretizing it as a polygon just trades one set of
//      corners for several.
//
//      That is what "fillet" as its OWN physical group is for.
//      The deck says
//
//        :exact_geometry => Dict("fillet" => (:circle, 0.65, 0.15, 0.05))
//
//      and src/kernel/mesh/exact_geometry.jl snaps the high-order
//      boundary nodes onto the true circle and blends the element
//      interiors (Kopriva, J. Sci. Comput. 26(3):301-327, 2006,
//      section 3 — linear blending transfinite map, which stays in
//      P^N and so preserves the discrete metric identities and the
//      free stream). The wall then sits on the circle to machine
//      precision, and the arc needs only a FEW linear segments:
//      at :nop => 4 three segments already give 13 boundary nodes
//      exactly on the circle. That is the whole benefit — do not
//      over-refine the arc to chase the geometry.
//
// The ramp's sharp leading edge is the same kind of point, so
// whatever this run shows transfers to rampCaoEtAl2021.
//
// Generate with:
//   gmsh -2 ffs_step_round.geo -o ffs_step_round.msh
// (the Mesh.* options below are part of the file, so this
//  reproduces the committed mesh exactly).
// ============================================================

// -------- Resolution --------
// Algorithm 11 (quasi-structured quad) subdivides, so lc = 0.05
// lands on an effective h of about 0.025 — the spacing of the
// sharp ffs_step_M7 mesh. ref = 2 halves it.
ref = 1;
lc  = 0.05 / ref;

// -------- Domain / step extents --------
xmin  = 0.0;
xmax  = 3.0;
ymin  = 0.0;
ymax  = 1.0;
xstep = 0.6;   // streamwise location of the step face
ystep = 0.2;   // step height

// -------- Fillet radius --------
// rfac = how many elements of the effective h = lc/2 the 90-degree
// turn is spread over. The exact-geometry snap means the arc does
// NOT need many segments; rfac controls the PHYSICAL radius, which
// is what the flow sees. Raise it to separate "the corner is
// singular" from "the corner is sharp".
rfac = 2;
r    = rfac * lc / 2;

// -------- Points --------
Point(1)  = {xmin,      ymin,      0, lc};   // floor, inflow corner
Point(2)  = {xstep,     ymin,      0, lc};   // floor, base of the step
Point(3)  = {xstep,     ystep - r, 0, lc};   // tangent point, step face
Point(4)  = {xstep + r, ystep,     0, lc};   // tangent point, step top
Point(5)  = {xmax,      ystep,     0, lc};   // outflow, bottom
Point(6)  = {xmax,      ymax,      0, lc};   // outflow, top
Point(7)  = {xmin,      ymax,      0, lc};   // top wall, inflow corner
Point(10) = {xstep + r, ystep - r, 0, lc};   // ARC CENTRE (not on the boundary)

// -------- Boundary curves --------
Line(1)   = {1, 2};        // floor              (wall)
Line(2)   = {2, 3};        // step vertical face (wall)
Circle(3) = {3, 10, 4};    // THE FILLET         (its own group: see above)
Line(4)   = {4, 5};        // step top           (wall)
Line(5)   = {5, 6};        // outflow
Line(6)   = {6, 7};        // top wall           (wall)
Line(7)   = {7, 1};        // inflow

Curve Loop(1)    = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};
Recombine Surface{1};

// -------- Meshing options --------
// Chosen by measurement over Algorithm 6/8/11 x RecombinationAlgorithm
// 2/3 x SubdivisionAlgorithm 0/1, scoring minSICN and the edge-length
// spread. This combination gives 4103 quads (the sharp mesh has 4032),
// minSICN 0.711, edge lengths 0.0189-0.0371 m. Algorithm 8 scores a
// slightly better minSICN (0.790) but its shortest edge is 0.0107,
// which would cut the explicit time step for no gain.
Mesh.Algorithm               = 11;   // quasi-structured quad
Mesh.RecombineAll            = 1;
Mesh.RecombinationAlgorithm  = 2;    // blossom
Mesh.Smoothing               = 100;
Mesh.ElementOrder            = 1;    // LINEAR: exact_geometry does the curving

// -------- Physical groups --------
// "fillet" is SEPARATE from "wall" for one reason only: it is the tag
// :exact_geometry names. user_bc.jl treats it as a free-slip wall like
// any other, which is what its `else` branch already does.
Physical Curve("inflow")   = {7};
Physical Curve("outflow")  = {5};
Physical Curve("wall")     = {1, 2, 4, 6};   // floor, step face, step top, roof
Physical Curve("fillet")   = {3};            // the arc — snapped to the circle
Physical Surface("domain") = {1};
