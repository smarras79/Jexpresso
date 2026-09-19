// ============================================================
// cylinder_M7: Mach-7 laminar flow over a circular cylinder.
//
// CompEuler/shock_circle's geometry — box [0,3] x [-1,1] with a
// cylinder of radius 0.2 centred at (1,0) — meshed as
// UNSTRUCTURED QUADS with the cell size clustered towards the
// cylinder so that the LAMINAR BOUNDARY LAYER and the surface
// heat flux are resolved.
//
// THE CIRCLE IS CURVED BY :exact_geometry, NOT BY THIS FILE.
// gmsh writes a LINEAR grid: the four Circle arcs below come back
// as straight segments whose endpoints merely happen to sit on
// the circle, so filling those elements with LGL nodes would put
// every high-order node on a CHORD and the wall the solver sees
// would be a polygon however large :nop is. On a no-slip wall
// that is fatal twice over: the polygon corners shed spurious
// vorticity, and the wall-normal direction the heat flux is
// computed along is wrong by O(h) at every node. The deck
// therefore carries
//
//   :exact_geometry => Dict("cylinder" => (:circle, 1.0, 0.0, 0.2))
//
// and src/kernel/mesh/exact_geometry.jl snaps the high-order
// boundary nodes onto the true circle and blends the element
// interiors (Kopriva, J. Sci. Comput. 26(3):301-327, 2006, §3).
// This is not optional on a curved wall — it is the default
// expectation for one. "cylinder" is its own Physical Curve for
// exactly that reason.
//
// RESOLUTION. With the free stream of initialize.jl (M = 7,
// T = 125 K, p = 5 Pa) the cylinder Reynolds number is
// Re_D = 1.0e4, so the laminar boundary layer is
// delta ~ R/sqrt(Re_R) = 2.8 mm on a 200 mm radius. The size
// field below puts 2.2 mm cells at the wall, which at :nop => 4
// is about six LGL nodes across delta — enough for a surface heat
// flux, and a MODEST cluster rather than a y+ = 1 mesh: the far
// field is 60-77 mm, a ratio of ~30, and there is no anisotropic
// stretching anywhere (these are near-square quads).
//
// Measured: 14296 quads, all quadrilateral, minSICN 0.609, no
// inverted cells, edge lengths 0.00223-0.0770 m.
//
// Generate with:
//   gmsh -2 cylinder_M7.geo -o cylinder_M7.msh
// (the Mesh.* options and the size field are part of this file,
//  so it reproduces the committed mesh exactly).
// ============================================================

xmin=0.0; xmax=3.0; ymin=-1.0; ymax=1.0;
xc=1.0; yc=0.0; radius=0.2;
lc_far=0.060; lc_wall=0.005;
Point(1)={xmin,ymin,0,lc_far}; Point(2)={xmax,ymin,0,lc_far};
Point(3)={xmax,ymax,0,lc_far}; Point(4)={xmin,ymax,0,lc_far};
Point(5)={xc,yc,0,lc_wall};
Point(6)={xc+radius,yc,0,lc_wall}; Point(7)={xc,yc+radius,0,lc_wall};
Point(8)={xc-radius,yc,0,lc_wall}; Point(9)={xc,yc-radius,0,lc_wall};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Circle(5)={6,5,7}; Circle(6)={7,5,8}; Circle(7)={8,5,9}; Circle(8)={9,5,6};
Curve Loop(1)={1,2,3,4}; Curve Loop(2)={5,6,7,8};
Plane Surface(1)={1,2};
Field[1]=Distance; Field[1].CurvesList={5,6,7,8}; Field[1].Sampling=400;
Field[2]=Threshold; Field[2].InField=1;
Field[2].SizeMin=lc_wall; Field[2].SizeMax=lc_far;
Field[2].DistMin=0.10; Field[2].DistMax=0.90;
Background Field=2;
Mesh.MeshSizeExtendFromBoundary=0; Mesh.MeshSizeFromPoints=0; Mesh.MeshSizeFromCurvature=0;
Recombine Surface{1};
Mesh.Algorithm=8; Mesh.RecombineAll=1; Mesh.RecombinationAlgorithm=3; Mesh.Smoothing=200;
Mesh.ElementOrder=1;
Physical Curve("inflow")={4}; Physical Curve("outflow")={2};
Physical Curve("top")={3}; Physical Curve("bottom")={1};
Physical Curve("cylinder")={5,6,7,8};
Physical Surface("domain")={1};
