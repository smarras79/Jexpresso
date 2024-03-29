// -----------------------------------------------------------------------------
//
//  Gmsh GEO tutorial 10
//
//  Mesh size fields
//
// -----------------------------------------------------------------------------

// In addition to specifying target mesh sizes at the points of the geometry
// (see `t1.geo') or using a background mesh (see `t7.geo'), you can use general
// mesh size "Fields".

// -----------------------------------------------------------------------------
// USER DEFINED QUANTITIES:
size_scaling = 100000;
xscaling = 100000;
yscaling = 100000;
zscaling = 1;

zlevels = 10;

lc = 0.1*size_scaling;
xmin = -1*xscaling; xmax = 1*xscaling;
ymin = -1*yscaling; ymax = 1*yscaling;
zmin = 0; zmax = 24000;

xc = (xmin + xmax)/2;
yc = (ymin + ymax)/2;
zc = 0.0;

//Refined box:
coarse_to_fine_ratio = 5; //--> how many times finer is the refined mesh w.r.t. coarsest mesh size
refined_extension = 50000;
xbox_min = (xc - refined_extension); xbox_max = (xc + refined_extension);
ybox_min = (yc - refined_extension); ybox_max = (yc + refined_extension);


// END USER DEFINED
// -----------------------------------------------------------------------------

Point(1) = {xmin, ymin, zmin, lc}; Point(2) = {xmax, ymin, zmin, lc};
Point(3) = {xmax, ymax, zmin,lc};  Point(4) = {xmin, ymax, zmin, lc};
//Point(5) = {xc, yc, zc, lc};

Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};

Curve Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

// Say we would like to obtain mesh elements with size lc/30 near curve 2 and
// point 5, and size lc elsewhere. To achieve this, we can use two fields:
// "Distance", and "Threshold". We first define a Distance field (`Field[1]') on
// points 5 and on curve 2. This field returns the distance to point 5 and to
// (100 equidistant points on) curve 2.
//Field[1] = Distance;
//Field[1].PointsList = {5};
//Field[1].CurvesList = {2};
//Field[1].Sampling = 100;

Recombine Surface {6}; //TRI --> QUAD

// We then define a `Threshold' field, which uses the return value of the
// `Distance' field 1 in order to define a simple change in element size
// depending on the computed distances
//
// SizeMax -                     /------------------
//                              /
//                             /
//                            /
// SizeMin -o----------------/
//          |                |    |
//        Point         DistMin  DistMax
/*Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc / 30;
Field[2].SizeMax = lc;
Field[2].DistMin = 0.15;
Field[2].DistMax = 0.5;
*/
// Say we want to modulate the mesh element sizes using a mathematical function
// of the spatial coordinates. We can do this with the MathEval field:
//Field[3] = MathEval;
//Field[3].F = "Cos(4*3.14*x) * Sin(4*3.14*y) / 10 + 0.101";


// We could also combine MathEval with values coming from other fields. For
// example, let's define a `Distance' field around point 1
Field[4] = Distance;
Field[4].PointsList = {1};

// We can then create a `MathEval' field with a function that depends on the
// return value of the `Distance' field 4, i.e., depending on the distance to
// point 1 (here using a cubic law, with minimum element size = lc / 100)
//Field[5] = MathEval;
//Field[5].F = Sprintf("F4^3 + %g", lc / 100);

// We could also use a `Box' field to impose a step change in element sizes
// inside a box
Field[6] = Box;
Field[6].VIn = lc / coarse_to_fine_ratio;
Field[6].VOut = lc;
Field[6].XMin = xbox_min;
Field[6].XMax = xbox_max;
Field[6].YMin = ybox_min;
Field[6].YMax = ybox_max;
Field[6].Thickness = 1*size_scaling; //larger value --> more extended refined region

// Many other types of fields are available: see the reference manual for a
// complete list. You can also create fields directly in the graphical user
// interface by selecting `Define->Size fields' in the `Mesh' module.

// Let's use the minimum of all the fields as the background mesh size field
Field[7] = Min;
Field[7].FieldsList = {2,  6};
Background Field = 7;


// To determine the size of mesh elements, Gmsh locally computes the minimum of
//
// 1) the size of the model bounding box;
// 2) if `Mesh.MeshSizeFromPoints' is set, the mesh size specified at
//    geometrical points;
// 3) if `Mesh.MeshSizeFromCurvature' is positive, the mesh size based on
//    curvature (the value specifying the number of elements per 2 * pi rad);
// 4) the background mesh size field;
// 5) any per-entity mesh size constraint.
//
// This value is then constrained in the interval [`Mesh.MeshSizeMin',
// `Mesh.MeshSizeMax'] and multiplied by `Mesh.MeshSizeFactor'. In addition,
// boundary mesh sizes are interpolated inside surfaces and/or volumes depending
// on the value of `Mesh.MeshSizeExtendFromBoundary' (which is set by default).
//
// When the element size is fully specified by a background mesh size field (as
// it is in this example), it is thus often desirable to set

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// This will prevent over-refinement due to small mesh sizes on the boundary.

// Finally, while the default "Frontal-Delaunay" 2D meshing algorithm
// (Mesh.Algorithm = 6) usually leads to the highest quality meshes, the
// "Delaunay" algorithm (Mesh.Algorithm = 5) will handle complex mesh size
// fields better - in particular size fields with large element size gradients:

Mesh.Algorithm = 5;

surfaceVector = Extrude {0, 0, zmax} {
  Surface{6};
  Layers{zlevels};
  Recombine;
};

//-------------------------------------------------------------------------------
//Boundary tagging
//-------------------------------------------------------------------------------
// At this level, Gmsh knows everything to display the rectangular surface 1 and
// to mesh it. An optional step is needed if we want to group elementary
// geometrical entities into more meaningful groups, e.g. to define some
// mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
// material ("steel", "carbon") properties.
//
// Such groups are called "Physical Groups" in Gmsh. By default, if physical
// groups are defined, Gmsh will export in output files only mesh elements that
// belong to at least one physical group. (To force Gmsh to save all elements,
// whether they belong to physical groups or not, set `Mesh.SaveAll=1;', or
// specify `-save_all' on the command line.) Physical groups are also identified
// by tags, i.e. strictly positive integers, that should be unique per dimension
// (0D, 1D, 2D or 3D). Physical groups can also be given names.
//
// Here we define a physical curve that groups the left, bottom and right curves
// in a single group (with prescribed tag 5); and a physical surface with name
// "My surface" (with an automatic tag) containing the geometrical surface 1:
//
Physical Volume("internal", 1) = 1;
Physical Surface("left", 2)    = {27};
Physical Surface("right", 3)   = {19};
Physical Surface("bottom", 4)  = {6};
Physical Surface("top", 5)     = {28};
Physical Surface("front", 6)   = {15};
Physical Surface("back", 7)    = {23};
//+
Show "*";
//+
Show "*";
