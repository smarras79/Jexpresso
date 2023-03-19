// -----------------------------------------------------------------------------
//
//  Gmsh GEO tutorial 1
//
//  Geometry basics, elementary entities, physical groups
//
// -----------------------------------------------------------------------------

// The simplest construction in Gmsh's scripting language is the
// `affectation'. The following command defines a new variable `lc':

lc = 1e-2;

gridsize_bottom = 0.2;
gridsize_top    = 0.2;

// This variable can then be used in the definition of Gmsh's simplest
// `elementary entity', a `Point'. A Point is uniquely identified by a tag (a
// strictly positive integer; here `1') and defined by a list of four numbers:
// three coordinates (X, Y and Z) and the target mesh size (lc) close to the
// point:

Point(1) = {-1, -1, 0, gridsize_bottom};
Point(2) = { 1, -1, 0, gridsize_bottom};
Point(3) = { 1,  1, 0, gridsize_top};
Point(4) = {-1,  1, 0, gridsize_top};

// Curves are Gmsh's second type of elementary entities, and, amongst curves,
// straight lines are the simplest. A straight line is identified by a tag and
// is defined by a list of two point tags. In the commands below, for example,
// the line 1 starts at point 1 and ends at point 2.
//
// Note that curve tags are separate from point tags - hence we can reuse tag
// `1' for our first curve. And as a general rule, elementary entity tags in
// Gmsh have to be unique per geometrical dimension.

Line(1) = {1, 2}; //bottom
Line(2) = {3, 2}; //right
Line(3) = {3, 4}; //top
Line(4) = {4, 1}; //left

// The third elementary entity is the surface. In order to define a simple
// rectangular surface from the four curves defined above, a curve loop has
// first to be defined. A curve loop is also identified by a tag (unique amongst
// curve loops) and defined by an ordered list of connected curves, a sign being
// associated with each curve (depending on the orientation of the curve to form
// a loop):
Curve Loop(1) = {4, 1, -2, 3};

// We can then define the surface as a list of curve loops (only one here,
// representing the external contour, since there are no holes--see `t4.geo' for
// an example of a surface with a hole):
Plane Surface(1) = {1};
Recombine Surface {1}; 

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
Physical Point("boundary",   1) = {1, 2, 3, 4};
Physical Curve("inflow",     2) = {4};
Physical Curve("outflow",    3) = {2};
Physical Curve("free_slip", 4) = {3};
Physical Curve("no_slip",   5) = {1};
Physical Surface("domain") = {1};

//
// You can save the mesh in older
// versions of the MSH format:
//
// - In the GUI: open `File->Export', enter your `filename.msh' and then pick
//   the version in the dropdown menu.
// - On the command line: use the `-format' option (e.g. `gmsh file.geo -format
//   msh2 -2').
// - In a `.geo' script: add `Mesh.MshFileVersion = x.y;' for any version
//   number `x.y'.
// - As an alternative method, you can also not specify the format explicitly,
//   and just choose a filename with the `.msh2' or `.msh4' extension.

// Note that starting with Gmsh 3.0, models can be built using other geometry
// kernels than the default built-in kernel. By specifying
//
//   SetFactory("OpenCASCADE");
//
// any subsequent command in the `.geo' file would be handled by the OpenCASCADE
// geometry kernel instead of the built-in kernel. Different geometry kernels
// have different features. With OpenCASCADE, instead of defining the surface by
// successively defining 4 points, 4 curves and 1 curve loop, one can define the
// rectangular surface directly with
//
//   Rectangle(2) = {.2, 0, 0, .1, .3};
//
// The underlying curves and points could be accessed with the `Boundary' or
// `CombinedBoundary' operators.
//
// See e.g. `t16.geo', `t18.geo', `t19.geo' or `t20.geo' for complete examples
// based on OpenCASCADE, and `demos/boolean' for more.
