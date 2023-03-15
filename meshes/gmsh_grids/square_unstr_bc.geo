// 
lc = 0.25;
//
Point(1) = {-1, -1, -1, lc};
Point(2) = {-1, 1, -1, lc};
Point(3) = {1, -1, 1, lc};
Point(4) = {1, 1, 1, lc};
//
Line(1) = {1, 2};
Line(2) = {4, 3};
Line(3) = {3, 1};
Line(4) = {2, 4};
//
Curve Loop(1) = {1, 4, 2, 3};
//
Plane Surface(1) = {1};
Recombine Surface{1};
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
Physical Point("solidPoints", 1) = {1, 2, 3, 4};
Physical Curve("boundary", 2) = {1, 2, 3};
Physical Curve("top",      3) = {4};
Physical Surface("domain") = {1};