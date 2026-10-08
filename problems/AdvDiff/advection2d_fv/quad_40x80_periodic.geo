// [-5,5]x[0,20], 40x80 quads; periodic via the "periodicx"/"periodicz" names, as ../advection2d_dg.
// julia --project=. -e 'using GridapGmsh: gmsh; gmsh.initialize(); gmsh.option.setNumber("Mesh.MshFileVersion", 4.1);
//   gmsh.open("problems/AdvDiff/advection2d_fv/quad_40x80_periodic.geo"); gmsh.model.mesh.generate(2);
//   gmsh.write("problems/AdvDiff/advection2d_fv/quad_40x80_periodic.msh"); gmsh.finalize()'
nx = 40;  ny = 80;
Point(1) = {-5,  0, 1};  Point(2) = {5,  0, 1};
Point(3) = { 5, 20, 1};  Point(4) = {-5, 20, 1};
Line(1) = {1, 2};  Line(2) = {2, 3};  Line(3) = {3, 4};  Line(4) = {4, 1};
Transfinite Curve{1, 3} = nx + 1;
Transfinite Curve{2, 4} = ny + 1;
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface{1};
Recombine Surface{1};
Physical Point("boundary")   = {1, 2, 3, 4};
Physical Curve("periodicx")  = {1, 3};
Physical Curve("periodicz")  = {2, 4};
Physical Surface("domain")   = {1};
