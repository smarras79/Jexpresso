// Gmsh project created on Fri Oct  1 15:50:24 2021
nelx  = 1;
nely  = 1;
nelz  = 1;

xmin =  0;
xmax = 4000;
ymin =  0;
ymax = 1000;
zmin =  0;
zmax = 2000;

lc1 = 1000;

Point(1) = {xmin, 0, zmin, lc1};
Point(2) = {xmax, 0, zmin, lc1};
Point(3) = {xmax, 0, zmax, lc1};
Point(4) = {xmin, 0, zmax, lc1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

npx = nelx + 1;
npz = nelz + 1;

Transfinite Line {1, -3} = npx Using Progression 1;
Transfinite Line {2, -4} = npz Using Progression 1;

Transfinite Surface {1};
Recombine Surface {1};

Extrude {0, ymax, 0} {
  Surface{1};
  Layers{nely};
  Recombine;
}

Physical Volume("volume") = {1};
Physical Surface("front") = {1};
Physical Surface("top") = {21};
Physical Surface("right") = {17};
Physical Surface("left") = {25};
Physical Surface("bottom") = {13};
Physical Surface("back") = {26};
