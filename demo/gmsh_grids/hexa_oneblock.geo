// Gmsh project created on Fri Oct  1 15:50:24 2021
nelx  = 10;
nely  = 10;
nelz  = 1;

xmin =  -1;
xmax =   1;
ymin =  -1;
ymax =   1;
zmin =   0;
zmax =  0.1;

lc1 = 0.1;

Point(1) = {xmin, ymin, 0, lc1};
Point(2) = {xmax, ymin, 0, lc1};
Point(3) = {xmax, ymax, 0, lc1};
Point(4) = {xmin, ymax, 0, lc1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

npx = nelx + 1;
npy = nely + 1;

Transfinite Line {1, -3} = npx Using Progression 1;
Transfinite Line {2, -4} = npy Using Progression 1;

Transfinite Surface {1};
Recombine Surface {1};

Extrude {0, 0, zmax} {
  Surface{1};
  Layers{nelz};
  Recombine;
}

Physical Volume("volume") = {1};
Physical Surface("front") = {1};
Physical Surface("top") = {21};
Physical Surface("right") = {17};
Physical Surface("left") = {25};
Physical Surface("bottom") = {13};
Physical Surface("back") = {26};
