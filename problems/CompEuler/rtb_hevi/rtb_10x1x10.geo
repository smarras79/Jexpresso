nelemx = 10;
nelemy =  1;
nelemz = 50;

xmin = -5000.0;
xmax =  5000.0;
ymin =     0.0;
ymax =  5000.0;
zmin =     0.0;
zmax = 10000.0;

gridsize = (xmax-xmin) / nelemx;

Point(1) = {xmin, ymin, zmin, gridsize};
Point(2) = {xmax, ymin, zmin, gridsize};
Point(3) = {xmax, ymin, zmax, gridsize};
Point(4) = {xmin, ymin, zmax, gridsize};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

npx = nelemx + 1;
npz = nelemz + 1;

// Horizontal sides
Transfinite Line {1, 3} = npx;
// Vertical sides
Transfinite Line {4, -2} = npz Using Progression 1.0;

Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {11};

Transfinite Surface {12};
Recombine Surface {12};

surfaceVector = Extrude {0,(ymax-ymin),0} {
  Surface{12};
  Layers{nelemy};
  Recombine;
};

    Physical Surface("front")  = {12};
    Physical Surface("back")   = {34};
    Physical Volume("internal") = {1};
    Physical Surface("top")    = {33};
    Physical Surface("bottom") = {25};
    Physical Surface("left")   = {21};
    Physical Surface("right")  = {29};
