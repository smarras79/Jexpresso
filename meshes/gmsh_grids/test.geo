nelemx = 20;
nelemy = 1;
nelemz = 12;

xmin =  -30000;
xmax =   30000;
ymin =  -6000;
ymax =   6000;
zmin =      0;
zmax =  20000;
gridsize = (xmax-xmin) / nelemx;

Point(1) = {xmin, ymin, zmin, gridsize};
Point(2) = {xmax, ymin, zmin, gridsize};
Point(3) = {xmax, ymin, zmax, gridsize};
Point(4) = {xmin, ymin, zmax, gridsize};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Number of points along each axis
npx = nelemx + 1;
npy = nelemy + 1;
npz = nelemz + 1;

// Define the splitting point for the vertical axis
z_half = zmin + (zmax - zmin) / 2;

// Define number of divisions in the lower and upper halves
npz_half = npz / 2;
npz_lower = Ceil(npz_half); // Uniform lower half
npz_upper = Ceil(npz_half); // Stretched upper half

// Uniform distribution in the lower half
Transfinite Line {4} = npz_lower Using Progression 1.0;

// Stretched distribution in the upper half
Transfinite Line {2} = npz_upper Using Progression 2.0; // Adjust Progression factor as needed

// Horizontal sides
Transfinite Line {1, 3} = npx; 

Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {11};

Transfinite Surface {12};
Recombine Surface {12};

surfaceVector = Extrude {0, (ymax - ymin), 0} {
  Surface{12};
  Layers{nelemy};
  Recombine;
};

// Physical groups for surfaces and volumes
Physical Surface("bottom") = {12};
Physical Volume("internal") = {1};
Physical Surface("back") = {25};
Physical Surface("front") = {33};
Physical Surface("periodicx") = {21, 29};
Physical Surface("top") = {34}; // from Plane Surface (6)

Show "*";
Show "*";
