nelemx = 64;
nelemy = 64;
nelemz = 36;

//
// Nelx = Lx/(nop*Dx)
//
xmin =  0;
xmax =	10240.0;
ymin =  0;
ymax =  10240.0;
//ymax =  10240.0;
zmin =   0;
zmax =  3500.0;
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
npy = nelemy + 1;
npz = nelemz + 1;

//Horizontal sides
Transfinite Line {1, 3} = npx; //Ceil((xmax-xmin)/gridsize) Using Progression 1;
//Vertical sides
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
//Coherence;

  /* surfaceVector contains in the following order:
     [0] - front surface (opposed to source surface)
     [1] - extruded volume
     [2] - bottom surface (belonging to 1st line in "Line Loop (6)")
     [3] - right surface (belonging to 2nd line in "Line Loop (6)")
     [4] - top surface (belonging to 3rd line in "Line Loop (6)")
     [5] - left surface (belonging to 4th line in "Line Loop (6)")
    */
    Physical Surface("periodicy") = {12,34};
    Physical Volume("internal") = {1};
    Physical Surface("top_wall") = {33};
    Physical Surface("MOST") = {25};
    Physical Surface("periodicx") = {21,29};
    // from Plane Surface (6) ...
  //+
Show "*";
//+
Show "*";
//+
Show "*";
//+
Show "*";
//+
Show "*";
