nelemx = 2;
nelemy =  2;
nelemz = 3;

xmin =      0;
xmax = 15000; //200000;
ymin =      0;
ymax =  2200;
zmin =      0;
zmax =  2400;
gridsize = xmax / nelemx;

gridsize_bottom = 1000;
gridsize_top    = 1750;

Point(1) = {xmin, ymin, zmin, gridsize_bottom};
Point(2) = {xmax, ymin, zmin, gridsize_bottom};
Point(3) = {xmax, ymin, zmax, gridsize_top};
Point(4) = {xmin, ymin, zmax, gridsize_top};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


npx = nelemx + 1;
npy = nelemy + 1;
npz = nelemz + 1;

//Horizontal sides
//Transfinite Line {1, 3} = npx; //Ceil((xmax-xmin)/gridsize) Using Progression 1;Line(4) = {4, 1};
//Vertical sides
//Transfinite Line {4, -2} = npz; // Ceil(R/gridsize) Using Progression 1.1;


Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {11};

//Transfinite Surface {12};
Recombine Surface {12};

surfaceVector = Extrude {0,(ymax-ymin),0} {
  Surface{12};
  Layers{1};
  Recombine;
};
//Coherence;

//    Physical Volume("internal") = 1;
/*    Physical Surface("left") = {21};
    Physical Surface("right") = {29};
    Physical Surface("bottom") = {25};
    Physical Surface("top") = {33};
    Physical Surface("front") = {12};
    Physical Surface("back") = {34};
    */
  //+
Show "*";
//+
Show "*";