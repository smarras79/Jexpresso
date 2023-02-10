nelemx = 20;
nelemy = 1;
nelemz = 1;

xmin = -1; //-nelemx;
xmax =	1; //nelemx;
ymin =  0; //-nelemy;
ymax =  0.1; //nelemy;
gridsize = (xmax-xmin) / nelemx;

Point(1) = {xmin, ymin, gridsize};
Point(2) = {xmax, ymin, gridsize};
Point(3) = {xmax, ymax, gridsize};
Point(4) = {xmin, ymax, gridsize};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

npx = nelemx + 1;
npy = nelemy + 1;

//Horizontal sides
Transfinite Line {1, 3} = npx; //Ceil((xmax-xmin)/gridsize) Using Progression 1;
//Vertical sides
Transfinite Line {4, -2} = npy Using Progression 1.0;


Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {11};

Transfinite Surface {12};
Recombine Surface {12};

/*surfaceVector = Extrude {0,(ymax-ymin),0} {
  Surface{12};
  Layers{nelemy};
  Recombine;
};*/
//Coherence;

  /* surfaceVector contains in the following order:
     [0] - front surface (opposed to source surface)
     [1] - extruded volume
     [2] - bottom surface (belonging to 1st line in "Line Loop (6)")
     [3] - right surface (belonging to 2nd line in "Line Loop (6)")
     [4] - top surface (belonging to 3rd line in "Line Loop (6)")
     [5] - left surface (belonging to 4th line in "Line Loop (6)")
    */
    Physical Surface("front") = {12};
    Physical Volume("internal") = {1};
    Physical Surface("bottom") = {25};
    Physical Surface("top") = {33};
    Physical Surface("left") = {21};
    Physical Surface("right") = {29};
    Physical Surface("back") = {34}; // from Plane Surface (6) ...
  //+
Show "*";
//+
Show "*";