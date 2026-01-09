nelemx = 32;
nelemy = 2;
nelemz = 15; //19;

h = 1.0;
xmin = 0.0;
xmax = 6.28;
ymin = 0.0;
ymax = 1.57;
zmin = -1.0;
zmax =  1.0;
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

// Calculate the progression ratio to achieve first/last points at ±0.994
// For domain [0,2] with 32 points, we want:
// - First interior point at z ≈ 0.994 (from bottom)
// - Last interior point at z ≈ 1.006 (from bottom), or equivalently z ≈ -0.994 from top
// This requires a geometric progression with ratio r such that:
// The first element length ≈ 0.994 and the sum of all elements = 2.0

// For your specific case, we need a progression ratio of approximately 1.12
progression_ratio = 0.15;

//Horizontal sides
Transfinite Line {1, 3} = npx; //Ceil((xmax-xmin)/gridsize) Using Progression 1;
//Vertical sides - modified to use geometric progression
Transfinite Curve {4, -2} = npz Using Bump progression_ratio;

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
    Physical Surface("wall_model_top") = {33};
    Physical Surface("wall_model_bottom") = {25};
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