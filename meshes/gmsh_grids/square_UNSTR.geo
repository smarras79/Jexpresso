nelemx = 20;
nelemy = 20;
nelemz = 1;

xmin = -5000.0;
xmax =  5000.0;
ymin =     0.0;
ymax = 10000.0;


// 
lc = (xmax - xmin)/nelemx;
//
Point(1) = {xmin, ymin, lc};
Point(2) = {xmax, ymin, lc};
Point(3) = {xmax, ymax, lc};
Point(4) = {xmin, ymax, lc};
//
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//
npx = nelemx + 1;
npy = nelemy + 1;

//Horizontal sides
Transfinite Line {1, 3} = npx; //Ceil((xmax-xmin)/gridsize) Using Progression 1;
//Vertical sides
Transfinite Line {4, -2} = npy Using Progression 1.0;
//
Curve Loop(1) = {4, 1, 2, 3};
//
Plane Surface(1) = {1};
Recombine Surface{1};



//-------------------------------------------------------------------------------
//Boundary tagging
//-------------------------------------------------------------------------------
Physical Point("boundary",   1) = {1, 2, 3, 4};
Physical Curve("free_slip",  2) = {1, 3, 2, 4};
Physical Surface("domain") = {1};
