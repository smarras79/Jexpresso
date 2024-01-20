pitch = 10000.0;
circle_r = pitch/2;
bottom_left_circle = pitch/2;

Point(1) = {0, 0, 0};
Point(2) = {pitch, 0, 0};
Point(3) = {pitch, pitch, 0};
Point(4) = {0, pitch, 0};
//
// Outer Box domain
//
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Points that make up the circle
Point(5) = {bottom_left_circle,bottom_left_circle,0.0};
Point(6) = {bottom_left_circle-circle_r*Cos(Pi/4.)/2.,bottom_left_circle+circle_r*Cos(Pi/4.)/2.,0.0};
Point(7) = {bottom_left_circle-circle_r*Cos(Pi/4.)/2.,bottom_left_circle-circle_r*Cos(Pi/4.)/2.,0.0};
Point(8) = {bottom_left_circle+circle_r*Cos(Pi/4.)/2.,bottom_left_circle-circle_r*Cos(Pi/4.)/2.,0.0};
Point(9) = {bottom_left_circle+circle_r*Cos(Pi/4.)/2.,bottom_left_circle+circle_r*Cos(Pi/4.)/2.,0.0};

// Curves that connect points defining circle
Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};
// specify closed loop to make surface
Curve Loop(1) = {-5,-8,-7,-6};
Surface(1) = {1};

// break area outside circle into four surfaces
// First define lines, the line loops, then the
// surfaces
Line(9) = {1,7};
Line(10) = {8,2};
Line(11) = {9,3};
Line(12) = {6,4};

Curve Loop(3) = {9,6,10,-1};
Curve Loop(4) = {-10,7,11,-2};
Curve Loop(5) = {-11,8,12,-3};
Curve Loop(6) = {-12,5,-9,-4};

Surface(2) = {3};
Surface(3) = {4};
Surface(4) = {5};
Surface(5) = {6};
Coherence;

// bottom left cell
Transfinite Line {1, 6} = 10 Using Progression 1;
Transfinite Line {9, 10} = 10 Using Progression 1;

Transfinite Line {2, 7} = 10 Using Progression 1;
Transfinite Line {10, 11} = 10 Using Progression 1;

Transfinite Line {3, 8} = 10 Using Progression 1;
Transfinite Line {11, 12} = 10 Using Progression 1;


Transfinite Line {4, 5} = 10 Using Progression 1;
Transfinite Line {12, 9} = 10 Using Progression 1;


Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {5};

Recombine Surface "*";


Physical Line("free_slip") = {1, 2, 3, 4};
Physical Surface("domain") = {1, 2, 3, 4, 5};
Physical Point("boundary",   1) = {1, 2, 3, 4};