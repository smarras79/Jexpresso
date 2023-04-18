D = 1.0;
R = D/2.0;
H = 4;
Ll = 4;
Lr = 6;

Lil = 2;
Lir = 2;
Hi = 2;

sz = 0.1;

Point(1) = {R, 0, 0, sz};
Point(2) = {R, R, 0, sz};
Point(3) = {-R, R, 0, sz};
Point(4) = {-R, 0, 0, sz};

Point(5) = {-Ll, 0, 0, sz};
Point(6) = {-Ll, H, 0, sz};
Point(7) = {Lr, H, 0, sz};
Point(8) = {Lr, 0, 0, sz};

Point(9)  = {Lir, 0, 0, sz};
Point(10) = {Lir, Hi, 0, sz};
Point(11) = {-Lil, Hi, 0, sz};
Point(12) = {-Lil, 0, 0, sz};

Point(13) = {Lir, R, 0, sz};
Point(14) = {R, Hi, 0, sz};
Point(15) = {-R, Hi, 0, sz};
Point(16) = {-Lil, R, 0, sz};

Point(17) = {Lir, H, 0, sz};
Point(18) = {-Lil, H, 0, sz};
Point(19) = {Lr, Hi, 0, sz};
Point(20) = {-Ll, Hi, 0, sz};
Point(21) = {R,H,0,sz};
Point(22) = {-R,H,0,sz};
Point(23) = {R,H,0,sz};
Point(24) = {Lr,R,0,sz};
Point(25) = {-Ll,R,0,sz};



//Line Loop(22) = {5, 6, 18, -11};
//Ruled Surface(23) = {22};
//Line Loop(24) = {18, 12, 13, -19};
//Ruled Surface(25) = {24};
//Line Loop(26) = {20, -14, -19, 7};
//Ruled Surface(27) = {26};
//Line Loop(28) = {16, -21, 20, 15};
//Ruled Surface(29) = {28};
//Line Loop(30) = {21, 17, -9, -8};
//Ruled Surface(31) = {30};



//Transfinite Line {8, 6, 17, 11} = RHO1 Using Progression 1;
//Transfinite Line {-16, 20, 19, 12} = RHO2 Using Progression 1.03;
//Transfinite Line {9, 21, 15} = RHO3 Using Progression 1.03;
//Transfinite Line {7, 14} = RHO4 Using Progression 1.03;
//Transfinite Line {-5, 18, -13} = RHO5 Using Progression 1.03;

//Transfinite Surface {23,25,27,29,31};
//Recombine Surface "*";

//Field[1] = Attractor;
//Field[1].NodesList = {6};
//Field[1].NNodesByEdge = 100;
//Field[1].EdgesList = {11,12,13,14,15};
//Background Field = 1;

//Extrude {0, 0, R} {
//  Surface {23,25,27,29,31,33};
//  Layers{1};
//  Recombine;
//}

Line(1) = {5, 25};
Line(2) = {12, 16};
Line(3) = {4, 3};
Line(4) = {1, 2};
Line(5) = {9, 13};
Line(6) = {8, 24};
Line(7) = {25, 20};
Line(8) = {16, 11};
Line(9) = {3, 15};
Line(10) = {2, 14};
Line(11) = {13, 10};
Line(12) = {24, 19};
Line(13) = {19, 7};
Line(14) = {10, 17};
Line(15) = {14, 21};
Line(16) = {15, 22};
Line(17) = {11, 18};
Line(18) = {20, 6};
Line(19) = {5, 12};
Line(20) = {25, 16};
Line(21) = {20, 11};
Line(22) = {6, 18};
Line(23) = {12, 4};
Line(24) = {16, 3};
Line(25) = {11, 15};
Line(26) = {18, 22};
Line(27) = {3, 2};
Line(28) = {15, 14};
Line(29) = {22, 21};
Line(30) = {1, 9};
Line(31) = {2, 13};
Line(32) = {14, 10};
Line(33) = {21, 17};
Line(34) = {9, 8};
Line(35) = {13, 24};
Line(36) = {10, 19};
Line(37) = {17, 7};

RHO1 = 20; 
RHO2 = 40;
RHO3 = 40;
RHO4 = 30; 
RHO5 = 30;

Line Loop(38) = {18, 22, -17, -21};
Ruled Surface(39) = {38};
Line Loop(40) = {17, 26, -16, -25};
Ruled Surface(41) = {40};
Line Loop(42) = {29, -15, -28, 16};
Ruled Surface(43) = {42};
Line Loop(44) = {33, -14, -32, 15};
Ruled Surface(45) = {44};
Line Loop(46) = {37, -13, -36, 14};
Ruled Surface(47) = {46};
Line Loop(48) = {7, 21, -8, -20};
Ruled Surface(49) = {48};
Line Loop(50) = {8, 25, -9, -24};
Ruled Surface(51) = {50};
Line Loop(52) = {9, 28, -10, -27};
Ruled Surface(53) = {52};
Line Loop(54) = {10, 32, -11, -31};
Ruled Surface(55) = {54};
Line Loop(56) = {11, 36, -12, -35};
Ruled Surface(57) = {56};
Line Loop(58) = {1, 20, -2, -19};
Ruled Surface(59) = {58};
Line Loop(60) = {2, 24, -3, -23};
Ruled Surface(61) = {60};
Line Loop(62) = {4, 31, -5, -30};
Ruled Surface(63) = {62};
Line Loop(64) = {5, 35, -6, -34};
Ruled Surface(65) = {64};
Transfinite Line {27, 28, 29} = 30 Using Progression 1;
Transfinite Line {-23, -24, -25,-26} = 40 Using Progression 1.02;
Transfinite Line {-19, -20, -21,-22} = 30 Using Progression 1.02;
Transfinite Line {30, 31, 32,33} = 30 Using Progression 1.05;
Transfinite Line {34, 35, 36,37} = 30 Using Progression 1;
Transfinite Line {1,2,3,4,5,6} = 20 Using Progression 1;
Transfinite Line {7,8,9,10,11,12} = 30 Using Progression 1.05;
Transfinite Line {13,14,15,16,17,18} = 20 Using Progression 1.02;

Transfinite Surface "*";
Recombine Surface "*";

surfaceVector[] = Extrude {0, 0, R} {
  Surface{39, 49, 59, 41, 51, 61, 43, 53, 45, 55, 63, 57, 47, 65};
  Layers{1};
  Recombine;
};


//Physical Surface("front") = surfaceVector[0];
//Physical Volume("internal") = surfaceVector[1];
//Physical Volume("") = surfaceVector[1];

//Printf("back 0 =  %g", surfaceVector[0]);
//Printf("back 0 =  %g", surfaceVector[1]);
//Printf("size =  %g", #surfaceVector[]);

//Mesh 3;