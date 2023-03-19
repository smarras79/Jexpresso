xmin_sq = -1;
ymin_sq = -1;
xmax_sq =  1;
ymax_sq =  1;

rxmin = -2;
rymin = -2;
rxmax =  2;
rymax =  2;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {xmax_sq, ymax_sq, 0, 1.0};
Point(3) = {rxmax, rymax, 0, 1.0};
Point(4) = {xmin_sq, ymax_sq, 0, 1.0};
Point(5) = {rxmin, rymax, 0, 1.0};
Point(6) = {xmin_sq, ymin_sq, 0, 1.0};
Point(7) = {rxmin, rymin, 0, 1.0};
Point(8) = {xmax_sq, ymin_sq, 0, 1.0};
Point(9) = {rxmax, rymin, 0, 1.0};

Line(1) = {2, 3};
Circle(2) = {3, 1, 5};
Line(3) = {5, 4};
Line(4) = {4, 2};
Line(5) = {6, 7};
Circle(6) = {7, 1, 9};
Line(7) = {9, 8};
Line(8) = {8, 6};
Circle(9) = {5, 1, 7};
Circle(10) = {9,1, 3};
Line(11) = {4, 6};
Line(12) = {8, 2};

// Lines of the square and arcs of the circle
Transfinite Line {2, 4, 8, 6} = 60;
Transfinite Line {9, 11, 12, 10} = 60;

// Lines connecting square and circle
Transfinite Line {1, 3, 5, 7} = 60;

Transfinite Surface(1) = {};
Recombine Surface "*";
