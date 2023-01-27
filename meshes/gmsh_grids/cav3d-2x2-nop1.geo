// ----------- BEGIN INPUT ------------

Mesh.ElementOrder = 1;
// ----------- END INPUT --------------

Point(1) = {0,  0, 0, 0.5};
Point(2) = {1,  0, 0, 0.5};
Point(3) = {1,  0, 1, 0.5};
Point(4) = {0,  0, 1, 0.5};
Point(5) = {0,  1, 0, 0.5};
Point(6) = {1,  1, 0, 0.5};
Point(7) = {1,  1, 1, 0.5};
Point(8) = {0,  1, 1, 0.5};

Line(1) = {5, 1};
Line(2) = {1, 4};
Line(3) = {4, 8};
Line(4) = {8, 5};
Line(5) = {1, 2};
Line(6) = {2, 6};
Line(7) = {6, 5};
Line(8) = {8, 7};
Line(9) = {7, 6};
Line(10) = {4, 3};
Line(11) = {3, 2};
Line(12) = {3, 7};

Line Loop(1) = {5, -11, -10, -2};
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface {1};

Line Loop(2) = {1, 5, 6, 7};
Plane Surface(2) = {2};
Transfinite Surface {2};
Recombine Surface {2};

Line Loop(3) = {1, 2, 3, 4};
Plane Surface(3) = {3};
Transfinite Surface {3};
Recombine Surface {3};

Line Loop(4) = {6, -9, -12, 11};
Plane Surface(4) = {4};
Transfinite Surface {4};
Recombine Surface {4};

Line Loop(5) = {7, -4, 8, 9};
Plane Surface(5) = {5};
Transfinite Surface {5};
Recombine Surface {5};

Line Loop(6) = {8, -12, -10, 3};
Plane Surface(6) = {6};
Transfinite Surface {6};
Recombine Surface {6};

Surface Loop(1) = {5, 2, 3, 1, 4, 6};
Volume(1) = {1};
Transfinite Volume {1};

Recombine Surface {11};
// Real entities
Physical Surface("top") = {5};
Physical Surface("wall") = {2, 3, 1, 6, 4};
Physical Volume  ("internalVolume") = {1};
Physical Line("edges_boundary") = {7, 4, 1, 2, 5, 6, 9, 8, 12, 11, 10, 3};
