// Gmsh project created on Thu Nov 27 03:32:31 2025
// Gmsh project created on Thu Nov 27 03:29:24 2025
// ============================================================================
// Hill2D Flow - Matches WRF hill flow configuration
// ============================================================================

nelemx = 107;
nelemy = 50;
xmin = -80000;
xmax =  80000;
ymin =  0;
ymax =  15000;

gridsize = (xmax - xmin) / nelemx;

Point(1) = {xmin, ymin, 0, gridsize};
Point(2) = {xmax, ymin, 0, gridsize};
Point(3) = {xmax, ymax, 0, gridsize};
Point(4) = {xmin, ymax, 0, gridsize};

Line(1) = {1, 2};  // Bottom
Line(2) = {2, 3};  // Right
Line(3) = {3, 4};  // Top
Line(4) = {4, 1};  // Left

npx = nelemx + 1;
npy = nelemy + 1;

Transfinite Line {1, 3} = npx;
Transfinite Line {4, -2} = npy Using Progression 1.0;

Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {11};
Transfinite Surface {12};
Recombine Surface {12};

// ============================================================================
// Physical Groups - Matches WRF Configuration
// ============================================================================

Physical Curve("bottom", 1) = {1};     // y-start: terrain surface (symmetric_ys)
Physical Curve("top", 2) = {3};        // y-end: free-slip (symmetric_ye)
Physical Curve("inflow", 3) = {4};     // x-start: open inflow (open_xs)
Physical Curve("outflow", 4) = {2};    // x-end: open outflow (open_xe)
Physical Surface("domain", 5) = {12};