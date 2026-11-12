// ============================================================
// 2.5D cross-section of the domain and grid of
//
//   J. E. Leake, M. G. Linton & S. K. Antiochos, "The role of reconnection in
//   the onset of solar eruptions", ApJ 934:10 (2022), Sec. 2.3,
//
// the x-z plane of problems/MHD/fluxEmergenceLeake2022/LK22.geo (units of
// L0 = 170 km): the paper's stretched cells,
//
//   d = d0 [1 + f (1 + tanh((|s0| - l)/w))],   s0 = the cell's outer face,
//   x: d0 = 0.5,  f = 6, l = 256/2.1, w = 256/10   (512 cells, x in [-277.2, 277.2])
//   z: d0 = 0.25, f = 2, l = 70,      w = 192/20   (768 cells, z in [-30, 530.5]),
//
// with the ELEMENT boundaries every kx-th face in x and every kz-th in z, so that
// with :nop => 4 the nodal spacing follows the paper's cells. The code's 2D coordinates are (x, y) with y
// the paper's z (up); the rope axis, the paper's y, is the invariant direction.
//
//   ref = 1  ->  32 x 48 elements   (k = 16)  LK22_2D_32x48.msh     (the 3D LK22_REF = 2 cross-section)
//   ref = 2  ->  64 x 96 elements   (k = 8)   LK22_2D_64x96.msh     (1/2 of the paper's resolution)
//   ref = 3  -> 128 x 192 elements  (k = 4)   LK22_2D_128x192.msh   (the paper's 512 x 768 cells at nop 4)
//   ref = 4  ->  32 x 96 elements   (kx = 16, kz = 8)  LK22_2D_32x96.msh  (ref 1 in x, ref 2 in z:
//                                    elements 2 tall at the bottom; README, "Vertical resolution";
//                                    the 3D LK22_REF = 5 cross-section)
//
// Boundary tags (the strings user_bc_dirichlet! receives):
//   "bottom" y = -30,  "top" y = Z,  "left" x = -X,  "right" x = X.
// In 2D the surface is named "domain".
//
// Generate with:
//   gmsh -2 LK22_2D.geo -setnumber ref 1 -o LK22_2D_32x48.msh
//   gmsh -2 LK22_2D.geo -setnumber ref 2 -o LK22_2D_64x96.msh
//   gmsh -2 LK22_2D.geo -setnumber ref 3 -o LK22_2D_128x192.msh
//   gmsh -2 LK22_2D.geo -setnumber ref 4 -o LK22_2D_32x96.msh
// ============================================================

DefineConstant[ ref = 1 ];
kx = 16; kz = 16;
If (ref == 2)
  kx = 8;  kz = 8;
EndIf
If (ref == 3)
  kx = 4;  kz = 4;
EndIf
If (ref == 4)
  kx = 16; kz = 8;
EndIf

// horizontal faces on [0, X]: 256 paper cells per half-width
lx = 256/2.1; wx = 256/10;
acc = 0; nh = 0;
xh[0] = 0;
For i In {1:256}
  s0 = 0.5*i;
  acc = acc + 0.5*(1 + 6*(1 + Tanh((s0 - lx)/wx)));
  If (Fmod(i, kx) == 0)
    nh = nh + 1;
    xh[nh] = acc;
  EndIf
EndFor
X = xh[nh];
// cumulative fractions of the full span [-X, X]: -xh[nh-1], ..., 0, ..., xh[nh]
For m In {1:2*nh}
  If (m <= nh)
    fx[m-1] = (X - xh[nh - m])/(2*X);
  Else
    fx[m-1] = (X + xh[m - nh])/(2*X);
  EndIf
  ox[m-1] = 1;
EndFor

// vertical faces from z = -30: 768 paper cells
acc = 0; nv = 0;
For i In {1:768}
  s0 = -30 + 0.25*i;
  acc = acc + 0.25*(1 + 2*(1 + Tanh((s0 - 70)/9.6)));
  If (Fmod(i, kz) == 0)
    nv = nv + 1;
    zv[nv] = acc;
  EndIf
EndFor
H = zv[nv];
Z = -30 + H;
For m In {1:nv}
  fz[m-1] = zv[m]/H;
  oz[m-1] = 1;
EndFor

Point(1) = {-X, -30, 0};
ex[] = Extrude {2*X, 0, 0} { Point{1}; Layers{ox[], fx[]}; };
ey[] = Extrude {0, H, 0} { Curve{ex[1]}; Layers{oz[], fz[]}; Recombine; };

e = 1e-6;
Physical Surface("domain") = {ey[1]};
Physical Curve("bottom") = Curve In BoundingBox{-X-e, -30-e, -e, X+e, -30+e, e};
Physical Curve("top")    = Curve In BoundingBox{-X-e, Z-e,   -e, X+e, Z+e,   e};
Physical Curve("left")   = Curve In BoundingBox{-X-e, -30-e, -e, -X+e, Z+e,  e};
Physical Curve("right")  = Curve In BoundingBox{X-e,  -30-e, -e, X+e,  Z+e,  e};

Mesh.ElementOrder = 1;
