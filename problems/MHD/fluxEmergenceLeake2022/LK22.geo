// ============================================================
// Domain and grid of the 3D flux-emergence parameter study of
//
//   J. E. Leake, M. G. Linton & S. K. Antiochos, "The role of reconnection in
//   the onset of solar eruptions", ApJ 934:10 (2022), Sec. 2.3,
//
// in units of L0 = 170 km. The paper starts from a uniform grid
//
//   x0, y0 in [-128, 128], z0 in [-30, 162],  512 x 512 x 768 cells,
//   (dx0, dy0, dz0) = (0.5, 0.5, 0.25),
//
// and gives every cell the size (printed as d = d0 + f tanh(...), Eqs. 19-20,
// which is negative near the centre; the form below reproduces every number
// the paper quotes: 0.085 -> 0.714 Mm horizontally, 0.0425 -> 0.2125 Mm
// vertically, X x Y x Z = [-47.1, 47.1]^2 x [-5.1, 90.2] Mm)
//
//   d = d0 [1 + f (1 + tanh((|s0| - l)/w))],   s0 = the cell's outer face,
//   x, y: d0 = 0.5,  f = 6, l = 256/2.1, w = 256/10
//   z:    d0 = 0.25, f = 2, l = 70,      w = 192/20,
//
// which places the faces at x, y in [-277.2, 277.2] and z in [-30, 530.5].
// The ELEMENT boundaries are every k-th face of that grid (every kx-th in x
// and y, kz-th in z), so with :nop => 4 the nodal spacing follows the paper's
// cells:
//
//   ref = 1  ->  16 x 16 x 24 elements   (k = 32)  LK22_16x16x24.msh    smoke test
//   ref = 2  ->  32 x 32 x 48 elements   (k = 16)  LK22_32x32x48.msh    1/4 of the paper's resolution
//   ref = 3  ->  64 x 64 x 96 elements   (k = 8)   LK22_64x64x96.msh    1/2
//   ref = 4  -> 128 x 128 x 192 elements (k = 4)   LK22_128x128x192.msh the paper's 512^2 x 768 at nop 4
//   ref = 5  ->  32 x 32 x 96 elements   (kx = 16, kz = 8)  LK22_32x32x96.msh  ref 2 in x and y,
//                ref 3 in z: elements 2 tall at the bottom (README, "Vertical resolution")
//
// Boundary tags (the strings user_bc_dirichlet! receives):
//   "bottom" z = -30,  "top" z = Z,  "left" x = -X,  "right" x = X,
//   "front" y = -X,    "back" y = X.
// In 3D the volume must be named "internal".
//
// Generate with:
//   gmsh -3 LK22.geo -setnumber ref 1 -o LK22_16x16x24.msh
//   gmsh -3 LK22.geo -setnumber ref 2 -o LK22_32x32x48.msh
//   gmsh -3 LK22.geo -setnumber ref 3 -o LK22_64x64x96.msh
//   gmsh -3 LK22.geo -setnumber ref 4 -o LK22_128x128x192.msh
//   gmsh -3 LK22.geo -setnumber ref 5 -o LK22_32x32x96.msh
// ============================================================

DefineConstant[ ref = 1 ];
kx = 32; kz = 32;
If (ref == 2)
  kx = 16; kz = 16;
EndIf
If (ref == 3)
  kx = 8;  kz = 8;
EndIf
If (ref == 4)
  kx = 4;  kz = 4;
EndIf
If (ref == 5)
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

Point(1) = {-X, -X, -30};
ex[] = Extrude {2*X, 0, 0} { Point{1}; Layers{ox[], fx[]}; };
ey[] = Extrude {0, 2*X, 0} { Curve{ex[1]}; Layers{ox[], fx[]}; Recombine; };
ez[] = Extrude {0, 0, H} { Surface{ey[1]}; Layers{oz[], fz[]}; Recombine; };

e = 1e-6;
Physical Volume("internal") = {ez[1]};
Physical Surface("bottom") = Surface In BoundingBox{-X-e, -X-e, -30-e, X+e, X+e, -30+e};
Physical Surface("top")    = Surface In BoundingBox{-X-e, -X-e, Z-e,   X+e, X+e, Z+e};
Physical Surface("left")   = Surface In BoundingBox{-X-e, -X-e, -30-e, -X+e, X+e, Z+e};
Physical Surface("right")  = Surface In BoundingBox{X-e,  -X-e, -30-e, X+e,  X+e, Z+e};
Physical Surface("front")  = Surface In BoundingBox{-X-e, -X-e, -30-e, X+e, -X+e, Z+e};
Physical Surface("back")   = Surface In BoundingBox{-X-e, X-e,  -30-e, X+e,  X+e, Z+e};

Mesh.ElementOrder = 1;
