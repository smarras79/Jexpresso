// ============================================================
// Domain and grid of the 3D flux-emergence simulations of
//
//   J. E. Leake, M. G. Linton & S. K. Antiochos, ApJ 787:46 (2014)  (Paper II),
//   J. E. Leake, M. G. Linton & T. Torok,         ApJ 778:99 (2013) (Paper I, Sec. 2.3),
//
// in units of L0 = 170 km:
//
//   x, y in [-X, X],   z in [-30, Z],   X = 127.6, Z = 209.9
//
// Paper I stretches its 304 cells per direction with (Eqs. 7-8)
//
//   x, y = +-[(1+f_h) chi + f_h w ln(cosh((chi-L_h)/w)/cosh(-L_h/w))],  chi in [0, 100]
//   z    = -30 + (1+f_v) chi + f_v w ln(cosh((chi-L_v)/w)/cosh(-L_v/w)),  chi in [0, 130]
//
// with f_h = 2.1, f_v = 1.83, L_h = 95, L_v = 100, w = 10 (the paper's
// x_max = 126.85 and z_max = 210.45 come from its cell-centred variant; the
// face formula gives 127.6 and 209.9). Here the ELEMENT boundaries follow the
// same map at nh elements per half-width and nv elements in z, so with
// :nop => 4 the nodal spacing follows the paper's cells:
//
//   ref = 1  ->  16 x 16 x 16 elements (nh = 8,  nv = 16)  LK_16x16x16.msh  smoke test
//   ref = 2  ->  38 x 38 x 38 elements (nh = 19, nv = 38)  LK_38x38x38.msh  half the paper's resolution
//   ref = 3  ->  76 x 76 x 76 elements (nh = 38, nv = 76)  LK_76x76x76.msh  the paper's 304^3 at nop 4
//
// The paper's centre spacing 0.658 is 100/152, its bottom spacing 0.428 is
// 130/304: at nop 4 each element spans four of the paper's cells.
//
// Boundary tags (the strings user_bc_dirichlet! receives):
//   "bottom" z = -30,  "top" z = Z,  "left" x = -X,  "right" x = X,
//   "front" y = -X,    "back" y = X.
// In 3D the volume must be named "internal": Geom.jl and mesh.jl drop that name
// from the boundary tags (2D meshes use "domain" for the surface instead).
//
// Generate with:
//   gmsh -3 LK.geo -setnumber ref 1 -o LK_16x16x16.msh
//   gmsh -3 LK.geo -setnumber ref 2 -o LK_38x38x38.msh
//   gmsh -3 LK.geo -setnumber ref 3 -o LK_76x76x76.msh
// ============================================================

DefineConstant[ ref = 1 ];
If (ref == 1)
  nh = 8;  nv = 16;
EndIf
If (ref == 2)
  nh = 19; nv = 38;
EndIf
If (ref == 3)
  nh = 38; nv = 76;
EndIf

fh = 2.1; fv = 1.83; Lh = 95; Lv = 100; w = 10;

// horizontal element boundaries on [0, X]
For k In {0:nh}
  chi = k*100/nh;
  xh[k] = (1 + fh)*chi + fh*w*Log(Cosh((chi - Lh)/w)/Cosh(-Lh/w));
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

// vertical element boundaries on [-30, Z]
For k In {0:nv}
  chi = k*130/nv;
  zv[k] = -30 + (1 + fv)*chi + fv*w*Log(Cosh((chi - Lv)/w)/Cosh(-Lv/w));
EndFor
Z = zv[nv];
For m In {1:nv}
  fz[m-1] = (zv[m] + 30)/(Z + 30);
  oz[m-1] = 1;
EndFor

Point(1) = {-X, -X, -30};
ex[] = Extrude {2*X, 0, 0} { Point{1}; Layers{ox[], fx[]}; };
ey[] = Extrude {0, 2*X, 0} { Curve{ex[1]}; Layers{ox[], fx[]}; Recombine; };
ez[] = Extrude {0, 0, Z + 30} { Surface{ey[1]}; Layers{oz[], fz[]}; Recombine; };

e = 1e-6;
Physical Volume("internal") = {ez[1]};
Physical Surface("bottom") = Surface In BoundingBox{-X-e, -X-e, -30-e, X+e, X+e, -30+e};
Physical Surface("top")    = Surface In BoundingBox{-X-e, -X-e, Z-e,   X+e, X+e, Z+e};
Physical Surface("left")   = Surface In BoundingBox{-X-e, -X-e, -30-e, -X+e, X+e, Z+e};
Physical Surface("right")  = Surface In BoundingBox{X-e,  -X-e, -30-e, X+e,  X+e, Z+e};
Physical Surface("front")  = Surface In BoundingBox{-X-e, -X-e, -30-e, X+e, -X+e, Z+e};
Physical Surface("back")   = Surface In BoundingBox{-X-e, X-e,  -30-e, X+e,  X+e, Z+e};

Mesh.ElementOrder = 1;
