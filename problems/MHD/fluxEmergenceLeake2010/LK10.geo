// ============================================================
// Domain and grid of the 2.5D flux-emergence simulations of
//
//   J. E. Leake, M. G. Linton & S. K. Antiochos, "Tests of dynamical flux
//   emergence as a mechanism for coronal mass ejection initiation",
//   ApJ 722:550 (2010), Sec. 2.1; Erratum ApJ 741:125 (2011),
//
// in units of L0 = 150 km (the normalizing length of the erratum):
//
//   x in [-X, X] = [-300, 300]   (-45 to 45 Mm)
//   y in [-20, Y] = [-20, 600]   (-3 to 90 Mm)
//
// The paper stretches its grid "in both x and y", from 0.0375 Mm (0.25 L0) at
// the surface to 0.15 Mm (1 L0) at the top, without giving the map. Here the
// ELEMENT boundaries follow a two-level tanh grading of the element size,
//
//   s(chi) = chi + a1 [chi + w1 ln(cosh((chi-c1)/w1)/cosh(c1/w1))]
//                + a2 [chi + w2 ln(cosh((chi-c2)/w2)/cosh(c2/w2))],
//
// whose slope ds/dchi is 1 in the core, 1 + 2 a1 = 2 past c1 and
// 1 + 2 a1 + 2 a2 = 10 past c2, sampled at n equal steps of chi in
// [0, chi_max] (chi_max solves s(chi_max) = half-width, resp. height):
//
//   x (|x| from 0):     a1 = 0.5, c1 = 40,  w1 = 6;  a2 = 4, c2 = 100, w2 = 15
//   y (from y = -20):   a1 = 0.5, c1 = 55,  w1 = 6;  a2 = 4, c2 = 160, w2 = 20
//
// so the core element size h0 holds over |x| < 35 and -20 < y < 30 (the
// convection zone, the tube, the photosphere, chromosphere and transition
// region), 2 h0 over the region of the paper's figures (|x| < 85, y < 130,
// i.e. 12.75 Mm and 19.5 Mm), and grows to ~9 h0 at the far boundaries.
//
//   ref = 1 -> h0 = 2.5  L0:  90 x  78 elements  LK10_ref1.msh  smoke test
//   ref = 2 -> h0 = 1.25 L0: 180 x 156 elements  LK10_ref2.msh  ~ the paper's
//              resolution at :nop => 4 (nodal spacing 0.31 L0 in the core)
//   ref = 3 -> h0 = 0.625 L0: 362 x 312 elements LK10_ref3.msh
//   ref = 4 -> h0 = 1.25 L0 over the WHOLE figure region (|x| < 85, y < 130)
//              and one grading level to ~10 h0 beyond it: 200 x 172 elements,
//              LK10_ref4.msh. Nodal spacing 0.33 L0 (49 km) where the dome and
//              the flux rope of the erratum's Fig. 3 evolve, as the paper's grid
//              there (37.5 km at the surface, 150 km at the top), against 0.63 L0
//              (94 km) for ref = 2; 1.23 times the elements of ref = 2, same Δt.
//
// Boundary tags (the strings user_bc_dirichlet! receives):
//   "bottom" y = -20,  "top" y = 600,  "left" x = -300,  "right" x = 300.
//
// Generate with:
//   gmsh -2 LK10.geo -setnumber ref 1 -format msh41 -o LK10_ref1.msh
//   gmsh -2 LK10.geo -setnumber ref 2 -format msh41 -o LK10_ref2.msh
//   gmsh -2 LK10.geo -setnumber ref 3 -format msh41 -o LK10_ref3.msh
//   gmsh -2 LK10.geo -setnumber ref 4 -format msh41 -o LK10_ref4.msh
// ============================================================

DefineConstant[ ref = 1 ];
If (ref == 1)
  nh = 45;  nv = 78;
EndIf
If (ref == 2)
  nh = 90;  nv = 156;
EndIf
If (ref == 3)
  nh = 181; nv = 312;
EndIf
If (ref == 4)
  nh = 100; nv = 172;
EndIf

X = 300.0;  Y0 = -20.0;  Y = 600.0;

// half-width map (|x| from 0 to X)
ax1 = 0.5; cx1 = 40.0;  wx1 = 6.0;  ax2 = 4.0; cx2 = 100.0; wx2 = 15.0; chix = 113.026993;
// height map (y - Y0 from 0 to Y - Y0)
ay1 = 0.5; cy1 = 55.0;  wy1 = 6.0;  ay2 = 4.0; cy2 = 160.0; wy2 = 20.0; chiy = 195.268206;
// ref = 4: one level, from h0 over the figure region to 1 + 2 a2 = 10 h0 beyond it
// (|x| = 105 and y - Y0 = 170 are 2.5 widths past the region's edges |x| = 85, y = 130)
If (ref == 4)
  ax1 = 0.0; ax2 = 4.5; cx2 = 105.0; wx2 = 8.0; chix = 124.472429;
  ay1 = 0.0; ay2 = 4.5; cy2 = 170.0; wy2 = 8.0; chiy = 214.999953;
EndIf

For k In {0:nh}
  chi = k*chix/nh;
  xh[k] = chi + ax1*(chi + wx1*Log(Cosh((chi - cx1)/wx1)/Cosh(cx1/wx1)))
              + ax2*(chi + wx2*Log(Cosh((chi - cx2)/wx2)/Cosh(cx2/wx2)));
EndFor
// cumulative fractions of the full span [-X, X] (exact end points)
For m In {1:2*nh}
  If (m <= nh)
    fx[m-1] = (xh[nh] - xh[nh - m])/(2*xh[nh]);
  Else
    fx[m-1] = (xh[nh] + xh[m - nh])/(2*xh[nh]);
  EndIf
  ox[m-1] = 1;
EndFor

For k In {0:nv}
  chi = k*chiy/nv;
  yv[k] = chi + ay1*(chi + wy1*Log(Cosh((chi - cy1)/wy1)/Cosh(cy1/wy1)))
              + ay2*(chi + wy2*Log(Cosh((chi - cy2)/wy2)/Cosh(cy2/wy2)));
EndFor
For m In {1:nv}
  fy[m-1] = yv[m]/yv[nv];
  oy[m-1] = 1;
EndFor

Point(1) = {-X, Y0, 0};
ex[] = Extrude {2*X, 0, 0} { Point{1}; Layers{ox[], fx[]}; };
ey[] = Extrude {0, Y - Y0, 0} { Curve{ex[1]}; Layers{oy[], fy[]}; Recombine; };

e = 1e-6;
Physical Surface("domain") = {ey[1]};
Physical Curve("bottom") = Curve In BoundingBox{-X-e, Y0-e, -e, X+e, Y0+e, e};
Physical Curve("top")    = Curve In BoundingBox{-X-e, Y-e,  -e, X+e, Y+e,  e};
Physical Curve("left")   = Curve In BoundingBox{-X-e, Y0-e, -e, -X+e, Y+e, e};
Physical Curve("right")  = Curve In BoundingBox{X-e,  Y0-e, -e, X+e,  Y+e, e};

Mesh.ElementOrder = 1;
