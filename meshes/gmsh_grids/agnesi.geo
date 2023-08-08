// User inputs (except for B.C. below)

//Domain:
xmin = -100.0; // um
xmax = 100;
zmax = 50.0;   // um

//Set nelx = -1 if you want to use dx, dz intead
nelx = -10;
nelz = 10;

dx = 1.2;
dz = 1.0;

stretchingx = 1.0;
stretchingz = 1.0;

//Mountain:
hm = 5.0;
a = 2.5;
xc = (xmin + xmax)/2.0;
// End user inputs (except for B.C. below)


ls = 0.0; //this is useless since we use transfinite with defined dx/dz or nelx/nelz
Point(1) = {xmin, 0, 0, ls};
Point(2) = {xmax, 0, 0, ls};
Point(3) = {xmax, 0, zmax, ls};
Point(4) = {xmin, 0, zmax, ls};

Line(1) = {3, 2};
Line(2) = {4, 3};
Line(3) = {1, 4};

pList[0] = 1; // First point label
nPoints = 401; // Number of discretization points (top-right point of the inlet region)
For i In {1 : nPoints}
    x = xmin + (xmax - xmin)*i/(nPoints + 1);
    pList[i] = newp;
    Point(pList[i]) = {x,
		       0,
		       hm*a*a/((x - xc)*(x - xc) + a*a), 
		       ls};
EndFor
pList[nPoints+1] = 2; // Last point label (top-left point of the outlet region)

Spline(newl) = pList[];

If (nelx < 0)
   Transfinite Line {2, 4}  = Ceil((xmax - xmin)/dx) Using Progression 1;
   Transfinite Line {-1, 3} = Ceil(zmax/dz) Using Progression stretchingz;
Else
	Transfinite Line {2, 4}  = nelx Using Progression 1;
	Transfinite Line {-1, 3} = nelz Using Progression stretchingz;
EndIf

Line Loop(11) = {1, 2, 3, -4};
Plane Surface(12) = {11};
Transfinite Surface {12};
Recombine Surface {12};


//BOUNDARY DATA:
Physical Point("boundary",   1) = {1, 2, 3, 4};
Physical Curve("periodic",   2) = {1, 3};
Physical Curve("free_slip",  3) = {4};
Physical Curve("Laguerre",   4) = {2};
//Physical Surface("domain")      = {1};
