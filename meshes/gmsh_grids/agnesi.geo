
        /*
	Profile of the axisymmetric stenosis following a cosine
	function dependent on the axial coordinate x [1].
 
	 f(x) = hm*a^2/( (x - xc)^2 + a^2 )
	 -- xm mountain height
	 -- xc mountain center

	 References:
 
	[1] Sridhar et al. (2022) GMD ClimateMachine LES paper. See Slack
	*/


	nelx = 50;
 	nelz = 10; 
 
	xmax = 100;
	xmin = 0;
	ymax = 50;
	zmax = 100;

	L = (xmax - xmin)/2;

	//Mountain parameters
	hm = 15;
	a  = 10;	
	xc = (xmax + xmin)/2;
 
	/*******************************************************************************************
	 * Do NOT touch below this point unless you want to change the mountain function at line ~49
	 *******************************************************************************************/
	ls = 5;
	Point(1) = {0,    0,    0, ls};
	Point(2) = {xmax, 0,    0, ls};
	Point(3) = {xmax, 0, zmax, ls};
	Point(4) = {0,    0, zmax, ls};
	
 
	pList[0] = 1; // First point label
	nPoints = 41; // Number of discretization points (top-right point of the inlet region)
	a2 = a*a;
	For i In {1 : nPoints}
	  x = xmax*i/(nPoints + 1);
	  pList[i] = newp;
	  Point(pList[i]) = {x, 0,
	  		  ( hm*a2/((x - xc)*(x - xc) + a2) ), //AGNESI MOUNTAIN
			  ls};
	                
	EndFor
	pList[nPoints+1] = 2; // Last point label (top-left point of the outlet region)	
	Spline(newl) = pList[];

	Line(2) = {2, 3};
	Line(3) = {3, 4};
	Line(4) = {4, 1};
 

 	Transfinite Line {3, -1} = nelx; //Using Progression 1.1;
	Transfinite Line {2, 4} = nelz; //Ceil(L/ls) Using Progression 1;

 
	Line Loop(11) = {4, 1, 2, 3};
	Plane Surface(12) = {11};
	Transfinite Surface {12};
	Recombine Surface {12};
 
	Extrude {0,ymax,0} {
	  Surface{12}; Layers{1}; Recombine;
	}
	Coherence;
 /*
	Physical Surface("symmetryLine") = {51, 37, 73};
	Physical Surface("frontAndBack") = {60, 38, 82, 16, 14, 12};
	Physical Surface("wall") = {59, 29, 81};
	Physical Surface("inlet") = {47};
	Physical Surface("outlet") = {77};
	Physical Volume("volume") = {2, 1, 3};
	*/