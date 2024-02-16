function BldTMat(f,d,r,InVal,Tot_Int_Pts,Tot_X_Pts)
    #BLDTMAT constructs Taylor polynomial coefficients 
    #   BldTMat constructs Taylor polynomial coefficients for a subsubinterval
    #   at each grid-point & stores in array mat.
    #
    #   INPUTS: f = d x [r+1] x Tot_Int_Pts array storing ODE driver function 
    #                   & its first r derivatives at each grid-point
    #                   evaluated at InVal 
    #           d = number of components of ODE driver function
    #           r = maximum number of derivatives of ODE driver function
    #           InVal = d x Tot_X_Pts array with initial condition for Taylor
    #                    polynomials in subinterval i
    #           Tot_Int_Pts = number of interior grid-points
    #           Tot_X_Pts = number of grid-points
    #
    #   OUTPUT:
    #           mat = d x [r+2] x Tot_Int_Pts array with Taylor polynomial
    #                       coefficients for a given subsubinterval. Each
    #                       page of mat stores the Taylor polynomial
    #                       coefficients at an interior grid-point.
    
    rmax = r + 1;       # degree of Taylor polynomials
    rmaxp1 = rmax + 1;  # number of terms in Taylor polynomials
    
    # for given subsubinterval mat will return the Taylor polynomial 
    #  coefficients for each component of approximation to ODE driver 
    #  function & interior grid-point
    
    mat = zeros(d, rmaxp1, Tot_Int_Pts)
    
    # evaluate Taylor polynomial coefficients by looping around the interior
    #       grid-points; then the d components; then the rmax coefficients
    
    
    for ll = 1:Tot_Int_Pts
       for k = 1:d
         for m = 1:rmax
            mat[k,m,ll] = f[k,(rmaxp1 - m), ll]/factorial(rmaxp1 - m)
         end
        
         # assign value of polynomial at start of subsubinterval. since InVal
         #     is defined at all grid-points; grid-point label is ll + 1 for
         #     InVal
    
         mat[k,rmaxp1, ll] = InVal[k,(ll+1)];   
       end 
    end
    
    return mat;
    
    end
    
    
    
