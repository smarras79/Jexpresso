function CalcFunc( F, J, Del_x, d, Tot_Int_Pts )
    #CALCFUNC evaluates the ODE driver function for 1D Nav.-Stokes dynamics
    #   CalcFunc evaluates ODE driver function for 1D compressible; inviscid
    #       Navier-Stokes flow through a nozzle.
    #
    #   INPUTS:
    #           F = d x Tot_Int_Pts array with flow fluxes at all interior
    #                   grid-points
    #           J = 1 x Tot_Int_Pts array with the source current for
    #                   Navier-Stokes momentum equation of motion at all()
    #                   interior grid-points
    #           Del_x = separation between grid-points
    #           d = number of components of ODE driver function
    #           Tot_Int_Pts = number of interior grid-points
    #
    #   OUTPUT:
    #           ff_vals = d x Tot_Int_Pts array with values of ODE driver
    #                       function at all interior grid-points
    
    ff_vals = zeros(d,Tot_Int_Pts)
    
    Fac1 = (-1)/(2*Del_x)
    
    # evaluate ODE driver function; loop over interior grid-points
    #   NOTE: 1. gridLabel gives grid-point label for each interior grid-point
    #         2. Since F is specified at all grid-points we must use gridLabel
    #               to get the values of F to left & right of an interior
    #               grid-point
    #         3. J is specified only at interior points so i is used to get()
    #               values of j at specified interior point
    
    for i = 1:Tot_Int_Pts
        gridLabel = i + 1
        
        global ff_vals[1,i] = Fac1*( F[1,gridLabel + 1] - F[1,gridLabel - 1])
        global ff_vals[2,i] = Fac1*( (F[2,gridLabel+1] - F[2,gridLabel-1]) - J[i] )
        global ff_vals[3,i] = Fac1*( F[3,gridLabel + 1] - F[3,gridLabel - 1] )
    end
    
    return ff_vals
    end
