function CalcfBvalspSW(U,Gamma,ff_vals,d,Tot_Int_Pts)
    #CALCFBVALSPSW assigns boundary values to ODE driver fcn - with shock-wave 
    #   CalcfBvalspSW uses flow boundary conditions to assign "values" to ODE
    #       driver function at nozzle entrance/exit when shock-wave present
    #       in nozzle. Actually calculating value of 
    #       dU[k,ll = 1, Tot_X_Pts]/dt at boundary points using boundary 
    #       conditions. Formulas below based on analysis in 7-23-2019 entry 
    #       in my research notebook.
    #
    #   INPUTS:
    #           U = d x Tot_X_Pts array storing flow variables at grid-points
    #           Gamma = ratio of specific heats
    #           ff_vals = d x Tot_Int_Pts array storing ODE driver function
    #                           values at interior grid-points
    #           d = number of components of ODE driver function
    #           Tot_Int_Pts = number of interior grid-points
    #
    #   OUTPUT:
    #           ff_Bvals = d x 2 array storing ODE driver function values at
    #                           nozzle entrance & exit()
    
    ff_Bvals = zeros(d,2)
    
    Fac1 = Gamma/2
    Fac3 = U[2,1]/U[1,1]
    
    # calculate ff_Bvals at nozzle entrance
    
    ff_Bvals[1,1] = 0
    ff_Bvals[2,1] = 2*ff_vals[2,1] - ff_vals[2,2]
    ff_Bvals[3,1] = Gamma*Fac3*ff_Bvals[2,1]
    
    # calculate use boundary conditions on flow variables to determine
    #   ff_Bvals at nozzle exit()
    
    for k = 1:(d-1)
        ff_Bvals[k,2] = 2*ff_vals[k,Tot_Int_Pts] - ff_vals[k,(Tot_Int_Pts-1)]
    end
    
    ff_Bvals[3,2] = Fac1*Fac3*( 2*ff_Bvals[2,2] - Fac3*ff_Bvals[1,2] )
    return ff_Bvals
    end
