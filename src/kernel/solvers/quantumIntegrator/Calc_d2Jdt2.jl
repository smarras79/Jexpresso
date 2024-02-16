function Calc_d2Jdt2(U,dffdt_vals,ff_vals,A,Gamma,d,Tot_Int_Pts)
    #CALC_D2JDT2 evaluates 2nd time derivative of source term at interior pts
    #   Calc_d2Jdt2 evaluates the second time derivative of the flow source
    #       term for 1D compressible; inviscid flow through nozzle at all()
    #       interior grid-points.
    #
    #   INPUTS:
    #           U = d x Tot_X_Pts array storing primary flow variables at all()
    #                   grid-points
    #           dffdt_vals = d x Tot_Int_Pts array storing first time
    #                           derivative of source term at interior
    #                           grid-points
    #           ff_vals = d x Tot_Int_Pts array storing ODE driver function
    #                           values at all interior grid-points
    #           A = 1 x Tot_X_Pts array storing nozzle area at all grid-points
    #           Gamma = ratio of specific heats
    #           d = number of primary flow variables at each grid-point
    #           Tot_Int_Pts = number of interior grid-points
    #
    #   OUTPUT:
    #           d2Jdt2 = 1 x Tot_Int_Pts array storing second time derivative
    #                       of flow source term at all interior grid-points
    
    d2Jdt2 = zeros(1,Tot_Int_Pts)
    
    TotIPtp1 = Tot_Int_Pts + 1
    
    fac1 = (Gamma - 1)/Gamma
    fac2 = Gamma/2
    
    # evaluate d2Jdt2 by looping over interior grid-points
    #           ( 2 <= ll <= TotIPtp1 )
    
    for ll = 2: TotIPtp1
        IPLabel = ll -1
        
        fac3 = U[2,ll]/U[1,ll]
        fac4 = log( A[ll+1]/A[ll-1] )
        
        fac5 = (fac1)^(2);          # (U[2,ll]/U[1,ll])^(2)
        fac6 = 1/U[1,ll];           # 1/U[1,ll]
        
        fac7 = (fac6)^(2);          # 1/(U[1,ll])^(2)
        fac8 = fac6*fac7;           # 1/(U[1,ll])^(3)
        
        fac9 = fac3*fac6;           # U[2,ll]/(U[1,ll])^(2)
        
        d2Jdt2[IPLabel] = fac4*fac1*( 
                            dffdt_vals[3,IPLabel] 
                             -fac2*( 2*fac3*dffdt_vals[2,IPLabel]
                              -fac5*dffdt_vals[1,IPLabel]
                               +2*fac6*(ff_vals[2,IPLabel])^(2)
                                +2*fac8*(U[2,ll]*ff_vals[1,IPLabel])^(2)
                            -4*fac9*ff_vals[1,IPLabel]*ff_vals[2,IPLabel] ) )
    end
    
    return d2Jdt2
    
    end
