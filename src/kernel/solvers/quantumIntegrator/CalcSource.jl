function CalcSource(U, A, Gamma, Tot_X_Pts)
    #CALCSOURCE evaluates the source current in the Navier-Stokes dynamics
    #   CalcSource evaluates the source current in the momentum equation of
    #       the compressible; inviscid Navier-Stokes dynamics at interior
    #       grid-points.
    #
    #       INPUTS:
    #               U = d x Tot_X_Pts array containing primary flow variables
    #               A = 1 x Tot_X_Pts array containing nozzle area at each
    #                       grid-point
    #               Gamma = ratio of specific heats
    #               Tot_X_Pts = number of grid-points
    #
    #       OUTPUTS:
    #               J = 1 x Tot_Int_Pts array containing source current
    
    Tot_Int_Pts = Tot_X_Pts -2
    
    J = zeros(1,Tot_Int_Pts)
    
    
    # calculate source current; loop over interior grid-points with
    #       2 <= ll <= Tot_Int_Pts +1
    
    for i = 2:(Tot_Int_Pts + 1)    # ll indexes the grid-points
        colPos = i -1;             # colPos indexes the interior grid-points
                                    #   it indexes the columns of J
                                    
        fac1 = ( (Gamma - 1)/Gamma )*log(A[i+1]/A[i-1])
        fac2 = Gamma/2
        
        term1 = U[1,i]*U[3,i]
        term2 = U[2,i]*U[2,i]
                                
        ratio1 = 1/U[1,i]
        
        global J[colPos] = fac1*ratio1*( term1 - fac2*term2 )
    end
        
    return J
end
