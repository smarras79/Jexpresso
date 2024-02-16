function CalcFlux(U, Gamma, d, Tot_X_Pts)
    #CALCFLUX evaluates flow fluxes at all grid-points for Nav-Stokes dynamics
    #   CalcFlux evaluates the flow fluxes for the compressible; inviscid
    #       Navier-Stokes dynamics at all grid-points.
    #
    #       INPUTS:
    #               U = d x Tot_X_Pts array containing primary flow variables
    #               Gamma = ratio of specific heats
    #               d = number of primary flow variables
    #               Tot_X_Pts = number of grid-points
    #
    #       OUTPUT:
    #               F = d x Tot_X_Pts array containing the flow fluxes
    
    F = zeros(d,Tot_X_Pts)
    
    fac1 = (3 - Gamma)/2
    fac2 = (Gamma - 1)/Gamma
    fac3 = (Gamma - 1)/2
    
    # calculate fluxes; loop over all grid-points 
    
    for i = 1:Tot_X_Pts 
                           
        F[1,i] = U[2,i]
        
        # introduce some useful definitions
        
        ratio1 = 1/U[1,i]
        ratio2 = (Gamma*U[2,i])/(U[1,i]*U[1,i])
        
        term1 = U[1,i]*U[3,i]; 
        term2 = U[2,i]*U[2,i]
        
        F[2,i] = ratio1*( fac2 * term1 + fac1*term2 )
        F[3,i] = ratio2*( term1 - fac3*term2 )
    end

    return F
 end
