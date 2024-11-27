function Calc_dffdt(dFdt, dJdt, Del_x, d, Tot_Int_Pts, params)
    #CALC_DFFDT evaluates first time derivative of ODE driver function
    #   Calc_dffdt evaluates first time derivative of ODE driver function at
    #       all interior grid-points.
    #
    #   INPUTS:
    #           dFdt = d x Tot_X_Pts array storing first time derivative of
    #                       flow flux at all grid-points
    #           dJdt = 1 x Tot_Int_Pts array storing first time derivative of
    #                       flow source term at all interior grid-points
    #           Del_x = distance between grid-points
    #           d = number of primary flow variables at a grid-point
    #           Tot_Int_pts = number of interior grid-points
    #
    #   OUTPUT:
    #           dffdt_vals = d x Tot_Int_Pts array storing first time
    #                           derivative of ODE driver function
    
    dffdt_vals = zeros(d,Tot_Int_Pts)
    
    TotIPtp1 = Tot_Int_Pts + 1
    # evaluate dffdt_vals at interior points [2 <= ll <= TotIPtp1]
    for i=1:params.mesh.nelem #over each element
        for j=1:params.mesh.ngl 
            #@info params.mesh.npoin
            ll = params.mesh.connijk[i, j, 1]
            #@info ll
            if (ll == 1)
                continue
            elseif (ll==params.mesh.npoin)
                continue
            elseif i==j
                continue
            end
            fac1 = (-1)/(2*params.mesh.Δx[j])
            IPLabel = ll -1
            
            dffdt_vals[1,IPLabel] = fac1*(dFdt[1,ll+1] - dFdt[1,ll-1])
            
            dffdt_vals[2,IPLabel] = fac1*(dFdt[2,ll+1] - dFdt[2,ll-1] 
                                            -dJdt[IPLabel])
            
            dffdt_vals[3,IPLabel] = fac1*(dFdt[3,ll+1] - dFdt[3,ll-1])
        end
    end


    # for ll = 2:TotIPtp1
    #     fac1 = (-1)/(2*params.mesh.Δx[ll])
    #     IPLabel = ll -1
        
    #     dffdt_vals[1,IPLabel] = fac1*(dFdt[1,ll+1] - dFdt[1,ll-1])
        
    #     dffdt_vals[2,IPLabel] = fac1*(dFdt[2,ll+1] - dFdt[2,ll-1] 
    #                                     -dJdt[IPLabel])
        
    #     dffdt_vals[3,IPLabel] = fac1*(dFdt[3,ll+1] - dFdt[3,ll-1])
    # end
    
    return dffdt_vals
    
    end
    
    
    
