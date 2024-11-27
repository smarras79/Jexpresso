function Calc_d2ffdt2(d2Fdt2,d2Jdt2,Del_x,d,Tot_Int_Pts, params)
    #CALC_D2FFDT2 evaluates second time derivative of ODE driver function
    #   Calc_d2ffdt2 evaluates second time derivative of ODE driver function
    #       for 1D compressible; inviscid flow through a nozzle at all()
    #       interior grid-points.
    #
    #   INPUTS:
    #           d2Fdt2 = d x Tot_X_Pts array storing second time derivative
    #                           values of flow flux at all grid-points
    #           d2Jdt2 = 1 x Tot_Int_Pts array storing second time derivative
    #                           values of flow source term at all interior
    #                           grid-points
    #           Del_x = distance between grid-points
    #           d = number of primary flow variables at a grid-point
    #           Tot_Int_Pts = number of interior grid-points
    #
    #   OUTPUT:
    #           d2ffdt2_vals = d x Tot_Int_Pts array storing second time
    #                               derivative values of ODE driver function.
    
    d2ffdt2_vals = zeros(d,Tot_Int_Pts)
    
    TotIPtp1 = Tot_Int_Pts + 1
    
    
    
    # evaluate d2ffdt2_vals at interior grid-points
    #           ( 2 <= ll <= TotIPtp1 )

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
            
            d2ffdt2_vals[1,IPLabel] = fac1*(d2Fdt2[1,ll+1] - d2Fdt2[1,ll-1])
            
            d2ffdt2_vals[2,IPLabel] = fac1*(d2Fdt2[2,ll+1] - d2Fdt2[2,ll-1] 
                                             -d2Jdt2[IPLabel])
                                          
            d2ffdt2_vals[3,IPLabel] = fac1*(d2Fdt2[3,ll+1] - d2Fdt2[3,ll-1]); 
         end
      end
    
   #  for ll = 2: TotIPtp1
   #    fac1 = (-1)/(2*params.mesh.Δx[ll])
   #     IPLabel = ll -1
       
   #     d2ffdt2_vals[1,IPLabel] = fac1*(d2Fdt2[1,ll+1] - d2Fdt2[1,ll-1])
       
   #     d2ffdt2_vals[2,IPLabel] = fac1*(d2Fdt2[2,ll+1] - d2Fdt2[2,ll-1] 
   #                                      -d2Jdt2[IPLabel])
                                    
   #     d2ffdt2_vals[3,IPLabel] = fac1*(d2Fdt2[3,ll+1] - d2Fdt2[3,ll-1]);    
   #  end
    
    return d2ffdt2_vals
    
    end
