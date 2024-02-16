function Calc_d2Fdt2(U,dffdt_vals,dffdt_Bvals,ff_vals,
    ff_Bvals,Gamma,d,Tot_X_Pts)
#CALC_D2FDT2 evaluates flow flux second time derivative at all grid-points
#   Calc_d2Fdt2 evaluates flow flux second time derivative at all()
#     grid-points for 1D compressible; inviscid flow through a nozzle.
#
#   INPUTS:
#           U = d xTot_X_Pts array storing primary flow variables at all()
#                       grid-points
#           dffdt_vals = d x Tot_Int_Pts array storing first time
#                           derivatives of ODE driver function at all()
#                           interior grid-points.
#           dffdt_Bvals = d x 2 array storing first time derivative of ODE
#                           driver function at nozzle entrance [exit] in
#                           cloumn 1 [2].
#           ff_vals = d x Tot_Int_Pts array storing ODE driver function
#                           values at interior grid-points.
#           ff_Bvals = d x 2 array storing ODE driver function values at
#                           nozzle entrance [exit] in cloumn 1 [2].
#           Gamma = ratio of specific heats
#           d = number of primary flow variables at a grid-point
#           Tot_X_Pts = number of grid-points
#
#   OUTPUT:
#           d2Fdt2 = d x Tot_X_Pts array storing second time derivative of
#                           flow fluxes at all grid-points

d2Fdt2 = zeros(d, Tot_X_Pts)

TotXPtm1 = Tot_X_Pts - 1

fac4 = Gamma - 1
fac5 = Gamma/2;         
fac6 = fac4*fac5;           # (Gamma - 1)*(Gamma/2)
fac7 = (3 - Gamma)/2
fac8 = fac4/Gamma;          # (Gamma - 1)/Gamma

# evaluate d2Fdt2 at boundaries

for ll = 1:2
     if (ll == 1) # at nozzle ENTRANCE
                     
               fac1 = U[2,ll]/U[1,ll]
               fac2 = U[3,ll]/U[1,ll]
               fac3 = fac1*fac2;       # U[2,ll]U[3,ll]/(U[1,ll])^2

               fac9 = (fac1)^(2);      # (U[2,ll]/U[1,ll])^(2)
               fac10 = 1/U[1,ll]
               fac11 = fac9*fac10;     # (U[2,ll])^(2)/(U[1,ll])^3

               fac12 = fac1*fac10;     # U[2,ll]/(U[1,ll]^2
               fac13 = fac2*fac10;     # U[3,ll]/(U[1,ll])^2
               fac14 = fac1*fac9;      # (U[2,ll]/U[1,ll])^3

               fac15 = fac3*fac10;     # (U[2,ll]U[3,ll])/(U[1,ll])^(3)
               fac16 = fac14*fac10;    # (U[2,ll])^(3)/(U[1,ll])^(4)

               d2Fdt2[1,ll] = dffdt_Bvals[2,ll]

               d2Fdt2[2,ll] = fac7*( fac1*(2*dffdt_Bvals[2,ll] 
                                             -fac1*dffdt_Bvals[1,ll])
                                   +2*fac10*(ff_Bvals[2,ll])^(2) 
                                   +2*fac11*(ff_Bvals[1,ll])^(2) 
                                   -4*fac12*ff_Bvals[1,ll]*ff_Bvals[2,ll])
                              +fac8*dffdt_Bvals[3,ll]
    
               d2Fdt2[3,ll] = 
                    Gamma*( fac1*dffdt_Bvals[3,ll] + fac2*dffdt_Bvals[2,ll] 
                              -fac3*dffdt_Bvals[1,ll] 
                              +2*(fac10*ff_Bvals[2,ll]*ff_Bvals[3,ll] 
                                   -fac12*ff_Bvals[1,ll]*ff_Bvals[3,ll] 
                                   -fac13*ff_Bvals[1,ll]*ff_Bvals[2,ll] 
                                   +fac15*(ff_Bvals[1,ll])^(2) ) ) 
                    -fac6*( 3*fac9*dffdt_Bvals[2,ll] 
               -2*fac14*dffdt_Bvals[1,ll] 
               +6*fac12*(ff_Bvals[2,ll])^(2) 
               +6*fac16*(ff_Bvals[1,ll])^(2) 
               -12*fac11*ff_Bvals[1,ll]*ff_Bvals[2,ll] ); 
     elseif(ll == 2)            # at nozzle EXIT
               fac1 = U[2,Tot_X_Pts]/U[1,Tot_X_Pts]
               fac2 = U[3,Tot_X_Pts]/U[1,Tot_X_Pts]
               fac3 = fac1*fac2;       #U[2,TotXP]U[3,TotXPt]/(U[1,TotXPt])^2

               fac9 = (fac1)^(2);      # (U[2,TotXPt]/U[1,TotXPt])^(2)
               fac10 = 1/U[1,Tot_X_Pts]
               fac11 = fac9*fac10;     # (U[2,TotXPt])^(2)/(U[1,TotXPt])^3

               fac12 = fac1*fac10;     # U[2,TotXPt]/(U[1,TotXPt]^2
               fac13 = fac2*fac10;     # U[3,TotXPt]/(U[1,TotXPt])^2
               fac14 = fac1*fac9;      # (U[2,TotXPt]/U[1,TotXPT])^3

               fac15 = fac3*fac10;     #(U[2,ll]U[3,TotXPT])/(U[1,TotXPt])^(3)
               fac16 = fac14*fac10;    # (U[2,TotXPt])^(3)/(U[1,TotXPt])^(4)

               d2Fdt2[1,Tot_X_Pts] = dffdt_Bvals[2,ll]

               d2Fdt2[2,Tot_X_Pts] = fac7*( fac1*(2*dffdt_Bvals[2,ll] 
                                             -fac1*dffdt_Bvals[1,ll])
                                   +2*fac10*(ff_Bvals[2,ll])^(2) 
                                   +2*fac11*(ff_Bvals[1,ll])^(2) 
                                   -4*fac12*ff_Bvals[1,ll]*ff_Bvals[2,ll])
                              +fac8*dffdt_Bvals[3,ll]
               d2Fdt2[3,Tot_X_Pts] = 
                    Gamma*( fac1*dffdt_Bvals[3,ll] + fac2*dffdt_Bvals[2,ll] 
                              -fac3*dffdt_Bvals[1,ll] 
                              +2*(fac10*ff_Bvals[2,ll]*ff_Bvals[3,ll] 
                                   -fac12*ff_Bvals[1,ll]*ff_Bvals[3,ll] 
                                   -fac13*ff_Bvals[1,ll]*ff_Bvals[2,ll] 
                                   +fac15*(ff_Bvals[1,ll])^(2) ) ) 
                    -fac6*( 3*fac9*dffdt_Bvals[2,ll] 
                         -2*fac14*dffdt_Bvals[1,ll] 
                         +6*fac12*(ff_Bvals[2,ll])^(2) 
                         +6*fac16*(ff_Bvals[1,ll])^(2) 
                         -12*fac11*ff_Bvals[1,ll]*ff_Bvals[2,ll] );             
     else
               disp( ["Unknown switch case label: " ll])
     end
end

# evaluate d2Fdt2 at interior points

for ll = 2:TotXPtm1
     IPLabel = ll - 1

     fac1 = U[2,ll]/U[1,ll]
     fac2 = U[3,ll]/U[1,ll]
     fac3 = fac1*fac2;       # U[2,ll]U[3,ll]/(U[1,ll])^2

     fac9 = (fac1)^(2);      # (U[2,ll]/U[1,ll])^(2)
     fac10 = 1/U[1,ll]
     fac11 = fac9*fac10;     # (U[2,ll])^(2)/(U[1,ll])^3

     fac12 = fac1*fac10;     # U[2,ll]/(U[1,ll]^2
     fac13 = fac2*fac10;     # U[3,ll]/(U[1,ll])^2
     fac14 = fac1*fac9;      # (U[2,ll]/U[1,ll])^3

     fac15 = fac3*fac10;     # (U[2,ll]U[3,ll])/(U[1,ll])^(3)
     fac16 = fac14*fac10;    # (U[2,ll])^(3)/(U[1,ll])^(4)

     d2Fdt2[1,ll] = dffdt_vals[2,IPLabel]

     d2Fdt2[2,ll] = fac7*( fac1*(2*dffdt_vals[2,IPLabel] 
                                        -fac1*dffdt_vals[1,IPLabel])
                              +2*fac10*(ff_vals[2,IPLabel])^(2) 
                              +2*fac11*(ff_vals[1,IPLabel])^(2) 
                         -4*fac12*ff_vals[1,IPLabel]*ff_vals[2,IPLabel])
                              +fac8*dffdt_vals[3,IPLabel]
     d2Fdt2[3,ll] = 
                    Gamma*( fac1*dffdt_vals[3,IPLabel]
                              + fac2*dffdt_vals[2,IPLabel] 
                              -fac3*dffdt_vals[1,IPLabel] 
                         +2*(fac10*ff_vals[2,IPLabel]*ff_vals[3,IPLabel] 
                              -fac12*ff_vals[1,IPLabel]*ff_vals[3,IPLabel] 
                              -fac13*ff_vals[1,IPLabel]*ff_vals[2,IPLabel] 
                              +fac15*(ff_vals[1,IPLabel])^(2) ) ) 
                    -fac6*( 3*fac9*dffdt_vals[2,IPLabel] 
                              -2*fac14*dffdt_vals[1,IPLabel] 
                              +6*fac12*(ff_vals[2,IPLabel])^(2) 
                              +6*fac16*(ff_vals[1,IPLabel])^(2) 
                         -12*fac11*ff_vals[1,IPLabel]*ff_vals[2,IPLabel] )
end


return d2Fdt2;

end