function Derivs( d,r,InitVal,Del_x,Gamma,Tot_Int_Pts,Tot_X_Pts,
    A,Shock_Flag, params, U)
#DERIVS evaluates ODE driver function f & its first r time derivatives
#   Derivs evaluates the ODE driver function f(z) & its first r time
#    derivatives at z = InitVal. 
#
#   INPUTS:
#           d = number of components of ODE driver function
#           r = maximum number of time derivatives to evaluate
#           InitVal = d x Tot_X_Pts array storing flow variables at
#                           all grid-points
#           Del_x = separation between grid-points
#           Gamma = ratio of specific heats
#           Tot_Int_Pts = number of interior grid-points
#                       = Tot_X_Pts - 2
#           Tot_X_Pts = number of grid-points
#           A = 1 x Tot_X_Pts array storing nozzle area at each grid-point
#           Shock_Flag = 0[1] if shock-wave absent [present] in nozzle
#
#   OUTPUT: 
#           ff = d x [r+1] x Tot_Int_Pts array where each page corresponds 
#                   to a grid-point; & on each page the column indexes 
#                   time derivatives; & row indexes the components of 
#                   f(z). 
#
# Support functions: CalcFlux; CalcSource; CalcFunc; CalcfBvalsmSW
#                       CalcfBvalspSW; Calc_dFdt; Calc_dJdt; Calc_dffdt; 
#                           CalcdfdtBvalsmSW; CalcdfdtBvalspSW
#                           Calc_d2Fdt2; Calc_d2Jdt2; Calc_d2ffdt2

rmax = r+1;     # max column label [store fcn & first r derivatives]

# F = d x Tot_X_Pts array - stores fluxes contributing to ff at grid-points

# J = 1 x Tot_Int_Pts array - stores source current contributing to ff
#       at interior grid-points

ff = zeros(d,rmax,Tot_Int_Pts); # store i-th derivative values in column i
      #  for each component of driver function &
      #  interior grid-point

     
# the following formulas were derived in my black research notebook entry
#   for 03/08/2019.

# Calculate ODE driver function ff at all interior grid-points
#
# 1.a Evaluate flow fluxes F which is a d x Tot_X_Pts array
#=
u = zeros(Float64, params.mesh.npoin*d)
for ip=1:Tot_X_Pts
    for ieq=0:d-1
        u[Tot_X_Pts*ieq + ip] = U[ieq+1, ip]
    end
end

u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)

for i=1:Tot_X_Pts
    user_flux!(@view(params.F[i,:]), @view(params.G[i,:]), NSD_1D(),
            @view(params.uaux[i,:]),
            @view(params.qp.qe[i,:]),       #pref
            params.mesh,
            params.CL,
            params.SOL_VARS_TYPE;
            neqs=params.neqs, ip=i)

    if params.inputs[:lsource]
        user_source!(@view(params.S[ip,:]),
                         @view(params.uaux[ip,:]),
                         @view(params.qp.qe[ip,:]), #ρref 
                         params.mesh.npoin,
                         params.CL,
                         params.SOL_VARS_TYPE;
                         neqs=params.neqs,
                         x=params.mesh.x[i],
                         y=params.mesh.y[i],
                         xmax=params.xmax,
                         xmin=params.xmin,
                         ymax=params.ymax)
    end
end

#F = CalcFlux(InitVal,Gamma,d,Tot_X_Pts, params)

# 1.b Evaluate source current J in momentum equation of motion which is()
#       a 1 x Tot_Int_Pts array

#J = CalcSource(InitVal,A,Gamma,Tot_X_Pts)

# 1.c Now calculate ff[d,1, Tot_Int_Pts] which stores driver function 
#       values. First calculate ff_vals = d x Tot_Int_Pts array which() 
#       stores ff values at all interior grid-point.

ff_vals = CalcFunc(Del_x,d,Tot_Int_Pts, params) =#

# 1.d Now store ff_vals in d x rmax x Tot_Int_Pts array ff

u = zeros(Float64, params.mesh.npoin*d)
for ip=1:Tot_X_Pts
    for ieq=0:d-2 #TODO: change back to d-1 after done with wave-fd
        u[Tot_X_Pts*ieq + ip] = U[ieq+1, ip]
    end
end

# replace with call to _build_rhs
#build_rhs!(params.RHS, u, params, 0.0)
inviscid_rhs_el!(u, params, params.inputs[:lsource], NSD_1D(), params.inputs[:AD])
# for m=1:d
#     U[m,1] = params.uaux[1,m]
#     U[m,Tot_X_Pts] = params.uaux[Tot_X_Pts, m]
# end
ff_vals = zeros(Float64, d, Tot_X_Pts)
for i=1:Tot_Int_Pts
    for j=1:d-1 #TODO: change back to d after wave1d
        ff_vals[j, i] = params.RHS[i+1, j]
    end
end

for ll = 1:Tot_Int_Pts
    for k = 1:d
        ff[k,1,ll] = ff_vals[k,ll+1]
    end
end

# determine first [time] derivatives of ODE driver function at all()
#   interior grid-points

#   2.a calculate "values" of ODE driver function at boundary points - I
#        really mean value of dU[k,1l =1,Tot_X_Pts ]/dt at boundary points
#        which are found from the boundary conditions. I abuse notation
#        and denote these derivatives by f[k,1] & f[k,Tot_X_Pts] which()
#        are stored in columns 1 & 2 of ff_Bvals.
#
#       ff_Bvals = d x 2 array storing "values" of ODE driver function
#                       at nozzle entrance [exit] in column 1 [2].

ff_Bvals = zeros(Float64, d, 2)
ff_Bvals[1, 1] = ff_vals[1, 1]
ff_Bvals[1, 2] = ff_vals[1, Tot_X_Pts]
ff_Bvals[2, 1] = ff_vals[2, 1]
ff_Bvals[2, 2] = ff_vals[2, Tot_X_Pts]
#= TODO: uncomment if necessary
if Shock_Flag .== 0
    ff_Bvals = CalcfBvalsmSW(InitVal,Gamma,ff_vals,d,Tot_Int_Pts)
elseif Shock_Flag .== 1
    ff_Bvals = CalcfBvalspSW(InitVal,Gamma,ff_vals,d,Tot_Int_Pts)
else
    disp(["Unkown Shock_Flag value: " int2str(Shock_Flag)])
end
=#

#if r>0
#   2.b calculate first time derivative of flow fluxes at all grid-points
#       
#       dFdt = d x Tot_X_Pts array storing the flow flux first time
#                               derivatives at all grid-points.

    dFdt = Calc_dFdt(InitVal,ff_vals,ff_Bvals,Gamma,d,Tot_X_Pts)

    #   2.c calculate first time derivative of flow source term at all()
    #           interior grid-points
    #
    #       dJdt = 1 x Tot_Int_Pts array storing the flow source term
    #                               first time derivative at all interior
    #                               grid-points

    dJdt = zeros(Float64, 1, Tot_Int_Pts)#Calc_dJdt(InitVal,ff_vals,A,Gamma,d,Tot_Int_Pts) #TODO: change back later

    #   2.d calculate dffdt_vals = d x Tot_Int_Pts array which() 
    #       stores values of first time derivative of ff at all interior 
    #       grid-points.

    dffdt_vals = Calc_dffdt(dFdt, dJdt, Del_x, d, Tot_Int_Pts, params)

    #   2.e Now store dffdt_vals in column 2 of d x rmax x Tot_Int_Pts array 
    #       ff

    for ll = 1:Tot_Int_Pts
        for k = 1:d
            ff[k,2,ll] = dffdt_vals[k,ll]
        end
    end

    # determine second [time] derivatives of ODE driver function at all()
    #   interior grid-points

    #   3.a calculate "values" first time derivative of ODE driver function 
    #        at boundary points - I really mean value of 
    #        d^(2)U[k,1l =1,Tot_X_Pts ]/dt^(2) at boundary points which are 
    #        found from the boundary conditions. I abuse notation & denote 
    #        these derivatives by df[k,1]/dt & df[k,Tot_X_Pts]/dt which are 
    #        stored in columns 1 & 2 of dffdt_Bvals.
    #
    #       dffdt_Bvals = d x 2 array storing "values" of first derivative of
    #                       ODE driver function at nozzle entrance [exit] in 
    #                       column 1 [2].



    if Shock_Flag .== 0
        dffdt_Bvals = CalcdfdtBvalsmSW(InitVal,Gamma,ff_Bvals,dffdt_vals,d,
                                        Tot_Int_Pts)
    elseif Shock_Flag .== 1
        dffdt_Bvals = CalcdfdtBvalspSW(InitVal,Gamma,ff_Bvals,dffdt_vals,d,
                                        Tot_Int_Pts)
    else()
        disp(["Unknown Shock_Flag value: " int2str(Shock_Flag)])
    end

    #   3.b calculate second time derivative of flow fluxes at all grid-points
    #       
    #       d2Fdt2 = d x Tot_X_Pts array storing the flow flux second time
    #                               derivatives at all grid-points.

    d2Fdt2 = Calc_d2Fdt2(InitVal,dffdt_vals,dffdt_Bvals,ff_vals,ff_Bvals,
                            Gamma,d,Tot_X_Pts)

    #   3.c calculate second time derivative of flow source term at all()
    #           interior grid-points
    #
    #       d2Jdt2 = 1 x Tot_Int_Pts array storing the flow source term
    #                               second time derivative at all interior
    #                               grid-points

    d2Jdt2 = zeros(Float64, 1, Tot_Int_Pts) #TODO: change back later #Calc_d2Jdt2(InitVal,dffdt_vals,ff_vals,A,Gamma,d,Tot_Int_Pts)

    #   3.d calculate d2ffdt2_vals = d x Tot_Int_Pts array which() 
    #       stores values of second time derivative of ff at all interior 
    #       grid-points.

    d2ffdt2_vals = Calc_d2ffdt2(d2Fdt2, d2Jdt2, Del_x, d, Tot_Int_Pts, params)

    #   3.e Now store d2ffdt2_vals in column 3 of d x rmax x Tot_Int_Pts array 
    #       ff

    for ll = 1:Tot_Int_Pts
        for k = 1:d
            ff[k,3,ll] = d2ffdt2_vals[k,ll]
        end
    end 
#end

return ff;

end
