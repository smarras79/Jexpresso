function fOrc(t, Start, TCoeffs, d, 
              rmaxp1, Tot_Int_Pts, Gamma, Del_x, A, 
              Shock_Flag, Exit_Pressure, params)
    #FORC evaluates ODE driver function f at l^{s}_[i](t) at interior grd-pts
    #
    #   fOrc evaluates the ODE driver function at l^{s}_[i](t) at the knot time
    #       t in subsubinterval j for each interior grid-point. ODE driver 
    #       function f is for 1D inviscid; compressible flow through a nozzle.
    #
    #   INPUTS: t = knot time value
    #           Start = starting time for subsubinterval j
    #           TCoeffs = d x rmaxp1 x Tot_Int_Pts array with Taylor 
    #                       Polynomial coefficients for l^{s}_[i](t) in 
    #                       subsubinterval j at each interior grid-point
    #           d = number of components of f & l^{s}_[i](t)
    #           rmaxp1 = number of terms/coefficients in a Taylor polynomial
    #           Tot_Int_Pts = number of interior grid-points
    #           Gamma = ratio of specific heats
    #           Del_x = distance between grid-points
    #           A = 1 x Tot_X_Pts array storing nozzle area at all grid-points
    #           Shock_Flag = 0 [1] if shock wave absent [present]
    #           Exit_Pressure = pressure at nozzle exit()
    #
    #   OUTPUT:
    #           f_Loc = d x Tot_Int_Pts array storing d components of ODE 
    #                   driver function f(U) at each interior grid-point at  
    #                   current value of primary flow variable U
    #
    #   Support functions: CalcBCmSW; CalcBCpSW; Calcf0

    # initialize parameter & array

    Tot_X_Pts = Tot_Int_Pts + 2; # number of grid-points
    #@info t Start
    delt = t - Start; #added +1 for bounds error, might have to change   # elapsed time from start of subsubinterval j
    lt = zeros(d,Tot_Int_Pts);   # stores l^{s}_[i](t) at each interior 
    #     grid-point

    PolyArr = zeros(1, rmaxp1);  # initialize to zero array storing Taylor 
    #    polynomial coefficients for l[t[j-1], i] 
    #    for given component & interior grid-point
    
    U = zeros(d, Tot_X_Pts);  # array to store primary flow variables

    # evaluate l^{s}_[i](t) at each interior grid-point, one component at time

    for ll = 1:Tot_Int_Pts    
        for m = 1: d           
            # load Taylor polynomial coefficients for m-th component at
            #   interior grid-point ll in array Poly

            for column = 1:rmaxp1
                #@info TCoeffs[m, column, ll]
                PolyArr[column] = TCoeffs[m, column, ll]
                #@info "poly col: " PolyArr[column]
            end
            #@info "poly is: " PolyArr
            newPoly = reverse(PolyArr)
            #@info "new poly is: " newPoly
            newPolyPoly = Polynomial(vec(newPoly))
            #@info "new poly polynomial ver is: " newPolyPoly

            # store value of m-th Taylor polynomial at elapsed time delt & 
            #       at interior grid-point ll

            #@info delt
            lt[m,ll] = newPolyPoly(delt)#Poly[Int(delt)+1]; 
            #@info lt[m, ll]
            #readline();
        end
    end

    # assign U at each interior grid-point

    for ll = 2: (Tot_X_Pts - 1)
        IP_Label = ll - 1

        for m = 1:d
            U[m,ll] = lt[m,IP_Label]
        end
    end

    # assign U at boundary points using flow boundary conditions

    # from here
    # if (Shock_Flag .== 0)
    #     U_Bvals = CalcBCmSW(U,A,Gamma,d,Tot_X_Pts)
    #     #@info U_Bvals
    #     #readline();
    # elseif (Shock_Flag .== 1)
    #     U_Bvals = CalcBCpSW(U,A,Gamma,d,Tot_X_Pts,Exit_Pressure)
    #     #@info U_Bvals
    #     #readline();
    # else
    #     disp(["Unknown Shock_Flag value: " int2str(Shock_Flag)])
    # end

    # for m = 1:d
    #     U[m,1] = U_Bvals[m,1]
    #     U[m,Tot_X_Pts] = U_Bvals[m,2]
    # end
    # # to here, replace with jexpresso bc's
    u = zeros(Float64, params.mesh.npoin*d)
    for ip=1:Tot_X_Pts
        for ieq=0:d-1
            u[Tot_X_Pts*ieq + ip] = U[ieq+1, ip]
        end
    end

    # replace with call to _build_rhs
    build_rhs!(params.RHS, u, params, 0.0)

    # for m=1:d
    #     U[m,1] = params.uaux[1,m]
    #     U[m,Tot_X_Pts] = params.uaux[Tot_X_Pts, m]
    # end

    rhs = zeros(Float64, d, Tot_Int_Pts)
    for i=1:Tot_Int_Pts
        for j=1:d
            rhs[j, i] = params.RHS[i+1, j]
        end
    end
    #@info rhs

    
    #f_Loc = rhs
    # evaluate f using Calcf0
    #f_Loc = Calcf0(d, Tot_X_Pts, Tot_Int_Pts, Gamma, Del_x, U, A, params)

    return rhs
end
