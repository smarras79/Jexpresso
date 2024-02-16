include("CalcBCmSW.jl");
include("CalcBCpSW.jl");

using Polynomials

function NextInCond(mmat, InitVal, hbar, d, r, A, Gamma,
    Tot_Int_Pts, Tot_X_Pts, Shock_Flag,
    Exit_Pressure)
    #NEXTINCOND determines initial condition for next subsubinterval 
    #   NextInCond evaluates the Taylor polynomials for a given subsubinterval
    #   at the subsubinterval's greatest value & uses this result as the
    #   initial value for the next subsubinterval.
    #
    #   INPUTS: mmat = d x [r+2] x Tot_Int_Pts array of Taylor polynomial 
    #                   coefficients for each driver function component &
    #                   interior grid-point
    #           InitVal = d x Tot_X_Pts array with initial values of primary
    #                       flow variables for current subsubinterval
    #           hbar = width of subsubinterval
    #           d = number of components of InitVal
    #           r + 1 = degree of Taylor polynomials
    #           A = 1 x Tot_X_Pts array storing nozzle area at all grid-points
    #           Gamma = ratio of specific heats
    #           Tot_Int_Pts = number of interior grid-points
    #           Tot_X_Pts = number of grid-points
    #           Shock_Flag = 0 [1] if shock wave present [absent] in divergent
    #                           part of nozzle
    #           Exit_Pressure = pressure at nozzle exit()
    #
    #   OUTPUT:
    #           NextInVal = d x Tot_X_Pts array with initial condition values 
    #                           for primary flow variables for next
    #                           subsubinterval
    #
    # Support functions: CalcBCmSW; CalcBCpSW

    rmax = r + 1       # degree of Taylor Polynomials
    rmaxp1 = rmax + 1  # number of terms in Taylor polynomials

    TylrPoly = zeros(1, rmaxp1) # array to store coefficients for a Taylor
    #   polynomial

    NextInVal = zeros(d, Tot_X_Pts)     # array to store initial conditions 
    #   for next subsubinterval 

    # calculate initial condition for next subsubinterval for the
    #       interior grid-points    [ 2 <= ll <= TotXPtm1 ]

    for ll = 1:Tot_Int_Pts
        GPLabel = ll + 1

        for pp = 1:d
            for mm = 1:rmax
                TylrPoly[mm] = mmat[pp, mm, ll]
            end

            TylrPoly[rmaxp1] = InitVal[pp, GPLabel]

            TylrPoly = reverse(TylrPoly);

            TylrPolyPoly = Polynomial(vec(TylrPoly))

            NextInVal[pp, GPLabel] = TylrPolyPoly(hbar)#Polynomials.polyval(TylrPoly,hbar)
        end
    end

    # use flow boundary conditions to determine new initial conditions for
    #       boundary grid-points

    if (Shock_Flag .== 0)
        U_Bvals = CalcBCmSW(InitVal, A, Gamma, d, Tot_X_Pts)
    elseif (Shock_Flag .== 1)
        U_Bvals = CalcBCpSW(InitVal, A, Gamma, d, Tot_X_Pts, Exit_Pressure)
    else
        disp([" Unknown Shock_Flag value: " Shock_Flag])
    end

    # store boundary values in NextInVal

    for p = 1:d
        NextInVal[p, 1] = U_Bvals[p, 1]      # values at nozzle entrance
        NextInVal[p, Tot_X_Pts] = U_Bvals[p, 2]  # values at nozzle exit()
    end

    return NextInVal

end
