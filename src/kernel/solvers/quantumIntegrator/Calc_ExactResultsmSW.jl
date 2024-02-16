using Polynomials

function Calc_ExactResultsmSW(Gamma, Tot_X_Pts, A)
    #CALC_EXACTRESULTSMSW assigns exact flow values at grd-pts - no shock-wave
    #   Calc_ExactResultsmSW assigns exact values of the flow variables at
    #       all grid-points for 1D inviscid compressible flow through a 
    #       symmetric nozzle when no shock-wave is present. All flow variables
    #       are dimensionless as in rest of simulation code. Formulas on which()
    #       this function is base appear in Anderson; "Computational Fluid 
    #       Dynamics"; Chapter 7; section 2.
    #
    #   INPUTS:
    #           Gamma = ratio of specific heats
    #           Tot_X_Pts = number of spatial grid-points
    #           A = 1 x Tot_X_Pts array storing nozzle area at all grid-points
    #
    #   OUTPUTS:
    #           M = 1 x Tot_X_Pts array storing Mach number at all grid-points
    #           Mrho = 1 x Tot_X_Pts array storing mass density at all grd-pts
    #           Press = 1 x Tot_X_Pts array storing pressure at all grd-pts
    #           Temp = 1 x Tot_X_Pts array storing temperature at all grd-pts
    #           Vel = 1 x Tot_X_Pts array storing velocity at all grd-pts

    # assign values for useful parameters
    fac1 = (Gamma - 1)/2
 
    ithroat = Int((Tot_X_Pts + 1)/2);   # grid-point index for nozzle throat
    istart = ithroat + 1;   # grid-point index for start of supersonic region
    istop = Tot_X_Pts;       # grid-point index for end of supersonic region

    # will use expon as a loop parameter so must be an integer
    expon = Int(round((Gamma + 1)/(Gamma - 1)))
    expo2 = -1/(Gamma - 1)
    expo1 = Gamma*expo2

    # Initialize array to store return arrays
    M = zeros(1, Tot_X_Pts)
    Mrho = zeros(1, Tot_X_Pts)
    Press = zeros(1, Tot_X_Pts)
    Temp = zeros(1, Tot_X_Pts)
    Vel = zeros(1, Tot_X_Pts)

    # throat is always at Mach 1 so ...
    M[ithroat] = 1

    # evaluate other flow variables at throat
    argum = 1 + fac1*(M[ithroat])^(2)
 
    Press[ithroat] = (argum)^(expo1)
    Mrho[ithroat] = (argum)^(expo2)
    Temp[ithroat] = 1/argum
    Vel[ithroat] = M[ithroat]*sqrt(Temp[ithroat])

    # begin loop over grid-points. Note since nozzle area is symmetric about
    # the nozzle throat only need to loop over the supersonic flow region
    for i = istart:istop
        # assign value to coefficient of z in polynomial resulting from
        #  area-Mach number [A-M] relation
        
        c1 = 18750 - 46656*(A[i])^(2)
        
        # specify A-M relation polynomial
        
        #p = Polynomial([1, 30, 375, 2500, 9375, c1, 15625])
        p = Polynomial([15625, c1, 9375, 2500, 375, 30, 1])
        # evaluate roots of p
        
        z = Polynomials.roots(p)
        # loop over roots: if a root is real & greater than 1 then assign
        #      its square root to M[i]; if root real & less than 1, assign
        #      its square root to M[isubson], where isubson = 2*ithroat - i
        
        for j = 1:expon
            if isreal(z[j])
                if real(z[j]) > 1
                    M[i] = sqrt(z[j])
                    # calculate other flow variables at i
                    
                    argum = 1 + fac1*z[j]
                    
                    Press[i] = (argum)^(expo1)
                    Mrho[i] = (argum)^(expo2)
                    Temp[i] = 1/argum
                    Vel[i] = M[i]*sqrt(Temp[i])
                elseif real(z[j]) < 1
                    isubson = 2*ithroat - i
                    
                    M[isubson] = sqrt(z[j])
                    
                    argum = 1 + fac1*z[j]
                    
                    Press[isubson] = (argum)^(expo1)
                    Mrho[isubson] = (argum)^(expo2)
                    Temp[isubson] = 1/argum
                    Vel[isubson] = M[isubson]*sqrt(Temp[isubson])
                end
            end
        end
    end
    
    return M, Mrho, Press, Temp, Vel
end

#Calc_ExactResultsmSW(0.2, 3, [28, 45.5, 67.5])