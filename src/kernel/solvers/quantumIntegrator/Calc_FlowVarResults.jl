function Calc_FlowVarResults(Gamma, Tot_X_Pts, A, InitVal)
    #CALC_FLOWVARRESULTS finds physical flow variables from simulation results
    #   Calc_FlowVarResults finds Mach number & dimensionless mass density
    #       pressure; temperature; & velocity at end of subinterval i at 
    #       all grid-points using simulation results.
    #
    #   INPUTS:
    #           Gamma = ratio of specific heats
    #           Tot_X_Pts = number of grid-points
    #           A = nozzle area at each grid-point
    #           InitVal = d x Tot_X_Pts array storing simulation flow variables
    #                       at end of subinterval i [here d = 3]
    #
    #   OUTPUTS:
    #           Mach_D = 1 x Tot_X_Pts array storing flow Mach number at each
    #                       grid-point
    #           Mrho_D = 1 x Tot_X_Pts array storing dimensionless mass density
    #                       at each grid-point
    #           Press_D = 1 x Tot_X_Pts array storing dimensionless pressure
    #                       at each grid-point
    #           Temp_D = 1 x Tot_X_Pts array storing dimensionless temperature
    #                       at each grid-point
    #           Vel_D = 1 x Tot_X_Pts array storing dimensionless velocity
    #                       at each grid-point

    # declare & initialize arrays to be returned

    Mach_D = zeros(1, Tot_X_Pts)
    Mrho_D = zeros(1, Tot_X_Pts)
    Press_D = zeros(1, Tot_X_Pts)
    Temp_D = zeros(1, Tot_X_Pts)
    Vel_D = zeros(1, Tot_X_Pts)

    # define useful factors

    fac1 = Gamma - 1
    fac2 = Gamma / 2

    # loop over all grid-points

    for i = 1:Tot_X_Pts
        Mrho_D[i] = InitVal[1, i] / A[i]
        Vel_D[i] = InitVal[2, i] / InitVal[1, i]
        Temp_D[i] = fac1 * ((InitVal[3, i] / InitVal[1, i])
                            -
                            fac2 * Vel_D[i] * Vel_D[i])
        Press_D[i] = Mrho_D[i] * Temp_D[i]
        Mach_D[i] = Vel_D[i] / sqrt(Temp_D[i])
    end

    return [Mach_D, Mrho_D, Press_D, Temp_D, Vel_D]
end
