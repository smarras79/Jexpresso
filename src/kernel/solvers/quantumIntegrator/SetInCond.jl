function SetInCond(Shock_Flag, In_Mass_Flow_E, Gamma, x, Del_x, A, d, Mrho_E, Temp_E, ICMFlowErrScale, ICrhoErrScale, ICtempErrScale)
    #SETINCOND assigns initial condition for flow variables U
    #   
    #   INPUT: 
    #           Shock_Flag = 0 [1] if shock wave absent [present]
    #           In_Mass_Flow_E = mass flow value for exact solution of 
    #                               subsonic to supersonic flow w/o shock-wave 
    #           Gamma = ratio of specific heats C_p/C_v
    #           x = 1 x Tot_Pts array with x grid-point locations
    #           Del_x = spatial separation between grid-points
    #           A = 1 x Tot_Pts array with nozzle area at grid-points
    #           d = number of flow variables at each grid-point
    #           Mrho_E = 1 x Tot_X_Pts array storing exact values of mass
    #                       density at each grid-point
    #           Temp_E = 1 x Tot_X_Pts array storing exact values of 
    #                       temperature at each grid-point
    #           ICMFlowErrScale = scale of random shift to initial mass flow() 
    #                               away from exact results 
    #           ICrhoErrScale = scale of random shift to initial mass density
    #                               away from exact results (no shift at nozzle  
    #                               throat | inlet)
    #           ICtempErrScale = scale of random shift to initial temperature
    #                               away from exact results (no shift at nozzle  
    #                               throat | inlet) 

    #
    #   OUTPUT:
    #           U = d x Tot_Pts array of flow variables
    #           Del_t = time-step size (determined by CFL condition)
    #           In_Mass_Flow = initial mass flow rate used in simulation:
    #                               includes random shift

    # determine number of grid-points
    Tot_Pts = length(x)

    # determine coefficients & set Tiny to small value for later use
    Coef1 = 1 /(Gamma -1)
    Coef2 = Gamma/2

    ithroat = (Tot_Pts + 1)/2

    # initialize arrays
    Mrho = zeros(1,Tot_Pts);    # will store shifted mass density
    Temp = zeros(1,Tot_Pts);    # will store shifted temperature
    U = zeros(d,Tot_Pts);       # will store flow variables
    V = zeros(1,Tot_Pts);       # will store velocity
    Loc_TSteps = zeros(1,(Tot_Pts - 2));    # will store local time-step

    # add random shift to In_Mass_Flow_E
    In_Mass_Flow = In_Mass_Flow_E * (1 + (-ICMFlowErrScale + (2 * ICMFlowErrScale * rand() )))

    # begin loop over grid-points - initialize mass density & temperature

    # Introduce a small shift
    #Tiny = 10^(-8)

    for i = 1:Tot_Pts
        if Shock_Flag .== 0
            if ( (i == 1) || (i .== ithroat) )
                Mrho[i] = Mrho_E[i]
                Temp[i] = Temp_E[i]
            elseif ( (i != 1) && (i != ithroat) )
                Mrho[i] = Mrho_E[i]*( 1 + (-ICrhoErrScale + 2*ICrhoErrScale*rand()) )
                
                if Mrho[i] .> 1
                    Mrho[i] = 1
                end
                
                Temp[i] = Temp_E[i]*( 1 + (-ICtempErrScale + 2*ICtempErrScale*rand()) )
                
                if Temp[i] .> 1
                    Temp[i] = 1
                end
            end
        elseif Shock_Flag .== 1
            #= 
            Anderson's initial condition
            
            if ( x[i] .< (0.5 + Tiny) )
                Mrho[i] = 1.0
                Temp[i] = 1.0
            elseif [ x[i] .< (x[ithroat] + Tiny) ]
                Mrho[i] = 1.0 - 0.366*(x[i] - 0.5)
                Temp[i] = 1.0 - 0.167*(x[i] - 0.5)
            elseif [ x[i] .< (2.1 + Tiny) ]
                Mrho[i] = 0.634 - 0.702*(x[i] - 1.5)
                Temp[i] = 0.833 - 0.4908*(x[i] - 1.5)
            elseif [ x[i] .> (2.1 + Tiny) ]
                Mrho[i] = 0.5892 + 0.10228*(x[i] - 2.1)
                Temp[i] = 0.93968 + 0.0622*(x[i] - 2.1)
            end
            =#
            if ( (i == 1) || (i .== ithroat) )
                Mrho[i] = Mrho_E[i]
                Temp[i] = Temp_E[i]
            elseif ( (i != 1) && (i != ithroat) )
                Mrho[i] = Mrho_E[i]*( 1 + (-ICrhoErrScale + 2*ICrhoErrScale*rand()) )
                
                if Mrho[i] .> 1
                    Mrho[i] = 1
                end
                
                Temp[i] = Temp_E[i]*( 1 + (-ICtempErrScale + 2*ICtempErrScale*rand()) )
                
                if Temp[i] .> 1
                    Temp[i] = 1
                end
            end
        end

        # assign initial condition to flow variables
        U[1,i] = Mrho[i]*A[i]
        U[2,i] = In_Mass_Flow
        U[3,i] = U[1,i]*( Coef1*Temp[i] + Coef2*(U[2,i]/U[1,i])^(2) )
    end

    # determine time-step Del_t using Courant-Friedrichs-Levy [CFL] stability
    #   condition with C = 0.5

    

    # Del_t is the min value in Loc_TSteps

    Delta_t = minimum(Loc_TSteps)

    return U, Delta_t, In_Mass_Flow
end