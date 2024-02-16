function InitParms( n, delta, err1, Tot_TSteps )
    #INITPARMS initializes parameters for quantum ODE algorithm
    #   InitParms initializes values for recursion level k; number of
    #   sub-subintervals N; & probability delta1 that approximate value 
    #   of integral violates its error upper bound. Fixes n so that quantum
    #   ODE algorithm time-step is less than the Courant_Friedrichs-Levy [CFL]
    #   time-step needed for numerical stability of Navier-Stokes simulation.
    #
    #   INPUTS: 
    #           n = initial number of primary subintervals
    #           delta = probability that quantum ODE algorithm solution 
    #                       violates its error upper bound 
    #           err1 = error upper bound in approximate value of integral 
    #                       output by quantum integration algorithm
    #           Tot_TSteps = total number of time-steps used in a MacCormack's
    #                           method solution of Nozzle flow problem. I
    #                           use it with Delta_t to set final time b in 
    #                           quantum ODE algorithm.
    #
    #   OUTPUTS:
    #
    #           N = number of subsubintervals
    #           delta1 = probability that approximate value of integral 
    #                       violates its error upper bound
    #           final_n = number of subintervals that insures hbar is()
    #                       less than CFL upper bound on time-step
    #
    #   REMARK: Fix final_n by requiring number of Taylor steps exceed
    #             Tot_TSteps - insures hbar less than CFL time-step Delta_t.


    # 3rd page of review contains these equations
    k = 1+ ceil(log(1/err1)/log(n))

    N = n^(k-1)
    
    Tylr_Steps = n*N
    
    ratio = Tot_TSteps/Tylr_Steps

    while ratio > 1
        n = n + 1
        
        k = 1 + ceil(log(1/err1)/log(n))
        N = n^(k-1)
        
        Tylr_Steps = n*N
        ratio = Tot_TSteps/Tylr_Steps
    end

    temp_n = (Tot_TSteps/ratio)^(1/k)

    final_n = ceil(temp_n)

    k = 1 + ceil(log(1/err1)/log(final_n))

    N = (final_n)^(k-1)

    delta1 = 1 - (1-delta)^(1/N)

    return  N, delta1, final_n, k 
end