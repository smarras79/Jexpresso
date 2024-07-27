function CalcBCmSW(U,A,Gamma,d,Tot_X_Pts, params)
    #CALCBCMSW evaluates flow variables at nozzle boundaries-shock-wave absent
    #   CalcBCmSW evaluates flow variables at boundaries of a nozzle. Flow is()
    #       1D compressible; inviscid flow with no shock-wave present.
    #   
    #   INPUTS:
    #           U = d x Tot_X_Pts array with current values of flow variables
    #                   at interior grid-points only. Values at boundary
    #                   grid-points are zero - will be updated in calling
    #                   program.
    #           A = 1 x Tot_X_Pts array with values of nozzle area at
    #                   grid-points.
    #           Gamma = ratio of specific heats
    #           d = number of flow variables at a grid-point
    #           Tot_X_Pts = number of grid-points
    #
    #   OUTPUT:
    #           U_Bvals = d x 2 array with values of flow variable at nozzle
    #                       entrance & exit.
    U_Bvals = zeros(d,2)
    
    # fac1 = 1/(Gamma -1)
    # fac2 = Gamma/2
    
    # GP_Nm1 = Tot_X_Pts -1
    # GP_Nm2 = Tot_X_Pts -2
    
    # # evaluate U at nozzle entrance:
    
    # U_Bvals[1,1] = A[1]
    # U_Bvals[2,1] = 2*U[2,2] - U[2,3]
    
    # fac3 = (U_Bvals[2,1]/U_Bvals[1,1])^2
    
    # U_Bvals[3,1] = U_Bvals[1,1]*(fac1 + (fac2*fac3))
    
    # # evaluate flow at nozzle exit:
    
    # for k = 1:d
    #     U_Bvals[k,2] = 2*U[k,GP_Nm1] - U[k,GP_Nm2]
    # end
    u = zeros(Float64, Tot_X_Pts*d)
    for ip=1:Tot_X_Pts
        for ieq=0:d-1
            u[Tot_X_Pts*ieq + ip] = U[ieq+1, ip]
        end
    end

    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    apply_boundary_conditions!(u, params.uaux, 0.0, params.qp.qe,
    params.mesh, params.metrics, params.basis,
    params.RHS, params.rhs_el, params.ubdy,
    params.ω, params.neqs, params.inputs, params.inputs[:AD], NSD_1D())

    for l=0:d-1
        U_Bvals[l+1, 1] = u[Tot_X_Pts*l + 1]
        U_Bvals[l+1, 2] = u[Tot_X_Pts*(l+1)]
    end
    
    return U_Bvals
    end
    
    