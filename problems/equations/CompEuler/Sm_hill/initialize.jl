function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
    Initialize fields for 2D mountain flow with microphysics
    """
    @info " Initialize fields for 2D CompEuler with microphysics ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #---------------------------------------------------------------------------------
    qvars = ("ρ", "ρu", "ρv", "hl", "ρqt", "ρqp")
    qoutvars = ["ρ", "ρu", "ρv", "hl", "ρqt", "ρqp", "T", "qn", "qc", "qi", "qr", "qs", "qg"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------
    
    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()
        
        if inputs[:lrestart] == true
            #
            # READ RESTART HDF5:
            #
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars))
            
            # Recalculate pressure from stored state
            for ip=1:mesh.npoin
                ρ  = q.qn[ip,1]
                ρhl = q.qn[ip,4]
                ρqt = q.qn[ip,5]
                ρqp = q.qn[ip,6]
                hl = ρhl/ρ
                qt = ρqt/ρ
                qp = ρqp/ρ
                P = perfectGasLaw_ρhltoP(PhysConst, ρ=ρ, hl=hl, qt=qt, qp=qp)
                q.qn[ip,end] = P
            
                ρe  = q.qe[ip,1]
                ρhle = q.qe[ip,4]
                ρqte = q.qe[ip,5]
                ρqpe = q.qe[ip,6]
                hle = ρhle/ρe
                qte = ρqte/ρe
                qpe = ρqpe/ρe
                Pe = perfectGasLaw_ρhltoP(PhysConst, ρ=ρe, hl=hle, qt=qte, qp=qpe)
                q.qe[ip,end] = Pe
            end
        
        else
            #
            # INITIAL STATE from scratch: Pure mountain flow (no thermal perturbation)
            #
            data = read_sounding(inputs[:sounding_file])
            background = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.y, data) 
            
            for ip = 1:mesh.npoin
                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]
                
                # Extract reference state from sounding
                θ_ref = background[ip,1]
                qv_ref = background[ip,2]/1000  # Convert g/kg to kg/kg
                u_ref = background[ip,3]
                v_ref = background[ip,4]
                pref = background[ip,5]
                
                # Calculate reference state
                Tref = θ_ref / (PhysConst.pref/pref)^(PhysConst.Rair/PhysConst.cp)
                Tv_ref = Tref*(1 + 0.61*qv_ref)
                ρref = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv_ref, Press=pref)
                hl_ref = PhysConst.cp*Tref + PhysConst.g*y
                
                # Initial state = reference state (no perturbations for mountain case)
                u = u_ref
                v = 0.0  # No vertical velocity initially
                
                # Moist pressure
                pref_m = ρref*Tv_ref*PhysConst.Rair
            
                if inputs[:SOL_VARS_TYPE] == PERT()
                    # Perturbation formulation: q = q_total - q_reference
                    q.qn[ip,1] = 0.0  # ρ - ρref (no density perturbation initially)
                    q.qn[ip,2] = 0.0  # ρu - ρref*u_ref (no momentum perturbation)
                    q.qn[ip,3] = 0.0  # ρv (no vertical momentum)
                    q.qn[ip,4] = 0.0  # ρhl - ρref*hl_ref (no enthalpy perturbation)
                    q.qn[ip,5] = 0.0  # ρqt - ρref*qv_ref (no moisture perturbation)
                    q.qn[ip,6] = 0.0  # ρqp (no precipitation)
                    q.qn[ip,end] = pref_m
                    
                    # Store reference state
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = ρref*u_ref
                    q.qe[ip,3] = 0.0
                    q.qe[ip,4] = ρref*hl_ref
                    q.qe[ip,5] = ρref*qv_ref
                    q.qe[ip,6] = 0.0
                    q.qe[ip,end] = pref_m
                else
                    # Total formulation: q = q_total
                    q.qn[ip,1] = ρref
                    q.qn[ip,2] = ρref*u
                    q.qn[ip,3] = 0.0  # ρv = 0
                    q.qn[ip,4] = ρref*hl_ref
                    q.qn[ip,5] = ρref*qv_ref
                    q.qn[ip,6] = 0.0  # No precipitation initially
                    q.qn[ip,end] = pref_m
                    
                    # Store reference state for analysis
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = ρref*u_ref
                    q.qe[ip,3] = 0.0
                    q.qe[ip,4] = ρref*hl_ref
                    q.qe[ip,5] = ρref*qv_ref
                    q.qe[ip,6] = 0.0
                    q.qe[ip,end] = pref_m
                end
            end
        end
    
        # Non-conservative form adjustments
        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[:,2] .= q.qn[:,2]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,3] .= q.qn[:,3]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,4] .= q.qn[:,4]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,5] .= q.qn[:,5]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,6] .= q.qn[:,6]./(q.qn[:,1] + q.qe[:,1])
                
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
                q.qe[:,6] .= q.qe[:,6]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]
                q.qn[:,5] .= q.qn[:,5]./q.qn[:,1]
                q.qn[:,6] .= q.qn[:,6]./q.qn[:,1]
                
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
                q.qe[:,6] .= q.qe[:,6]./q.qe[:,1]
            end
        end

    else
        # GPU initialization
        if (inputs[:SOL_VARS_TYPE] == PERT())
            lpert = true
        else
            lpert = false
        end
        data = read_sounding(inputs[:sounding_file])
        background = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.y, data)
        PhysConst = PhysicalConst{TFloat}()
        
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, background, mesh.x, mesh.y, PhysConst, lpert; ndrange = (mesh.npoin))
    end
    
    @info "Pressure range: ", maximum(q.qe[:,end]), minimum(q.qe[:,end])
    @info " Initialize fields for 2D CompEuler with microphysics ........................ DONE "
    return q
end

@kernel function initialize_gpu!(qn, qe, background, x, y, PhysConst, lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    x1 = x[ip]
    y1 = y[ip]

    # Extract reference state from sounding
    θ_ref = background[ip,1]
    qv_ref = T(background[ip,2]/T(1000))  # g/kg to kg/kg
    u_ref = background[ip,3]
    v_ref = background[ip,4]
    pref = background[ip,5]
    
    # Calculate reference state
    Tref = T(θ_ref / (PhysConst.pref/pref)^(PhysConst.Rair/PhysConst.cp))
    Tv_ref = Tref*(T(1) + T(0.61)*qv_ref)
    ρref = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv_ref, Press=pref)
    hl_ref = PhysConst.cp*Tref + PhysConst.g*y1
    
    # Initial velocities
    u = u_ref
    v = T(0.0)
    
    # Moist pressure
    pref_m = ρref*Tv_ref*PhysConst.Rair

    if (lpert)
        # Perturbation formulation
        qn[ip,1] = T(0.0)  # No density perturbation
        qn[ip,2] = T(0.0)  # No u-momentum perturbation
        qn[ip,3] = T(0.0)  # No v-momentum
        qn[ip,4] = T(0.0)  # No enthalpy perturbation
        qn[ip,5] = T(0.0)  # No moisture perturbation
        qn[ip,6] = T(0.0)  # No precipitation
        qn[ip,end] = pref_m
    else
        # Total formulation
        qn[ip,1] = ρref
        qn[ip,2] = ρref*u
        qn[ip,3] = T(0.0)  # ρv = 0
        qn[ip,4] = ρref*hl_ref
        qn[ip,5] = ρref*qv_ref
        qn[ip,6] = T(0.0)  # No precipitation
        qn[ip,end] = pref_m
    end

    # Store reference state
    qe[ip,1] = ρref
    qe[ip,2] = ρref*u
    qe[ip,3] = T(0.0)
    qe[ip,4] = ρref*hl_ref
    qe[ip,5] = ρref*qv_ref
    qe[ip,6] = T(0.0)
    qe[ip,end] = pref_m
end