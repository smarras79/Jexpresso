function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize analytic neutral ABL ........................ "
    end
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρw", "ρθ"]
    qoutvars = ["ρ", "u", "v", "w", "θ", "p"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------

    if inputs[:lrestart] == true
        #
        # READ RESTART HDF5:
        #
        q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars))

        PhysConst = PhysicalConst{Float64}()
        
        
        for ip=1:mesh.npoin
            ρ  = q.qn[ip,1]
            ρθ = q.qn[ip,5]
            θ  = ρθ/ρ
            P = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
            q.qn[ip,end] = P
        
            ρe  = q.qe[ip,1]
            ρθe = q.qe[ip,5]
            θe  = ρθe/ρ
            Pe = perfectGasLaw_ρθtoP(PhysConst, ρ=ρe, θ=θe)
            q.qe[ip,end] = Pe
        end
        
    else
        lanalytic = false
        if lanalytic == true
            if (inputs[:backend] == CPU())    
                PhysConst = PhysicalConst{Float64}()
                
                #
                # INITIAL STATE from scratch:
                #    
                zi  = 840.0 #m
                θ1  = 289.0 #K
                θ2  = 297.5 #K
                Ts  = 290.4
                H   = PhysConst.Rair*Ts/PhysConst.g
                amp = 0.25
                for ip = 1:mesh.npoin
                    
                    z = mesh.z[ip]

                    if z <= zi
                        θ = θ1
                    else
                        θ = θ2 + cbrt(z - zi)
                    end
                    p    = PhysConst.pref*exp(-z/H)
                    pref = p

                    randnoise = 0.0
                    if z < 800.0
                        randnoise = 2*amp*(rand() - 1.0)
                    end
                    θ    = θ + randnoise
                    
                    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)    #kg/m³
                    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=pref) #kg/m³

                    u = 10.0
                    v = 0.0
                    w = 0.0

                    if inputs[:SOL_VARS_TYPE] == PERT()
                        q.qn[ip,1] = ρ - ρref
                        q.qn[ip,2] = ρ*u - ρref*u
                        q.qn[ip,3] = ρ*v - ρref*v
                        q.qn[ip,4] = ρ*w - ρref*w
                        q.qn[ip,5] = ρ*θ - ρref*θ
                        q.qn[ip,end] = p
                        
                        #Store initial background state for plotting and analysis of pertuebations
                        q.qe[ip,1] = 0.0*ρref
                        q.qe[ip,2] = 0.0*u
                        q.qe[ip,3] = 0.0*v
                        q.qe[ip,4] = 0.0*w
                        q.qe[ip,5] = 0.0*ρref*θ
                        q.qe[ip,end] = 0.0*pref
                    else
                        q.qn[ip,1] = ρ
                        q.qn[ip,2] = ρ*u
                        q.qn[ip,3] = ρ*v
                        q.qn[ip,4] = ρ*w
                        q.qn[ip,5] = ρ*θ
                        q.qn[ip,end] = p

                        #Store initial background state for plotting and analysis of pertuebations
                        q.qe[ip,1] = q.qn[ip,1]
                        q.qe[ip,2] = q.qn[ip,2]
                        q.qe[ip,3] = q.qn[ip,3]
                        q.qe[ip,4] = q.qn[ip,4]
                        q.qe[ip,5] = q.qn[ip,5]
                        q.qe[ip,end] = q.qn[ip,end]
                    end
                end
            end    
        else
            PhysConst = PhysicalConst{Float64}()
            
            #
            # INITIAL STATE from scratch:
            #
            data = read_sounding(inputs[:sounding_file])
            
            sounding_nlevels = size(data)[1]
            sounding_nvars   = size(data)[2]
            data_with_p      = zeros(Float64, sounding_nlevels, sounding_nvars+1)
            
            # Get surface pressure from sounding header (first line: pressure, theta, qv)
            # This should be read from the surface values line in your sounding file
            p_surface = 100000.0  # 1000 hPa converted to Pa, adjust based on your sounding reader
            
            # Copy original data (height, theta, qv, u, v)
            data_with_p[:,1:sounding_nvars] = copy(data[:,:])
            
            # Calculate pressure from hydrostatic equilibrium with iterative balancing
            # Assuming data structure: [height, theta, qv, u, v]
            
            # Initialize pressure array
            pressure = zeros(Float64, sounding_nlevels)

            # Initial pressure calculation
            for i = 1:sounding_nlevels
                z_current  = data[i,1] # height (m)
                θ_current  = data[i,2] # potential temperature (K)
                qv_current = data[i,3] # water vapor mixing ratio (g/kg)
                
                # Convert mixing ratio from g/kg to kg/kg
                qv_kg_per_kg = qv_current / 1000.0
                
                if i == 1
                    # For first level, use hydrostatic relation from surface
                    # Calculate virtual potential temperature
                    # θ_v = θ * (1 + 0.61 * qv) for mixing ratio
                    θ_v = θ_current * (1.0 + 0.61 * qv_kg_per_kg)
                    
                    # Virtual temperature: T_v = θ_v * (p/p0)^(R/cp)
                    T_v = θ_v * (p_surface/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                    
                    # Hydrostatic equation with virtual temperature
                    pressure[i] = p_surface * exp(-PhysConst.g * z_current / (PhysConst.Rair * T_v))
                else
                    # For subsequent levels, integrate hydrostatic equation
                    z_prev  = data[i-1,1]
                    θ_prev  = data[i-1,2]
                    qv_prev = data[i-1,3] / 1000.0  # Convert to kg/kg
                    p_prev  = pressure[i-1]
                    
                    # Average values for integration
                    z_avg = 0.5 * (z_current + z_prev)
                    θ_avg = 0.5 * (θ_current + θ_prev)
                    qv_avg = 0.5 * (qv_kg_per_kg + qv_prev)
                    dz = z_current - z_prev
                    
                    # Virtual potential temperature
                    θ_v_avg = θ_avg * (1.0 + 0.61 * qv_avg)
                    
                    # Virtual temperature for hydrostatic calculation
                    T_v_avg = θ_v_avg * (p_prev/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                    
                    # Hydrostatic equation: dp/dz = -ρg = -pg/(R_v*T_v)
                    # where R_v is the gas constant for moist air
                    pressure[i] = p_prev * exp(-PhysConst.g * dz / (PhysConst.Rair * T_v_avg))
                end
            end

            # Iterative hydrostatic balance correction
            max_iterations = 10
            tolerance      = 1.0  # Pa
            converged      = false
            
            if rank == 0
                @info "Starting iterative hydrostatic balance correction..."
            end

            for iter = 1:max_iterations
                pressure_old   = copy(pressure)
                max_correction = 0.0
                
                # Forward pass: recalculate pressure using updated values
                for i = 2:sounding_nlevels
                    z_current  = data[i,1]
                    z_prev     = data[i-1,1]
                    θ_current  = data[i,2]
                    θ_prev     = data[i-1,2]
                    qv_current = data[i,3]/1000.0
                    qv_prev    = data[i-1,3]/1000.0
                    
                    dz = z_current - z_prev
                    
                    # Use updated pressure from previous level
                    p_prev = pressure[i-1]
                    
                    # Calculate virtual temperature at both levels
                    θ_v_prev = θ_prev * (1.0 + 0.61 * qv_prev)
                    θ_v_current = θ_current * (1.0 + 0.61 * qv_current)
                    
                    T_v_prev = θ_v_prev * (p_prev/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                    
                    # More accurate integration using trapezoidal rule
                    # First estimate for current level
                    p_est = p_prev * exp(-PhysConst.g * dz / (PhysConst.Rair * T_v_prev))
                    
                    # Calculate virtual temperature at current level with estimated pressure
                    T_v_current = θ_v_current * (p_est/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                    
                    # Average virtual temperature for better integration
                    T_v_avg = 0.5 * (T_v_prev + T_v_current)
                    
                    # Refined pressure calculation
                    pressure[i] = p_prev * exp(-PhysConst.g * dz / (PhysConst.Rair * T_v_avg))
                    
                    # Track maximum correction
                    correction     = abs(pressure[i] - pressure_old[i])
                    max_correction = max(max_correction, correction)
                end
                
                # Check for convergence
                if max_correction < tolerance
                    converged = true
                    if rank == 0
                        @info "Hydrostatic balance converged after $iter iterations (max correction: $(max_correction) Pa)"
                    end
                    break
                end
                
                if rank == 0 && iter % 3 == 0
                    @info "Iteration $iter: max pressure correction = $(max_correction) Pa"
                end
            end

            if !converged && rank == 0
                @warn "Hydrostatic balance did not converge after $max_iterations iterations (max correction: $(max_correction) Pa)"
            end
            
            # Verify hydrostatic balance
            max_imbalance = 0.0
            for i = 2:sounding_nlevels
                z_current   = data[i,1]
                z_prev      = data[i-1,1]
                θ_current   = data[i,2]
                θ_prev      = data[i-1,2]
                qv_current  = data[i,3] / 1000.0
                qv_prev     = data[i-1,3] / 1000.0
                
                dz          = z_current - z_prev
                p_current   = pressure[i]
                p_prev      = pressure[i-1]
                
                # Calculate average virtual temperature
                θ_v_prev    = θ_prev * (1.0 + 0.61 * qv_prev)
                θ_v_current = θ_current * (1.0 + 0.61 * qv_current)
                T_v_prev    = θ_v_prev * (p_prev/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                T_v_current = θ_v_current * (p_current/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                T_v_avg     = 0.5 * (T_v_prev + T_v_current)
                
                # Calculate expected pressure from hydrostatic equation
                p_expected  = p_prev * exp(-PhysConst.g * dz / (PhysConst.Rair * T_v_avg))
                
                # Calculate imbalance
                imbalance = abs(p_current - p_expected)
                max_imbalance = max(max_imbalance, imbalance)
            end
            
            if rank == 0
                @info "Final hydrostatic balance verification: max imbalance = $(max_imbalance) Pa"
                relative_imbalance = max_imbalance / p_surface * 100.0
                @info "Relative imbalance: $(relative_imbalance)%"
            end
            
            # Store calculated pressures
            data_with_p[:, sounding_nvars+1] .= pressure
            
            #Interpolate
            data_interpolate = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.z, data_with_p)

            amp = 0.25
            for ip = 1:mesh.npoin
                randnoise = 0.0
                if mesh.z[ip] < 800.0
                    randnoise = 2*amp*(rand() - 1.0)
                end
                θ     = data_interpolate[ip,1] + randnoise  # theta from column 2
                qv    = data_interpolate[ip,2] / 1000.0     # qv from column 3, convert g/kg to kg/kg
                Press = data_interpolate[ip,5]              # pressure from column 6 (newly calculated)
                ρ     = perfectGasLaw_θPtoρ(PhysConst; Press=Press, θ=θ)           
                u     = data_interpolate[ip,3]  # u from column 4
                v     = data_interpolate[ip,4]  # v from column 5
                w     = 0.0
                
                q.qn[ip,1] = ρ
                q.qn[ip,2] = ρ*u
                q.qn[ip,3] = ρ*v
                q.qn[ip,4] = ρ*w
                q.qn[ip,5] = ρ*θ
                q.qn[ip,6] = Press
                
                #Store initial data_interpolate state for plotting and analysis of perturbations
                q.qe[ip,1] = q.qn[ip,1]
                q.qe[ip,2] = q.qn[ip,2]
                q.qe[ip,3] = q.qn[ip,3]
                q.qe[ip,4] = q.qn[ip,4]
                q.qe[ip,5] = q.qn[ip,5]
                q.qe[ip,6] = q.qn[ip,6]
            end
        end
    end
    
    if rank == 0
        @info " Initialize analytic neutral ABL ........................ DONE "
    end
    
    return q
end


function user_get_adapt_flags(inputs, old_ad_lvl, q, qe, connijk, nelem, ngl)
    adapt_flags = KernelAbstractions.zeros(CPU(), TInt, Int64(nelem))
    ips         = KernelAbstractions.zeros(CPU(), TInt, ngl * ngl)
    tol         = 1.0
    max_level   = inputs[:amr_max_level] 
    
    for iel = 1:nelem
        m = 1
        for i = 1:ngl
            for j = 1:ngl
                ips[m] = connijk[iel, i, j]
                m += 1
            end
        end
        # @info q[ips,4] - qe[ips,4]
        theta      = q[ips, 4] ./ q[ips, 1]
        theta_ref  = qe[ips, 4] ./ qe[ips, 1]
        dtheta     = theta - theta_ref
        # @info dtheta
        if any(dtheta .> tol) && (old_ad_lvl[iel] < max_level)
            adapt_flags[iel] = refine_flag
        end
        if all(dtheta .< tol)
            adapt_flags[iel] = coarsen_flag
        end
    end
    return adapt_flags
end
