function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
    Initialize 2D thermal rising bubble with background state from sounding
    """
    @info " Initialize 2D thermal bubble with sounding ........................ "
    @info "=== CHECKING SOL_VARS_TYPE ==="
    @info "SOL_VARS_TYPE = $(inputs[:SOL_VARS_TYPE])"
    @info "Is PERT? $(inputs[:SOL_VARS_TYPE] == PERT())"
    @info "Is TOTAL? $(inputs[:SOL_VARS_TYPE] == TOTAL())"
    @info "=============================="
    
    #---------------------------------------------------------------------------------
    # Solution variables
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρθ", "hl", "qtr"]
    qoutvars = ["ρ", "u", "v", "θ", "hl", "qtr"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------
    
    if (inputs[:backend] == CPU())    
        PhysConst = PhysicalConst{Float64}()
        
        if (inputs[:case] === "rtb")
            comm = MPI.COMM_WORLD
            rank = MPI.Comm_rank(comm)
            
            #==========================================================================
            STEP 1: Read and process sounding data
            ==========================================================================#
            if rank == 0
                @info "Reading sounding from: $(inputs[:sounding_file])"
            end
            
            data = read_sounding(inputs[:sounding_file])
            
            sounding_nlevels = size(data)[1]
            sounding_nvars   = size(data)[2]
            data_with_p      = zeros(Float64, sounding_nlevels, sounding_nvars+1)
            
            # Surface pressure (adjust based on your sounding)
            p_surface = 100000.0  # 1000 hPa in Pa
            
            # Copy original data (height, theta, qv, u, v)
            # NOTE: For dry case, qv column should be zeros or very small
            data_with_p[:,1:sounding_nvars] = copy(data[:,:])
            
            #==========================================================================
            STEP 2: Calculate pressure hydrostatically
            ==========================================================================#
            pressure = zeros(Float64, sounding_nlevels)
            
            # Initial pressure calculation
            for i = 1:sounding_nlevels
                z_current  = data[i,1]  # height (m)
                θ_current  = data[i,2]  # potential temperature (K)
                qv_current = data[i,3]  # water vapor mixing ratio (g/kg) - should be ~0 for dry
                
                # Convert mixing ratio from g/kg to kg/kg
                qv_kg_per_kg = qv_current / 1000.0
                
                if i == 1
                    # First level: use surface pressure
                    θ_v = θ_current * (1.0 + 0.61 * qv_kg_per_kg)
                    T_v = θ_v * (p_surface/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                    pressure[i] = p_surface * exp(-PhysConst.g * z_current / (PhysConst.Rair * T_v))
                else
                    # Subsequent levels: integrate hydrostatically
                    z_prev  = data[i-1,1]
                    θ_prev  = data[i-1,2]
                    qv_prev = data[i-1,3] / 1000.0
                    p_prev  = pressure[i-1]
                    
                    θ_avg  = 0.5 * (θ_current + θ_prev)
                    qv_avg = 0.5 * (qv_kg_per_kg + qv_prev)
                    dz     = z_current - z_prev
                    
                    θ_v_avg = θ_avg * (1.0 + 0.61 * qv_avg)
                    T_v_avg = θ_v_avg * (p_prev/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                    
                    pressure[i] = p_prev * exp(-PhysConst.g * dz / (PhysConst.Rair * T_v_avg))
                end
            end
            
            #==========================================================================
            STEP 3: Iterative hydrostatic balance correction
            ==========================================================================#
            max_iterations = 10
            tolerance      = 1.0  # Pa
            converged      = false
            
            if rank == 0
                @info "Starting iterative hydrostatic balance correction..."
            end
            
            for iter = 1:max_iterations
                pressure_old   = copy(pressure)
                max_correction = 0.0
                
                for i = 2:sounding_nlevels
                    z_current  = data[i,1]
                    z_prev     = data[i-1,1]
                    θ_current  = data[i,2]
                    θ_prev     = data[i-1,2]
                    qv_current = data[i,3]/1000.0
                    qv_prev    = data[i-1,3]/1000.0
                    
                    dz = z_current - z_prev
                    p_prev = pressure[i-1]
                    
                    θ_v_prev    = θ_prev * (1.0 + 0.61 * qv_prev)
                    θ_v_current = θ_current * (1.0 + 0.61 * qv_current)
                    
                    T_v_prev = θ_v_prev * (p_prev/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                    p_est    = p_prev * exp(-PhysConst.g * dz / (PhysConst.Rair * T_v_prev))
                    
                    T_v_current = θ_v_current * (p_est/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                    T_v_avg     = 0.5 * (T_v_prev + T_v_current)
                    
                    pressure[i] = p_prev * exp(-PhysConst.g * dz / (PhysConst.Rair * T_v_avg))
                    
                    correction     = abs(pressure[i] - pressure_old[i])
                    max_correction = max(max_correction, correction)
                end
                
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
                @warn "Hydrostatic balance did not converge after $max_iterations iterations"
            end
            
            # Store calculated pressures
            data_with_p[:, sounding_nvars+1] .= pressure
            
            if rank == 0
                @info "Pressure range: $(minimum(pressure)) to $(maximum(pressure)) Pa"
            end
            
            #==========================================================================
            STEP 4: Interpolate sounding to mesh points
            ==========================================================================#
            data_interpolate = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.y, data_with_p)
            
            if rank == 0
                @info "Interpolated sounding to $(mesh.npoin) grid points"
            end
            
            #==========================================================================
            STEP 5: Define thermal bubble parameters
            ==========================================================================#
            # Bubble center and dimensions
            xc = 0.0        # Center at x=0
            yc = 2000.0     # 2 km height
            rx = 10000.0    # Elliptical x-radius (10 km)
            rz = 1500.0     # Elliptical z-radius (1.5 km)
            
            # Bubble amplitude
            θc = 3.0        # K (potential temperature perturbation)
            
            # Tracer parameters (optional)
            xc_tracer = -2000.0
            yc_tracer = 4000.0
            r0_tracer = 1000.0
            qtrc = 1.0
            
            if rank == 0
                @info "Bubble parameters:"
                @info "  Center: (x,y) = ($xc, $yc) m"
                @info "  Radii: (rx,rz) = ($rx, $rz) m"
                @info "  θ perturbation: $θc K"
                @info "🔍 DOMAIN CHECK:"  # ← ADD THIS
                @info "  Domain X: [$(minimum(mesh.x)), $(maximum(mesh.x))] m"  # ← ADD THIS
                @info "  Domain Y: [$(minimum(mesh.y)), $(maximum(mesh.y))] m"  # ← ADD THIS
            end
            
            #==========================================================================
            STEP 6: Initialize each grid point
            ==========================================================================#
            for ip = 1:mesh.npoin
                x = mesh.x[ip]
                y = mesh.y[ip]
                
                # Get background state from interpolated sounding
                θ_bg   = data_interpolate[ip,1]  # Potential temperature
                qv_bg  = data_interpolate[ip,2] / 1000.0  # Convert g/kg to kg/kg (should be ~0)
                u_bg   = data_interpolate[ip,3]  # u velocity
                v_bg   = data_interpolate[ip,4]  # v velocity (this is actually w in 2D x-z)
                Press_bg = data_interpolate[ip,5]  # Pressure
                
                # Calculate distance from bubble center (elliptical)
                r = sqrt((x - xc)^2/rx^2 + (y - yc)^2/rz^2)
                
                # Add thermal perturbation
                Δθ = 0.0
                if r <= 1.0
                    Δθ = θc * cospi(r/2)^2
                end
                
                # Add tracer perturbation (optional)
                r_tracer = sqrt((x - xc_tracer)^2 + (y - yc_tracer)^2)
                Δqtr = 0.0
                if r_tracer <= r0_tracer
                    Δqtr = qtrc * cospi(r_tracer/(2*r0_tracer))^2
                end
                
                # Total state
                θ = θ_bg + Δθ
                u = u_bg
                v = v_bg
                qtr = Δqtr
                
                # Calculate density from equation of state
                # For background state
                ρ_bg = perfectGasLaw_θPtoρ(PhysConst; θ=θ_bg, Press=Press_bg)
                
                # For perturbed state - recalculate pressure if needed
                # For simplicity, we can use the same pressure (small perturbation assumption)
                # Or recalculate: this is more accurate
                Press = Press_bg  # Simple approach
                ρ = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=Press)
                
                # Calculate enthalpy
                T_bg = θ_bg * (Press_bg/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                T    = θ * (Press/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                hl_bg = PhysConst.cp*T_bg + PhysConst.g*y
                hl    = PhysConst.cp*T + PhysConst.g*y
                
                #======================================================================#
                # Store solution based on perturbation or total variable formulation
                #======================================================================#
                if inputs[:SOL_VARS_TYPE] == PERT()
                    # Store perturbations
                    q.qn[ip,1] = ρ - ρ_bg
                    q.qn[ip,2] = ρ*u - ρ_bg*u_bg
                    q.qn[ip,3] = ρ*v - ρ_bg*v_bg
                    q.qn[ip,4] = ρ*θ - ρ_bg*θ_bg
                    q.qn[ip,5] = ρ*hl - ρ_bg*hl_bg
                    q.qn[ip,6] = qtr
                    q.qn[ip,end] = Press_bg
                    
                    # Store background state
                    q.qe[ip,1] = ρ_bg
                    q.qe[ip,2] = u_bg
                    q.qe[ip,3] = v_bg
                    q.qe[ip,4] = ρ_bg*θ_bg
                    q.qe[ip,5] = ρ_bg*hl_bg
                    q.qe[ip,6] = 0.0
                    q.qe[ip,end] = Press_bg
                else
                    # Store total state
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*θ
                    q.qn[ip,5] = ρ*hl
                    q.qn[ip,6] = qtr
                    q.qn[ip,end] = Press
                    
                    # Store background state for reference
                    q.qe[ip,1] = ρ_bg
                    q.qe[ip,2] = u_bg
                    q.qe[ip,3] = v_bg
                    q.qe[ip,4] = ρ_bg*θ_bg
                    q.qe[ip,5] = ρ_bg*hl_bg
                    q.qe[ip,6] = 0.0
                    q.qe[ip,end] = Press_bg
                end
            end
            
            if rank == 0
                @info "Initialization complete!"
                @info "  ρ range: $(minimum(q.qn[:,1])) to $(maximum(q.qn[:,1]))"
                @info "  θ range: $(minimum(q.qn[:,4]./q.qn[:,1])) to $(maximum(q.qn[:,4]./q.qn[:,1]))"
            end
            
            #==========================================================================#
            # ADD BUBBLE DIAGNOSTICS HERE ← ← ← ← ← ← ← ← ← ← ← ← ← ← ← ← ← ← ← ←#
            #==========================================================================#
            if rank == 0
                bubble_points = 0
                max_theta_pert = 0.0
                
                for ip = 1:mesh.npoin
                    x = mesh.x[ip]
                    y = mesh.y[ip]
                    r = sqrt((x - xc)^2/rx^2 + (y - yc)^2/rz^2)
                    
                    if r <= 1.0
                        bubble_points += 1
                        
                        if inputs[:SOL_VARS_TYPE] == PERT()
                            # In PERT mode: q.qn has perturbations
                            ρ_pert = q.qn[ip,1]
                            ρθ_pert = q.qn[ip,4]
                            ρ_bg = q.qe[ip,1]
                            ρθ_bg = q.qe[ip,4]
                            
                            # Total values
                            ρ_total = ρ_pert + ρ_bg
                            ρθ_total = ρθ_pert + ρθ_bg
                            
                            θ_total = ρθ_total / ρ_total
                            θ_bg = ρθ_bg / ρ_bg
                            
                            θ_pert = θ_total - θ_bg
                            max_theta_pert = max(max_theta_pert, θ_pert)
                        else
                            # In TOTAL mode: q.qn has total values
                            ρ_total = q.qn[ip,1]
                            ρθ_total = q.qn[ip,4]
                            ρ_bg = q.qe[ip,1]
                            ρθ_bg = q.qe[ip,4]
                            
                            θ_total = ρθ_total / ρ_total
                            θ_bg = ρθ_bg / ρ_bg
                            
                            θ_pert = θ_total - θ_bg
                            max_theta_pert = max(max_theta_pert, θ_pert)
                        end
                    end
                end
                
                @info "=========================================="
                @info "BUBBLE DIAGNOSTICS:"
                @info "  Domain x range: [$(minimum(mesh.x)), $(maximum(mesh.x))] m"
                @info "  Domain y range: [$(minimum(mesh.y)), $(maximum(mesh.y))] m"
                @info "  Bubble center: (x,y) = ($xc, $yc) m"
                @info "  Bubble radii: (rx,rz) = ($rx, $rz) m"
                @info "  Points inside bubble (r≤1): $bubble_points"
                @info "  Max θ perturbation: $max_theta_pert K (expected: ~$θc K)"
                @info "=========================================="
            end
            
        else
            error("ERROR: CompEuler: initialize.jl:\n assign value to inputs[:case]")
        end
        
        #==========================================================================#
        # Handle non-conservative form if needed
        #==========================================================================# 
        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[:,2] .= q.qn[:,2]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,3] .= q.qn[:,3]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,4] .= q.qn[:,4]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,5] .= q.qn[:,5]./(q.qn[:,1] + q.qe[:,1])
                
                q.qe[:,4] .= q.qe[:,4]./q.qe[:,1]
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]
                q.qn[:,5] .= q.qn[:,5]./q.qn[:,1]
                
                q.qe[:,4] .= q.qe[:,4]./q.qe[:,1]
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            end
        end
        
    else
        # GPU implementation would go here
        error("GPU backend not yet implemented for sounding-based initialization")
    end
    
    @info " Initialize 2D thermal bubble with sounding ........................ DONE "
    
    return q
end