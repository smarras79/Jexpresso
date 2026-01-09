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
    # Solution variables - CHANGED: No enthalpy, use direct θ equation
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρθ", "qtr", "qtr2"]
    qoutvars = ["ρ", "u", "v", "θ", "qtr", "qtr2"]
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
            
            # Surface pressure
            p_surface = 100000.0  # 1000 hPa in Pa
            
            # Copy original data (height, theta, qv, u, v)
            data_with_p[:,1:sounding_nvars] = copy(data[:,:])
            
            #==========================================================================
            STEP 2: Calculate pressure hydrostatically
            ==========================================================================#
            pressure = zeros(Float64, sounding_nlevels)
            
            # Initial pressure calculation
            for i = 1:sounding_nlevels
                z_current  = data[i,1]
                θ_current  = data[i,2]
                qv_current = data[i,3]
                
                qv_kg_per_kg = qv_current / 1000.0
                
                if i == 1
                    θ_v = θ_current * (1.0 + 0.61 * qv_kg_per_kg)
                    T_v = θ_v * (p_surface/PhysConst.pref)^(PhysConst.Rair/PhysConst.cp)
                    pressure[i] = p_surface * exp(-PhysConst.g * z_current / (PhysConst.Rair * T_v))
                else
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
            tolerance      = 1.0
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
            xc = 0.0
            yc = 2500.0
            r0 = 2000.0
            rx = 10000.0
            rz = 1500.0
            θc = 3.0
            
            # Tracer parameters
            xc_tracer = -2000.0
            yc_tracer = 4000.0
            r0_tracer = 1000.0
            qtrc = 1.0
            
            if rank == 0
                @info "Bubble parameters:"
                @info "  Center: (x,y) = ($xc, $yc) m"
                @info "  Radius: $r0 m"
                @info "  θ perturbation: $θc K"
            end
    
            #==========================================================================
            STEP 6: Initialize each grid point
            ==========================================================================#
            for ip = 1:mesh.npoin
                x = mesh.x[ip]
                y = mesh.y[ip]
                
                # Get background state from interpolated sounding
                θ_bg     = data_interpolate[ip,1]
                qv_bg    = data_interpolate[ip,2] / 1000.0
                u_bg     = data_interpolate[ip,3]
                v_bg     = data_interpolate[ip,4]
                Press_bg = data_interpolate[ip,5]
                
                # Calculate distance from bubble center
               # r = sqrt((x - xc)^2 + (y - yc)^2)
               r = sqrt((x - xc)^2/(rx^2) + (y - yc)^2/(rz^2))
                
                # Add thermal perturbation
                Δθ = 0.0
                if r <= 1.0
                   # Δθ = θc * (1.0 - r/r0)
                    Δθ = θc * cospi(r/2)^2
                end
                
                
                # Add tracer perturbation
                r_tracer = sqrt((x - xc_tracer)^2 + (y - yc_tracer)^2)
                Δqtr = 0.0
                if r_tracer < r0_tracer
                    Δqtr = qtrc * (1.0 - r_tracer/r0_tracer)
                end
                
                # Total state
                θ = θ_bg + Δθ
                u = u_bg
                v = v_bg
                qtr = Δqtr
                qtr2 = 0.0  # Second tracer (unused)

                
                
                # Calculate density from equation of state
                ρ_bg = perfectGasLaw_θPtoρ(PhysConst; θ=θ_bg, Press=Press_bg)
                Press = Press_bg  # Small perturbation assumption
                ρ = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=Press)
                
                #======================================================================#
                # Store solution - NO ENTHALPY, just ρθ directly
                #======================================================================#
                if inputs[:SOL_VARS_TYPE] == PERT()
                    # Store perturbations
                    q.qn[ip,1] = ρ - ρ_bg
                    q.qn[ip,2] = ρ*u - ρ_bg*u_bg
                    q.qn[ip,3] = ρ*v - ρ_bg*v_bg
                    q.qn[ip,4] = ρ*θ - ρ_bg*θ_bg
                    q.qn[ip,5] = qtr
                    q.qn[ip,6] = qtr2
                    q.qn[ip,end] = Press_bg
                    
                    # Store background state
                    q.qe[ip,1] = ρ_bg
                    q.qe[ip,2] = u_bg
                    q.qe[ip,3] = v_bg
                    q.qe[ip,4] = ρ_bg*θ_bg
                    q.qe[ip,5] = 0.0
                    q.qe[ip,6] = 0.0
                    q.qe[ip,end] = Press_bg
                else
                    # Store total state
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*θ
                    q.qn[ip,5] = qtr
                    q.qn[ip,6] = qtr2
                    q.qn[ip,end] = Press
                    
                    # Store background state for reference
                    q.qe[ip,1] = ρ_bg
                    q.qe[ip,2] = u_bg
                    q.qe[ip,3] = v_bg
                    q.qe[ip,4] = ρ_bg*θ_bg
                    q.qe[ip,5] = 0.0
                    q.qe[ip,6] = 0.0
                    q.qe[ip,end] = Press_bg
                end
            end
            
            if rank == 0
                @info "Initialization complete!"
                @info "  ρ range: $(minimum(q.qn[:,1])) to $(maximum(q.qn[:,1]))"
                @info "  θ range: $(minimum(q.qn[:,4]./q.qn[:,1])) to $(maximum(q.qn[:,4]./q.qn[:,1]))"
            end

            #==========================================================================
            BUBBLE DIAGNOSTIC - moved outside previous if block
            ==========================================================================#
            # Get global domain bounds
            global_xmin = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)
            global_xmax = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
            global_ymin = MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm)
            global_ymax = MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm)
            
            if rank == 0
                @info "=== BUBBLE DIAGNOSTIC ==="
                @info "GLOBAL Domain X range: $global_xmin to $global_xmax"
                @info "GLOBAL Domain Y range: $global_ymin to $global_ymax"
                @info "Bubble center: x=0, y=2500"
            end
            
            # Check on ALL ranks if bubble exists
            has_bubble_local = 0
            for ip = 1:mesh.npoin
                x = mesh.x[ip]
                y = mesh.y[ip]
                if abs(x) < 500.0 && abs(y - 2500.0) < 500.0
                    has_bubble_local = 1
                    break
                end
            end
            
            has_bubble_global = MPI.Allreduce(has_bubble_local, MPI.MAX, comm)
            
            if rank == 0
                if has_bubble_global > 0
                    @info "✅ Bubble region EXISTS in the mesh"
                else
                    @error "❌ NO MESH POINTS NEAR BUBBLE CENTER!"
                end
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
                
                q.qe[:,4] .= q.qe[:,4]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]

                q.qe[:,4] .= q.qe[:,4]./q.qe[:,1]
            end
        end
        
    else
        error("GPU backend not yet implemented for sounding-based initialization")
    end
    
    @info " Initialize 2D thermal bubble with sounding ........................ DONE "
    
    return q
end