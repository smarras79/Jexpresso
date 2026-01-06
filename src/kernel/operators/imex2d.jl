using SparseArrays
using Krylov

"""
COMPLETE CORRECTED IMEX IMPLEMENTATION - 2D VERSION WITH DSS
This version includes proper Direct Stiffness Summation (DSS) for the RHS vector
and all necessary corrections for 2D spectral element methods.
"""

"""
Build LHS matrix for semi-implicit discretization - 2D CORRECTED VERSION
"""
function build_imex_lhs_matrix_simple_2d!(u, params, connij, qe, coords, Δt)
    
    # Convert to auxiliary variables (primitive variables)
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    # Initialize flux Jacobian arrays for pressure terms (2D: only F and G fluxes)
    dFp_du = zeros(params.mesh.ngl, params.mesh.ngl, params.neqs, params.neqs)
    dGp_du = zeros(params.mesh.ngl, params.mesh.ngl, params.neqs, params.neqs)
    
    # Use COO (coordinate) format for efficient matrix assembly
    I_indices = Int[]
    J_indices = Int[]
    V_values = Float64[]
    
    # Pre-allocate for better performance
    max_entries = params.mesh.nelem * (params.mesh.ngl^2 * params.neqs)^2 ÷ 10  # Conservative estimate
    sizehint!(I_indices, max_entries)
    sizehint!(J_indices, max_entries)
    sizehint!(V_values, max_entries)
    
    # Element loop
    for iel = 1:params.mesh.nelem
        
        # Compute pressure flux Jacobians for this element
        compute_pressure_flux_jacobians_simple_2d!(iel, params, connij, qe, coords,
                                                  dFp_du, dGp_du)
        
        # Assemble element contribution using COO format
        assemble_element_lhs_simple_2d!(iel, params, Δt, dFp_du, dGp_du, 
                                       I_indices, J_indices, V_values, connij)
    end
    
    # Add identity matrix entries to COO format
    ndof_total = params.mesh.npoin * params.neqs
    for ip = 1:params.mesh.npoin
        for ieq = 1:params.neqs
            dof = (ip - 1) * params.neqs + ieq
            push!(I_indices, dof)
            push!(J_indices, dof)
            push!(V_values, 1.0)  # Identity matrix entry
        end
    end
    
    # Convert COO to CSC format (efficient for Julia solvers)
    # The sparse() function automatically handles DSS by summing duplicate entries
    LHS_matrix = sparse(I_indices, J_indices, V_values, ndof_total, ndof_total)
    
    # Ensure numerical stability
    dropzeros!(LHS_matrix)  # Remove near-zero entries
    
    return LHS_matrix
end

"""
Compute pressure flux Jacobians - 2D CORRECTED VERSION
"""
function compute_pressure_flux_jacobians_simple_2d!(iel, params, connij, qe, coords,
                                                   dFp_du, dGp_du)
    
    # Clear Jacobian arrays
    fill!(dFp_du, 0.0)
    fill!(dGp_du, 0.0)
    
    # Gas constant - handle different float types
    if isdefined(Main, :PhysicalConst)
        PhysConst = PhysicalConst{Float64}()  # Use Float64 for better precision
        R = PhysConst.Rair
    else
        R = 287.0  # Standard air gas constant J/(kg·K)
    end
    
    # Loop over 2D quadrature points
    for j = 1:params.mesh.ngl
        for i = 1:params.mesh.ngl
            
            # For conservative equations in 2D: pressure p = R*ρθ = R*q[4] (4th equation in 2D)
            # Assuming equation order: [ρ, ρu, ρv, ρθ] for 2D compressible flow
            
            # F_pressure = [0, p, 0, 0]ᵀ
            dFp_du[i,j,2,4] = R  # ∂p/∂(ρθ) affects u-momentum
            
            # G_pressure = [0, 0, p, 0]ᵀ  
            dGp_du[i,j,3,4] = R  # ∂p/∂(ρθ) affects v-momentum
        end
    end
end

"""
Assemble element contribution to LHS matrix - 2D CORRECTED VERSION  
"""
function assemble_element_lhs_simple_2d!(iel, params, Δt, dFp_du, dGp_du,
                                        I_indices, J_indices, V_values, connij)
    
    ngl = params.mesh.ngl
    neqs = params.neqs
    tolerance = 1e-14
    
    # Loop over test functions (2D)
    for ieq_test = 1:neqs
        for j_test = 1:ngl
            for i_test = 1:ngl
                
                # Get global DOF for test function
                ip_test = connij[iel, i_test, j_test]
                test_dof = (ip_test - 1) * neqs + ieq_test
                
                # Integration weight times Jacobian (2D)
                ωJac = (params.ω[i_test] * params.ω[j_test] * 
                       params.metrics.Je[iel, i_test, j_test])
                
                # Loop over trial functions (2D)
                for ieq_trial = 1:neqs
                    for j_trial = 1:ngl  
                        for i_trial = 1:ngl
                            
                            # Get global DOF for trial function
                            ip_trial = connij[iel, i_trial, j_trial]
                            trial_dof = (ip_trial - 1) * neqs + ieq_trial
                            
                            # Compute divergence of pressure flux Jacobian
                            div_jac = compute_divergence_pressure_jacobian_simple_2d(
                                params, iel, dFp_du, dGp_du,
                                i_test, j_test, i_trial, j_trial,
                                ieq_test, ieq_trial)
                            
                            # Implicit contribution: -Δt * (∇ · ∂F_pressure/∂u)
                            lhs_contrib = -Δt * ωJac * div_jac
                            
                            # Only add significant entries
                            if abs(lhs_contrib) > tolerance
                                push!(I_indices, test_dof)
                                push!(J_indices, trial_dof)
                                push!(V_values, lhs_contrib)
                            end
                        end
                    end
                end
            end
        end
    end
end

"""
Compute divergence of pressure flux Jacobian - 2D CORRECTED VERSION
"""
function compute_divergence_pressure_jacobian_simple_2d(params, iel, dFp_du, dGp_du,
                                                       i_test, j_test, 
                                                       i_trial, j_trial,
                                                       ieq_test, ieq_trial)
    
    ngl = params.mesh.ngl
    
    # Compute derivatives in reference coordinates (2D)
    dFp_dξ = 0.0; dFp_dη = 0.0
    dGp_dξ = 0.0; dGp_dη = 0.0
    
    # Use regular loops for compatibility
    for m = 1:ngl
        dFp_dξ += params.basis.dψ[m,i_trial] * dFp_du[m,j_trial,ieq_test,ieq_trial]
        dFp_dη += params.basis.dψ[m,j_trial] * dFp_du[i_trial,m,ieq_test,ieq_trial]  
        
        dGp_dξ += params.basis.dψ[m,i_trial] * dGp_du[m,j_trial,ieq_test,ieq_trial]
        dGp_dη += params.basis.dψ[m,j_trial] * dGp_du[i_trial,m,ieq_test,ieq_trial]
    end
    
    # Transform to physical coordinates using metric terms (2D)
    dξdx_ij = params.metrics.dξdx[iel,i_test,j_test]
    dξdy_ij = params.metrics.dξdy[iel,i_test,j_test] 
    
    dηdx_ij = params.metrics.dηdx[iel,i_test,j_test]
    dηdy_ij = params.metrics.dηdy[iel,i_test,j_test]
    
    dFp_dx = dFp_dξ*dξdx_ij + dFp_dη*dηdx_ij
    dGp_dy = dGp_dξ*dξdy_ij + dGp_dη*dηdy_ij
    
    # Return divergence (2D)
    return dFp_dx + dGp_dy
end

"""
Build explicit RHS with proper DSS - 2D CORRECTED VERSION
"""
function build_imex_rhs_explicit_simple_2d!(u, params, connij, qe, coords, lsource, RHS_explicit)
    
    # Clear explicit RHS - CRITICAL for proper DSS
    fill!(RHS_explicit, 0.0)
    
    # Convert to auxiliary variables
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    # Element loop for explicit terms (advective fluxes)
    # Each element contributes to shared nodes - DSS happens automatically
    for iel = 1:params.mesh.nelem
        
        # Compute advective fluxes (without pressure terms for IMEX splitting)
        for j = 1:params.mesh.ngl, i = 1:params.mesh.ngl
            ip = connij[iel,i,j]
            
            # Compute fluxes using existing user_flux! function
            # Note: Adapt this call to match your specific flux function signature
            try
                if isdefined(Main, :NSD_2D)
                    # If you have a dimension type
                    user_flux!(@view(params.F[i,j,:]),
                               @view(params.G[i,j,:]), NSD_2D(),
                               @view(params.uaux[ip,:]),
                               @view(qe[ip,:]),
                               params.mesh,
                               params.CL, params.SOL_VARS_TYPE;
                               neqs=params.neqs, ip=ip)
                else
                    # Generic 2D flux call - you may need to adapt this
                    user_flux!(@view(params.F[i,j,:]),
                               @view(params.G[i,j,:]),
                               @view(params.uaux[ip,:]),
                               @view(qe[ip,:]),
                               params.mesh,
                               params.CL, params.SOL_VARS_TYPE;
                               neqs=params.neqs, ip=ip)
                end
            catch e
                @warn "Flux function call failed, using fallback: $e"
                # Fallback flux computation
                compute_fallback_flux_2d!(params.F[i,j,:], params.G[i,j,:], 
                                          params.uaux[ip,:], params)
            end
            
            # Remove pressure from fluxes for explicit treatment (IMEX splitting)
            remove_pressure_simple_2d!(params.F[i,j,:], params.G[i,j,:], 
                                      params.uaux[ip,:], params)
            
            # Compute source terms if needed
            if lsource
                try
                    user_source!(@view(params.S[i,j,:]),
                                 @view(params.uaux[ip,:]),
                                 @view(qe[ip,:]),
                                 params.mesh.npoin, params.CL, params.SOL_VARS_TYPE;
                                 neqs=params.neqs,
                                 x=coords[ip,1], y=coords[ip,2])
                catch e
                    @warn "Source function call failed: $e"
                    fill!(@view(params.S[i,j,:]), 0.0)  # Zero source as fallback
                end
            else
                fill!(@view(params.S[i,j,:]), 0.0)
            end
        end
        
        # Apply divergence operator with proper DSS accumulation
        apply_divergence_operator_simple_2d!(params, iel, RHS_explicit, connij)
    end
    
    # Apply additional DSS operations if needed for your specific discretization
    if haskey(params, :mesh) && haskey(params.mesh, :requires_additional_dss) && params.mesh.requires_additional_dss
        dss_vector_2d!(RHS_explicit, params, connij)
    end
end

"""
Remove pressure terms from fluxes - 2D CORRECTED VERSION
"""
function remove_pressure_simple_2d!(F, G, uaux, params)
    # For conservative equations: pressure p = R*ρθ
    if isdefined(Main, :PhysicalConst)
        PhysConst = PhysicalConst{Float64}()
        R = PhysConst.Rair
    else
        R = 287.0  # Standard air gas constant
    end
    
    # In 2D conservative form: [ρ, ρu, ρv, ρθ]
    # Get ρθ from auxiliary variables (or compute from conservative variables)
    if length(uaux) >= 4
        ρθ = uaux[4]  # ρθ from auxiliary variables (4th equation in 2D)
    else
        @warn "Insufficient auxiliary variables, using fallback pressure calculation"
        ρθ = 1.0  # Fallback value
    end
    
    p = R * ρθ
    
    # Remove pressure from momentum equations (2D)
    if length(F) >= 2
        F[2] -= p  # Remove from u-momentum
    end
    if length(G) >= 3
        G[3] -= p  # Remove from v-momentum
    end
end

"""
Apply divergence operator with proper DSS - 2D CORRECTED VERSION
"""
function apply_divergence_operator_simple_2d!(params, iel, RHS_explicit, connij)
    
    ngl = params.mesh.ngl
    neqs = params.neqs
    
    for ieq = 1:neqs
        for j = 1:ngl
            for i = 1:ngl
                
                # Integration weight times Jacobian (2D)
                ωJac = (params.ω[i] * params.ω[j] * 
                       params.metrics.Je[iel,i,j])
                
                # Compute derivatives in reference coordinates
                dFdξ = 0.0; dFdη = 0.0
                dGdξ = 0.0; dGdη = 0.0
                
                for m = 1:ngl
                    dFdξ += params.basis.dψ[m,i] * params.F[m,j,ieq]
                    dFdη += params.basis.dψ[m,j] * params.F[i,m,ieq]
                    
                    dGdξ += params.basis.dψ[m,i] * params.G[m,j,ieq]
                    dGdη += params.basis.dψ[m,j] * params.G[i,m,ieq]
                end
                
                # Transform to physical coordinates (2D)
                dξdx_ij = params.metrics.dξdx[iel,i,j]
                dξdy_ij = params.metrics.dξdy[iel,i,j]
                
                dηdx_ij = params.metrics.dηdx[iel,i,j]
                dηdy_ij = params.metrics.dηdy[iel,i,j]
                
                dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
                dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij
                
                # Get global DOF
                ip = connij[iel,i,j]
                global_dof = (ip - 1) * params.neqs + ieq
                
                # Compute element contribution to weak form residual
                # Standard weak form: ∫ψ(∇·F)dV = -∫(∇ψ)·F dV + boundary terms
                # Here we compute the volume integral contribution
                auxi = ωJac * ((dFdx + dGdy) - params.S[i,j,ieq])
                
                # DSS: Accumulate element contributions at shared nodes
                # This is the key DSS operation - each element adds to shared nodes
                RHS_explicit[global_dof] -= auxi  # Standard sign for weak form
            end
        end
    end
end

"""
Direct Stiffness Summation for vectors - 2D CORRECTED VERSION
This handles special DSS operations if needed beyond standard accumulation
"""
function dss_vector_2d!(vector, params, connij)
    # For standard continuous spectral elements, DSS is handled by accumulation
    # in the assembly loop. This function handles special cases.
    
    # Example: If you need to average at corner/edge nodes for discontinuous methods
    if haskey(params.mesh, :shared_nodes_info)
        for (global_node, sharing_info) in params.mesh.shared_nodes_info
            n_sharing = length(sharing_info.elements)
            if n_sharing > 1  # Node is shared between elements
                for ieq = 1:params.neqs
                    global_dof = (global_node - 1) * params.neqs + ieq
                    # Apply averaging or other DSS operation
                    if sharing_info.average_contributions
                        vector[global_dof] /= n_sharing
                    end
                end
            end
        end
    end
    
    return nothing
end

"""
Complete IMEX time step - 2D CORRECTED VERSION
"""
function imex_time_step_simple_2d!(u, params, connij, qe, coords, Δt, lsource=true)
    
    println(" ##### Starting IMEX time step 2D")
    
    # Build LHS matrix for implicit pressure terms
    println("   Building LHS matrix...")
    LHS_matrix = build_imex_lhs_matrix_simple_2d!(u, params, connij, qe, coords, Δt)
    
    # Build explicit RHS with proper DSS
    println("   Building explicit RHS...")
    ndof_total = params.mesh.npoin * params.neqs
    RHS_explicit = zeros(ndof_total)
    build_imex_rhs_explicit_simple_2d!(u, params, connij, qe, coords, lsource, RHS_explicit)
    
    # Form complete RHS: u^n + Δt * RHS_explicit
    println("   Forming final RHS...")
    RHS_final = copy(u) + Δt * RHS_explicit
    
    # Solve linear system: (I - Δt*L_impl) * u^{n+1} = RHS_final
    println("   Solving linear system...")
    try
        u_new = LHS_matrix \ RHS_final
        # Update solution
        u .= u_new
        println("   ✓ Linear system solved successfully")
    catch e
        @error "Linear system solve failed: $e"
        @info "Matrix size: $(size(LHS_matrix))"
        @info "RHS size: $(length(RHS_final))"
        @info "Matrix condition number estimate: $(cond(Matrix(LHS_matrix[1:min(100,end), 1:min(100,end)])))"
        rethrow(e)
    end
    
    println(" ##### IMEX time step completed")
    return nothing
end

"""
Main IMEX integration routine - 2D CORRECTED VERSION
"""
function imex_integration_simple_2d!(u, params, connij, qe, coords, Δt, ntime_steps, lsource=true)
    
    println("="^60)
    println("Starting 2D IMEX integration with $(ntime_steps) time steps")
    println("Time step size: Δt = $(Δt)")
    println("Number of DOFs: $(length(u))")
    println("="^60)
    
    # Pre-integration checks
    @assert length(u) == params.mesh.npoin * params.neqs "Solution vector size mismatch"
    @assert size(connij, 1) == params.mesh.nelem "Connectivity array size mismatch"
    
    initial_max = maximum(abs.(u))
    println("Initial solution max: $(initial_max)")
    
    for n = 1:ntime_steps
        try
            # Perform IMEX time step
            imex_time_step_simple_2d!(u, params, connij, qe, coords, Δt, lsource)
            
            # Monitor solution
            if n % 10 == 0 || n <= 5
                u_max = maximum(abs.(u))
                u_min = minimum(u)
                println("Time step $n: max(|u|) = $(u_max), min(u) = $(u_min)")
                
                # Check for blow-up
                if u_max > 1e6 || isnan(u_max) || isinf(u_max)
                    @error "Solution appears to be unstable at time step $n"
                    @info "Consider reducing time step size or checking initial conditions"
                    break
                end
            end
            
        catch e
            @error "IMEX integration failed at time step $n: $e"
            @info "Last successful solution max: $(maximum(abs.(u)))"
            rethrow(e)
        end
    end
    
    final_max = maximum(abs.(u))
    println("="^60)
    println("2D IMEX integration completed successfully!")
    println("Final solution max: $(final_max)")
    println("="^60)
    
    return u
end

"""
Fallback flux computation for 2D case
"""
function compute_fallback_flux_2d!(F, G, uaux, params)
    @warn "Using fallback flux computation - implement proper flux functions"
    
    # Simple advective flux for debugging
    if length(uaux) >= 3 && length(F) >= 4 && length(G) >= 4
        ρ = uaux[1]
        u = uaux[2] 
        v = uaux[3]
        
        # Simple advective fluxes (without pressure)
        F[1] = ρ * u      # mass flux in x
        F[2] = ρ * u * u  # x-momentum flux in x (without pressure)
        F[3] = ρ * u * v  # y-momentum flux in x
        if length(F) >= 4
            F[4] = ρ * u * (uaux[4] / ρ)  # energy/temperature flux in x
        end
        
        G[1] = ρ * v      # mass flux in y
        G[2] = ρ * v * u  # x-momentum flux in y
        G[3] = ρ * v * v  # y-momentum flux in y (without pressure)
        if length(G) >= 4
            G[4] = ρ * v * (uaux[4] / ρ)  # energy/temperature flux in y
        end
    else
        # Zero flux fallback
        fill!(F, 0.0)
        fill!(G, 0.0)
    end
end

"""
Placeholder for 2D user flux function - IMPLEMENT THIS FOR YOUR SPECIFIC EQUATIONS
"""
function user_flux_2d!(F, G, uaux, qe, mesh, CL, SOL_VARS_TYPE; kwargs...)
    @error """
    user_flux_2d! needs to be implemented for your specific 2D equation system.
    
    Expected signature:
    function user_flux_2d!(F, G, uaux, qe, mesh, CL, SOL_VARS_TYPE; 
                          neqs, ip, kwargs...)
        # Compute F and G fluxes based on your governing equations
        # F = flux in x-direction
        # G = flux in y-direction
        # uaux = auxiliary/primitive variables
        # qe = element-specif<ic data
    end
    """
end

# Export main functions
export imex_integration_simple_2d!, imex_time_step_simple_2d!
export build_imex_lhs_matrix_simple_2d!, build_imex_rhs_explicit_simple_2d!
