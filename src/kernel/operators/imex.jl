using SparseArrays
using Krylov

"""
SIMPLIFIED IMEX IMPLEMENTATION
This version focuses on the core matrix-based approach without problematic matrix-free functions
"""

"""
Build LHS matrix for semi-implicit discretization - SIMPLIFIED VERSION
"""
function build_imex_lhs_matrix_simple!(u, params, connijk, qe, coords, Δt)
    
    # Convert to auxiliary variables (primitive variables)
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    # Initialize flux Jacobian arrays for pressure terms
    dFp_du = zeros(params.mesh.ngl, params.mesh.ngl, params.mesh.ngl, params.neqs, params.neqs)
    dGp_du = zeros(params.mesh.ngl, params.mesh.ngl, params.mesh.ngl, params.neqs, params.neqs)
    dHp_du = zeros(params.mesh.ngl, params.mesh.ngl, params.mesh.ngl, params.neqs, params.neqs)
    
    # Use COO (coordinate) format for efficient matrix assembly
    I_indices = Int[]
    J_indices = Int[]
    V_values = Float64[]
    
    # Pre-allocate for better performance
    max_entries = params.mesh.nelem * (params.mesh.ngl^3 * params.neqs)^2 ÷ 10  # Conservative estimate
    sizehint!(I_indices, max_entries)
    sizehint!(J_indices, max_entries)
    sizehint!(V_values, max_entries)
    
    # Element loop
    for iel = 1:params.mesh.nelem
        
        # Compute pressure flux Jacobians for this element
        compute_pressure_flux_jacobians_simple!(iel, params, connijk, qe, coords,
                                               dFp_du, dGp_du, dHp_du)
        
        # Assemble element contribution using COO format
        assemble_element_lhs_simple!(iel, params, Δt, dFp_du, dGp_du, dHp_du, 
                                    I_indices, J_indices, V_values, connijk)
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
    LHS_matrix = sparse(I_indices, J_indices, V_values, ndof_total, ndof_total)
    
    # Ensure numerical stability
    dropzeros!(LHS_matrix)  # Remove near-zero entries
    
    return LHS_matrix
end

"""
Compute pressure flux Jacobians - SIMPLIFIED VERSION
"""
function compute_pressure_flux_jacobians_simple!(iel, params, connijk, qe, coords,
                                                dFp_du, dGp_du, dHp_du)
    
    # Clear Jacobian arrays
    fill!(dFp_du, 0.0)
    fill!(dGp_du, 0.0) 
    fill!(dHp_du, 0.0)
    
    # Gas constant
    PhysConst = PhysicalConst{Float32}()
    R = PhysConst.Rair
    
    # Loop over quadrature points
    for k = 1:params.mesh.ngl
        for j = 1:params.mesh.ngl
            for i = 1:params.mesh.ngl
                
                # For conservative equations: pressure p = R*ρθ = R*q[5]
                # F_pressure = [0, p, 0, 0, 0]ᵀ
                dFp_du[i,j,k,2,5] = R  # ∂p/∂(ρθ) affects u-momentum
                
                # G_pressure = [0, 0, p, 0, 0]ᵀ  
                dGp_du[i,j,k,3,5] = R  # ∂p/∂(ρθ) affects v-momentum
                
                # H_pressure = [0, 0, 0, p, 0]ᵀ
                dHp_du[i,j,k,4,5] = R  # ∂p/∂(ρθ) affects w-momentum
            end
        end
    end
end

"""
Assemble element contribution - SIMPLIFIED VERSION  
"""
function assemble_element_lhs_simple!(iel, params, Δt, dFp_du, dGp_du, dHp_du,
                                     I_indices, J_indices, V_values, connijk)
    
    ngl = params.mesh.ngl
    neqs = params.neqs
    tolerance = 1e-14
    
    # Loop over test functions
    for ieq_test = 1:neqs
        for k_test = 1:ngl
            for j_test = 1:ngl
                for i_test = 1:ngl
                    
                    # Get global DOF for test function
                    ip_test = connijk[iel, i_test, j_test, k_test]
                    test_dof = (ip_test - 1) * neqs + ieq_test
                    
                    # Integration weight times Jacobian
                    ωJac = (params.ω[i_test] * params.ω[j_test] * params.ω[k_test] * 
                           params.metrics.Je[iel, i_test, j_test, k_test])
                    
                    # Loop over trial functions  
                    for ieq_trial = 1:neqs
                        for k_trial = 1:ngl
                            for j_trial = 1:ngl  
                                for i_trial = 1:ngl
                                    
                                    # Get global DOF for trial function
                                    ip_trial = connijk[iel, i_trial, j_trial, k_trial]
                                    trial_dof = (ip_trial - 1) * neqs + ieq_trial
                                    
                                    # Compute divergence of pressure flux Jacobian
                                    div_jac = compute_divergence_pressure_jacobian_simple(
                                        params, iel, dFp_du, dGp_du, dHp_du,
                                        i_test, j_test, k_test, i_trial, j_trial, k_trial,
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
    end
end

"""
Compute divergence of pressure flux Jacobian - SIMPLIFIED VERSION
"""
function compute_divergence_pressure_jacobian_simple(params, iel, dFp_du, dGp_du, dHp_du,
                                                    i_test, j_test, k_test, 
                                                    i_trial, j_trial, k_trial,
                                                    ieq_test, ieq_trial)
    
    ngl = params.mesh.ngl
    
    # Compute derivatives in reference coordinates
    dFp_dξ = 0.0; dFp_dη = 0.0; dFp_dζ = 0.0
    dGp_dξ = 0.0; dGp_dη = 0.0; dGp_dζ = 0.0  
    dHp_dξ = 0.0; dHp_dη = 0.0; dHp_dζ = 0.0
    
    # Use regular loops (no @turbo to avoid issues)
    for m = 1:ngl
        dFp_dξ += params.basis.dψ[m,i_trial] * dFp_du[m,j_trial,k_trial,ieq_test,ieq_trial]
        dFp_dη += params.basis.dψ[m,j_trial] * dFp_du[i_trial,m,k_trial,ieq_test,ieq_trial]  
        dFp_dζ += params.basis.dψ[m,k_trial] * dFp_du[i_trial,j_trial,m,ieq_test,ieq_trial]
        
        dGp_dξ += params.basis.dψ[m,i_trial] * dGp_du[m,j_trial,k_trial,ieq_test,ieq_trial]
        dGp_dη += params.basis.dψ[m,j_trial] * dGp_du[i_trial,m,k_trial,ieq_test,ieq_trial]
        dGp_dζ += params.basis.dψ[m,k_trial] * dGp_du[i_trial,j_trial,m,ieq_test,ieq_trial]
        
        dHp_dξ += params.basis.dψ[m,i_trial] * dHp_du[m,j_trial,k_trial,ieq_test,ieq_trial] 
        dHp_dη += params.basis.dψ[m,j_trial] * dHp_du[i_trial,m,k_trial,ieq_test,ieq_trial]
        dHp_dζ += params.basis.dψ[m,k_trial] * dHp_du[i_trial,j_trial,m,ieq_test,ieq_trial]
    end
    
    # Transform to physical coordinates using metric terms
    dξdx_ij = params.metrics.dξdx[iel,i_test,j_test,k_test]
    dξdy_ij = params.metrics.dξdy[iel,i_test,j_test,k_test] 
    dξdz_ij = params.metrics.dξdz[iel,i_test,j_test,k_test]
    
    dηdx_ij = params.metrics.dηdx[iel,i_test,j_test,k_test]
    dηdy_ij = params.metrics.dηdy[iel,i_test,j_test,k_test]
    dηdz_ij = params.metrics.dηdz[iel,i_test,j_test,k_test]
    
    dζdx_ij = params.metrics.dζdx[iel,i_test,j_test,k_test]
    dζdy_ij = params.metrics.dζdy[iel,i_test,j_test,k_test] 
    dζdz_ij = params.metrics.dζdz[iel,i_test,j_test,k_test]
    
    dFp_dx = dFp_dξ*dξdx_ij + dFp_dη*dηdx_ij + dFp_dζ*dζdx_ij
    dGp_dy = dGp_dξ*dξdy_ij + dGp_dη*dηdy_ij + dGp_dζ*dζdy_ij  
    dHp_dz = dHp_dξ*dξdz_ij + dHp_dη*dηdz_ij + dHp_dζ*dζdz_ij
    
    # Return divergence
    return dFp_dx + dGp_dy + dHp_dz
end

"""
Build explicit RHS - SIMPLIFIED VERSION
"""
function build_imex_rhs_explicit_simple!(u, params, connijk, qe, coords, lsource, RHS_explicit)
    
    # Clear explicit RHS
    fill!(RHS_explicit, 0.0)
    
    # Convert to auxiliary variables
    u2uaux!(@view(params.uaux[:,:]), u, params.neqs, params.mesh.npoin)
    
    # Element loop for explicit terms (advective fluxes)
    for iel = 1:params.mesh.nelem
        
        # Compute advective fluxes (same as original RHS but without pressure terms)
        for k = 1:params.mesh.ngl, j = 1:params.mesh.ngl, i = 1:params.mesh.ngl
            ip = connijk[iel,i,j,k]
            
            # Use your existing flux computation but remove pressure
            if !(params.inputs[:lsaturation])
                user_flux!(@view(params.F[i,j,k,:]),
                          @view(params.G[i,j,k,:]),
                          @view(params.H[i,j,k,:]),
                          @view(params.uaux[ip,:]),
                          @view(qe[ip,:]),
                          params.mesh,
                          params.CL, params.SOL_VARS_TYPE;
                          neqs=params.neqs, ip=ip)
            else
                user_flux!(@view(params.F[i,j,k,:]),
                          @view(params.G[i,j,k,:]),
                          @view(params.H[i,j,k,:]),
                          @view(params.uaux[ip,:]),
                          @view(qe[ip,:]),
                          params.mesh, params.thermo_params,
                          params.CL, params.SOL_VARS_TYPE;
                          neqs=params.neqs, ip=ip,
                          x=coords[ip,1], y=coords[ip,2], z=coords[ip,3])
            end
            
            # Remove pressure from fluxes for explicit treatment
            remove_pressure_simple!(params.F[i,j,k,:], params.G[i,j,k,:], params.H[i,j,k,:], 
                                   params.uaux[ip,:], params)
            
            if lsource
                # Add source terms
                user_source!(@view(params.S[i,j,k,:]),
                            @view(params.uaux[ip,:]),
                            @view(qe[ip,:]),
                            params.mesh.npoin,
                            params.CL, params.SOL_VARS_TYPE;
                            neqs=params.neqs,
                            x=coords[ip,1], y=coords[ip,2], z=coords[ip,3],
                            xmax=params.xmax, xmin=params.xmin, zmax=params.zmax)
            end
        end
        
        # Apply divergence operator
        apply_divergence_operator_simple!(params, iel, RHS_explicit, connijk)
    end
end

"""
Remove pressure terms from fluxes - SIMPLIFIED VERSION
"""
function remove_pressure_simple!(F, G, H, uaux, params)
    # For conservative equations: pressure p = R*ρθ
    PhysConst = PhysicalConst{TFloat}()
    R = PhysConst.Rair
    ρθ = uaux[5]  # ρθ from auxiliary variables
    p = R * ρθ
    
    # Remove pressure from momentum equations
    F[2] -= p  # Remove from u-momentum
    G[3] -= p  # Remove from v-momentum
    H[4] -= p  # Remove from w-momentum
end

"""
Apply divergence operator - SIMPLIFIED VERSION
"""
function apply_divergence_operator_simple!(params, iel, RHS_explicit, connijk)
    
    println(" ##### apply_divergence_operator_simple!!  ")
    
    ngl = params.mesh.ngl
    neqs = params.neqs
    
    for ieq = 1:neqs
        for k = 1:ngl
            for j = 1:ngl
                for i = 1:ngl
                    
                    ωJac = (params.ω[i] * params.ω[j] * params.ω[k] * 
                           params.metrics.Je[iel,i,j,k])
                    
                    # Compute derivatives in reference coordinates (no @turbo)
                    dFdξ = 0.0; dFdη = 0.0; dFdζ = 0.0
                    dGdξ = 0.0; dGdη = 0.0; dGdζ = 0.0
                    dHdξ = 0.0; dHdη = 0.0; dHdζ = 0.0
                    
                    for m = 1:ngl
                        dFdξ += params.basis.dψ[m,i] * params.F[m,j,k,ieq]
                        dFdη += params.basis.dψ[m,j] * params.F[i,m,k,ieq]
                        dFdζ += params.basis.dψ[m,k] * params.F[i,j,m,ieq]
                        
                        dGdξ += params.basis.dψ[m,i] * params.G[m,j,k,ieq]
                        dGdη += params.basis.dψ[m,j] * params.G[i,m,k,ieq]
                        dGdζ += params.basis.dψ[m,k] * params.G[i,j,m,ieq]
                        
                        dHdξ += params.basis.dψ[m,i] * params.H[m,j,k,ieq]
                        dHdη += params.basis.dψ[m,j] * params.H[i,m,k,ieq]
                        dHdζ += params.basis.dψ[m,k] * params.H[i,j,m,ieq]
                    end
                    
                    # Transform to physical coordinates
                    dξdx_ij = params.metrics.dξdx[iel,i,j,k]
                    dξdy_ij = params.metrics.dξdy[iel,i,j,k]
                    dξdz_ij = params.metrics.dξdz[iel,i,j,k]
                    
                    dηdx_ij = params.metrics.dηdx[iel,i,j,k]
                    dηdy_ij = params.metrics.dηdy[iel,i,j,k]
                    dηdz_ij = params.metrics.dηdz[iel,i,j,k]
                    
                    dζdx_ij = params.metrics.dζdx[iel,i,j,k]
                    dζdy_ij = params.metrics.dζdy[iel,i,j,k]
                    dζdz_ij = params.metrics.dζdz[iel,i,j,k]
                    
                    dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij + dFdζ*dζdx_ij
                    dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij + dGdζ*dζdy_ij
                    dHdz = dHdξ*dξdz_ij + dHdη*dηdz_ij + dHdζ*dζdz_ij
                    
                    # Get global DOF
                    ip = connijk[iel,i,j,k]
                    global_dof = (ip - 1) * params.neqs + ieq
                    
                    # Add to explicit RHS
                    auxi = ωJac * ((dFdx + dGdy + dHdz) - params.S[i,j,k,ieq])
                    RHS_explicit[global_dof] -= auxi
                end
            end
        end
    end

    println(" ##### apply_divergence_operator_simple!!  END ")
end

"""
Complete IMEX time step - SIMPLIFIED VERSION
"""
function imex_time_step_simple!(u, params, connijk, qe, coords, Δt, lsource=true)

    println(" ##### imex_time_step_simple!  ")
     println(" a")
    # Build LHS matrix for implicit pressure terms
    LHS_matrix = build_imex_lhs_matrix_simple!(u, params, connijk, qe, coords, Δt)
         println(" b")
    # Build explicit RHS
    ndof_total = params.mesh.npoin * params.neqs
    RHS_explicit = zeros(ndof_total)
    build_imex_rhs_explicit_simple!(u, params, connijk, qe, coords, lsource, RHS_explicit)
         println(" c")
    # Form complete RHS: u^n + Δt * RHS_explicit
    RHS_final = copy(u) + Δt * RHS_explicit
         println(" d")
    # Solve linear system: (I - Δt*L_impl) * u^{n+1} = RHS_final
    u_new = LHS_matrix \ RHS_final
         println(" e")
    # Update solution
    u .= u_new
         println(" f")
    
    println(" ##### imex_time_step_simple!  end ")
    return nothing
end

"""
Main IMEX integration routine - SIMPLIFIED VERSION
"""
function imex_integration_simple!(u, params, connijk, qe, coords, Δt, ntime_steps, lsource=true)
    
    println("Starting IMEX integration with $(ntime_steps) time steps...")
    
    for n = 1:ntime_steps
        # Perform IMEX time step
        imex_time_step_simple!(u, params, connijk, qe, coords, Δt, lsource)
        
        if n % 100 == 0
            u_max = maximum(abs.(u))
            println("Time step $n: max(|u|) = $(u_max)")
        end
    end
    
    println("IMEX integration completed.")
    return u
end
