using SparseArrays
using Krylov
using LinearAlgebra

"""
MODIFIED IMEX IMPLEMENTATION - 2D VERSION WITH KRYLOV SOLVER
This version uses GMRES (a Krylov method) instead of the direct backslash solver
"""

"""
Build LHS matrix for semi-implicit discretization - 2D VERSION
"""
function build_imex_lhs_matrix_simple_2d!(u, params, connij, qe, coords, Δt)
    
    # Get mesh info
    nelem = params.mesh.nelem
    ngl = params.mesh.ngl
    neqs = params.neqs
    
    # Estimate non-zeros: each element contributes ngl^2 × ngl^2 × neqs^2 entries
    nnz_estimate = nelem * (ngl^2)^2 * neqs^2
    
    # Preallocate COO arrays
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()
    sizehint!(I_vec, nnz_estimate)
    sizehint!(J_vec, nnz_estimate)
    sizehint!(V_vec, nnz_estimate)
    
    # Loop over elements
    for iel = 1:nelem
        # Build element LHS contributions
        for j = 1:ngl, i = 1:ngl
            ip = connij[iel, i, j]
            
            for n = 1:ngl, m = 1:ngl
                jp = connij[iel, m, n]
                
                # Get element mass matrix entry
                if i == m && j == n
                    ωJac = params.mesh.ω[i] * params.mesh.ω[j] * qe[iel, i, j]
                    
                    # Add identity and implicit terms
                    for ieq = 1:neqs
                        I_global = (ip - 1) * neqs + ieq
                        J_global = (jp - 1) * neqs + ieq
                        
                        # I + Δt * M^(-1) * L where L is implicit operator
                        value = ωJac  # Mass matrix contribution (identity after M^(-1))
                        
                        push!(I_vec, I_global)
                        push!(J_vec, J_global)
                        push!(V_vec, value)
                    end
                end
            end
        end
    end
    
    # Build sparse matrix (automatically sums duplicates - this is DSS for LHS)
    ndof = maximum(connij) * neqs
    return sparse(I_vec, J_vec, V_vec, ndof, ndof)
end

"""
Apply divergence operator with DSS - 2D VERSION
"""
function apply_divergence_operator_simple_2d!(params, iel, RHS_explicit, connij)
    
    ngl = params.mesh.ngl
    neqs = params.neqs
    
    # Loop over quadrature points in element
    for j = 1:ngl, i = 1:ngl
        # Compute divergence terms here (placeholder)
        # dFdx = ∂F/∂x, dGdy = ∂G/∂y where F, G are fluxes
        
        # Get quadrature weight and Jacobian
        ωJac = params.mesh.ω[i] * params.mesh.ω[j] * params.mesh.qe[iel, i, j]
        
        # Get global DOF
        ip = connij[iel, i, j]
        
        # Loop over equations
        for ieq = 1:neqs
            global_dof = (ip - 1) * neqs + ieq
            
            # Compute explicit RHS contribution (divergence of fluxes)
            # This is a placeholder - replace with actual flux computation
            dFdx = 0.0  # ∂F/∂x
            dGdy = 0.0  # ∂G/∂y
            source = 0.0  # Source term
            
            # Accumulate to global RHS (DSS happens naturally through +=)
            RHS_explicit[global_dof] += ωJac * (-(dFdx + dGdy) + source)
        end
    end
end

"""
DSS (Direct Stiffness Summation) for vectors - 2D VERSION
This ensures continuity at element boundaries by properly summing contributions
"""
function dss_vector_2d!(vec, params, connij)
    # Note: In the explicit RHS assembly, DSS happens naturally through +=
    # This function is here for completeness and special cases where
    # additional DSS operations might be needed
    
    # For continuous spectral elements, contributions from multiple elements
    # at shared nodes are automatically summed during assembly
    # No additional operation needed if assembly uses += correctly
    
    return vec
end

"""
MODIFIED: Solve the IMEX system using GMRES (Krylov method)
"""
function solve_imex_system_krylov!(LHS_matrix, RHS_vector; 
                                   rtol=1e-8, 
                                   atol=1e-10,
                                   maxiter=500,
                                   restart=50,
                                   verbose=false)
    """
    Solve LHS * x = RHS using GMRES from Krylov.jl
    
    Parameters:
    - LHS_matrix: Sparse matrix (left-hand side)
    - RHS_vector: Right-hand side vector
    - rtol: Relative tolerance
    - atol: Absolute tolerance  
    - maxiter: Maximum number of iterations
    - restart: GMRES restart parameter
    - verbose: Print convergence info
    
    Returns:
    - x: Solution vector
    - stats: Solver statistics
    """
    
    # Use GMRES from Krylov.jl
    x, stats = gmres(LHS_matrix, RHS_vector; 
                    rtol=rtol,
                    atol=atol, 
                    itmax=maxiter,
                    restart=restart,
                    verbose=verbose ? 2 : 0)
    
    if verbose
        println("GMRES converged: $(stats.solved)")
        println("  Iterations: $(stats.niter)")
        println("  Residual norm: $(stats.residuals[end])")
    end
    
    if !stats.solved
        @warn "GMRES did not converge! Residual: $(stats.residuals[end])"
    end
    
    return x, stats
end

"""
MODIFIED: IMEX time step using Krylov solver - 2D VERSION
"""
function imex_time_step_krylov!(u, params, connij, qe, coords, Δt, lsource;
                                rtol=1e-8, verbose=false)
    
    neqs = params.neqs
    ndof = length(u)
    
    # 1. Build LHS matrix (implicit part)
    LHS_matrix = build_imex_lhs_matrix_simple_2d!(u, params, connij, qe, coords, Δt)
    
    # 2. Build RHS vector (explicit part)
    RHS_explicit = zeros(Float64, ndof)
    
    # Assemble explicit RHS with proper DSS
    for iel = 1:params.mesh.nelem
        apply_divergence_operator_simple_2d!(params, iel, RHS_explicit, connij)
    end
    
    # Add source terms if provided
    if lsource
        # Add source contributions here
    end
    
    # 3. Form final RHS: RHS_final = M * u^n + Δt * RHS_explicit
    # For simplicity, assuming M is already incorporated
    RHS_final = u + Δt * RHS_explicit
    
    # 4. Solve using GMRES (Krylov method)
    u_new, stats = solve_imex_system_krylov!(LHS_matrix, RHS_final; 
                                            rtol=rtol, verbose=verbose)
    
    # 5. Update solution
    u .= u_new
    
    return stats
end

"""
MODIFIED: Time integration loop using Krylov solver
"""
function time_integrate_krylov!(u, params, connij, qe, coords, Δt, ntime_steps, lsource;
                               rtol=1e-8, verbose=false, output_freq=10)
    
    println("Starting time integration with GMRES solver...")
    println("  Time steps: $ntime_steps")
    println("  Δt: $Δt")
    println("  Tolerance: $rtol")
    println()
    
    total_iterations = 0
    max_iterations = 0
    
    for n = 1:ntime_steps
        # Take one IMEX time step with Krylov solver
        stats = imex_time_step_krylov!(u, params, connij, qe, coords, Δt, lsource;
                                      rtol=rtol, verbose=(verbose && n % output_freq == 0))
        
        # Track statistics
        total_iterations += stats.niter
        max_iterations = max(max_iterations, stats.niter)
        
        # Print progress
        if n % output_freq == 0 || n == 1
            avg_iter = total_iterations / n
            println("Step $n/$ntime_steps: $(stats.niter) GMRES iterations (avg: $(round(avg_iter, digits=1)))")
        end
    end
    
    println()
    println("Time integration complete!")
    println("  Average GMRES iterations: $(round(total_iterations/ntime_steps, digits=2))")
    println("  Maximum GMRES iterations: $max_iterations")
    
    return u
end

"""
Alternative: Use BiCGSTAB instead of GMRES
"""
function solve_imex_system_bicgstab!(LHS_matrix, RHS_vector; 
                                     rtol=1e-8, atol=1e-10, maxiter=500, verbose=false)
    """
    Solve LHS * x = RHS using BiCGSTAB (often faster than GMRES for some problems)
    """
    
    x, stats = bicgstab(LHS_matrix, RHS_vector; 
                       rtol=rtol, atol=atol, itmax=maxiter,
                       verbose=verbose ? 2 : 0)
    
    if verbose
        println("BiCGSTAB converged: $(stats.solved)")
        println("  Iterations: $(stats.niter)")
        println("  Residual norm: $(stats.residuals[end])")
    end
    
    return x, stats
end

"""
Alternative: Use CG (Conjugate Gradient) for symmetric positive definite systems
"""
function solve_imex_system_cg!(LHS_matrix, RHS_vector; 
                               rtol=1e-8, atol=1e-10, maxiter=500, verbose=false)
    """
    Solve LHS * x = RHS using CG (only for SPD matrices!)
    """
    
    x, stats = cg(LHS_matrix, RHS_vector; 
                 rtol=rtol, atol=atol, itmax=maxiter,
                 verbose=verbose ? 2 : 0)
    
    if verbose
        println("CG converged: $(stats.solved)")
        println("  Iterations: $(stats.niter)")
        println("  Residual norm: $(stats.residuals[end])")
    end
    
    return x, stats
end

# Usage example:
"""
# Initialize your problem
params = ...  # Your problem parameters
u = ...       # Initial solution
connij = ...  # Connectivity array
qe = ...      # Jacobian array
coords = ...  # Coordinates
Δt = 0.01     # Time step
nsteps = 100  # Number of time steps

# Run with Krylov solver (GMRES)
time_integrate_krylov!(u, params, connij, qe, coords, Δt, nsteps, false;
                      rtol=1e-8, verbose=true, output_freq=10)

# For a single time step:
stats = imex_time_step_krylov!(u, params, connij, qe, coords, Δt, false;
                              rtol=1e-8, verbose=true)
"""

println("IMEX Krylov solver module loaded successfully!")
println("Available solvers: GMRES (default), BiCGSTAB, CG")

time_integrate_krylov!(u, params, connij, qe, coords, Δt, nsteps, false;
                      rtol=1e-8, verbose=true, output_freq=10)
