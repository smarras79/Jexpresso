# imex.jl - IMEX functionality for Jexpresso

using DifferentialEquations
using LinearAlgebra
using SparseArrays

"""
        global_index(iel, i, j, ieq, ngl, neqs)
        
    Convert element-local indices to global DOF index
    """
function global_index(iel, i, j, ieq, ngl, neqs)
    return (iel - 1) * ngl * ngl * neqs + (j - 1) * ngl * neqs + (i - 1) * neqs + ieq
end

"""
        create_jacobian_prototype(p)
        
    Create sparse matrix prototype for Jacobian structure
    This helps the Newton-Krylov solver understand sparsity pattern
    """
function create_jacobian_prototype(p)
    neqs = p.qp.neqs
    ngl = p.mesh.ngl
    nelem = p.mesh.nelem
    n = neqs * ngl * ngl * nelem
    
    # Create sparsity pattern based on element connectivity
    I = Int[]
    J = Int[]
    
    for iel = 1:nelem
        for ieq = 1:neqs
            for j = 1:ngl, i = 1:ngl
                # Global index for this DOF
                idx_global = global_index(iel, i, j, ieq, ngl, neqs)
                
                # Local coupling within element
                for jeq = 1:neqs
                    for jj = 1:ngl, ii = 1:ngl
                        jdx_global = global_index(iel, ii, jj, jeq, ngl, neqs)
                        push!(I, idx_global)
                        push!(J, jdx_global)
                    end
                end
            end
        end
    end
    
    # Create sparse matrix prototype
    V = ones(length(I))
    return sparse(I, J, V, n, n)
end

"""
        get_imex_algorithm(inputs, jac_prototype)
        
    Select appropriate IMEX algorithm based on user inputs
    """
function get_imex_algorithm(inputs, jac_prototype)
    
    imex_scheme = get(inputs, :imex_scheme, :KenCarp4)
    linear_solver = get(inputs, :imex_linear_solver, :GMRES)

    # Configure linear solver - use DifferentialEquations.jl defaults to avoid conflicts
    linsolve_options = Dict{Symbol, Any}()
    #=    
    if haskey(inputs, :gmres_rtol)
        linsolve_options[:rtol] = inputs[:gmres_rtol]
    end
    if haskey(inputs, :gmres_atol) 
        linsolve_options[:atol] = inputs[:gmres_atol]
    end
    if haskey(inputs, :gmres_maxiter)
        linsolve_options[:maxiter] = inputs[:gmres_maxiter] 
    end
    
    # Store solver options for later use
    solver_options = (
        rtol = get(inputs, Symbol("$(lowercase(string(linear_solver)))_rtol"), 1e-6),
        atol = get(inputs, Symbol("$(lowercase(string(linear_solver)))_atol"), 1e-12),
        maxiter = get(inputs, Symbol("$(lowercase(string(linear_solver)))_maxiter"), 100)
    )
    
    # Select IMEX scheme - jac_prototype goes to problem, not algorithm
    if imex_scheme == :KenCarp4
        algorithm = KenCarp4(autodiff=false, concrete_jac=nothing)
    elseif imex_scheme == :KenCarp3
        algorithm = KenCarp3(autodiff=false, concrete_jac=nothing)
    elseif imex_scheme == :ARKODE
        algorithm = ARKODE()
    else
        @warn "Unknown IMEX scheme $imex_scheme, defaulting to KenCarp4"
        algorithm = KenCarp4(autodiff=false, concrete_jac=nothing)
    end=#

    algorithm = SSPRK54()
    return algorithm
end

function get_imex_algorithm_or(inputs, jac_prototype)
    
    imex_scheme = get(inputs, :imex_scheme, :KenCarp4)
    linear_solver = get(inputs, :imex_linear_solver, :GMRES)
    
    # Configure linear solver - use DifferentialEquations.jl defaults to avoid conflicts
    linsolve_options = Dict{Symbol, Any}()
    
    if haskey(inputs, :gmres_rtol)
        linsolve_options[:rtol] = inputs[:gmres_rtol]
    end
    if haskey(inputs, :gmres_atol) 
        linsolve_options[:atol] = inputs[:gmres_atol]
    end
    if haskey(inputs, :gmres_maxiter)
        linsolve_options[:maxiter] = inputs[:gmres_maxiter] 
    end
    
    # Store solver options for later use
    solver_options = (
        rtol = get(inputs, Symbol("$(lowercase(string(linear_solver)))_rtol"), 1e-6),
        atol = get(inputs, Symbol("$(lowercase(string(linear_solver)))_atol"), 1e-12),
        maxiter = get(inputs, Symbol("$(lowercase(string(linear_solver)))_maxiter"), 100)
    )
    
    # Select IMEX scheme - jac_prototype goes to problem, not algorithm
    if imex_scheme == :KenCarp4
        algorithm = KenCarp4(autodiff=false, concrete_jac=nothing)
    elseif imex_scheme == :KenCarp3
        algorithm = KenCarp3(autodiff=false, concrete_jac=nothing)
    elseif imex_scheme == :ARKODE
        algorithm = ARKODE()
    else
        @warn "Unknown IMEX scheme $imex_scheme, defaulting to KenCarp4"
        algorithm = KenCarp4(autodiff=false, concrete_jac=nothing)
    end
    
    return algorithm
end

"""
        create_imex_functions(p, acoustic_implicit=true)
        
    Create IMEX RHS functions for stratified atmospheric flow
    acoustic_implicit: if true, treat acoustic waves implicitly
    """
function _build_rhs_explicit!(RHS, u, p, time,
                             acoustic_implicit, F_ex, G_ex, F_im, G_im, S_cache)

    T       = Float64
    SD      = p.SD
    VT      = p.VT
    QT      = p.QT
    CL      = p.CL
    AD      = p.AD
    neqs    = p.neqs
    ngl     = p.mesh.ngl
    nelem   = p.mesh.nelem
    npoin   = p.mesh.npoin
    lsource = p.inputs[:lsource]
    xmin    = p.mesh.xmin
    xmax    = p.mesh.xmax
    ymin    = p.mesh.ymin
    ymax    = p.mesh.ymax
    zmin    = p.mesh.zmin
    zmax    = p.mesh.zmax    

    if SD == NSD_1D()
        comm = MPI.COMM_WORLD
    else
        comm = p.mesh.parts.comm
    end
    mpisize = MPI.Comm_size(comm)

    #-----------------------------------------------------------------------------------
    # Inviscid rhs:
    #-----------------------------------------------------------------------------------
    resetRHSToZero_inviscid!(p)
    
    u2uaux!(@view(p.uaux[:,:]), u, p.neqs, p.mesh.npoin)

    resetbdyfluxToZero!(p)  
    
    apply_boundary_conditions_dirichlet!(u, p.uaux, time, p.qp.qe,
                                         p.mesh.x, p.mesh.y, p.mesh.z, 
                                         p.metrics.nx, p.metrics.ny, p.metrics.nz, p.mesh.npoin, p.mesh.npoin_linear, 
                                         p.mesh.poin_in_bdy_edge, p.mesh.poin_in_bdy_face, p.mesh.nedges_bdy, p.mesh.nfaces_bdy, p.mesh.ngl, 
                                         p.mesh.ngr, p.mesh.nelem_semi_inf, p.basis.ψ, p.basis.dψ,
                                         xmax, ymax, zmax, xmin, ymin, zmin, p.RHS, p.rhs_el, p.ubdy,
                                         p.mesh.connijk_lag, p.mesh.bdy_edge_in_elem, 
                                         p.mesh.bdy_edge_type, p.mesh.bdy_face_in_elem, p.mesh.bdy_face_type,
                                         p.mesh.connijk, p.metrics.Jef, p.S_face, 
                                         p.S_flux, p.F_surf, p.M_surf_inv, p.M_edge_inv, p.Minv,
                                         p.mp.Tabs, p.mp.qn,
                                         p.ω, neqs, p.inputs, AD, SD)

    
   
    #= Compute explicit terms
    for iel = 1:nelem
    compute_split_fluxes!(F_ex, G_ex, F_im, G_im, S_cache, 
    iel, p, t, acoustic_implicit)
    
    # Apply explicit expansion
    expansion_split!(u, neqs, ngl,
    p.basis.dψ, p.ω,
    F_ex, G_ex, S_cache, true,  # explicit=true
    p.metrics.Je,
    p.metrics.dξdx, p.metrics.dξdy,
    p.metrics.dηdx, p.metrics.dηdy,
    rhs_el, iel, p.CL, p.QT, p.SD, p.AD)
    
    end=#

    imex_inviscid_rhs_el!(u, p, p.mesh.connijk, p.qp.qe,
                          acoustic_implicit, F_ex, G_ex, F_im, G_im, S_cache,
                          p.mesh.x, p.mesh.y, p.mesh.z,
                          lsource, 
                          p.mp.S_micro,
                          p.mp.qn,
                          p.mp.flux_lw, p.mp.flux_sw,
                          SD) #QUI

    outfile = "u_IMEX.txt"
    open(outfile, "w") do f
        for i in u
            @printf(f, "%16.8f\n", i)
        end
    end
    #@mystop("STOP IMEX")
    
    #RHStoDU!(du, @view(p.RHS[:,:]), p.neqs, p.mesh.npoin)
    
    DSS_rhs!(p.RHS, p.rhs_el, p.mesh.connijk, nelem, ngl, neqs, SD, AD)

    apply_boundary_conditions_neumann!(u, p.uaux, time, p.qp.qe,
                                       p.mesh.x, p.mesh.y, p.mesh.z, 
                                       p.metrics.nx, p.metrics.ny, p.metrics.nz, p.mesh.npoin, p.mesh.npoin_linear,
                                       p.mesh.poin_in_bdy_edge, p.mesh.poin_in_bdy_face, p.mesh.nedges_bdy, p.mesh.nfaces_bdy, p.mesh.ngl,
                                       p.mesh.ngr, p.mesh.nelem_semi_inf, p.basis.ψ, p.basis.dψ,
                                       xmax, ymax, zmax, xmin, ymin, zmin, p.RHS, p.rhs_el, p.ubdy,
                                       p.mesh.connijk_lag, p.mesh.bdy_edge_in_elem, 
                                       p.mesh.bdy_edge_type, p.mesh.bdy_face_in_elem, p.mesh.bdy_face_type,
                                       p.mesh.connijk, p.metrics.Jef, p.S_face, 
                                       p.S_flux, p.F_surf, p.M_surf_inv, p.M_edge_inv, p.Minv,
                                       p.WM.τ_f, p.WM.wθ,
                                       p.mp.Tabs, p.mp.qn,
                                       p.ω, neqs, p.inputs, AD, SD) 
    
    DSS_global_RHS!(@view(p.RHS[:,:]), p.pM, p.neqs)

    for ieq=1:neqs
        divide_by_mass_matrix!(@view(p.RHS[:,ieq]), p.vaux, p.Minv, neqs, npoin, AD)
    end
    
    
    outfile = "RHS_IMEX_DSS.txt"
    open(outfile, "w") do f
        for i in p.RHS
            @printf(f, "%16.8f\n", i)
        end
    end

    outfile = "u_IMEX_DSS.txt"
    open(outfile, "w") do f
        for i in u
            @printf(f, "%16.8f\n", i)
        end
    end
   # @mystop("STOP IMEX")

    
end

function imex_inviscid_rhs_el!(u, p, connijk::Array{Int64,4}, qe::Matrix{Float64},
                               acoustic_implicit, F_ex, G_ex, F_im, G_im, S,
                               x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, 
                               lsource, S_micro_vec, qn_vec, flux_lw_vec,
                               flux_sw_vec, SD::NSD_2D)

    xmin = p.xmin; xmax = p.xmax; ymax = p.ymax

    ngl = p.mesh.ngl
    neqs = p.qp.neqs
    
    for iel = 1:p.mesh.nelem
        fill!(F_ex, 0.0); fill!(G_ex, 0.0)
        fill!(F_im, 0.0); fill!(G_im, 0.0)
        fill!(S, 0.0)
        
        for j = 1:p.mesh.ngl, i=1:p.mesh.ngl
            ip = connijk[iel,i,j,1]

            # Compute total flux (reuse existing user_flux! function)
            user_flux!(@view(p.F[i,j,:]), @view(p.G[i,j,:]), p.SD,
                       @view(p.uaux[ip,:]),
                       @view(p.qp.qe[ip,:]),
                       p.mesh,
                       p.CL, p.SOL_VARS_TYPE;
                       neqs=neqs, ip=ip)
            
            # Split fluxes based on physics
            split_atmospheric_fluxes!(@view(F_ex[i,j,:]), @view(G_ex[i,j,:]),
                                      @view(F_im[i,j,:]), @view(G_im[i,j,:]),
                                      @view(p.F[i,j,:]), @view(p.G[i,j,:]),
                                      @view(p.uaux[ip,:]),
                                      acoustic_implicit, p.SOL_VARS_TYPE)
            
            # Compute source terms (explicit treatment)
            if lsource
                user_source!(@view(S[i,j,:]),
                             @view(p.uaux[ip,:]),
                             @view(p.qp.qe[ip,:]),
                             p.mesh.npoin, p.CL, p.SOL_VARS_TYPE;
                             neqs=neqs, 
                             x=0.0,#x[ip], 
                             y=0.0,#y[ip], 
                             xmax=p.xmax, 
                             xmin=p.xmin, 
                             ymax=p.ymax)
            end
            
            # Apply explicit expansion
            expansion_split!(u, neqs, ngl,
                             p.basis.dψ, p.ω,
                             F_ex, G_ex, S, true, #explicit=true, 
                             p.metrics.Je,
                             p.metrics.dξdx, p.metrics.dξdy,
                             p.metrics.dηdx, p.metrics.dηdy,
                             p.rhs_el, iel, p.CL, p.QT, p.SD, p.AD)
        end
    end
    
end

function create_imex_functions(p, acoustic_implicit=true)
    
    # Get problem dimensions
    neqs  = p.qp.neqs
    ngl   = p.mesh.ngl  
    nelem = p.mesh.nelem
    npoin = p.mesh.npoin
    
    resetRHSToZero_inviscid!(p)
    
    # Preallocate working arrays
    F_ex = zeros(ngl, ngl, neqs)  # Explicit flux in ξ
    G_ex = zeros(ngl, ngl, neqs)  # Explicit flux in η
    F_im = zeros(ngl, ngl, neqs)  # Implicit flux in ξ  
    G_im = zeros(ngl, ngl, neqs)  # Implicit flux in η
    S_cache = zeros(ngl, ngl, neqs)  # Source terms
    
    # Create Jacobian prototype
    jac_prototype = create_jacobian_prototype(p)
    
    """
        Explicit RHS function
    """
    function rhs_explicit!(du, u, p, t) #analogue of create_rhs! in the fully explicit case (rhs.jl)
        #fill!(du, 0.0)

        _build_rhs_explicit!(@view(p.RHS[:,:]), u, p, t,
                            acoustic_implicit, F_ex, G_ex, F_im, G_im, S_cache)

        
        RHStoDU!(du, @view(p.RHS[:,:]), p.neqs, p.mesh.npoin)
    
        #= Compute explicit terms
        for iel = 1:nelem
        compute_split_fluxes!(F_ex, G_ex, F_im, G_im, S_cache, 
        iel, p, t, acoustic_implicit)
        
        # Apply explicit expansion
        expansion_split!(u, neqs, ngl,
        p.basis.dψ, p.ω,
        F_ex, G_ex, S_cache, true,  # explicit=true
        p.metrics.Je,
        p.metrics.dξdx, p.metrics.dξdy,
        p.metrics.dηdx, p.metrics.dηdy,
        rhs_el, iel, p.CL, p.QT, p.SD, p.AD)
        
        end
        RHStoDU!(du, @view(p.RHS[:,:]), p.neqs, p.mesh.npoin)
        =#
        
        return nothing
    end
    
    """
        Implicit RHS function
            """  
    function rhs_implicit!(du, u, p, t)
        return nothing
        #fill!(du, 0.0)

        #_build_rhs_explicit!(@view(p.RHS[:,:]), u, p, t,
        #                    acoustic_implicit,
        #                    F_ex, G_ex, F_im, G_im, S_cache)

        #RHStoDU!(du, @view(p.RHS[:,:]), p.neqs, p.mesh.npoin)
        
        # Compute implicit terms
        #= for iel = 1:nelem
        compute_split_fluxes!(F_ex, G_ex, F_im, G_im, S_cache,
        iel, p, t, acoustic_implicit)
        
        # Apply implicit expansion (no source terms)
        expansion_split!(u, neqs, ngl,
        p.basis.dψ, p.ω,
        F_im, G_im, zeros(ngl, ngl, neqs), false,  # explicit=false
        p.metrics.Je,
        p.metrics.dξdx, p.metrics.dξdy,
        p.metrics.dηdx, p.metrics.dηdy,
        du, iel, p.CL, p.QT, p.SD, p.AD)
        end
        =#
        return nothing
    end
    
    return rhs_explicit!, rhs_implicit!, jac_prototype
end

"""
        compute_split_fluxes!(F_ex, G_ex, F_im, G_im, S, u_el, iel, p, t, acoustic_implicit)
        
    Compute flux splitting for IMEX integration
    """
function compute_split_fluxes!(F_ex, G_ex, F_im, G_im, S,
                               iel, p, t, acoustic_implicit)

    ###
    ### !!!!!
    #### NOT SUPPOSED TO BE USED ANYMORE. I moved this to part of _build_rhs_explicit!"
    ### !!!!!
    ###

    #=
    ngl = p.mesh.ngl
    neqs = p.qp.neqs
    
    fill!(F_ex, 0.0); fill!(G_ex, 0.0)
    fill!(F_im, 0.0); fill!(G_im, 0.0)
    fill!(S, 0.0)
    
    for j = 1:ngl, i = 1:ngl
        ip = p.mesh.connijk[iel, i, j, 1]
        
        # Compute total flux (reuse existing user_flux! function)
        user_flux!(@view(p.F[i,j,:]), @view(p.G[i,j,:]), p.SD,
                   @view(p.uaux[ip,:]),
                   @view(p.qp.qe[ip,:]),
                   p.mesh,
                   p.CL, p.SOL_VARS_TYPE;
                   neqs=neqs, ip=ip)
        
        # Split fluxes based on physics
        split_atmospheric_fluxes!(@view(F_ex[i,j,:]), @view(G_ex[i,j,:]),
                                  @view(F_im[i,j,:]), @view(G_im[i,j,:]),
                                  @view(p.F[i,j,:]), @view(p.G[i,j,:]),
                                  @view(p.uaux[ip,:]),
                                  acoustic_implicit, p.SOL_VARS_TYPE)
        
        # Compute source terms (explicit treatment)
        lsource = true
        if lsource
            user_source!(@view(S[i,j,:]),
                         @view(p.uaux[ip,:]),
                         @view(p.qp.qe[ip,:]),
                         p.mesh.npoin, p.CL, p.SOL_VARS_TYPE;
                         neqs=neqs, 
                         x=0.0, #p.x[ip], 
                         y=0.0, #p.y[ip], 
                         xmax=p.xmax, 
                         xmin=p.xmin, 
                         ymax=p.ymax)
        end
    end
    =#
    return nothing
end

"""
        split_atmospheric_fluxes!(F_ex, G_ex, F_im, G_im, F_total, G_total, u_prim, 
                                 acoustic_implicit, SOL_VARS_TYPE)
                                 
    Split total fluxes into explicit (advective) and implicit (acoustic) parts
    """
function split_atmospheric_fluxes!(F_ex, G_ex, F_im, G_im, F_total, G_total, u_prim,
                                   acoustic_implicit, SOL_VARS_TYPE)

    if acoustic_implicit
        # Acoustic terms implicit, advection explicit
        ρ = u_prim[1]
        u = u_prim[2]
        v = u_prim[3] 
        p = u_prim[5]  # pressure
        
        # Explicit fluxes (advective nonlinear terms)
        F_ex[1] = 0.0               # Mass continuity (acoustic part implicit)
        F_ex[2] = ρ * u * u         # Nonlinear advection ρu²
        F_ex[3] = ρ * u * v         # Cross advection ρuv
        F_ex[4] = F_total[4]        # Full potential temperature flux
        
        G_ex[1] = 0.0               # Mass continuity (acoustic part implicit)  
        G_ex[2] = ρ * v * u         # Cross advection ρvu
        G_ex[3] = ρ * v * v         # Nonlinear advection ρv²
        G_ex[4] = G_total[4]        # Full potential temperature flux
        
        # Implicit fluxes (linear acoustic terms)
        F_im[1] = ρ * u             # Mass flux (linear)
        F_im[2] = p                 # Pressure gradient
        F_im[3] = 0.0               # No cross-pressure term
        F_im[4] = 0.0               # No acoustic θ flux
        
        G_im[1] = ρ * v             # Mass flux (linear)
        G_im[2] = 0.0               # No cross-pressure term  
        G_im[3] = p                 # Pressure gradient
        G_im[4] = 0.0               # No acoustic θ flux
        
    else
        # Full explicit (original behavior)
        F_ex[:] = F_total[:]
        G_ex[:] = G_total[:]
        fill!(F_im, 0.0)
        fill!(G_im, 0.0)
    end
end

"""
        expansion_split!(u, neqs, ngl, dψ, ω, F, G, S, explicit, Je, dξdx, dξdy, dηdx, dηdy,
                        rhs_el, iel, CL, QT, SD, AD)
                        
    Spectral element expansion for split (IMEX) terms
    """
function expansion_split!(u, neqs, ngl, dψ, ω, F, G, S, explicit,
                          Je, dξdx, dξdy, dηdx, dηdy,
                          rhs_el, iel,  ::CL, QT::Inexact, SD::NSD_2D, AD::ContGal)
    
    for ieq = 1:neqs
        for j = 1:ngl
            for i = 1:ngl
                ωJac = ω[i] * ω[j] * Je[iel, i, j]

                dFdξ = 0.0; dFdη = 0.0
                dGdξ = 0.0; dGdη = 0.0
                
                @turbo for k = 1:ngl
                    dFdξ += dψ[k, i] * F[k, j, ieq]
                    dFdη += dψ[k, j] * F[i, k, ieq]
                    
                    dGdξ += dψ[k, i] * G[k, j, ieq]
                    dGdη += dψ[k, j] * G[i, k, ieq]
                end
                
                dξdx_ij = dξdx[iel, i, j]; dξdy_ij = dξdy[iel, i, j]
                dηdx_ij = dηdx[iel, i, j]; dηdy_ij = dηdy[iel, i, j]

                dFdx = dFdξ * dξdx_ij + dFdη * dηdx_ij
                dGdy = dGdξ * dξdy_ij + dGdη * dηdy_ij
                
                if explicit
                    # Include source terms for explicit part
                    auxi = ωJac * ((dFdx + dGdy) - S[i, j, ieq])
                else
                    # No source terms for implicit part
                    auxi = ωJac * (dFdx + dGdy)
                end
                
                rhs_el[iel, i, j, ieq] -= auxi
            end
        end
    end
end

"""
        solve_imex_problem(prob, algorithm, inputs, callback, dosetimes, rank)
        
    Solve IMEX problem with appropriate time stepping and tolerances
    """
function solve_imex_problem(prob, algorithm, inputs, callback, dosetimes, rank)
    
    # IMEX-specific time stepping parameters
    dt_initial = get(inputs, :imex_dt_initial, inputs[:Δt])
    reltol = get(inputs, :imex_reltol, 1e-6)
    abstol = get(inputs, :imex_abstol, 1e-8)
    adaptive = get(inputs, :imex_adaptive, true)
    
    println(" #   IMEX dt_initial = ", dt_initial)
    println(" #   IMEX tolerances: rel=", reltol, ", abs=", abstol)
    
    # Solve with IMEX algorithm - only use supported solve() parameters
    solution = solve(prob, algorithm;
                     dt = dt_initial,
                     callback = CallbackSet(callback),
                     tstops = dosetimes,
                     save_everystep = false,
                     adaptive = adaptive,
                     reltol = reltol,
                     abstol = abstol,
                     saveat = range(inputs[:tinit], inputs[:tend], 
                                    length=inputs[:ndiagnostics_outputs]),
                     maxiters = get(inputs, :imex_maxiters, 1e6))
    
    return solution
end

"""
        amr_strategy_imex!(inputs, p, u, t, partitioned_model)
        
    AMR strategy adapted for IMEX problems
    """
function amr_strategy_imex!(inputs, p, u, t, partitioned_model)
    # Call existing AMR but create SplitODEProblem instead of ODEProblem
    # This is a placeholder - implement based on your existing AMR code
    
    new_p, new_partitioned_model = amr_strategy!(inputs, p, u, t, partitioned_model)
    
    # Create new IMEX functions
    f_ex!, f_im!, jac_prototype = create_imex_functions(new_p, 
                                                        get(inputs, :imex_acoustic_implicit, true))
    
    # Create new SplitODEProblem
    new_prob = SplitODEProblem(f_ex!, f_im!, u, (t, inputs[:tend]), new_p)
    
    return new_prob, new_partitioned_model
end
