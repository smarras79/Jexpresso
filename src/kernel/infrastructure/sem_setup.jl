include("../mesh/restructure_for_periodicity.jl")

# ─── SEM preprocess cache helpers ─────────────────────────────────────────────
# Only `metrics` and `matrix` are cached — pure numerical arrays with no MPI
# state.  `mesh`, `partitioned_model`, and `mesh.parts` are intentionally
# excluded: they contain PartitionedArrays.MPIArray objects whose embedded MPI
# communicator handles are process-local integers that become invalid after
# JLD2 deserialisation, causing MPI_ERR_COMM on the next use.
#
# rank is obtained via MPI inside each helper so these functions are not
# sensitive to the caller passing LinearAlgebra.rank (a Function) instead of
# an integer — which happens because `drivers.jl` has no local `rank` binding
# and `using LinearAlgebra` exports a symbol of the same name.
function _preprocess_cache_path(inputs::Dict, Nξ::Int, Qξ::Int, nparts::Int)
    rank = MPI.Comm_rank(get_mpi_comm())
    dir  = let d = dirname(inputs[:gmsh_filename]); isempty(d) ? "." : d end
    suffix = nparts > 1 ? "_rank$(rank)" : ""
    return joinpath(dir, "PREPROCESS_nop$(Nξ)$(suffix).jld2")
end

function _try_load_sem_cache(path::String)
    rank = MPI.Comm_rank(get_mpi_comm())
    isfile(path) || return (nothing, nothing)
    try
        d = JLD2.load(path)
        return (d["metrics"], d["matrix"])
    catch e
        rank == 0 && @warn "Ignoring unreadable SEM cache $path" exception=(e, catch_backtrace())
        return (nothing, nothing)
    end
end

function _save_sem_cache(path::String, metrics, matrix)
    rank = MPI.Comm_rank(get_mpi_comm())
    try
        JLD2.jldsave(path; metrics, matrix)
        rank == 0 && @info "Saved SEM preprocess cache: $path"
    catch e
        rank == 0 && @warn "Failed to save SEM cache $path" exception=(e, catch_backtrace())
    end
end
# ──────────────────────────────────────────────────────────────────────────────

function sem_setup(inputs::Dict, nparts, distribute, rank, args...)
    
    adapt_flags, partitioned_model_coarse, omesh = _handle_optional_args4amr(args...)

    if inputs[:lwarmup] == true && rank == 0
        inputs[:lwarmup] = false #set it false to prevent it from entering again.

        println(BLUE_FG(string(" # JIT pre-compilation of large problem ...")))
        
        input_mesh             = inputs[:gmsh_filename]
        inputs[:gmsh_filename] = inputs[:gmsh_filename_c]
        sem_dummy              = sem_setup(inputs, nparts, distribute, rank)
        inputs[:gmsh_filename] = input_mesh
        
        # --- MEMORY CLEANUP ---
        sem_dummy = nothing 
        
        println(BLUE_FG(string(" # JIT pre-compilation of large problem ... END")))
               
    end
    
    fx        = zeros(Float64,1,1)
    fy        = zeros(Float64,1,1)
    fz        = zeros(Float64,1,1)
    fy_lag    = zeros(Float64,1,1)
    phys_grid = zeros(Float64,1,1)
    Nξ        = inputs[:nop]
    AD        = inputs[:AD]
    CL        = inputs[:CL]
    
    lexact_integration = inputs[:lexact_integration]
    SOL_VARS_TYPE      = inputs[:SOL_VARS_TYPE]
    volume_flux        = inputs[:volume_flux]
    
    connijk_original          = zeros(TInt,1,1,1,1)
    poin_in_bdy_face_original = zeros(TInt,1,1,1)
    
    x_original = zeros(1,1)
    y_original = zeros(1,1)
    z_original = zeros(1,1)

    #--------------------------------------------------------
    # Create/read mesh
    # return mesh::St_mesh
    # and Build interpolation nodes
    #             the user decides among LGL, GL, etc. 
    # Return:
    # ξ = ND.ξ.ξ
    # ω = ND.ξ.ω
    #--------------------------------------------------------
    if isnothing(adapt_flags)
        mesh, partitioned_model = mod_mesh_mesh_driver(inputs, nparts, distribute)
    else
        mesh, partitioned_model, uaux_new = mod_mesh_mesh_driver(inputs, nparts, distribute, args...)
    end
    if (inputs[:xscale] != 1.0 && inputs[:xdisp] != 0.0)
        mesh.x .= (@view(mesh.x[:]) .+ TFloat(inputs[:xdisp])) .*TFloat(inputs[:xscale]*0.5)
    elseif (inputs[:xscale] != 1.0)
        mesh.x[:] = @view(mesh.x[:])*TFloat(inputs[:xscale]*0.5)
    elseif (inputs[:xdisp] != 0.0)
        mesh.x[:] .= (@view(mesh.x[:]) .+ TFloat(inputs[:xdisp]))
    end
    mesh.xmin = minimum(mesh.coords[:,1])
    mesh.xmax = maximum(mesh.coords[:,1])
    if (inputs[:yscale] != 1.0 && inputs[:ydisp] != 0.0)
        mesh.y[:] .= (mesh.y[:] .+ inputs[:ydisp]) .*inputs[:yscale] * 0.5
    elseif(inputs[:yscale] != 1.0)
        mesh.y[:] .= (mesh.y[:]) .*inputs[:yscale]*0.5
    elseif(inputs[:ydisp] != 0.0)
        mesh.y[:] .= (mesh.y[:] .+ inputs[:ydisp])
    end
    if mesh.nsd == 2
        mesh.ymin = minimum(mesh.coords[:,2])
        mesh.ymax = maximum(mesh.coords[:,2])
    end
    
    #--------------------------------------------------------
    # Build interpolation and quadrature points/weights
    # (reuse the nodes already computed during mesh construction)
    #--------------------------------------------------------
    ξω  = mesh.ξω
    if length(args) > 3
        interp  = args[4]
        project = args[5]
    else
        interp, project = build_projection_1d(ξω.ξ)
    end
    
    ξ,ω = ξω.ξ, ξω.ω
    if lexact_integration
        #
        # Exact quadrature:
        # Quadrature order (Q = N+1) ≠ polynomial order (N)
        #
        QT  = Exact() #Quadrature Type
        QT_String = "Exact"
        Qξ  = Nξ + 1
        
        ξωQ   = basis_structs_ξ_ω!(inputs[:quadrature_nodes], Qξ)
        ξq, ω = ξωQ.ξ, ξωQ.ω
    else  
        #
        # Inexact quadrature:
        # Quadrature and interpolation orders coincide (Q = N)
        #
        QT  = Inexact() #Quadrature Type
        QT_String = "Inexact"
        Qξ  = Nξ
        ξωq = ξω
        ξq  = ξ
        ω   = ξω.ω
        ωb  = barycentric_weights(ξ)
    end
    SD = mesh.SD

    # ── Preprocess cache setup ────────────────────────────────────────────────
    preprocess_cache = _preprocess_cache_path(inputs, Nξ, Qξ, nparts)
    MPI.Comm_rank(get_mpi_comm()) == 0 && @info preprocess_cache
    cached_metrics, cached_matrix = _try_load_sem_cache(preprocess_cache)
    loaded_from_cache = !isnothing(cached_metrics)
    if loaded_from_cache
        metrics = cached_metrics
        matrix  = cached_matrix
        MPI.Comm_rank(get_mpi_comm()) == 0 && @info " Loaded SEM preprocess cache — skipping metric terms and matrix build"
    end
    # ─────────────────────────────────────────────────────────────────────────

    #--------------------------------------------------------
    # Build Lagrange polynomials:
    #
    # Return:
    # ψ     = basis.ψ[N+1, Q+1]
    # dψ/dξ = basis.dψ[N+1, Q+1]
    #--------------------------------------------------------
    if (mesh.nsd > 1)
        
        #
        # 2D/3D grids (from GMSH)
        #
        if (mesh.lLaguerre)
            basis1 = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])
            ξω2 = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta], inputs[:backend])
            ξ2,ω2 = ξω2.ξ, ξω2.ω
            basis2 = build_Interpolation_basis!(ScaledLaguerreBasis(), ξ2, ξ2, inputs[:laguerre_beta], TFloat, inputs[:backend])
            basis = (basis1, basis2)
            ω1 = ω
            ω = (ω1,ω2)
            if (inputs[:lfilter])
                ξω3 = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.ngr-1, inputs[:backend])
                ξ3,ω3 = ξω3.ξ, ξω3.ω
                if (inputs[:backend] == CPU())
                    fx = init_filter(mesh.ngl-1,ξ,inputs[:mu_x],mesh,inputs, rank)
                    fy = init_filter(mesh.ngl-1,ξ,inputs[:mu_y],mesh,inputs, rank)
                    fy_lag = init_filter(mesh.ngr-1,ξ2,inputs[:mu_y],mesh,inputs, rank)
                else
                    ξ_gl = KernelAbstractions.zeros(CPU(), Float64, Int64(mesh.ngl))
                    KernelAbstractions.copyto!(CPU(),ξ_gl,ξ)
                    ξ_gr = KernelAbstractions.zeros(CPU(), Float64, Int64(mesh.ngr))
                    KernelAbstractions.copyto!(CPU(),ξ_gr,ξ2)
                    fx_1 = init_filter(mesh.ngl-1,ξ_gl,inputs[:mu_x],mesh,inputs, rank)
                    fy_1 = init_filter(mesh.ngl-1,ξ_gl,inputs[:mu_y],mesh,inputs, rank)
                    fy_lag_1 = init_filter(mesh.ngr-1,ξ_gr,inputs[:mu_y],mesh,inputs, rank)
                    fx = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                    fy = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                    fy_lag = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngr), Int64(mesh.ngr))
                    KernelAbstractions.copyto!(inputs[:backend], fx, fx_1)
                    KernelAbstractions.copyto!(inputs[:backend], fy, fy_1)
                    KernelAbstractions.copyto!(inputs[:backend], fy_lag, fy_lag_1)
                end
            end
            
            #if (inputs[:lwarp]) warp_mesh!(mesh,inputs) end
            
            if !loaded_from_cache
                if (rank == 0) @info " Build metrics ......" end
                metrics1 = allocate_metrics(SD, mesh.nelem, mesh.nedges_bdy, Qξ, TFloat, inputs[:backend])
                @time build_metric_terms!(metrics1, mesh, basis1, Nξ, Qξ, ξ, ω1, TFloat, COVAR(), SD; backend = inputs[:backend])

                metrics2 = allocate_metrics_laguerre(SD, mesh.nelem_semi_inf, mesh.nedges_bdy, Qξ, mesh.ngr, TFloat, inputs[:backend])
                build_metric_terms!(metrics2, mesh, basis1, basis2, Nξ, Qξ, mesh.ngr, mesh.ngr, ξ, ω1, ω2, TFloat, COVAR(), SD; backend = inputs[:backend])

                metrics = (metrics1, metrics2)
                if (rank == 0) @info " Build metrics ...... DONE" end

                matrix = matrix_wrapper_laguerre(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat;
                                                 ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend], interp)
                _save_sem_cache(preprocess_cache, metrics, matrix)
            end
            
        else
            if (rank == 0) @info " Build interpolation bases ......" end
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])
            if (rank == 0) @info " Build interpolation bases ...... END" end
            ω1 = ω
            ω = ω1
            if (inputs[:lfilter])
                if (inputs[:backend] == CPU())
                    fx = init_filter(mesh.ngl-1,ξ,inputs[:mu_x],mesh,inputs, rank)
                    fy = init_filter(mesh.ngl-1,ξ,inputs[:mu_y],mesh,inputs, rank)
                    if (mesh.nsd >2)
                        fz = init_filter(mesh.ngl-1,ξ,inputs[:mu_z],mesh,inputs, rank)
                    end
                else
                    ξ_temp = KernelAbstractions.zeros(CPU(), Float64, Int64(mesh.ngl))
                    KernelAbstractions.copyto!(CPU(),ξ_temp,ξ)
                    fx_1 = init_filter(mesh.ngl-1,ξ_temp,inputs[:mu_x],mesh,inputs, rank)
                    fy_1 = init_filter(mesh.ngl-1,ξ_temp,inputs[:mu_y],mesh,inputs, rank)
                    fx = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                    fy = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                    KernelAbstractions.copyto!(inputs[:backend], fx, fx_1)
                    KernelAbstractions.copyto!(inputs[:backend], fy, fy_1)
                    if (mesh.nsd > 2)
                        fz_1 = init_filter(mesh.ngl-1,ξ_temp,inputs[:mu_z],mesh,inputs, rank)
                        fz = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                        KernelAbstractions.copyto!(inputs[:backend], fz, fz_1)
                    end
                end
            end
            if (inputs[:lphysics_grid])
                phys_grid = init_phys_grid(mesh, inputs,inputs[:nlay_pg],inputs[:nx_pg],inputs[:ny_pg],mesh.xmin,mesh.xmax,mesh.ymin,mesh.ymax,mesh.zmin,mesh.zmax,inputs[:backend])
            end

            if !loaded_from_cache
                #--------------------------------------------------------
                # Build metric terms
                #--------------------------------------------------------
                if (rank == 0) @info " Build metrics ......" end
                metrics = allocate_metrics(SD, mesh.nelem, mesh.nedges_bdy, Qξ, TFloat, inputs[:backend])
                @time build_metric_terms!(metrics, mesh, basis, Nξ, Qξ, ξ, ω, TFloat, COVAR(), SD; backend = inputs[:backend])
                if (rank == 0) @info " Build metrics ...... END" end

                if (rank == 0) @info " Matrix wrapper ......" end
                matrix = matrix_wrapper(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace],
                            ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend], interp)
                if (rank == 0) @info " Matrix wrapper ...... END" end
                _save_sem_cache(preprocess_cache, metrics, matrix)
            end
        end
    else
        #
        # 1D grids (native)
        #
        if(inputs[:llaguerre_1d_right] || inputs[:llaguerre_1d_left])

            basis1 = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])
            ξω2    = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta], inputs[:backend])
            ξ2,ω2  = ξω2.ξ, ξω2.ω
            basis2 = build_Interpolation_basis!(ScaledLaguerreBasis(), ξ2, ξ2, inputs[:laguerre_beta], TFloat, inputs[:backend])
            basis  = (basis1, basis2)
            ω1     = ω
            ω      = (ω1,ω2)
            if !loaded_from_cache
                #--------------------------------------------------------
                # Build metric terms
                #--------------------------------------------------------
                if (rank == 0) @info " Build metrics ......" end
                metrics1 = allocate_metrics(SD, mesh.nelem, mesh.nedges_bdy, Qξ, TFloat, inputs[:backend])
                build_metric_terms!(metrics1, mesh, basis[1], Nξ, Qξ, ξ, ω, TFloat, COVAR(), SD; backend = inputs[:backend])
                metrics2 = allocate_metrics(SD, mesh.nelem_semi_inf, mesh.nedges_bdy, mesh.ngr, TFloat, inputs[:backend])
                build_metric_terms_1D_Laguerre!(metrics2, mesh, basis[2], mesh.ngr, mesh.ngr, ξ2, ω2, inputs, TFloat, COVAR(), SD; backend = inputs[:backend])

                metrics = (metrics1, metrics2)
                if (rank == 0) @info " Build metrics ...... DONE" end
                matrix = matrix_wrapper_laguerre(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend], interp)
                _save_sem_cache(preprocess_cache, metrics, matrix)
            end
        else
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])

            ω1 = ω
            ω = ω1

            if (inputs[:lperiodic_1d])
                @time restructure4periodicity_1D!(mesh, mesh.coords,
                                                  mesh.xmax, mesh.xmin,mesh.ymax,mesh.ymin,mesh.zmax,mesh.zmin,
                                                  mesh.poin_in_bdy_face,mesh.poin_in_bdy_edge,
                                                  mesh.ngl,mesh.ngr,mesh.nelem,mesh.npoin,mesh.nsd,mesh.bdy_edge_type,
                                                  mesh.bdy_face_type,mesh.bdy_face_in_elem,mesh.bdy_edge_in_elem,
                                                  mesh.connijk,mesh.connijk_lag,mesh.npoin_linear,mesh.nelem_semi_inf,
                                                  inputs,inputs[:backend])
            end

            if !loaded_from_cache
                #--------------------------------------------------------
                # Build metric terms
                #--------------------------------------------------------
                metrics = allocate_metrics(SD, mesh.nelem, mesh.nedges_bdy, Qξ, TFloat, inputs[:backend])
                @time build_metric_terms!(metrics, mesh, basis, Nξ, Qξ, ξ, ω, TFloat, COVAR(), SD; backend = inputs[:backend])

                matrix = matrix_wrapper(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat;
                                        ldss_laplace=inputs[:ldss_laplace],
                                        ldss_differentiation=inputs[:ldss_differentiation],
                                        backend = inputs[:backend], interp)
                _save_sem_cache(preprocess_cache, metrics, matrix)
            end
        end
    end

    if inputs[:lelementLearning]
        #
        # Element learning
        #
        new_matrix      = remove_arrays!(mesh.poin_in_edge, mesh.poin_in_bdy_edge)
        new_no_bdy_poin = replace_shared_values!(new_matrix, mesh.poin_in_bdy_edge)
        
        mesh.∂O       = unroll_positive_unique(new_no_bdy_poin)
        mesh.Γ        = unroll_positive_unique(mesh.poin_in_bdy_edge)
        mesh.∂τ       = vcat(mesh.Γ, mesh.∂O)
        mesh.Io       = view(mesh.internal_poin_in_elem, :)
        mesh.lengthIo = length(mesh.Io)
        mesh.O        = vcat(mesh.Io, mesh.∂O)
        mesh.length∂O = length(mesh.∂O)
        mesh.length∂τ = length(mesh.∂τ)
        mesh.lengthΓ  = length(mesh.Γ)
    end
    
    #--------------------------------------------------------
    # Build matrices
    #--------------------------------------------------------
    if isnothing(adapt_flags)
        return (; QT, CL, AD, SOL_VARS_TYPE, volume_flux, mesh, metrics, basis, ωb, ω, ξ, matrix, fx, fy, fy_lag, fz, phys_grid, 
                connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original, interp, project, nparts, distribute), partitioned_model
    else
        return (; QT, CL, AD, SOL_VARS_TYPE, volume_flux, mesh, metrics, basis, ωb, ω, ξ, matrix, fx, fy, fy_lag, fz, phys_grid, 
                connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original, interp, project, nparts, distribute), partitioned_model, uaux_new
    end
    @info " SEM _SETUP COMPLETE!!!"
end
