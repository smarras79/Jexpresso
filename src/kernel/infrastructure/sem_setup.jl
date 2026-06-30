include("../mesh/restructure_for_periodicity.jl")

# ─── SEM preprocess cache helpers ─────────────────────────────────────────────
# Only `metrics` and `matrix` are cached — pure numerical arrays with no MPI
# state.  Cache file path is computed by _preprocess_cache_path (see
# coupling/couplingStructs.jl) and is per-rank in parallel runs, so each rank
# reads/writes its own slice without contention.  JEXPRESSO_PREFETCHED_SEM_CACHE
# (populated by je_prefetch_caches! before the with_mpi block) is honoured
# when Alya coupling is active, keeping JLD2 JIT + disk I/O off the
# Alya-blocking path.  rank is obtained via MPI inside each helper so callers
# don't need to thread it through.
function _try_load_sem_cache(path::String; gmsh_path::String="",
                              inputs=nothing, nparts::Int=1)
    rank = MPI.Comm_rank(get_mpi_comm())
    # Prefer the prefetched payload before any disk I/O.
    if JEXPRESSO_PREFETCHED_SEM_CACHE[] !== nothing
        return JEXPRESSO_PREFETCHED_SEM_CACHE[]
    end
    # Pre-load validity check via fingerprint-only read. See
    # _check_cache_validity (couplingStructs.jl) for the rationale -
    # this dodges JLD2 reconstruct failures from custom-struct shape
    # changes and auto-deletes the stale file.
    if inputs !== nothing
        valid, _ = _check_cache_validity(path, inputs, nparts; gmsh_path=gmsh_path)
        valid || return (nothing, nothing)
    else
        isfile(path) || return (nothing, nothing)
        if !isempty(gmsh_path) && _cache_is_stale(path, gmsh_path)
            rank == 0 && println(" # SEM cache $path is older than $gmsh_path — discarding stale cache")
            return (nothing, nothing)
        end
    end
    try
        d = JLD2.load(path)
        if inputs !== nothing
            if !haskey(d, "fingerprint")
                rank == 0 && println(" # SEM cache $path has no fingerprint — discarding and rebuilding")
                return (nothing, nothing)
            end
            saved_fp = d["fingerprint"]
            if !(saved_fp isa Dict) || !_cache_fingerprint_matches(saved_fp, inputs, nparts)
                rank == 0 && println(" # SEM cache $path fingerprint mismatch — discarding and rebuilding")
                return (nothing, nothing)
            end
        end
        # Re-inject g_dss_cache=nothing into the loaded NamedTuple.
        # We strip it on save (it holds MPI handles - MPI.Comm, MultiRequest -
        # that are tied to the current MPI session and aren't safe to
        # serialize). It's also dead in `matrix`: params_setup.jl rebuilds
        # its own g_dss_cache via setup_assembler() and stores it on
        # `params`. The field is kept in the loaded NamedTuple as
        # `nothing` purely to preserve the original return shape for any
        # downstream code that uses propertynames() or similar reflection.
        loaded_matrix = d["matrix"]
        if loaded_matrix isa NamedTuple && !haskey(loaded_matrix, :g_dss_cache)
            loaded_matrix = merge(loaded_matrix, (; g_dss_cache = nothing))
        end
        return (d["metrics"], loaded_matrix)
    catch e
        rank == 0 && @warn "Ignoring unreadable SEM cache $path" exception=(e, catch_backtrace())
        # Auto-delete the unreadable file so the next save writes a
        # clean replacement. Common after struct-shape changes that JLD2
        # cannot reconstruct. Per-rank file in parallel runs, so no
        # race - each rank owns and deletes its own slice.
        try
            isfile(path) && rm(path; force=true)
        catch _
            # best-effort
        end
        return (nothing, nothing)
    end
end

# Strip non-serializable / MPI-bound entries from `matrix` before saving.
# Currently: g_dss_cache (holds MPI.Comm + MPI.MultiRequest, unique per
# MPI session). Add new MPI-bound fields here as they appear.
function _matrix_for_cache(matrix)
    matrix isa NamedTuple || return matrix
    keep = filter(k -> k !== :g_dss_cache, propertynames(matrix))
    return NamedTuple{keep}(map(k -> getproperty(matrix, k), keep))
end

function _save_sem_cache(path::String, metrics, matrix; inputs=nothing, nparts::Int=1)
    isempty(path) && return
    rank = MPI.Comm_rank(get_mpi_comm())
    try
        _ensure_cache_dir(path)
        fp = inputs === nothing ? Dict{String,Any}() : _cache_fingerprint(inputs, nparts)
        slim_matrix = _matrix_for_cache(matrix)
        JLD2.jldsave(path; metrics, matrix = slim_matrix, fingerprint = fp)
        rank == 0 && println(" # Saved SEM preprocess cache: $path")
    catch e
        rank == 0 && @warn "Failed to save SEM cache $path" exception=(e, catch_backtrace())
    end
end
# ──────────────────────────────────────────────────────────────────────────────

function sem_setup(inputs::Dict, nparts, distribute, args...)
    
    comm = distribute.comm
    rank = MPI.Comm_rank(comm)
    adapt_flags, partitioned_model_coarse, omesh = _handle_optional_args4amr(args...)
    
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
        # AMR restart: replace coarse/preadapt mesh with the forest saved at checkpoint
        if get(inputs, :lrestart_amr, false)
            # Resolve checkpoint iteration index: explicit or auto-detected from simulation.pvd
            iout = if haskey(inputs, :restart_vtk_iout) && inputs[:restart_vtk_iout] > 0
                inputs[:restart_vtk_iout]
            else
                pvd_path = joinpath(inputs[:output_dir], "simulation.pvd")
                _, i = read_pvd_last_entry(pvd_path)
                inputs[:restart_vtk_iout] = i   # cache for initialize()
                i
            end
            forest_file = joinpath(inputs[:output_dir], "iter_$(iout)", "iter_$(iout).p4est")
            if rank == 0 
                println(" # AMR restart: loading p4est forest from $forest_file")
            end

            # Build OctreeDistributedDiscreteModel from the loaded forest
            @outputrootonly loaded_model = load_p4est_checkpoint_model(partitioned_model, forest_file)

            # Build Jexpresso mesh struct from loaded model using the AMR-adapt path
            # with no-op flags (all nothing_flag = no further adaptation).
            # Build interp/project from the coarse mesh polynomial order first.
            ξω_base = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.nop, inputs[:backend])
            interp_base, project_base = build_projection_1d(ξω_base.ξ)
            nothing_flags = KernelAbstractions.zeros(CPU(), TInt, 0)   # empty → all nothing_flag
            dummy_uaux    = KernelAbstractions.zeros(CPU(), TFloat, (mesh.npoin, 1))
            mesh, partitioned_model, _ = mod_mesh_mesh_driver(
                inputs, nparts, distribute,
                nothing_flags, loaded_model, mesh,
                interp_base, project_base, dummy_uaux)
            # Restore ad_lvl from the loaded p4est forest.  The glue-propagation
            # path inside mod_mesh_mesh_driver starts from the coarse mesh (ad_lvl=0)
            # and cannot recover the true refinement levels of the checkpoint.
            # mesh.ad_lvl = read_ad_lvl_from_p4est(partitioned_model.ptr_pXest)
        end
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
    # mesh.xmin = minimum(mesh.coords[:,1])
    # mesh.xmax = maximum(mesh.coords[:,1])
    mesh.xmax = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
    mesh.xmin = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)
    if (inputs[:yscale] != 1.0 && inputs[:ydisp] != 0.0)
        mesh.y[:] .= (mesh.y[:] .+ inputs[:ydisp]) .*inputs[:yscale] * 0.5
    elseif(inputs[:yscale] != 1.0)
        mesh.y[:] .= (mesh.y[:]) .*inputs[:yscale]*0.5
    elseif(inputs[:ydisp] != 0.0)
        mesh.y[:] .= (mesh.y[:] .+ inputs[:ydisp])
    end
    if mesh.nsd == 2
        mesh.ymax = MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm)
        mesh.ymin = MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm)
    end

    # ── Keep mesh.coords in sync with the (possibly scaled/displaced) nodes ──────
    # mesh.coords was filled from mesh.x/mesh.y when the mesh was read, BEFORE the
    # affine xscale/yscale/xdisp/ydisp transform above. The Dirichlet-BC routines
    # (BCs.jl) read node positions from mesh.coords, while the source / initial /
    # exact fields read mesh.x/mesh.y, and the boundary detection compares
    # coords[:,1] against the (scaled) mesh.xmin/xmax. Without this re-sync a
    # scaled mesh evaluates boundary data at the UNSCALED coordinates ⇒
    # inconsistent BCs and a blown-up manufactured-solution error. No-op when no
    # scaling/displacement was applied (coords already equal x/y from the read).
    if (inputs[:xscale] != 1.0 || inputs[:xdisp] != 0.0 ||
        inputs[:yscale] != 1.0 || inputs[:ydisp] != 0.0) && !isempty(mesh.coords)
        np = min(size(mesh.coords, 1), length(mesh.x))
        @views mesh.coords[1:np, 1] .= mesh.x[1:np]
        if size(mesh.coords, 2) >= 2
            @views mesh.coords[1:np, 2] .= mesh.y[1:np]
        end
    end

    #--------------------------------------------------------
    # Build interpolation and quadrature points/weights
    #--------------------------------------------------------
    ξω  = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.nop, inputs[:backend])   
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
    end
    SD = mesh.SD

    # ── SEM preprocess cache load ────────────────────────────────────────────
    # Reuse metrics + matrix from a previous run with the same fingerprint.
    # Path is per-rank in parallel runs, so no file race. JEXPRESSO_PREFETCHED_SEM_CACHE
    # (populated by je_prefetch_caches! before the with_mpi block) is honoured
    # for the Alya coupling path. When the cache is unusable (stale, missing,
    # fingerprint mismatch, or :luse_mesh_cache=false) we fall through and
    # rebuild as normal, then save at the end of each metric-build branch.
    #
    # Cross-rank consistency: same pattern as the mesh cache. If any rank
    # fails to load (missing/stale/mismatched file), every rank rebuilds.
    # Prevents silent inconsistency where one rank uses cached metrics tied
    # to one partition and another rank builds fresh metrics for a different
    # one.
    preprocess_cache = _preprocess_cache_path(inputs, Nξ, Qξ, nparts)
    cached_metrics, cached_matrix = _try_load_sem_cache(preprocess_cache;
                                                        gmsh_path=get(inputs, :gmsh_filename, ""),
                                                        inputs=inputs, nparts=nparts)
    local_loaded = !isnothing(cached_metrics)
    loaded_from_cache = nparts > 1 ?
        (MPI.Allreduce(local_loaded ? 1 : 0, MPI.MIN, comm) == 1) :
        local_loaded
    if loaded_from_cache
        rank == 0 && println(" # Loaded SEM preprocess cache — skipping metric terms and matrix build: $preprocess_cache")
    elseif local_loaded
        # We had a usable local cache but some peer didn't — drop ours so
        # downstream code doesn't accidentally use it.
        cached_metrics, cached_matrix = nothing, nothing
        rank == 0 && println(" # SEM cache: some ranks failed to load — discarding all and rebuilding")
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
            
            if loaded_from_cache
                metrics = cached_metrics
                matrix  = cached_matrix
            else
                if (rank == 0) println(" # Build metrics ......") end
                metrics1 = allocate_metrics(SD, mesh.nelem, mesh.nedges_bdy, Qξ, TFloat, inputs[:backend])
                @time build_metric_terms!(metrics1, mesh, basis1, Nξ, Qξ, ξ, ω1, TFloat, COVAR(), SD; backend = inputs[:backend])

                metrics2 = allocate_metrics_laguerre(SD, mesh.nelem_semi_inf, mesh.nedges_bdy, Qξ, mesh.ngr, TFloat, inputs[:backend])
                build_metric_terms!(metrics2, mesh, basis1, basis2, Nξ, Qξ, mesh.ngr, mesh.ngr, ξ, ω1, ω2, TFloat, COVAR(), SD; backend = inputs[:backend])

                metrics = (metrics1, metrics2)
                if (rank == 0) println(" # Build metrics ...... DONE") end

                matrix = matrix_wrapper_laguerre(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat;
                                                 ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend], interp)
                _save_sem_cache(preprocess_cache, metrics, matrix; inputs=inputs, nparts=nparts)
            end

        else
            if (rank == 0) println(" # Build interpolation bases ......") end
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])
            if (rank == 0) println(" # Build interpolation bases ...... END") end
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
            #--------------------------------------------------------
            # Build metric terms
            #--------------------------------------------------------
            #if (mesh.nsd > 2)
            #    if (inputs[:lwarp]) warp_mesh_3D!(mesh,inputs) end
            #else
            #    if (inputs[:lwarp]) warp_mesh!(mesh,inputs) end
            #end
            if loaded_from_cache
                metrics = cached_metrics
                matrix  = cached_matrix
                # phys_grid still needs to be built even from cache — it is
                # not part of the cached payload (it depends on runtime
                # decomposition state).
                if (inputs[:lphysics_grid])
                    phys_grid = init_phys_grid(mesh, inputs,inputs[:nlay_pg],inputs[:nx_pg],inputs[:ny_pg],mesh.xmin,mesh.xmax,mesh.ymin,mesh.ymax,mesh.zmin,mesh.zmax,inputs[:backend])
                end
            else
                if (rank == 0) println(" # Build metrics ......") end
                metrics = allocate_metrics(SD, mesh.nelem, mesh.nedges_bdy, Qξ, TFloat, inputs[:backend])
                @mpi_time build_metric_terms!(metrics, mesh, basis, Nξ, Qξ, ξ, ω, TFloat, COVAR(), SD; backend = inputs[:backend])
                if (rank == 0) println(" # Build metrics ...... END") end

                if (inputs[:lphysics_grid])
                    phys_grid = init_phys_grid(mesh, inputs,inputs[:nlay_pg],inputs[:nx_pg],inputs[:ny_pg],mesh.xmin,mesh.xmax,mesh.ymin,mesh.ymax,mesh.zmin,mesh.zmax,inputs[:backend])
                end
                if (rank == 0) println(" # Build periodicity infrastructure ......") end

                if (rank == 0) println(" # Matrix wrapper ......") end
                matrix = matrix_wrapper(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace],
                            ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend], interp)
                if (rank == 0)  println(" # Matrix wrapper ...... END") end
                _save_sem_cache(preprocess_cache, metrics, matrix; inputs=inputs, nparts=nparts)
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
            #--------------------------------------------------------
            # Build metric terms
            #--------------------------------------------------------
            if loaded_from_cache
                metrics = cached_metrics
                matrix  = cached_matrix
            else
                if (rank == 0) println(" # Build metrics ......") end
                metrics1 = allocate_metrics(SD, mesh.nelem, mesh.nedges_bdy, Qξ, TFloat, inputs[:backend])
                build_metric_terms!(metrics1, mesh, basis[1], Nξ, Qξ, ξ, ω, TFloat, COVAR(), SD; backend = inputs[:backend])
                metrics2 = allocate_metrics(SD, mesh.nelem_semi_inf, mesh.nedges_bdy, mesh.ngr, TFloat, inputs[:backend])
                build_metric_terms_1D_Laguerre!(metrics2, mesh, basis[2], mesh.ngr, mesh.ngr, ξ2, ω2, inputs, TFloat, COVAR(), SD;backend = inputs[:backend])

                metrics = (metrics1, metrics2)
                if (rank == 0) println(" # Build metrics ...... DONE") end
                matrix = matrix_wrapper_laguerre(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend], interp)
                _save_sem_cache(preprocess_cache, metrics, matrix; inputs=inputs, nparts=nparts)
            end
        else
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])

            ω1 = ω
            ω = ω1
            #--------------------------------------------------------
            # Build metric terms
            #--------------------------------------------------------
            if loaded_from_cache
                metrics = cached_metrics
                matrix  = cached_matrix
                # 1D periodicity restructuring mutates `mesh.connijk`
                # etc. and is not stored in the cache — re-run after a
                # cache load to keep mesh state consistent with the
                # cached metrics/matrix (which were built post-restructure).
                if (inputs[:lperiodic_1d])
                    @time restructure4periodicity_1D!(mesh, mesh.coords,
                                                      mesh.xmax, mesh.xmin,mesh.ymax,mesh.ymin,mesh.zmax,mesh.zmin,
                                                      mesh.poin_in_bdy_face,mesh.poin_in_bdy_edge,
                                                      mesh.ngl,mesh.ngr,mesh.nelem,mesh.npoin,mesh.nsd,mesh.bdy_edge_type,
                                                      mesh.bdy_face_type,mesh.bdy_face_in_elem,mesh.bdy_edge_in_elem,
                                                      mesh.connijk,mesh.connijk_lag,mesh.npoin_linear,mesh.nelem_semi_inf,
                                                      inputs,inputs[:backend])
                end
            else
                metrics = allocate_metrics(SD, mesh.nelem, mesh.nedges_bdy, Qξ, TFloat, inputs[:backend])
                @time build_metric_terms!(metrics, mesh, basis, Nξ, Qξ, ξ, ω, TFloat, COVAR(), SD; backend = inputs[:backend])

                if (inputs[:lperiodic_1d])
                    @time restructure4periodicity_1D!(mesh, mesh.coords,
                                                      mesh.xmax, mesh.xmin,mesh.ymax,mesh.ymin,mesh.zmax,mesh.zmin,
                                                      mesh.poin_in_bdy_face,mesh.poin_in_bdy_edge,
                                                      mesh.ngl,mesh.ngr,mesh.nelem,mesh.npoin,mesh.nsd,mesh.bdy_edge_type,
                                                      mesh.bdy_face_type,mesh.bdy_face_in_elem,mesh.bdy_edge_in_elem,
                                                      mesh.connijk,mesh.connijk_lag,mesh.npoin_linear,mesh.nelem_semi_inf,
                                                      inputs,inputs[:backend])
                end
                matrix = matrix_wrapper(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat;
                                        ldss_laplace=inputs[:ldss_laplace],
                                        ldss_differentiation=inputs[:ldss_differentiation],
                                        backend = inputs[:backend], interp)
                _save_sem_cache(preprocess_cache, metrics, matrix; inputs=inputs, nparts=nparts)
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
        return (; QT, CL, AD, SOL_VARS_TYPE, volume_flux, mesh, metrics, basis, ξ, ω, matrix, fx, fy, fy_lag, fz, phys_grid,
                connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original, interp, project, nparts, distribute), partitioned_model
    else
        return (; QT, CL, AD, SOL_VARS_TYPE, volume_flux, mesh, metrics, basis, ξ, ω, matrix, fx, fy, fy_lag, fz, phys_grid,
                connijk_original, poin_in_bdy_face_original, x_original, y_original, z_original, interp, project, nparts, distribute), partitioned_model, uaux_new
    end
    
end
