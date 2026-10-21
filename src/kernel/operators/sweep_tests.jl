function test_stress_combined_ncf_unified_exponential(inputs, nparts, distribute, sem, partitioned_model;
                                                        seed::Int=42, frac::Float64=0.01, max_e_ext::Int=4,
                                                        κ_val::Float64=0.0001, I0::Float64=4.2)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()

    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(sem.mesh.gnelem; seed=seed, frac=frac, max_e_ext=max_e_ext)
    rank == 0 && println("Generated $(length(spatial_target_gels)) spatial targets, " *
                          "$(length(angular_target_gels)) angular targets")

    angular_targets = Tuple{Int,Int}[(angular_target_gels[i], angular_target_e_exts[i])
                                      for i = 1:length(angular_target_gels)]

    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=collect(Int, spatial_target_gels),
                                    angular_targets=angular_targets)

    err = run_unified_exponential_check(setup, comm; κ_val=κ_val, I0=I0)
    rank == 0 && println("STRESS TEST unified sweep exponential — max abs error: $err")
    return err
end

function test_stress_combined_ncf_unified_manufactured(inputs, nparts, distribute, sem, partitioned_model;
                                                        seed::Int=42, frac::Float64=0.01, max_e_ext::Int=4,
                                                        I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                        κ_val=0.001, σ_val=0.0005,
                                                        max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()

    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(sem.mesh.gnelem; seed=seed, frac=frac, max_e_ext=max_e_ext)
    rank == 0 && println("Generated $(length(spatial_target_gels)) spatial targets, " *
                          "$(length(angular_target_gels)) angular targets")

    angular_targets = Tuple{Int,Int}[(angular_target_gels[i], angular_target_e_exts[i])
                                      for i = 1:length(angular_target_gels)]

    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=collect(Int, spatial_target_gels),
                                    angular_targets=angular_targets)

    err = run_unified_manufactured_check(setup, comm; I0=I0, A=A, b=b, c=c, d=d,
                                          κ_val=κ_val, σ_val=σ_val, max_iter=max_iter, tol=tol)
    rank == 0 && println("STRESS TEST unified sweep manufactured — max abs error: $err")
    return err
end

function run_unified_exponential_check(setup, comm; κ_val::Float64=0.0001, I0::Float64=4.2)
    mesh = setup.mesh; metrics = setup.metrics
    ψ, dψ, ω = setup.sem.basis.ψ, setup.sem.basis.dψ, setup.sem.ω
    ngl = mesh.ngl; n_dirs = setup.n_dirs

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg);
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    bc_fn_all = Dict{Tuple{Int,Int}, Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωxd = sin(θ)*cos(ϕ); Ωyd = sin(θ)*sin(ϕ); Ωzd = cos(θ)
            bc_fn_all[(iel,d)] = (x,y,z,ip,nx,ny,nz) -> I0*exp(-κ_val*(Ωxd*(x-r_bc[1])+Ωyd*(y-r_bc[2])+Ωzd*(z-r_bc[3])))
        end
    end

    x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
        setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence, setup.periodic_face_correspondence,
        setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group, setup.neighbor_is_periodic,
        setup.same_rank_kind, setup.same_rank_partner, setup.same_rank_is_periodic,
        setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
        setup.remote_gid_coords,
        setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh,
        setup.gid_exch, setup.periodic_gid_exch,
        κ, bc_fn_all, comm)


    # before the existing error loop:
    sum_local = 0.0; sum2_local = 0.0; cnt_local = 0
    max_err_local = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωxd = sin(θ)*cos(ϕ); Ωyd = sin(θ)*sin(ϕ); Ωzd = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωxd*(mesh.x[gip]-r_bc[1])+Ωyd*(mesh.y[gip]-r_bc[2])+Ωzd*(mesh.z[gip]-r_bc[3])))
                e = abs(x_dg[(iel,d)][lidx(i,j,k,ngl)] - exact)
		max_err_local = max(max_err_local, e); sum_local += e; sum2_local += e*e; cnt_local += 1
            end
        end
    end
    max_err_global =  MPI.Allreduce(max_err_local, MPI.MAX, comm)
    st = MPI.Allreduce(Float64[sum_local, sum2_local, cnt_local], +, comm)
    fp = solution_fingerprint(mesh, x_dg, comm)            # collective: every rank calls it
    if MPI.Comm_rank(comm) == 0
    	println(@sprintf("exponential field error: mean |e| = %.6e, rms = %.6e, over %d values", st[1]/st[3], sqrt(st[2]/st[3]), Int(st[3])))
    	foreach(println, fp)
    end
    return max_err_global
end

function run_unified_manufactured_check(setup, comm;
                                         I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                         κ_val=0.001, σ_val=0.0005, max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(comm)
    mesh = setup.mesh; metrics = setup.metrics
    ψ, dψ, ω = setup.sem.basis.ψ, setup.sem.basis.dψ, setup.sem.ω
    ngl = mesh.ngl; n_dirs = setup.n_dirs

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 4π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg);

    groups, representative_of = group_elements_by_angular_mesh(mesh)
    coeff_all = Dict{Int, Vector{Float64}}(); C_all = Dict{Int, Matrix{Float64}}()
    Ωm_all = Dict{Int, Vector{NTuple{3,Float64}}}()
    for rep in keys(groups)
        em = mesh.extra_mesh[rep]
        w_ext = zeros(Float64, em.extra_npoin)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
                w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
            end
        end
        W = sum(w_ext)
        C_all[rep] = [σ_val*w_ext[m]/W for kk=1:em.extra_npoin, m=1:em.extra_npoin]
        Ωm_all[rep] = [(sin(em.extra_coords[1,m])*cos(em.extra_coords[2,m]),
                        sin(em.extra_coords[1,m])*sin(em.extra_coords[2,m]),
                        cos(em.extra_coords[1,m])) for m = 1:em.extra_npoin]
        coeff_all[rep] = σ_val .* w_ext ./ W
    end

    q_ext_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    bc_fn_all = Dict{Tuple{Int,Int}, Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        rep = representative_of[iel]
        C_iel = C_all[rep]; Ω_m = Ωm_all[rep]
        for dk = 1:em.extra_npoin
            θk = em.extra_coords[1,dk]; ϕk = em.extra_coords[2,dk]
            Ωxk = sin(θk)*cos(ϕk); Ωyk = sin(θk)*sin(ϕk); Ωzk = cos(θk)
            q_vec = zeros(npoin_cg)
            for k=1:ngl,j=1:ngl,i=1:ngl
                g = mesh.connijk[iel,i,j,k]
                x_,y_,z_ = mesh.x[g], mesh.y[g], mesh.z[g]
                grad = dI_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                Ω_dot_grad = Ωxk*grad[1] + Ωyk*grad[2] + Ωzk*grad[3]
                I_true = I_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                scatter_true = sum(C_iel[dk,m]*I_manuf(x_,y_,z_, Ω_m[m]...) for m = 1:em.extra_npoin)
                q_vec[g] = Ω_dot_grad + κ_val*I_true - scatter_true
            end
            q_ext_all[(iel,dk)] = q_vec
            bc_fn_all[(iel,dk)] = (x,y,z,ip,nx,ny,nz) -> I_manuf(x,y,z,Ωxk,Ωyk,Ωzk)
	end
    end

    I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
        I_cg_all[(iel,dd)] = zeros(ngl^3)
    end

    for it = 1:max_iter
        extra_source_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
        for iel = 1:mesh.nelem
            em = mesh.extra_mesh[iel]
            coeff = coeff_all[representative_of[iel]]
            local_gips = local_gips_lidx(mesh, iel)
	    for dk = 1:em.extra_npoin
                q_k = copy(q_ext_all[(iel,dk)])
                q_k_local = @view q_k[local_gips] #@view q_k[[mesh.connijk[iel,i,j,k] for k=1:ngl,j=1:ngl,i=1:ngl][:]]
                for m = 1:em.extra_npoin
                    coeff[m] == 0.0 && continue
                    q_k_local .+= coeff[m] .* I_cg_all[(iel,m)]
                end
                extra_source_all[(iel,dk)] = q_k
            end
        end

        x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
            setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence, setup.periodic_face_correspondence,
            setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group, setup.neighbor_is_periodic,
            setup.same_rank_kind, setup.same_rank_partner, setup.same_rank_is_periodic,
            setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
            setup.remote_gid_coords,
            setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh,
            setup.gid_exch, setup.periodic_gid_exch,
            κ, bc_fn_all, comm;
            extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all)

        max_delta_local = 0.0
        for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
            max_delta_local = max(max_delta_local, maximum(abs.(x_dg[(iel,dd)] .- I_cg_all[(iel,dd)])))
            I_cg_all[(iel,dd)] = copy(x_dg[(iel,dd)])
        end
        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        rank == 0 && println("  iter $it: max_delta=$max_delta_global")
        max_delta_global < tol && (rank==0 && println("Converged after $it iterations"); break)
    end

    max_err_local = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for dd = 1:n_dirs[iel]
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                g = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[g],mesh.y[g],mesh.z[g], Ωx,Ωy,Ωz)
                max_err_local = max(max_err_local, abs(I_cg_all[(iel,dd)][lidx(i,j,k,ngl)] - I_true))
            end
        end
    end
    return MPI.Allreduce(max_err_local, MPI.MAX, comm)
end

"""
    test_unified_component(inputs, nparts, distribute, sem, partitioned_model;
                            mode, problem, spatial_gel=315, angular_gel=242, angular_e_ext=4)

mode:    :none | :spatial | :angular | :combined
problem: :exponential | :manufactured

Targets are GLOBAL element ids, so a given case selects the same physical
elements at any rank count, and the error should match to all printed digits.
"""
function test_unified_component(inputs, nparts, distribute, sem, partitioned_model;
                                 mode::Symbol, problem::Symbol,
                                 spatial_gel::Int=315, angular_gel::Int=242, angular_e_ext::Int=4)
    mode in (:none, :spatial, :angular, :combined) || error("bad mode $mode")
    problem in (:exponential, :manufactured) || error("bad problem $problem")
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()

    spatial_targets = (mode == :spatial || mode == :combined) ? Int[spatial_gel] : Int[]
    angular_targets = (mode == :angular || mode == :combined) ? Tuple{Int,Int}[(angular_gel, angular_e_ext)] :
                                                                Tuple{Int,Int}[]

    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=spatial_targets, angular_targets=angular_targets)

    err = problem == :exponential ? run_unified_exponential_check(setup, comm) :
                                    run_unified_manufactured_check(setup, comm)
    rank == 0 && println("COMPONENT TEST [mode=$mode, problem=$problem, nranks=$(MPI.Comm_size(comm))] " *
                          "— max abs error: $err")
    return err
end

function test_bidirectional_disambiguation(mesh, metrics, neighbor_elem, neighbor_rank, neighbor_gel,
                                            exchange, value_cache, Ωx, Ωy, Ωz, comm)
    rank = MPI.Comm_rank(comm)
    ngl  = mesh.ngl
    nelem = mesh.nelem

    is_inflow_all = Vector{Vector{Bool}}(undef, nelem)
    is_active_all = Vector{Vector{Bool}}(undef, nelem)
    for iel = 1:nelem
        ii = Vector{Bool}(undef,6); ia = Vector{Bool}(undef,6)
        classify_element_faces!(ii, ia, iel, Ωx, Ωy, Ωz, ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)
        is_inflow_all[iel] = ii; is_active_all[iel] = ia
    end

    # ── Send a FABRICATED, self-describing payload: each rank sends its OWN
    # element's global id as the "value" (not a real DG solve) — so we can
    # check, on receipt, whether what arrives at each of MY faces matches
    # what SHOULD have arrived for that SPECIFIC face, using neighbor_gel
    # as ground truth. This isolates the exchange/matching mechanism
    # completely from the actual sweep/solve logic. ─────────────────────────
    send_vals     = Vector{Float64}(undef, length(exchange.send_keys))
    send_gel_vals = Vector{Int}(undef, length(exchange.send_keys))
    for (idx,(iel,fid,pos)) in enumerate(exchange.send_keys)
        send_vals[idx]     = Float64(mesh.el2gel[iel])   # fabricated: my own element's gel, AS the value
        send_gel_vals[idx] = Int(mesh.el2gel[iel])
    end
    send_and_receive!(value_cache, send_vals)
    send_and_receive!(exchange.cache_gel, send_gel_vals)

    remote_value = Dict{Tuple{Int,Int}, Vector{Float64}}()
    remote_gel   = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ((iel,fid),_) in neighbor_rank
        n = length(exchange.face_local_order[(TInt(iel),TInt(fid))])
        remote_value[(iel,fid)] = fill(NaN, n)
        remote_gel[(iel,fid)]   = fill(-1, n)
    end
    for m = 1:length(value_cache.combined_recv_data)
        val = value_cache.combined_recv_data[m]
        recv_gel = exchange.cache_gel.combined_recv_data[m]
        for pos in exchange.recv_perm[m]
            iel, fid, p = exchange.send_keys[pos]
            remote_value[(iel,fid)][p] = val
            remote_gel[(iel,fid)][p]   = recv_gel
        end
    end
    # ── Verify: for EVERY (iel,fid) in neighbor_rank — inflow AND outflow,
    # since this test checks the EXCHANGE, not just what the sweep happens
    # to consume — the received value must equal the received gel (since
    # we fabricated value==gel), AND the received gel must equal
    # neighbor_gel[(iel,fid)] (the INDEPENDENTLY established correct
    # remote element, from resolve_rank_cut_neighbors). ─────────────────────
    n_checked = 0
    n_value_gel_mismatch = 0     # received value != received gel (internal inconsistency)
    n_wrong_neighbor      = 0     # received gel != expected neighbor_gel (WRONG face matched)
    n_missing             = 0     # never received anything at all

    for ((iel,fid), expected_gel) in neighbor_rank
        rv = remote_value[(iel,fid)]
        rg = remote_gel[(iel,fid)]
        for p = 1:length(rv)
            n_checked += 1
            if isnan(rv[p]) || rg[p] == -1
                n_missing += 1
                continue
            end
            Int(round(rv[p])) != rg[p] && (n_value_gel_mismatch += 1)
            rg[p] != Int(neighbor_gel[(TInt(iel),TInt(fid))]) && (n_wrong_neighbor += 1)
        end
    end

    println("Rank $rank, Ω=($Ωx,$Ωy,$Ωz): bidirectional disambiguation check — " *
             "checked=$n_checked, missing=$n_missing, " *
             "value/gel internal mismatch=$n_value_gel_mismatch, " *
             "wrong-neighbor matches=$n_wrong_neighbor " *
             (n_missing==0 && n_value_gel_mismatch==0 && n_wrong_neighbor==0 ? "✓ PASS" : "✗ FAIL"))

    return (n_checked=n_checked, n_missing=n_missing,
            n_value_gel_mismatch=n_value_gel_mismatch, n_wrong_neighbor=n_wrong_neighbor)
end

function run_unified_manufactured_check_phase(setup, comm, phase_fn::Function;
                                               I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                               κ_val=0.001, σ_val=0.0005, g::Float64=0.0,
                                               max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(comm)
    mesh = setup.mesh; metrics = setup.metrics
    ψ, dψ, ω = setup.sem.basis.ψ, setup.sem.basis.dψ, setup.sem.ω
    ngl = mesh.ngl; n_dirs = setup.n_dirs

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg);

    groups, representative_of = group_elements_by_angular_mesh(mesh)
    C_all = Dict{Int, Matrix{Float64}}()
    Ωm_all = Dict{Int, Vector{NTuple{3,Float64}}}()
    for rep in keys(groups)
        em = mesh.extra_mesh[rep]
        w_ext = zeros(Float64, em.extra_npoin)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
                w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
            end
        end
        thetas = em.extra_coords[1,:]; phis = em.extra_coords[2,:]
        npd = em.extra_npoin
        # C[dk,m] = σ · Φ_HG(Ω_dk,Ω_m;g) · w_ext[m]  — NON-ISOTROPIC
        C = zeros(npd, npd)
        for dk = 1:npd, m = 1:npd
            C[dk,m] = σ_val * phase_fn(thetas[dk], thetas[m], phis[dk], phis[m], g) * w_ext[m]
        end
        C_all[rep] = C
        Ωm_all[rep] = [(sin(thetas[m])*cos(phis[m]), sin(thetas[m])*sin(phis[m]), cos(thetas[m])) for m = 1:npd]
    end

    q_ext_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    bc_fn_all = Dict{Tuple{Int,Int}, Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        rep = representative_of[iel]
        C_iel = C_all[rep]; Ω_m = Ωm_all[rep]
        for dk = 1:em.extra_npoin
            θk = em.extra_coords[1,dk]; ϕk = em.extra_coords[2,dk]
            Ωxk = sin(θk)*cos(ϕk); Ωyk = sin(θk)*sin(ϕk); Ωzk = cos(θk)
            q_vec = zeros(npoin_cg)
            for k=1:ngl,j=1:ngl,i=1:ngl
                g_ = mesh.connijk[iel,i,j,k]
                x_,y_,z_ = mesh.x[g_], mesh.y[g_], mesh.z[g_]
                grad = dI_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                Ω_dot_grad = Ωxk*grad[1] + Ωyk*grad[2] + Ωzk*grad[3]
                I_true = I_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                scatter_true = sum(C_iel[dk,m]*I_manuf(x_,y_,z_, Ω_m[m]...) for m = 1:em.extra_npoin)
                q_vec[g_] = Ω_dot_grad + κ_val*I_true - scatter_true
            end
            q_ext_all[(iel,dk)] = q_vec
            bc_fn_all[(iel,dk)] = (x,y,z,ip,nx,ny,nz) -> I_manuf(x,y,z,Ωxk,Ωyk,Ωzk)
	end
    end

    I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
        I_cg_all[(iel,dd)] = zeros(ngl^3)
    end

    for it = 1:max_iter
        extra_source_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
        for iel = 1:mesh.nelem
            em = mesh.extra_mesh[iel]
            C_iel = C_all[representative_of[iel]]
            local_gips = local_gips_lidx(mesh, iel)
	    for dk = 1:em.extra_npoin
                q_k = copy(q_ext_all[(iel,dk)])
                q_k_local = @view q_k[local_gips] #@view q_k[[mesh.connijk[iel,i,j,k] for k=1:ngl,j=1:ngl,i=1:ngl][:]]
                for m = 1:em.extra_npoin
                    C_iel[dk,m] == 0.0 && continue
                    q_k_local .+= C_iel[dk,m] .* I_cg_all[(iel,m)]
                end
                extra_source_all[(iel,dk)] = q_k
            end
        end

        x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
            setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence, setup.periodic_face_correspondence,
            setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group, setup.neighbor_is_periodic,
            setup.same_rank_kind, setup.same_rank_partner, setup.same_rank_is_periodic,
            setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
            setup.remote_gid_coords,
            setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh,
            setup.gid_exch, setup.periodic_gid_exch,
            κ, bc_fn_all, comm;
            extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all)

        max_delta_local = 0.0
        for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
            max_delta_local = max(max_delta_local, maximum(abs.(x_dg[(iel,dd)] .- I_cg_all[(iel,dd)])))
            I_cg_all[(iel,dd)] = copy(x_dg[(iel,dd)])
        end
        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        rank == 0 && println("  iter $it: max_delta=$max_delta_global")
        max_delta_global < tol && (rank==0 && println("Converged after $it iterations"); break)
    end

    max_err_local = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for dd = 1:n_dirs[iel]
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gg = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[gg],mesh.y[gg],mesh.z[gg],Ωx,Ωy,Ωz)
                max_err_local = max(max_err_local, abs(I_cg_all[(iel,dd)][lidx(i,j,k,ngl)] - I_true))
            end
        end
    end
    return MPI.Allreduce(max_err_local, MPI.MAX, comm)
end

function test_stress_combined_ncf_unified_manufactured_hg(inputs, nparts, distribute, sem, partitioned_model;
                                                           seed::Int=42, frac::Float64=0.01, max_e_ext::Int=4,
                                                           g::Float64=0.6, κ_val=0.001, σ_val=0.0005,
                                                           max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(sem.mesh.gnelem; seed=seed, frac=frac, max_e_ext=max_e_ext)
    rank == 0 && println("Generated $(length(spatial_target_gels)) spatial, $(length(angular_target_gels)) angular targets")
    angular_targets = Tuple{Int,Int}[(angular_target_gels[i], angular_target_e_exts[i]) for i = 1:length(angular_target_gels)]

    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=collect(Int, spatial_target_gels), angular_targets=angular_targets)

    err = run_unified_manufactured_check_phase(setup, comm, phase_function_hg;
                                                κ_val=κ_val, σ_val=σ_val, g=g, max_iter=max_iter, tol=tol)
    rank == 0 && println("STRESS TEST unified sweep HG scattering (g=$g) — max abs error: $err")
    return err
end

function test_stress_combined_ncf_unified_manufactured_mie(inputs, nparts, distribute, sem, partitioned_model;
                                                            seed::Int=42, frac::Float64=0.01, max_e_ext::Int=4,
                                                            ripple_N::Float64=8.0, κ_val=0.001, σ_val=0.0005,
                                                            max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(sem.mesh.gnelem; seed=seed, frac=frac, max_e_ext=max_e_ext)
    rank == 0 && println("Generated $(length(spatial_target_gels)) spatial, $(length(angular_target_gels)) angular targets")
    angular_targets = Tuple{Int,Int}[(angular_target_gels[i], angular_target_e_exts[i]) for i = 1:length(angular_target_gels)]

    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=collect(Int, spatial_target_gels), angular_targets=angular_targets)

    err = run_unified_manufactured_check_phase(setup, comm, phase_function_mie_like;
                                                κ_val=κ_val, σ_val=σ_val, g=ripple_N, max_iter=max_iter, tol=tol)
    rank == 0 && println("STRESS TEST unified sweep Mie-like scattering (ripple_N=$ripple_N) — max abs error: $err")
    return err
end

"""
    run_unified_manufactured_check_variable(setup, comm, phase_fn;
                                             κ_fn, σ_fn, g_fn, I0, A, b, c, d, ...)

Manufactured-solution check with SPATIALLY VARYING σ(x,y,z) and g(x,y,z)
(element-averaged for the phase matrix, exact per-node for the σ
multiplication) — the discrete RHS uses the SAME per-element C and
per-node σ as the solve, so exactness is preserved despite the
generalization.
"""
function run_unified_manufactured_check_variable(setup, comm, phase_fn::Function;
                                                   κ_fn::Function, σ_fn::Function, g_fn::Function,
                                                   I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                   max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(comm)
    mesh = setup.mesh; metrics = setup.metrics
    ψ, dψ, ω = setup.sem.basis.ψ, setup.sem.basis.dψ, setup.sem.ω
    ngl = mesh.ngl; n_dirs = setup.n_dirs



    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    npoin_cg = length(mesh.x)
    κ = zeros(npoin_cg); σ = zeros(npoin_cg)
    for ip = 1:npoin_cg
        κ[ip] = κ_fn(mesh.x[ip], mesh.y[ip], mesh.z[ip])
        σ[ip] = σ_fn(mesh.x[ip], mesh.y[ip], mesh.z[ip])
    end

    # ── Per-element C (pure phase-function weights, NO σ baked in),
    # using the ELEMENT-AVERAGED g. ──────────────────────────────────────────
    C_iel_all = Dict{Int, Matrix{Float64}}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        npd = em.extra_npoin
        w_ext = zeros(Float64, npd)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
                w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
            end
        end
        thetas = em.extra_coords[1,:]; phis = em.extra_coords[2,:]
        gsum = 0.0
        for k=1:ngl,j=1:ngl,i=1:ngl
            gip = mesh.connijk[iel,i,j,k]
            gsum += g_fn(mesh.x[gip], mesh.y[gip], mesh.z[gip])
        end
        g_iel = gsum / ngl^3
        C = zeros(npd, npd)
        for dk = 1:npd, m = 1:npd
            C[dk,m] = phase_fn(thetas[dk], thetas[m], phis[dk], phis[m], g_iel) * w_ext[m]
        end
        C_iel_all[iel] = C
    end

    # ── RHS: uses the SAME element-averaged g and per-node σ as the solve,
    # so the discrete scattering term is EXACTLY what the coupled loop
    # will compute — required for the manufactured solution to be exact. ────
    q_ext_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    bc_fn_all = Dict{Tuple{Int,Int}, Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        C_iel = C_iel_all[iel]
        thetas = em.extra_coords[1,:]; phis = em.extra_coords[2,:]
        Ωm = [(sin(thetas[m])*cos(phis[m]), sin(thetas[m])*sin(phis[m]), cos(thetas[m])) for m = 1:em.extra_npoin]
        for dk = 1:em.extra_npoin
            θk = thetas[dk]; ϕk = phis[dk]
            Ωxk = sin(θk)*cos(ϕk); Ωyk = sin(θk)*sin(ϕk); Ωzk = cos(θk)
            q_vec = zeros(npoin_cg)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                x_,y_,z_ = mesh.x[gip], mesh.y[gip], mesh.z[gip]
                grad = dI_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                Ω_dot_grad = Ωxk*grad[1] + Ωyk*grad[2] + Ωzk*grad[3]
                I_true = I_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                scatter_true = 0.0
                for m = 1:em.extra_npoin
                    C_iel[dk,m] == 0.0 && continue
                    scatter_true += C_iel[dk,m] * σ[gip] * I_manuf(x_,y_,z_, Ωm[m]...)
                end
                q_vec[gip] = Ω_dot_grad + κ[gip]*I_true - scatter_true
            end
            q_ext_all[(iel,dk)] = q_vec
            bc_fn_all[(iel,dk)] = (x,y,z,ip,nx,ny,nz) -> I_manuf(x,y,z,Ωxk,Ωyk,Ωzk)
        end
    end

    I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
        I_cg_all[(iel,dd)] = zeros(ngl^3)
    end

    for it = 1:max_iter
        extra_source_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
        for iel = 1:mesh.nelem
            em = mesh.extra_mesh[iel]
            C_iel = C_iel_all[iel]
            #local_gips = [mesh.connijk[iel,i,j,k] for k=1:ngl,j=1:ngl,i=1:ngl][:]
            local_gips = local_gips_lidx(mesh, iel)
	    for dk = 1:em.extra_npoin
                q_k = copy(q_ext_all[(iel,dk)])
                q_k_local = @view q_k[local_gips]
                for m = 1:em.extra_npoin
                    C_iel[dk,m] == 0.0 && continue
                    I_m = I_cg_all[(iel,m)]
                    for (local_idx, gip) in enumerate(local_gips)
                        q_k_local[local_idx] += C_iel[dk,m] * σ[gip] * I_m[local_idx]
                    end
                end
                extra_source_all[(iel,dk)] = q_k
            end
        end

        x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
            setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence, setup.periodic_face_correspondence,
            setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group, setup.neighbor_is_periodic,
            setup.same_rank_kind, setup.same_rank_partner, setup.same_rank_is_periodic,
            setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
            setup.remote_gid_coords,
            setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh,
            setup.gid_exch, setup.periodic_gid_exch,
            κ, bc_fn_all, comm;
            extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all)

        max_delta_local = 0.0
        for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
            max_delta_local = max(max_delta_local, maximum(abs.(x_dg[(iel,dd)] .- I_cg_all[(iel,dd)])))
            I_cg_all[(iel,dd)] = copy(x_dg[(iel,dd)])
        end
        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        rank == 0 && println("  iter $it: max_delta=$max_delta_global")
        max_delta_global < tol && (rank==0 && println("Converged after $it iterations"); break)
    end

    max_err_local = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for dd = 1:n_dirs[iel]
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[gip],mesh.y[gip],mesh.z[gip],Ωx,Ωy,Ωz)
                max_err_local = max(max_err_local, abs(I_cg_all[(iel,dd)][lidx(i,j,k,ngl)] - I_true))
            end
        end
    end
    return MPI.Allreduce(max_err_local, MPI.MAX, comm)
end

function test_stress_combined_ncf_unified_manufactured_variable(inputs, nparts, distribute, sem, partitioned_model;
                                                                  seed::Int=42, frac::Float64=0.01, max_e_ext::Int=4,
                                                                  κ_val::Float64=0.001, max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(sem.mesh.gnelem; seed=seed, frac=frac, max_e_ext=max_e_ext)
    angular_targets = Tuple{Int,Int}[(angular_target_gels[i], angular_target_e_exts[i]) for i = 1:length(angular_target_gels)]
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=collect(Int, spatial_target_gels), angular_targets=angular_targets)

    mesh = setup.mesh
    Lx = mesh.xmax - mesh.xmin; Ly = mesh.ymax - mesh.ymin
    xmid = mesh.xmin + 0.5*Lx

    κ_fn(x,y,z) = κ_val
    # σ varies smoothly across the domain — larger on one side, smaller on the other.
    σ_fn(x,y,z) = 0.0003 + 0.0007*(0.5 + 0.5*sin(2π*(x-mesh.xmin)/Lx))
    # g varies too, e.g. mimicking clear-sky (g≈0) vs cloudy (g≈0.7) regions.
    g_fn(x,y,z) = x < xmid ? 0.05 : 0.7

    err = run_unified_manufactured_check_variable(setup, comm, phase_function_hg;
                                                   κ_fn=κ_fn, σ_fn=σ_fn, g_fn=g_fn, max_iter=max_iter, tol=tol)
    rank == 0 && println("STRESS TEST unified sweep variable σ,g HG scattering — max abs error: $err")
    return err
end

function test_unified_case_manufactured(inputs, nparts, distribute, sem, partitioned_model;
                                         spatial_target_iel::Int=0, angular_target_iel::Int=0, angular_e_ext::Int=4,
                                         I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                         κ_val=0.001, σ_val=0.0005, g_val=0.0, max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    spatial_targets = spatial_target_iel == 0 ? Int[] : Int[spatial_target_iel]
    angular_targets = angular_target_iel == 0 ? Tuple{Int,Int}[] : Tuple{Int,Int}[(angular_target_iel, angular_e_ext)]

    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=spatial_targets, angular_targets=angular_targets)
    mesh = setup.mesh; metrics = setup.metrics


    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    # ── RadiativeCase, built so its source_fn/phase_fn/kappa_fn/sigma_fn
    # produce EXACTLY the same discrete equation the manufactured field
    # solves — same principle as every prior manufactured test. ─────────────
    npoin_cg = length(mesh.x)
    thetas_by_iel = Dict{Int,Vector{Float64}}(); phis_by_iel = Dict{Int,Vector{Float64}}()
    wext_by_iel = Dict{Int,Vector{Float64}}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        thetas_by_iel[iel] = em.extra_coords[1,:]; phis_by_iel[iel] = em.extra_coords[2,:]
        w = zeros(em.extra_npoin)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
                w[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ]*em.ωθ[iθ]*em.ωθ[iϕ]
            end
        end
	@info sum(w)
        wext_by_iel[iel] = w
    end

    source_fn = (x,y,z,ip,θ,ϕ) -> begin
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        grad = dI_manuf(x,y,z,Ωx,Ωy,Ωz)
        Ω_dot_grad = Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3]
        I_true = I_manuf(x,y,z,Ωx,Ωy,Ωz)

        iel_found = 0
        for iel = 1:mesh.nelem, k=1:mesh.ngl,j=1:mesh.ngl,i=1:mesh.ngl
            mesh.connijk[iel,i,j,k] == ip && (iel_found = iel; break)
        end
        thetas = thetas_by_iel[iel_found]; phis = phis_by_iel[iel_found]; w = wext_by_iel[iel_found]
        scatter_true = 0.0
        for m = 1:length(w)
            Ωxm = sin(thetas[m])*cos(phis[m]); Ωym = sin(thetas[m])*sin(phis[m]); Ωzm = cos(thetas[m])
            scatter_true += phase_function_hg(θ, thetas[m], ϕ, phis[m], g_val) * w[m] * σ_val * I_manuf(x,y,z,Ωxm,Ωym,Ωzm)
        end
        Ω_dot_grad + κ_val*I_true - scatter_true
    end

    case = RadiativeCase(
        kappa_fn = (x,y,z,ip) -> κ_val,
        sigma_fn = (x,y,z,ip) -> σ_val,
        phase_fn = phase_function_hg, g_const = g_val,
        bc_fn = (x,y,z,ip,nx,ny,nz,θ,ϕ) -> I_manuf(x,y,z,sin(θ)*cos(ϕ),sin(θ)*sin(ϕ),cos(θ)),
        source_fn = source_fn)

    I_cg_all, κ, σ = run_unified_case_solve(setup, metrics, comm, case; max_iter=max_iter, tol=tol)

    max_err_local = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for dd = 1:em.extra_npoin
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:mesh.ngl,j=1:mesh.ngl,i=1:mesh.ngl
                gip = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[gip],mesh.y[gip],mesh.z[gip],Ωx,Ωy,Ωz)
                max_err_local = max(max_err_local, abs(I_cg_all[(iel,dd)][lidx(i,j,k,mesh.ngl)] - I_true))
            end
        end
    end
    
    debug_find_worst_error(mesh, I_cg_all, I_manuf)

    max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)
    rank == 0 && println("CASE-FRAMEWORK manufactured (HG) — max abs error: $max_err_global")
    return max_err_global
end

"""
Standalone check: builds C_iel EXACTLY the way run_unified_case_solve does
(via build_phase_coupling_matrices_elementwise) and compares it against
the RHS's own inline scattering computation, using the EXACT manufactured
solution at every direction m (not I_cg_all) — isolates whether C_iel's
construction itself matches what the RHS assumes, independent of the
outer iteration's convergence state.
"""
function debug_compare_scatter_construction(mesh, case::RadiativeCase, iel_target::Int, dk_target::Int, ip_target::Int,
                                             I_manuf::Function, g_val::Float64, σ_val::Float64)
    rank = MPI.Comm_rank(get_mpi_comm())

    C_iel_all = build_phase_coupling_matrices_elementwise(mesh, case)
    C_iel = C_iel_all[iel_target]

    em = mesh.extra_mesh[iel_target]
    thetas = em.extra_coords[1,:]; phis = em.extra_coords[2,:]
    w_ext = zeros(em.extra_npoin)
    for e_ext = 1:em.extra_nelem
        nop = em.extra_nop[e_ext]
        for iθ=1:nop+1, iϕ=1:nop+1
            ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
            w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
        end
    end

    x = mesh.x[ip_target]; y = mesh.y[ip_target]; z = mesh.z[ip_target]
    θdk = thetas[dk_target]; ϕdk = phis[dk_target]

    println("Rank $rank: === comparing scatter construction for (iel=$iel_target, dk=$dk_target, ip=$ip_target) ===")
    println("Rank $rank: em.extra_npoin=$(em.extra_npoin), length(thetas)=$(length(thetas)), length(w_ext)=$(length(w_ext))")
    println("Rank $rank: g_val=$g_val, σ_val=$σ_val, (x,y,z)=($x,$y,$z), θdk=$θdk, ϕdk=$ϕdk")

    # ── Path A: via C_iel (what run_unified_case_solve actually uses) ────────
    scatter_via_Ciel = 0.0
    for m = 1:em.extra_npoin
        Ωxm = sin(thetas[m])*cos(phis[m]); Ωym = sin(thetas[m])*sin(phis[m]); Ωzm = cos(thetas[m])
        scatter_via_Ciel += C_iel[dk_target,m] * σ_val * I_manuf(x,y,z,Ωxm,Ωym,Ωzm)
    end
    println("Rank $rank: scatter_true via C_iel[dk,:] * σ * I_manuf(m) = $scatter_via_Ciel")

    # ── Path B: RHS's own direct/inline computation ──────────────────────────
    scatter_via_direct = 0.0
    for m = 1:em.extra_npoin
        Ωxm = sin(thetas[m])*cos(phis[m]); Ωym = sin(thetas[m])*sin(phis[m]); Ωzm = cos(thetas[m])
        scatter_via_direct += phase_function_hg(θdk, thetas[m], ϕdk, phis[m], g_val) * w_ext[m] * σ_val * I_manuf(x,y,z,Ωxm,Ωym,Ωzm)
    end
    println("Rank $rank: scatter_true via DIRECT phase_function_hg*w_ext*σ*I_manuf(m) = $scatter_via_direct")

    println("Rank $rank: DIFFERENCE = $(scatter_via_Ciel - scatter_via_direct)  " *
            (isapprox(scatter_via_Ciel, scatter_via_direct; atol=1e-9) ? "✓ MATCH" : "✗ MISMATCH"))

    # ── Path C: check C_iel[dk,m] itself against phase_fn(...)*w_ext[m] one
    # entry at a time, in case ANY single entry (not just the sum) differs. ──
    n_entry_mismatch = 0
    for m = 1:em.extra_npoin
        expected = phase_function_hg(θdk, thetas[m], ϕdk, phis[m], g_val) * w_ext[m]
        actual = C_iel[dk_target, m]
        if !isapprox(expected, actual; atol=1e-12, rtol=1e-9)
            n_entry_mismatch += 1
            n_entry_mismatch <= 5 && println("Rank $rank: ✗ entry m=$m: C_iel=$actual, expected=$expected")
        end
    end
    println("Rank $rank: entry-by-entry check: $n_entry_mismatch/$(em.extra_npoin) mismatches  " *
            (n_entry_mismatch == 0 ? "✓" : "✗"))

    return scatter_via_Ciel, scatter_via_direct
end

function debug_check_qext_against_source_fn(mesh, case::RadiativeCase, q_ext_all, iel_target::Int, dk_target::Int, ip_target::Int,
                                             I_manuf::Function, dI_manuf::Function, κ_val::Float64, g_val::Float64, σ_val::Float64)
    rank = MPI.Comm_rank(get_mpi_comm())

    em = mesh.extra_mesh[iel_target]
    thetas = em.extra_coords[1,:]; phis = em.extra_coords[2,:]
    w_ext = zeros(em.extra_npoin)
    for e_ext = 1:em.extra_nelem
        nop = em.extra_nop[e_ext]
        for iθ=1:nop+1, iϕ=1:nop+1
            ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
            w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
        end
    end

    x = mesh.x[ip_target]; y = mesh.y[ip_target]; z = mesh.z[ip_target]
    θk = thetas[dk_target]; ϕk = phis[dk_target]
    Ωxk = sin(θk)*cos(ϕk); Ωyk = sin(θk)*sin(ϕk); Ωzk = cos(θk)

    grad = dI_manuf(x,y,z,Ωxk,Ωyk,Ωzk)
    Ω_dot_grad = Ωxk*grad[1] + Ωyk*grad[2] + Ωzk*grad[3]
    I_true = I_manuf(x,y,z,Ωxk,Ωyk,Ωzk)
    scatter_true = 0.0
    for m = 1:em.extra_npoin
        Ωxm = sin(thetas[m])*cos(phis[m]); Ωym = sin(thetas[m])*sin(phis[m]); Ωzm = cos(thetas[m])
        scatter_true += phase_function_hg(θk, thetas[m], ϕk, phis[m], g_val) * w_ext[m] * σ_val * I_manuf(x,y,z,Ωxm,Ωym,Ωzm)
    end
    q_expected = Ω_dot_grad + κ_val*I_true - scatter_true

    q_actual = q_ext_all[(iel_target,dk_target)][ip_target]

    println("Rank $rank: === checking q_ext_all[($iel_target,$dk_target)][$ip_target] vs hand-computed expectation ===")
    println("Rank $rank: q_expected (hand, direct) = $q_expected")
    println("Rank $rank: q_actual   (via case.source_fn, through build_case_fields) = $q_actual")
    println("Rank $rank: DIFFERENCE = $(q_actual - q_expected)  " *
            (isapprox(q_actual, q_expected; atol=1e-9) ? "✓ MATCH" : "✗ MISMATCH"))
end

function debug_find_worst_error(mesh, I_cg_all, I_manuf::Function)
    rank = MPI.Comm_rank(get_mpi_comm())
    worst_err = 0.0
    worst_loc = nothing
    ngl = mesh.ngl
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for dd = 1:em.extra_npoin
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            haskey(I_cg_all, (iel,dd)) || continue
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[gip],mesh.y[gip],mesh.z[gip],Ωx,Ωy,Ωz)
                val = I_cg_all[(iel,dd)][lidx(i,j,k,ngl)]
                err = abs(val - I_true)
                if err > worst_err
                    worst_err = err
                    worst_loc = (iel=iel, d=dd, i=i, j=j, k=k, gip=gip,
                                 x=mesh.x[gip], y=mesh.y[gip], z=mesh.z[gip],
                                 θ=θ, ϕ=ϕ, val=val, I_true=I_true)
                end
            end
        end
    end
    println("Rank $rank: WORST error = $worst_err at $worst_loc")
    return worst_err, worst_loc
end

function debug_check_ordering_match(mesh, iel_target::Int)
    ngl = mesh.ngl
    #local_gips = [mesh.connijk[iel_target,i,j,k] for k=1:ngl,j=1:ngl,i=1:ngl][:]
    local_gips = local_gips_lidx(mesh, iel)
    n_mismatch = 0
    idx = 0
    for k=1:ngl, j=1:ngl, i=1:ngl
        idx += 1
        expected_lidx = i + (j-1)*ngl + (k-1)*ngl^2
        if idx != expected_lidx
            n_mismatch += 1
            n_mismatch <= 5 && println("MISMATCH at (i=$i,j=$j,k=$k): comprehension idx=$idx, lidx formula=$expected_lidx")
        end
    end
    println("Checked $idx entries, $n_mismatch order mismatches between comprehension flattening and lidx formula")
end

function test_longwave_case_isothermal(inputs, nparts, distribute, sem, partitioned_model; T_iso::Float64=288.0, κ_val::Float64=0.01)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model)
    mesh = setup.mesh; metrics = setup.metrics
    bdy = make_boundary_predicates(mesh)
    σ_SB = 5.670374419e-8
    B_T = σ_SB * T_iso^4 / π

    atmos_data = Atmosphere_State{Float64,mesh.npoin}(
        t_current=fill(T_iso, mesh.npoin), t_back=fill(T_iso, mesh.npoin))
    κ_lw = fill(κ_val, mesh.npoin)
    σ_lw = zeros(mesh.npoin)

    case = build_longwave_case(atmos_data, κ_lw, σ_lw, bdy; ε_surface=1.0, T_space=T_iso)
    I_cg_all, κ, σ = run_unified_case_solve(setup, metrics, comm, case; max_iter=200, tol=1e-10)

    ngl = mesh.ngl
    max_err_local = 0.0
    for iel = 1:mesh.nelem, dd = 1:setup.n_dirs[iel], k=1:ngl,j=1:ngl,i=1:ngl
        max_err_local = max(max_err_local, abs(I_cg_all[(iel,dd)][lidx(i,j,k,ngl)] - B_T))
    end
    max_err = MPI.Allreduce(max_err_local, MPI.MAX, comm)
    rank == 0 && println("LW isothermal equilibrium check (B(T)=$B_T) — max abs error: $max_err")
    return max_err
end

function test_shortwave_case_zero_scatter(inputs, nparts, distribute, sem, partitioned_model;
                                            S0::Float64=1361.0, μ0::Float64=0.5, κ_val::Float64=0.001)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model)
    mesh = setup.mesh; metrics = setup.metrics
    bdy = make_boundary_predicates(mesh)

    κ_sw = fill(κ_val, mesh.npoin)
    σ_sw = zeros(mesh.npoin)
    g_eff = zeros(mesh.npoin)
    κ_ext_sw = κ_sw .+ σ_sw
    z_prof, τ_from_TOA = build_sw_lateral_bc_profile(mesh, κ_ext_sw, mesh.ngl)
    sw = SWParams(S0, μ0, 0.0, 0.35, z_prof, τ_from_TOA)
    F_dir, G_dir, Q_dir = compute_direct_beam(mesh, κ_sw, σ_sw, sw)

    case = build_shortwave_diffuse_case(κ_sw, σ_sw, g_eff, sw, F_dir, bdy; surface_albedo=0.0)
    I_cg_all, κ, σ = run_unified_case_solve(setup, metrics, comm, case; max_iter=200, tol=1e-10)

    ngl = mesh.ngl
    max_diffuse = 0.0
    for iel = 1:mesh.nelem, dd = 1:setup.n_dirs[iel], k=1:ngl,j=1:ngl,i=1:ngl
        max_diffuse = max(max_diffuse, abs(I_cg_all[(iel,dd)][lidx(i,j,k,ngl)]))
    end
    max_diffuse = MPI.Allreduce(max_diffuse, MPI.MAX, comm)
    rank == 0 && println("SW zero-scatter check — max |diffuse I| (should be ≈0): $max_diffuse")
    return max_diffuse
end

"""
Sweeps GMRES configurations (restart length, max restarts) against the
periodic manufactured case, reporting final error for each config, plus
a plain-Jacobi baseline for reference. No scattering, small κ — same
regime that exposed the slow-convergence problem.
"""
function test_unified_case_gmres_sweep(inputs, nparts, distribute, sem, partitioned_model;
                                       κ_vals::Vector{Float64}=[1e-4, 1e-5, 1e-6],
                                       restart_values::Vector{Int}=[10, 20, 40],
                                       max_restarts_values::Vector{Int}=[5, 10, 20],
                                       jacobi_max_iter::Int=30,
                                       tol::Float64=1e-9,
                                       I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model)
    mesh = setup.mesh; metrics = setup.metrics
    ngl = mesh.ngl

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    function max_err_of(I_cg_all, κ_val)
        max_err_local = 0.0
        for iel = 1:mesh.nelem, dd = 1:setup.n_dirs[iel]
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[gip],mesh.y[gip],mesh.z[gip],Ωx,Ωy,Ωz)
                max_err_local = max(max_err_local, abs(I_cg_all[(iel,dd)][lidx(i,j,k,ngl)] - I_true))
            end
        end
        return MPI.Allreduce(max_err_local, MPI.MAX, comm)
    end

    for κ_val in κ_vals
        rank == 0 && println("\n################ κ_val = $κ_val ################")

        source_fn = (x,y,z,ip,θ,ϕ) -> begin
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            grad = dI_manuf(x,y,z,Ωx,Ωy,Ωz)
            Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3] + κ_val*I_manuf(x,y,z,Ωx,Ωy,Ωz)
        end
        case = RadiativeCase(
            kappa_fn = (x,y,z,ip) -> κ_val,
            sigma_fn = (x,y,z,ip) -> 0.0,
            bc_fn = (x,y,z,ip,nx,ny,nz,θ,ϕ) -> I_manuf(x,y,z,sin(θ)*cos(ϕ),sin(θ)*sin(ϕ),cos(θ)),
            source_fn = source_fn)

        # ── Jacobi baseline ───────────────────────────────────────────────────
        rank == 0 && println("--- Jacobi baseline ---")
        I_cg_jac, _, _ = run_unified_case_solve(setup, metrics, comm, case;
                                                 max_iter=jacobi_max_iter, tol=tol, solver=:jacobi)
        err_jac = max_err_of(I_cg_jac, κ_val)
        rank == 0 && println("  [Jacobi] final max abs error: $err_jac")

        # ── GMRES configuration sweep ──────────────────────────────────────────
        for restart in restart_values, max_restarts in max_restarts_values
            rank == 0 && println("--- GMRES: restart=$restart, max_restarts=$max_restarts " *
                                  "(max matvecs = $(restart*max_restarts)) ---")
            I_cg_gmres, _, _ = run_unified_case_solve(setup, metrics, comm, case;
                                                       tol=tol, solver=:gmres,
                                                       gmres_restart=restart, gmres_max_restarts=max_restarts)
            err_gmres = max_err_of(I_cg_gmres, κ_val)
            rank == 0 && println("  [GMRES restart=$restart,max_restarts=$max_restarts] final max abs error: $err_gmres")
        end
    end
end

function test_unified_case_gmres_precond_sweep(inputs, nparts, distribute, sem, partitioned_model;
                                                κ_vals::Vector{Float64}=[1e-5],
                                                restart_values::Vector{Int}=[10, 20, 40],
                                                max_restarts_values::Vector{Int}=[5, 10, 20],
                                                jacobi_max_iter::Int=40,
                                                tol::Float64=1e-9,
                                                I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model)
    mesh = setup.mesh; metrics = setup.metrics
    ngl = mesh.ngl

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    function max_err_of(I_cg_all)
        max_err_local = 0.0
        for iel = 1:mesh.nelem, dd = 1:setup.n_dirs[iel]
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[gip],mesh.y[gip],mesh.z[gip],Ωx,Ωy,Ωz)
                max_err_local = max(max_err_local, abs(I_cg_all[(iel,dd)][lidx(i,j,k,ngl)] - I_true))
            end
        end
        return MPI.Allreduce(max_err_local, MPI.MAX, comm)   # COLLECTIVE — must run on EVERY rank
    end

    for κ_val in κ_vals
        rank == 0 && println("\n################ κ_val = $κ_val ################")
        source_fn = (x,y,z,ip,θ,ϕ) -> begin
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            grad = dI_manuf(x,y,z,Ωx,Ωy,Ωz)
            Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3] + κ_val*I_manuf(x,y,z,Ωx,Ωy,Ωz)
        end
        case = RadiativeCase(
            kappa_fn = (x,y,z,ip) -> κ_val,
            sigma_fn = (x,y,z,ip) -> 0.0,
            bc_fn = (x,y,z,ip,nx,ny,nz,θ,ϕ) -> I_manuf(x,y,z,sin(θ)*cos(ϕ),sin(θ)*sin(ϕ),cos(θ)),
            source_fn = source_fn)

        rank == 0 && println("--- Jacobi baseline ---")
        I_cg_jac, _, _ = run_unified_case_solve(setup, metrics, comm, case;
                                                 max_iter=jacobi_max_iter, tol=tol, solver=:jacobi)
        err_jac = max_err_of(I_cg_jac)
        rank == 0 && println("  [Jacobi] final max abs error: $err_jac")

        for restart in restart_values, max_restarts in max_restarts_values
            rank == 0 && println("--- GMRES (plain): restart=$restart, max_restarts=$max_restarts ---")
            I_cg_g, _, _ = run_unified_case_solve(setup, metrics, comm, case;
                                                   tol=tol, solver=:gmres,
                                                   gmres_restart=restart, gmres_max_restarts=max_restarts)
            err_g = max_err_of(I_cg_g)
            rank == 0 && println("  [GMRES] final max abs error: $err_g")

            rank == 0 && println("--- GMRES (preconditioned): restart=$restart, max_restarts=$max_restarts ---")
            I_cg_gp, _, _ = run_unified_case_solve(setup, metrics, comm, case;
                                                    tol=tol, solver=:gmres_precond,
                                                    gmres_restart=restart, gmres_max_restarts=max_restarts)
            err_gp = max_err_of(I_cg_gp)
            rank == 0 && println("  [GMRES-precond] final max abs error: $err_gp")
        end
    end
end

"""
Sweeps isotropic scattering strength (as fractions of κ) against the
periodic manufactured case, using GMRES(restart=40,max_restarts=10) —
the configuration that was sufficient even for the hardest κ=1e-6,
σ=0 case. Checks whether weak scattering reduces the matvec count
needed, as the coupling-collapses-slow-modes argument predicts.
"""
function test_unified_case_gmres_scattering_sweep(inputs, nparts, distribute, sem, partitioned_model;
                                                   κ_val::Float64=1e-6,
                                                   σ_fracs::Vector{Float64}=[1.0, 0.1, 0.01, 0.001],
                                                   gmres_restart::Int=40, gmres_max_restarts::Int=10,
                                                   tol::Float64=1e-9,
                                                   I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model)
    mesh = setup.mesh; metrics = setup.metrics
    ngl = mesh.ngl

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    function max_err_of(I_cg_all)
        max_err_local = 0.0
        for iel = 1:mesh.nelem, dd = 1:setup.n_dirs[iel]
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[gip],mesh.y[gip],mesh.z[gip],Ωx,Ωy,Ωz)
                max_err_local = max(max_err_local, abs(I_cg_all[(iel,dd)][lidx(i,j,k,ngl)] - I_true))
            end
        end
        return MPI.Allreduce(max_err_local, MPI.MAX, comm)   # COLLECTIVE — unconditional on every rank
    end

    # ── Shared angular-mesh weights for the manufactured scattering RHS
    # (isotropic phase function: Φ ≡ 1/(4π)) — built ONCE, reused for
    # every σ_frac since the angular mesh doesn't change across the sweep. ──
    w_ext_by_iel = Dict{Int,Vector{Float64}}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        w = zeros(em.extra_npoin)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
                w[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
            end
        end
        w_ext_by_iel[iel] = w
    end

    for σ_frac in σ_fracs
        σ_val = σ_frac * κ_val
        rank == 0 && println("\n################ κ_val=$κ_val, σ_val=$σ_val (σ/κ=$σ_frac) ################")

        source_fn = (x,y,z,ip,θ,ϕ) -> begin
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            grad = dI_manuf(x,y,z,Ωx,Ωy,Ωz)
            Ω_dot_grad = Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3]
            I_true = I_manuf(x,y,z,Ωx,Ωy,Ωz)

            # Isotropic scattering term: σ * (1/4π) * ∫I dΩ' = σ * G_true(x)/(4π),
            # where G_true(x) is the EXACT mean intensity of I_manuf at this
            # point — using the SAME base angular mesh (no refinement in this
            # test) across every element, so w_ext_by_iel[1] is representative.
            iel_rep = 1
            thetas = mesh.extra_mesh[iel_rep].extra_coords[1,:]
            phis   = mesh.extra_mesh[iel_rep].extra_coords[2,:]
            w = w_ext_by_iel[iel_rep]
            scatter_true = 0.0
            for m = 1:length(w)
                Ωxm = sin(thetas[m])*cos(phis[m]); Ωym = sin(thetas[m])*sin(phis[m]); Ωzm = cos(thetas[m])
                scatter_true += (1/(4π)) * w[m] * σ_val * I_manuf(x,y,z,Ωxm,Ωym,Ωzm)
            end
            Ω_dot_grad + κ_val*I_true - scatter_true
        end

        case = RadiativeCase(
            kappa_fn = (x,y,z,ip) -> κ_val,
            sigma_fn = (x,y,z,ip) -> σ_val,
            phase_fn = (θ,θ1,ϕ,ϕ1,g) -> 1/(4π),   # isotropic
            bc_fn = (x,y,z,ip,nx,ny,nz,θ,ϕ) -> I_manuf(x,y,z,sin(θ)*cos(ϕ),sin(θ)*sin(ϕ),cos(θ)),
            source_fn = source_fn)

        I_cg_gmres, _, _ = run_unified_case_solve(setup, metrics, comm, case;
                                                   tol=tol, solver=:gmres,
                                                   gmres_restart=gmres_restart, gmres_max_restarts=gmres_max_restarts)
        err_gmres = max_err_of(I_cg_gmres)
        rank == 0 && println("  [σ/κ=$σ_frac] final max abs error: $err_gmres")
    end
end

function test_unified_case_extrap_sweep(inputs, nparts, distribute, sem, partitioned_model;
                                         κ_val::Float64=1e-4, tol::Float64=1e-5, max_iter::Int=50,
                                         I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model)
    mesh = setup.mesh; metrics = setup.metrics; ngl = mesh.ngl

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    source_fn = (x,y,z,ip,θ,ϕ) -> begin
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        grad = dI_manuf(x,y,z,Ωx,Ωy,Ωz)
        Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3] + κ_val*I_manuf(x,y,z,Ωx,Ωy,Ωz)
    end
    case = RadiativeCase(kappa_fn=(x,y,z,ip)->κ_val, sigma_fn=(x,y,z,ip)->0.0,
        bc_fn=(x,y,z,ip,nx,ny,nz,θ,ϕ)->I_manuf(x,y,z,sin(θ)*cos(ϕ),sin(θ)*sin(ϕ),cos(θ)),
        source_fn=source_fn)

    function max_err_of(I_cg_all)
        e = 0.0
        for iel = 1:mesh.nelem, dd = 1:setup.n_dirs[iel]
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                e = max(e, abs(I_cg_all[(iel,dd)][lidx(i,j,k,ngl)] - I_manuf(mesh.x[gip],mesh.y[gip],mesh.z[gip],Ωx,Ωy,Ωz)))
            end
        end
        return MPI.Allreduce(e, MPI.MAX, comm)
    end

    for solver in (:jacobi_extrap, :jacobi_extrap_block)
        rank == 0 && println("\n=== solver = $solver (κ=$κ_val) ===")
        I_cg, _, _ = run_unified_case_solve(setup, metrics, comm, case;
                        max_iter=max_iter, tol=tol, solver=solver,
                        gmres_restart=40, gmres_max_restarts=10)
        err = max_err_of(I_cg)
        rank == 0 && println("  [$solver] final max abs error: $err")
    end
end

# --- reset hook, copied from your construction code ---
function make_angular_mesh_builder(inputs, backend)
    inputs[:extra_dimensions] == 2 || error("angular reset implemented for extra_dimensions == 2 only")
    return function ()
        ξω  = basis_structs_ξ_ω!(inputs[:interpolation_nodes], inputs[:extra_dimensions_order], inputs[:backend])
        basis = build_Interpolation_basis!(LagrangeBasis(), ξω.ξ, ξω.ξ, TFloat, inputs[:backend])
        return make_extra_mesh_2D(inputs[:extra_dimensions_nelemx], inputs[:extra_dimensions_nelemy],
                                  inputs[:extra_dimensions_order],
                                  inputs[:extra_dimensions_xmin], inputs[:extra_dimensions_xmax],
                                  inputs[:extra_dimensions_ymin], inputs[:extra_dimensions_ymax],
                                  basis, backend, inputs, true)
    end
end
reset_angular_mesh!(mesh, builder) = (for iel = 1:mesh.nelem; mesh.extra_mesh[iel] = builder(); end)

# --- geometry / levels ---
function ang_nodes(em, e)
    n = em.extra_nop[e] + 1
    θs = Float64[]; ϕs = Float64[]
    for a = 1:n, b = 1:n
        ip = em.extra_connijk[e, a, b]
        push!(θs, em.extra_coords[1,ip]); push!(ϕs, em.extra_coords[2,ip])
    end
    (any(<(1e-6), ϕs) && any(>(π + 1e-9), ϕs)) && (ϕs = [ϕ < 1e-6 ? 2π : ϕ for ϕ in ϕs])
    return θs, ϕs
end
ang_theta_extent(θs) = maximum(θs) - minimum(θs)
base_theta_extent(em) = maximum(ang_theta_extent(ang_nodes(em, e)[1]) for e = 1:em.extra_nelem)   # call right after a reset
ang_level_of(em, e, dθ0) = round(Int, log2(dθ0 / ang_theta_extent(ang_nodes(em, e)[1])))
ang_diameter(θs, ϕs) = hypot(ang_theta_extent(θs), (maximum(ϕs)-minimum(ϕs))*sin(sum(θs)/length(θs)))


# --- quadrature criterion ---
function ang_quad_weights(em)
    w = zeros(em.extra_npoin)
    for e = 1:em.extra_nelem
        nop = em.extra_nop[e]
        for iθ = 1:nop+1, iϕ = 1:nop+1
            w[em.extra_connijk[e,iθ,iϕ]] += em.extra_metrics.Je[e,iθ,iϕ]*em.ωθ[iθ]*em.ωθ[iϕ]
        end
    end
    return w
end

function phase_integral_ref(phase_fn, g; n::Int=20001)          # 2π ∫ Φ(Θ) sinΘ dΘ, composite Simpson
    h = π/(n-1); s = 0.0
    for k = 0:n-1
        f = phase_fn(π/2, π/2, k*h, 0.0, g) * sin(k*h)
        s += (k == 0 || k == n-1) ? f : (isodd(k) ? 4f : 2f)
    end
    return 2π*s*h/3
end

function sun_row_quad_error(em, phase_fn, g, θs, ϕs, ref)
    w = ang_quad_weights(em); s = 0.0
    for m = 1:em.extra_npoin
        s += phase_fn(θs, em.extra_coords[1,m], ϕs, em.extra_coords[2,m], g) * w[m]
    end
    return abs(s/ref - 1)
end

function levels_from_table(table, g_bucket, tol)
    errs = table[g_bucket]
    for L = 0:length(errs)-1
        errs[L+1] <= tol && return L
    end
    return length(errs) - 1
end

function run_lobe_manufactured(infra, metrics, comm, g_fn::Function;
                                κ_val=1e-3, σ_val=5e-4, A=1.0, ρ=0.8, θs=2.4, ϕs=0.7,
                                near_angle=0.25, normalize_rows::Bool=false,
                                max_iter::Int=300, tol::Float64=1e-9, verbose_iter::Bool=false)
    rank = MPI.Comm_rank(comm)
    mesh = infra.mesh
    ψ, dψ, ω = infra.sem.basis.ψ, infra.sem.basis.dψ, infra.sem.ω
    ngl = mesh.ngl; n_dirs = infra.n_dirs; npoin = mesh.npoin
    Ωs = (sin(θs)*cos(ϕs), sin(θs)*sin(ϕs), cos(θs))
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz

    # lobe_r has mean 1 over the sphere. Convolving it with HG(g) gives lobe_{g·r} exactly.
    lobe(θ, ϕ, r) = 4π * phase_function_hg(θ, θs, ϕ, ϕs, r)
    I_exact(x,y,z,θ,ϕ) = A*sin(kx*x+ky*y+kz*z)*lobe(θ,ϕ,ρ)

    κ = fill(κ_val, npoin)
    C_all = Dict{Int,Matrix{Float64}}()
    q_ext_all = Dict{Tuple{Int,Int},Vector{Float64}}()
    bc_fn_all = Dict{Tuple{Int,Int},Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]; npd = em.extra_npoin
        w = ang_quad_weights(em)
        g_iel = 0.0                                              # element-averaged g, as in the production coupling
        for k=1:ngl, j=1:ngl, i=1:ngl
            gip = mesh.connijk[iel,i,j,k]
            g_iel += g_fn(mesh.x[gip], mesh.y[gip], mesh.z[gip]) / ngl^3
        end
        θ = em.extra_coords[1,1:npd]; ϕ = em.extra_coords[2,1:npd]
        C = zeros(npd, npd)
        for dk = 1:npd, m = 1:npd
            C[dk,m] = phase_function_hg(θ[dk], θ[m], ϕ[dk], ϕ[m], g_iel) * w[m]
        end
        normalize_rows && (C ./= sum(C; dims=2))                 # every row sums to exactly 1
        C_all[iel] = C
        for dk = 1:npd
            Ωx = sin(θ[dk])*cos(ϕ[dk]); Ωy = sin(θ[dk])*sin(ϕ[dk]); Ωz = cos(θ[dk])
            lob   = lobe(θ[dk], ϕ[dk], ρ)
            lob_s = lobe(θ[dk], ϕ[dk], g_iel*ρ)                  # exact continuous scattering integral
            q_vec = zeros(npoin)
            for k=1:ngl, j=1:ngl, i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                arg = kx*mesh.x[gip] + ky*mesh.y[gip] + kz*mesh.z[gip]
                q_vec[gip] = A*cos(arg)*(Ωx*kx+Ωy*ky+Ωz*kz)*lob + κ_val*A*sin(arg)*lob - σ_val*A*sin(arg)*lob_s
            end
            q_ext_all[(iel,dk)] = q_vec
            θk = θ[dk]; ϕk = ϕ[dk]
            bc_fn_all[(iel,dk)] = (x,y,z,ip,nx,ny,nz) -> I_exact(x,y,z,θk,ϕk)
        end
    end

    I_cg_all = Dict{Tuple{Int,Int},Vector{Float64}}()
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        I_cg_all[(iel,d)] = zeros(ngl^3)
    end

    iters = 0; dg = Inf
    for it = 1:max_iter
        iters = it
        extra_source_all = Dict{Tuple{Int,Int},Vector{Float64}}()
        for iel = 1:mesh.nelem
            em = mesh.extra_mesh[iel]; C = C_all[iel]
            #local_gips = [mesh.connijk[iel,i,j,k] for k=1:ngl, j=1:ngl, i=1:ngl][:]
            local_gips = local_gips_lidx(mesh, iel)
	    for dk = 1:em.extra_npoin
                q_k = copy(q_ext_all[(iel,dk)])
                q_loc = @view q_k[local_gips]
                for m = 1:em.extra_npoin
                    q_loc .+= (σ_val*C[dk,m]) .* I_cg_all[(iel,m)]
                end
                extra_source_all[(iel,dk)] = q_k
            end
        end
        x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
            infra.neighbor_elem, infra.neighbor_fid, infra.face_correspondence, infra.periodic_face_correspondence,
            infra.neighbor_rank, infra.neighbor_kind, infra.neighbor_gel, infra.neighbor_gel_group, infra.neighbor_is_periodic,
            infra.same_rank_kind, infra.same_rank_partner, infra.same_rank_is_periodic,
            infra.ncf_as_child, infra.ncf_L1L2, infra.ncf_as_parent,
            infra.remote_gid_coords,
            infra.angular_deps, infra.angular_deps_cr, infra.remote_angular_mesh,
            infra.gid_exch, infra.periodic_gid_exch,
            κ, bc_fn_all, comm;
            extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all)
        dl = 0.0
        for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
            dl = max(dl, maximum(abs.(x_dg[(iel,dd)] .- I_cg_all[(iel,dd)])))
            I_cg_all[(iel,dd)] = copy(x_dg[(iel,dd)])
        end
        dg = MPI.Allreduce(dl, MPI.MAX, comm)                    # collective, unconditional
        verbose_iter && rank == 0 && println("    iter $it: max_delta=$dg")
        dg < tol && break
    end

    err_all = 0.0; err_near = 0.0; err_far = 0.0
    for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
        em = mesh.extra_mesh[iel]; θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
        near = ang_angle(Ωs, θ, ϕ) <= near_angle
        for k=1:ngl, j=1:ngl, i=1:ngl
            gip = mesh.connijk[iel,i,j,k]
            e = abs(I_cg_all[(iel,dd)][lidx(i,j,k,ngl)] - I_exact(mesh.x[gip],mesh.y[gip],mesh.z[gip],θ,ϕ))
            err_all = max(err_all, e)
            near ? (err_near = max(err_near, e)) : (err_far = max(err_far, e))
        end
    end
    err_all  = MPI.Allreduce(err_all,  MPI.MAX, comm)
    err_near = MPI.Allreduce(err_near, MPI.MAX, comm)
    err_far  = MPI.Allreduce(err_far,  MPI.MAX, comm)
    ndirs    = MPI.Allreduce(sum(n_dirs), MPI.SUM, comm)
    return (err_all=err_all, err_near=err_near, err_far=err_far, ndirs=ndirs,
            iters=iters, converged=(dg < tol), peak=A*(1+ρ)/(1-ρ)^2)
end

function test_adaptive_angular_checks(inputs, nparts, distribute, sem, partitioned_model;
        κ_val=1e-3, σ_val=5e-4, ρ=0.5, θs=2.4, ϕs=0.7,
        tol_quad=1e-3, max_level=3, include_harsh::Bool=false, uniform_level::Int=0,
        max_iter::Int=300, tol::Float64=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm()); comm = get_mpi_comm()
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model)   # no refinement
    mesh = setup.mesh; metrics = setup.metrics
    builder = make_angular_mesh_builder(inputs, inputs[:backend])
    npoin = mesh.npoin
    σ_node = fill(σ_val, npoin); F_node = ones(npoin)
    xmid = mesh.xmin + 0.5*(mesh.xmax-mesh.xmin); wtr = (mesh.xmax-mesh.xmin)/8

    scenarios = Any[(label="constant g=0.6, ρ=$ρ",        g_lo=0.6, g_hi=0.6, ρ=ρ),
                    (label="varying g 0.2->0.6, ρ=$ρ",    g_lo=0.2, g_hi=0.6, ρ=ρ)]
    include_harsh && push!(scenarios, (label="harsh: g 0.2->0.85, ρ=0.8", g_lo=0.2, g_hi=0.85, ρ=0.8))

    rows = NamedTuple[]
    function run_stage(infra, g_fn, sc, near_angle, stage, normalize)
        r = run_lobe_manufactured(infra, metrics, comm, g_fn; κ_val=κ_val, σ_val=σ_val, ρ=sc.ρ, θs=θs, ϕs=ϕs,
                                  near_angle=near_angle, normalize_rows=normalize, max_iter=max_iter, tol=tol)
        push!(rows, (scenario=sc.label, stage=stage, normalize=normalize, r=r))
        rank == 0 && println("  [$stage, normalize=$normalize] dirs=$(r.ndirs) iters=$(r.iters) converged=$(r.converged)  " *
            "err_near=$(r.err_near) err_far=$(r.err_far) err_all=$(r.err_all)  (near/peak=$(r.err_near/r.peak))")
        return r
    end

    rank == 0 && println("\n=== CHECK 0: quadrature row sums on the base angular mesh ===")
    phase_row_sum_report(mesh, phase_function_hg, [0.2, 0.6, 0.85], comm)

    for (si, sc) in enumerate(scenarios)
        rank == 0 && println("\n################ scenario: $(sc.label) ################")
        g_fn(x,y,z) = sc.g_lo + (sc.g_hi-sc.g_lo)*0.5*(1 + tanh((x-xmid)/wtr))
        g_node = [g_fn(mesh.x[ip], mesh.y[ip], mesh.z[ip]) for ip = 1:npoin]
        near_angle = 2*phase_lobe_halfwidth(phase_function_hg, sc.ρ)      # twice the solution lobe's half-width
        rank == 0 && println("  near_angle = $near_angle rad")

        infraA = if si == 1
            setup                                                         # mesh is still the base mesh here
        else
            reset_angular_mesh!(mesh, builder)
            build_unified_ncf_infrastructure(mesh, setup.sem, comm)
        end
        rank == 0 && println("--- A: base angular mesh ---")
        resA = run_stage(infraA, g_fn, sc, near_angle, :A_base, false)
        run_stage(infraA, g_fn, sc, near_angle, :A_base, true)

        if si == 1                                                        # CHECK 1: max_level=0 must reproduce A exactly
            rank == 0 && println("--- CHECK 1: adaptive path with max_level=0 must reproduce A ---")
            infra0, _, _ = adaptive_angular_regrid!(mesh, setup.sem, infraA, comm, phase_function_hg,
                                                    g_node, σ_node, F_node, θs, ϕs, builder; tol_quad=tol_quad, max_level=0)
            res0 = run_stage(infra0, g_fn, sc, near_angle, :B_level0, false)
            same = res0.ndirs == resA.ndirs && isapprox(res0.err_all, resA.err_all; rtol=1e-10)
            rank == 0 && println("  CHECK 1: " * (same ? "✓ PASS" : "✗ FAIL (dirs $(res0.ndirs) vs $(resA.ndirs), err $(res0.err_all) vs $(resA.err_all)) — bug in reset/rebuild path"))
        end

        rank == 0 && println("--- B: adaptive (max_level=$max_level, tol_quad=$tol_quad) ---")
        infraB, _, info = adaptive_angular_regrid!(mesh, setup.sem, infraA, comm, phase_function_hg,
                                                   g_node, σ_node, F_node, θs, ϕs, builder; tol_quad=tol_quad, max_level=max_level)
        println("  rank $rank: dirs $(info.n_dirs_before) -> $(info.n_dirs_after); sun-row quad error vs level: " *
                join(["g=$g: $(round.(e; sigdigits=3))" for (g,e) in sort(collect(info.table))], "; "))
        rank == 0 && println("  levels used: min=$(minimum(values(info.levels))) max=$(maximum(values(info.levels)))")
        phase_row_sum_report(mesh, phase_function_hg, unique([sc.g_lo, sc.g_hi]), comm)
        run_stage(infraB, g_fn, sc, near_angle, :B_adaptive, false)
        run_stage(infraB, g_fn, sc, near_angle, :B_adaptive, true)

        if uniform_level > 0                                              # optional; cost grows like 4^level
            rank == 0 && println("--- C: uniform level $uniform_level ---")
            reset_angular_mesh!(mesh, builder)
            dθ0 = mesh.nelem > 0 ? base_theta_extent(mesh.extra_mesh[1]) : 0.0
            for iel = 1:mesh.nelem
                refine_lobe_to_level!(mesh, iel, θs, ϕs, 4.0, uniform_level, dθ0)    # R=4 > π covers the sphere
            end
            infraC = build_unified_ncf_infrastructure(mesh, setup.sem, comm)
            run_stage(infraC, g_fn, sc, near_angle, :C_uniform, false)
            run_stage(infraC, g_fn, sc, near_angle, :C_uniform, true)
        end
    end

    if rank == 0
        println("\n================ SUMMARY (near-sun error, B vs A, ratio < 1 means refinement helped) ================")
        for sc in scenarios, nrm in (false, true)
            a = first(r for r in rows if r.scenario == sc.label && r.stage == :A_base && r.normalize == nrm).r
            b = first(r for r in rows if r.scenario == sc.label && r.stage == :B_adaptive && r.normalize == nrm).r
            println("$(sc.label), normalize=$nrm: near $(b.err_near/a.err_near)  far $(b.err_far/a.err_far)  " *
                    "dirs $(b.ndirs/a.ndirs)x  iters $(a.iters) -> $(b.iters)")
        end
    end
    return rows
end


function phase_lobe_halfwidth(phase_fn::Function, g::Float64; nscan::Int=4000)
    f0 = phase_fn(π/2, π/2, 0.0, 0.0, g)               # θ=θ1=π/2, so cosΘ = cos(ϕ-ϕ1)
    for k = 1:nscan
        Θ = π*k/nscan
        phase_fn(π/2, π/2, Θ, 0.0, g) < 0.5*f0 && return Θ
    end
    return Float64(π)
end

function balance_angular_levels(L_local::Dict{Int,Int}, global_adj, gel_to_class, forced_zero::Set{Int}, comm)
    flat = Float64[]
    for (k,v) in L_local; push!(flat, Float64(k)); push!(flat, Float64(v)); end
    all_n = MPI.Allgather(Int32(length(L_local)), comm)
    all_flat = MPI.Allgatherv(flat, Cint.(all_n .* 2), comm)
    L = Dict{Int,Int}()
    for m = 1:(length(all_flat)÷2)
        L[Int(round(all_flat[2m-1]))] = Int(round(all_flat[2m]))
    end
    for g in forced_zero; haskey(L,g) && (L[g] = 0); end
    changed = true
    while changed                                          # lowering only, so it converges and is identical on every rank
        changed = false
        for g in keys(L)
            new = L[g]
            for n in get(global_adj, g, ()); haskey(L, Int(n)) && (new = min(new, L[Int(n)] + 1)); end
            for c in get(gel_to_class, g, ()); haskey(L, Int(c)) && (new = min(new, L[Int(c)])); end
            new < L[g] && (L[g] = new; changed = true)
        end
    end
    return L
end

sun_direction(μ0, φ0) = (θs = π - acos(clamp(μ0,0.0,1.0)); (sin(θs)*cos(φ0), sin(θs)*sin(φ0), cos(θs)))
ang_angle(Ωs, θ, ϕ) = acos(clamp(Ωs[1]*sin(θ)*cos(ϕ) + Ωs[2]*sin(θ)*sin(ϕ) + Ωs[3]*cos(θ), -1.0, 1.0))
ang_elem_size(θs, ϕs) = max(maximum(θs)-minimum(θs), (maximum(ϕs)-minimum(ϕs))*sin(sum(θs)/length(θs)))
ang_min_dist(Ωs, θs, ϕs) = minimum(ang_angle(Ωs, θs[k], ϕs[k]) for k in eachindex(θs))

function phase_row_sum_report(mesh, phase_fn, g_list, comm)
    rank = MPI.Comm_rank(comm)
    for g in g_list
        worst = 0.0; seen = Set{Int}()
        for iel = 1:mesh.nelem
            em = mesh.extra_mesh[iel]; npd = em.extra_npoin
            npd in seen && continue; push!(seen, npd)
            w = ang_quad_weights(em)
            θ = em.extra_coords[1,1:npd]; ϕ = em.extra_coords[2,1:npd]
            for k = 1:npd
                s = 0.0
                for m = 1:npd; s += phase_fn(θ[k], θ[m], ϕ[k], ϕ[m], g) * w[m]; end
                worst = max(worst, abs(s - 1))
            end
        end
        worst = MPI.Allreduce(worst, MPI.MAX, comm)
        rank == 0 && println("  g=$g: max |row sum - 1| over all rows = $worst")
    end
end

function ang_region_dist(θsun, ϕsun, θlo, θhi, ϕlo, ϕhi)      # angle from the sun to the element's (θ,ϕ) rectangle
    θc = clamp(θsun, θlo, θhi)
    ϕc = if ϕlo <= ϕsun <= ϕhi
        ϕsun
    else
        dlo = min(abs(ϕsun-ϕlo), 2π-abs(ϕsun-ϕlo)); dhi = min(abs(ϕsun-ϕhi), 2π-abs(ϕsun-ϕhi))
        dlo <= dhi ? ϕlo : ϕhi
    end
    return acos(clamp(sin(θsun)*sin(θc)*cos(ϕsun-ϕc) + cos(θsun)*cos(θc), -1.0, 1.0))
end

function refine_lobe_to_level!(mesh, iel, θsun, ϕsun, R, L, dθ0; max_passes::Int=500)
    for lvl = 1:L
        for _ = 1:max_passes
            em = mesh.extra_mesh[iel]; hit = 0
            for e = 1:em.extra_nelem
                ang_level_of(em, e, dθ0) >= lvl && continue
                θs, ϕs = ang_nodes(em, e)
                ang_region_dist(θsun, ϕsun, minimum(θs), maximum(θs), minimum(ϕs), maximum(ϕs)) <= R || continue
                hit = e; break
            end
            hit == 0 && break
            n_before = em.extra_nelem
            adapt_single_angular_element_conforming!(mesh, iel, hit)
            verify_angular_phi_seam_conformity(mesh, iel)
            mesh.extra_mesh[iel].extra_nelem > n_before || error("angular refinement made no progress (iel=$iel, e_ext=$hit)")
        end
    end
end

function sun_row_error_table!(mesh, builder, phase_fn, g_buckets, θsun, ϕsun, lobe_k, max_level, dθ0)
    table = Dict{Float64,Vector{Float64}}()
    for g in g_buckets
        ref = phase_integral_ref(phase_fn, g)
        R = lobe_k * phase_lobe_halfwidth(phase_fn, g)
        mesh.extra_mesh[1] = builder()                          # element 1 is the scratch mesh
        errs = zeros(max_level+1)
        for L = 0:max_level
            L > 0 && refine_lobe_to_level!(mesh, 1, θsun, ϕsun, R, L, dθ0)
            errs[L+1] = sun_row_quad_error(mesh.extra_mesh[1], phase_fn, g, θsun, ϕsun, ref)
        end
        table[g] = errs
    end
    mesh.extra_mesh[1] = builder()
    return table
end

function adaptive_angular_regrid!(mesh, sem, infra, comm, phase_fn, g_node, σ_node, F_node,
                                   θsun::Float64, ϕsun::Float64, builder::Function;
                                   tol_quad::Float64=1e-3, max_level::Int=3,
                                   strength_tol::Float64=0.05, lobe_k::Float64=2.0)
    ngl = mesh.ngl
    snapshots = Dict{Int,NamedTuple}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        snapshots[iel] = (coords = copy(em.extra_coords[:, 1:em.extra_npoin]), npoin = em.extra_npoin)
    end
    n_before = MPI.Allreduce(sum((s.npoin for s in values(snapshots)); init=0), +, comm)

    s_el = zeros(mesh.nelem); g_el = zeros(mesh.nelem)
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        gip = mesh.connijk[iel,i,j,k]
        s_el[iel] += σ_node[gip]*F_node[gip]/ngl^3; g_el[iel] += g_node[gip]/ngl^3
    end
    s_ref = MPI.Allreduce(maximum(s_el; init=0.0), MPI.MAX, comm)
    active = [s_ref > 0 && s_el[iel] > strength_tol*s_ref for iel = 1:mesh.nelem]
    bucket(g) = round(g; digits=2)
    g_buckets = sort(unique(Float64[bucket(g_el[iel]) for iel = 1:mesh.nelem if active[iel]]))
    hw = Dict(g => phase_lobe_halfwidth(phase_fn, g) for g in g_buckets)

    reset_angular_mesh!(mesh, builder)
    dθ0 = mesh.nelem > 0 ? base_theta_extent(mesh.extra_mesh[1]) : 0.0
    table = (mesh.nelem > 0 && !isempty(g_buckets)) ?
            sun_row_error_table!(mesh, builder, phase_fn, g_buckets, θsun, ϕsun, lobe_k, max_level, dθ0) :
            Dict{Float64,Vector{Float64}}()
    reset_angular_mesh!(mesh, builder)

    L_loc = Dict{Int,Int}()
    for iel = 1:mesh.nelem
        L_loc[Int(mesh.el2gel[iel])] = active[iel] ? levels_from_table(table, bucket(g_el[iel]), tol_quad) : 0
    end
    sp = Set{Int}()
    for (iel,_) in keys(infra.ncf_as_child);  push!(sp, iel); end
    for (iel,_) in keys(infra.ncf_as_parent); push!(sp, iel); end
    gel_to_class = build_periodic_equivalence_classes(mesh, infra.same_rank_kind, infra.same_rank_partner,
        infra.same_rank_is_periodic, infra.neighbor_kind, infra.neighbor_gel_group, infra.neighbor_is_periodic, comm)
    global_adj, global_refined = build_global_element_adjacency_and_refined_set(mesh, infra.neighbor_elem, sp,
        infra.neighbor_rank, infra.neighbor_gel_group, gel_to_class, comm)
    forced_zero = Set{Int}(Int(g) for g in global_refined)
    for g in global_refined, n in get(global_adj, g, ()); push!(forced_zero, Int(n)); end
    L = balance_angular_levels(L_loc, global_adj, gel_to_class, forced_zero, comm)

    for iel = 1:mesh.nelem
        lvl = L[Int(mesh.el2gel[iel])]; lvl == 0 && continue
        refine_lobe_to_level!(mesh, iel, θsun, ϕsun, lobe_k*hw[bucket(g_el[iel])], lvl, dθ0)
    end
    infra_new = build_unified_ncf_infrastructure(mesh, sem, comm)
    n_after = MPI.Allreduce(sum((mesh.extra_mesh[i].extra_npoin for i = 1:mesh.nelem); init=0), +, comm)
    return infra_new, snapshots, (levels=L, table=table, n_dirs_before=n_before, n_dirs_after=n_after)
end

function test_ncf_local_vs_remote(setup, comm; ε::Float64=1.0, n_dirs_sample::Int=3, top::Int=6)
    MPI.Comm_size(comm) == 1 || error("run on 1 rank, so every NCF interface is same-rank")
    mesh = setup.mesh; ngl = mesh.ngl; Lx = mesh.xmax - mesh.xmin
    mid = isodd(ngl) ? (ngl + 1) ÷ 2 : 0
    ds = unique(round.(Int, range(1, minimum(setup.n_dirs); length=n_dirs_sample)))

    function field(iel, d)
        gel = Int(mesh.el2gel[iel]); v = zeros(ngl^3)
        for k=1:ngl, j=1:ngl, i=1:ngl
            ip = mesh.connijk[iel,i,j,k]
            v[lidx(i,j,k,ngl)] = sin(2π*(mesh.x[ip] + 0.7*mesh.y[ip] + 0.4*mesh.z[ip])/Lx + 0.3*d) + ε*((gel*7919 + d*31) % 13)/13
        end
        return v
    end
    x_dg = Dict{Tuple{Int,Int},Vector{Float64}}((iel,d) => field(iel,d) for iel in 1:mesh.nelem for d in ds)

    function face_data!(rgc, gr, iel, fid)            # what the cross-rank exchange would deliver for this face
        gel = TInt(mesh.el2gel[iel]); irange, jrange, krange = face_ip_range(fid, ngl)
        haskey(rgc, gel) || (rgc[gel] = Dict{Int,NTuple{3,Float64}}())
        for k in krange, j in jrange, i in irange
            ip = mesh.connijk[iel,i,j,k]; gid = Int(mesh.ip2gip[ip])
            rgc[gel][gid] = (mesh.x[ip], mesh.y[ip], mesh.z[ip])
            for d in ds; gr[(gel,d,gid)] = x_dg[(iel,d)][lidx(i,j,k,ngl)]; end
        end
    end
    newrgc() = Dict{TInt,Dict{Int,NTuple{3,Float64}}}()
    newgr()  = Dict{Tuple{TInt,Int,Int},Float64}()
    clear_face_plans!()

    # ---- parent -> child (child is the fine element) ----
    rgc = newrgc(); gr = newgr(); nf = 0; worst_fine = 0.0
    for ((c_iel, c_fid), idx) in setup.ncf_as_child; face_data!(rgc, gr, Int(mesh.pip[idx]), OPP_FACE[c_fid]); end
    for ((c_iel, c_fid), idx) in setup.ncf_as_child
        p_iel = Int(mesh.pip[idx]); L1, L2 = setup.ncf_L1L2[idx]
        pl = fine_plan!(mesh, rgc, c_iel, c_fid, TInt(mesh.el2gel[p_iel]), ngl)
        for d in ds
            loc = build_dg_ncf_parent_to_child_local(mesh, x_dg, c_iel, p_iel, c_fid, d, mesh.IPc_list[:,idx], mesh.IPp_list[:,idx], L1, L2, ngl)
            out = Vector{Float64}(undef, ngl*ngl)
            eval_fine_plan!(out, pl, d, gr, ngl) || error("fine plan blocked with complete data")
            worst_fine = max(worst_fine, maximum(abs.(loc .- out))); nf += 1
        end
    end

    # ---- child -> parent (parent is the coarse element) ----
    rgc2 = newrgc(); gr2 = newgr(); rows = Tuple[]; nc = 0
    for ((p_iel, p_fid), group) in setup.ncf_as_parent, idx in group
        face_data!(rgc2, gr2, Int(mesh.cip[idx]), Int(mesh.lfid[idx]))
    end
    for ((p_iel, p_fid), group) in setup.ncf_as_parent
        child_gels = TInt[TInt(mesh.el2gel[Int(mesh.cip[idx])]) for idx in group]
        entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx]) for idx in group]
        pl = coarse_plan!(mesh, rgc2, p_iel, p_fid, child_gels, ngl)
        for d in ds
            loc = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, p_iel, p_fid, d, entries, ngl)
            out = Matrix{Float64}(undef, ngl, ngl)
            eval_coarse_plan!(out, pl, d, gr2, ngl) || error("coarse plan blocked with complete data")
            Δ = abs.(loc .- out); nc += 1
            nd = count(>(1e-12), Δ); on_mid = 0; off_mid = 0
            for p1 in 1:ngl, p2 in 1:ngl
                Δ[p1,p2] > 1e-12 || continue
                (p1 == mid || p2 == mid) ? (on_mid += 1) : (off_mid += 1)
            end
            push!(rows, (maximum(Δ), Int(mesh.el2gel[p_iel]), p_fid, d, nd, on_mid, off_mid))
        end
    end
    sort!(rows; by = r -> -r[1])
    println("NCF same-rank vs remote-plan check: parent->child $nf evaluations, worst |Δ| = $worst_fine")
    println("child->parent $nc evaluations, $(count(r -> r[1] > 1e-12, rows)) differ, worst |Δ| = $(isempty(rows) ? 0.0 : rows[1][1]) (ε=$ε)")
    println("differing coarse-face nodes: $(sum(r -> r[6], rows; init=0)) on child-child midlines, $(sum(r -> r[7], rows; init=0)) elsewhere")
    for r in rows[1:min(top, length(rows))]
        println("  parent gel $(r[2]) fid $(r[3]) d=$(r[4]): max|Δ|=$(r[1]), $(r[5]) nodes differ ($(r[6]) on midlines, $(r[7]) off)")
    end
    clear_face_plans!()
end

function ncf_equivalence_job(inputs, nparts, distribute, sem, partitioned_model; seed=42, frac=0.01, max_e_ext=4)
    spatial, _, _ = generate_stress_test_targets(sem.mesh.gnelem; seed=seed, frac=frac, max_e_ext=max_e_ext)
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=collect(Int, spatial), angular_targets=Tuple{Int,Int}[])
    test_ncf_local_vs_remote(setup, get_mpi_comm())
end
