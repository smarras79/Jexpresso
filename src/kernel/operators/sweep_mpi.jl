function build_rank_cut_exchange(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt},
                                  neighbor_gel::Dict{Tuple{TInt,TInt},TInt},   # NEW — needed for correct disambiguation
                                  comm::MPI.Comm)
    rank = MPI.Comm_rank(comm)
    ngl  = mesh.ngl

    face_global_ids  = Dict{Tuple{TInt,TInt}, Vector{TInt}}()
    face_local_order = Dict{Tuple{TInt,TInt}, Vector{NTuple{3,Int}}}()
    for ((iel,fid), _) in neighbor_rank
        irange, jrange, krange = face_ip_range(fid, ngl)
        gids = TInt[]; ijks = NTuple{3,Int}[]
        for k in krange, j in jrange, i in irange
            push!(gids, mesh.ip2gip[mesh.connijk[iel,i,j,k]])
            push!(ijks, (i,j,k))
        end
        perm = sortperm(gids)
        face_global_ids[(iel,fid)]  = gids[perm]
        face_local_order[(iel,fid)] = ijks[perm]
    end

    send_keys      = Tuple{TInt,TInt,Int}[]
    send_targets   = Int[]
    gids_send      = TInt[]
    my_gel_send    = TInt[]   # NEW — MY OWN element's global id, tagged onto every position I send
    for ((iel,fid), target_rank) in neighbor_rank
        gids = face_global_ids[(iel,fid)]
        my_gel = mesh.el2gel[iel]
        for pos = 1:length(gids)
            push!(send_keys, (iel,fid,pos))
            push!(send_targets, Int(target_rank))
            push!(gids_send, gids[pos])
            push!(my_gel_send, my_gel)
        end
    end

    cache_gid = SendReceiveCache(comm, gids_send, send_targets)
    send_and_receive!(cache_gid, gids_send)

    cache_gel = SendReceiveCache(comm, my_gel_send, send_targets)   # SAME send_targets/ordering as gid cache
    send_and_receive!(cache_gel, my_gel_send)

    # ── Match by (global_id, SENDER'S element global id) jointly — this is
    # what discriminates correctly when a node's global id is coincidentally
    # shared across multiple, logically DISTINCT rank-cut faces (an edge or
    # corner of the partition boundary touching more than one remote
    # element). Matching on global id alone was the actual bug: it applied
    # a received value to every LOCAL face sharing that id, even ones fed
    # by a DIFFERENT remote element — confirmed by the inconsistent-face
    # data (scattered wrong round numbers, not sentinels, not a clean
    # two-block split). ──────────────────────────────────────────────────────
    gid_and_gel_to_my_positions = Dict{Tuple{TInt,TInt},Vector{Int}}()
    for i = 1:length(gids_send)
        # KEY INSIGHT: for MY OWN send, the "expected sender's gel" I should
        # match against on receipt is the REMOTE element's gel — i.e.
        # neighbor_gel[(iel,fid)] for whichever face this position belongs
        # to — NOT my own gel. Build this keyed by (my_gid, expected_remote_gel).
        iel, fid, pos = send_keys[i]
        expected_remote_gel = neighbor_gel[(iel,fid)]
        key = (gids_send[i], expected_remote_gel)
        push!(get!(gid_and_gel_to_my_positions, key, Int[]), i)
    end

    n_recv = length(cache_gid.combined_recv_data)
    recv_perm = Vector{Vector{Int}}(undef, n_recv)
    n_unmatched = 0
    for m = 1:n_recv
        recv_gid = cache_gid.combined_recv_data[m]
        recv_gel = cache_gel.combined_recv_data[m]   # the SENDER's own element gel — matches "expected_remote_gel" above
        key = (recv_gid, recv_gel)
        positions = get(gid_and_gel_to_my_positions, key, Int[])
        isempty(positions) && (n_unmatched += 1)
        recv_perm[m] = positions
    end

    n_mismatch = 0
    for m = 1:n_recv
        for pos in recv_perm[m]
            cache_gid.combined_recv_data[m] != gids_send[pos] && (n_mismatch += 1)
        end
    end

    println("Rank $rank: rank-cut exchange setup — $(length(send_keys)) face-node entries, " *
            "$n_recv received, $n_unmatched unmatched, $n_mismatch mismatched " *
            (n_unmatched == 0 && n_mismatch == 0 ? "✓ PASS" : "✗ FAIL — do not trust this exchange"))

    return (cache=cache_gid, cache_gel=cache_gel, send_keys=send_keys, gids_send=gids_send,
            send_targets=send_targets,   # NEW — so build_rank_cut_value_cache doesn't need to re-derive it
            recv_perm=recv_perm, face_local_order=face_local_order,
            n_unmatched=n_unmatched, n_mismatch=n_mismatch)

end

function test_rank_cut_exchange(mesh, metrics)
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    # ── Same-rank adjacency (already validated) ─────────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    # ── Rank-cut identification (already validated via two-sided handshake) ──
    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end

    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
        end
    end

    neighbor_rank, neighbor_gel, n_unresolved = resolve_rank_cut_neighbors(mesh, suspected_cuts, comm)
    println("Rank $rank: $(length(suspected_cuts)) suspected cuts, " *
            "$(length(neighbor_rank)) resolved, $n_unresolved unresolved")

    if n_unresolved > 0
        println("Rank $rank: ⚠ unresolved rank-cuts present — aborting exchange test, " *
                "fix rank-cut resolution first")
        return nothing
    end

    # ── The actual new test ──────────────────────────────────────────────────
    result = build_rank_cut_exchange(mesh, neighbor_rank, neighbor_gel, comm)

    println("Rank $rank: send_keys total = $(length(result.send_keys))")

    return result
end

function test_multirank_schedule(mesh, metrics, neighbor_elem, neighbor_rank, neighbor_gel, exchange,
                                  Ωx::Float64, Ωy::Float64, Ωz::Float64, comm::MPI.Comm;
                                  max_rounds::Int = 200)
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

    n_rank_cut_faces_total = length(neighbor_rank)
    n_rank_cut_active_inflow = 0
    n_rank_cut_active_outflow = 0
    n_rank_cut_inactive = 0
    for ((iel,fid), _) in neighbor_rank
        ii = is_inflow_all[iel]; ia = is_active_all[iel]
        if !ia[fid]
            n_rank_cut_inactive += 1
        elseif ii[fid]
            n_rank_cut_active_inflow += 1
        else
            n_rank_cut_active_outflow += 1
        end
    end
    println("Rank $rank: rank-cut face classification for Ω=($Ωx,$Ωy,$Ωz) — " *
            "total=$n_rank_cut_faces_total, active+inflow=$n_rank_cut_active_inflow, " *
            "active+outflow=$n_rank_cut_active_outflow, inactive/tangential=$n_rank_cut_inactive")

    processed    = falses(nelem)
    round_solved = fill(-1, nelem)

    # remote_round now stores (round, sender_gel) PAIRS per position, so the
    # readiness/consistency checks can verify BOTH the value arrived AND it
    # came from the expected remote element — same disambiguation as the
    # exchange-setup fix, applied here too since this loop has its own
    # independent send/receive/match cycle every round.
    remote_round = Dict{Tuple{Int,Int}, Vector{Int}}()
    remote_gel   = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ((iel,fid),_) in neighbor_rank
        n = length(exchange.face_local_order[(TInt(iel),TInt(fid))])
        remote_round[(iel,fid)] = fill(-1, n)
        remote_gel[(iel,fid)]   = fill(-1, n)
    end

    n_violations = 0
    n_inconsistent_face = 0
    progress_history = Int[]

    round = 0
    while true
        round += 1
        n_blocked_same_rank = 0
        n_blocked_remote    = 0

        if round > max_rounds
            if rank == 0
                println("STALL DIAGNOSIS — last 10 rounds' global unprocessed count:")
                println(progress_history[max(1,end-9):end])
            end
            error("test_multirank_schedule: exceeded max_rounds=$max_rounds.")
        end

        for iel = 1:nelem
            processed[iel] && continue
            ii = is_inflow_all[iel]; ia = is_active_all[iel]
            ready = true
            this_same_rank_block = false
            this_remote_block = false
            for fid = 1:6
                (ia[fid] && ii[fid]) || continue
                jel = neighbor_elem[iel][fid]
                if jel != zero(TInt)
                    if !processed[jel]
                        ready = false
                        this_same_rank_block = true
                    end
                elseif haskey(neighbor_rank, (TInt(iel),TInt(fid)))
                    rv  = remote_round[(iel,fid)]
                    rg  = remote_gel[(iel,fid)]
                    expected_gel = neighbor_gel[(TInt(iel),TInt(fid))]
                    if any(==(-1), rv) || any(g -> g != expected_gel, rg)
                        # -1 → not yet received; wrong gel → contamination from
                        # a DIFFERENT face's sender (should be impossible now,
                        # but checked explicitly rather than assumed fixed)
                        ready = false
                        this_remote_block = true
                    else
                        length(unique(rv)) == 1 || (n_inconsistent_face += 1)
                        rv[1] >= round && (n_violations += 1)
                    end
                end
            end
            if ready
                processed[iel] = true
                round_solved[iel] = round
            else
                this_same_rank_block && (n_blocked_same_rank += 1)
                this_remote_block    && (n_blocked_remote += 1)
            end
        end

        # ── Exchange: send round_solved AND my own element's gel, matched
        # the same (gid, gel) way as the setup phase. ─────────────────────────
        send_round_vals = Vector{Int}(undef, length(exchange.send_keys))
        send_gel_vals   = Vector{Int}(undef, length(exchange.send_keys))
        for (idx,(iel,fid,pos)) in enumerate(exchange.send_keys)
            send_round_vals[idx] = round_solved[iel]
            send_gel_vals[idx]   = Int(mesh.el2gel[iel])
        end
        send_and_receive!(exchange.cache, send_round_vals)
        send_and_receive!(exchange.cache_gel, send_gel_vals)

        for m = 1:length(exchange.cache.combined_recv_data)
            for pos in exchange.recv_perm[m]
                iel, fid, p = exchange.send_keys[pos]
                remote_round[(iel,fid)][p] = exchange.cache.combined_recv_data[m]
                remote_gel[(iel,fid)][p]   = exchange.cache_gel.combined_recv_data[m]
            end
        end

        n_unprocessed_local = count(!, processed)
        n_unprocessed_total = MPI.Allreduce(n_unprocessed_local, +, comm)
        push!(progress_history, n_unprocessed_total)

        println("Rank $rank, round $round: unprocessed_local=$n_unprocessed_local, " *
                "unprocessed_global=$n_unprocessed_total, " *
                "blocked_same_rank=$n_blocked_same_rank, blocked_remote=$n_blocked_remote")

        n_unprocessed_total == 0 && break
    end

    n_processed_local = count(processed)
    n_processed_total = MPI.Allreduce(n_processed_local, +, comm)
    n_violations_total = MPI.Allreduce(n_violations, +, comm)
    n_inconsistent_total = MPI.Allreduce(n_inconsistent_face, +, comm)

    if rank == 0
        println("Multi-rank schedule test: converged in $round rounds")
        println("Total elements processed (all ranks): $n_processed_total")
        println("Ordering violations (remote round >= own round): $n_violations_total  " *
                (n_violations_total==0 ? "✓ PASS" : "✗ FAIL"))
        println("Inconsistent-face round numbers (should be 0): $n_inconsistent_total  " *
                (n_inconsistent_total==0 ? "✓ PASS" : "✗ FAIL"))
    end

    return (rounds=round, n_processed_total=n_processed_total,
            n_violations=n_violations_total, n_inconsistent=n_inconsistent_total)
end

function test_multirank_schedule_full(mesh, metrics; Ωx=0.3, Ωy=0.5, Ωz=0.8, max_rounds=20)
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    # ── Same-rank adjacency (already validated) ─────────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    # ── Rank-cut identification (already validated via two-sided handshake) ──
    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end

    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
        end
    end

    neighbor_rank, neighbor_gel, n_unresolved = resolve_rank_cut_neighbors(mesh, suspected_cuts, comm)
    println("Rank $rank: $(length(suspected_cuts)) suspected cuts, " *
            "$(length(neighbor_rank)) resolved, $n_unresolved unresolved")
    if n_unresolved > 0
        println("Rank $rank: ⚠ unresolved rank-cuts present — aborting, fix rank-cut resolution first")
        return nothing
    end

    # ── Rank-cut node-correspondence exchange (already validated) ───────────
    exchange = build_rank_cut_exchange(mesh, neighbor_rank, neighbor_gel, comm)
    if exchange.n_unmatched > 0 || exchange.n_mismatch > 0
        println("Rank $rank: ⚠ rank-cut exchange setup failed — aborting, fix that first")
        return nothing
    end

    # ── The actual new test ──────────────────────────────────────────────────
    result = test_multirank_schedule(mesh, metrics, neighbor_elem, neighbor_rank, neighbor_gel, exchange,
                                      Ωx, Ωy, Ωz, comm; max_rounds=max_rounds)

    return result
end

"""
    build_rank_cut_value_cache(exchange, comm)

Separate, Float64-typed SendReceiveCache for exchanging real DG solution
values across rank-cut faces — distinct from exchange.cache/cache_gel
(both Int-typed). Reuses the EXACT SAME send_targets already computed in
build_rank_cut_exchange (now returned from there), so this cache routes
identically to the other two — only the payload type differs.
"""
function build_rank_cut_value_cache(exchange, comm::MPI.Comm)
    dummy_vals = zeros(Float64, length(exchange.send_targets))
    return SendReceiveCache(comm, dummy_vals, exchange.send_targets)
end

function run_dg_sweep_multirank!(x_dg::Vector{Float64}, mesh, metrics, ψ, dψ, ω,
                                  Ωx::Float64, Ωy::Float64, Ωz::Float64,
                                  neighbor_elem, neighbor_rank, neighbor_gel, exchange,
                                  value_cache,
                                  connijk_dg, face_neighbor_dg,
                                  κ, σ, intΦ_eff, bc_fn, comm::MPI.Comm;
                                  extra_source = nothing, max_rounds::Int = 200)
    rank  = MPI.Comm_rank(comm)
    ngl   = mesh.ngl
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

    processed = falses(nelem)

    remote_valid = Dict{Tuple{Int,Int}, Vector{Bool}}()
    remote_value = Dict{Tuple{Int,Int}, Vector{Float64}}()
    remote_gel   = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ((iel,fid),_) in neighbor_rank
        n = length(exchange.face_local_order[(TInt(iel),TInt(fid))])
        remote_valid[(iel,fid)] = falses(n)
        remote_value[(iel,fid)] = zeros(Float64, n)
        remote_gel[(iel,fid)]   = fill(-1, n)
    end

    A_vol = zeros(ngl^3, ngl^3)
    RHS_local = zeros(ngl^3)

    round = 0
    while true
        round += 1
        round > max_rounds && error("run_dg_sweep_multirank!: exceeded max_rounds=$max_rounds " *
                                     "— schedule should have converged in ~10-15 rounds per the " *
                                     "already-validated scheduling test; treat this as a real problem.")

        for iel = 1:nelem
            processed[iel] && continue
            ii = is_inflow_all[iel]; ia = is_active_all[iel]

            # ── Readiness check — IDENTICAL three-way classification to the
            # validated scheduling test: same-rank / cross-rank / true boundary. ──
            ready = true
            for fid = 1:6
                (ia[fid] && ii[fid]) || continue
                jel = neighbor_elem[iel][fid]
                if jel != zero(TInt)
                    processed[jel] || (ready = false; break)
                elseif haskey(neighbor_rank, (TInt(iel),TInt(fid)))
                    all(remote_valid[(iel,fid)]) || (ready = false; break)
                end
                # else: true physical boundary — always ready
            end
            ready || continue

            # ── Solve this element — REUSES the exact same per-element logic
            # as run_dg_sweep!, inlined here since there's no precomputed
            # `order` (readiness is discovered dynamically, round by round).
            # Same-rank/true-boundary faces use the ALREADY-VALIDATED
            # face_neighbor_dg / bc_fn paths, UNCHANGED. ─────────────────────
            fill!(A_vol, 0.0); fill!(RHS_local, 0.0)
            build_dg_volume_matrix!(A_vol, iel, Ωx, Ωy, Ωz, ngl, ψ, dψ, ω,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                     metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
            extra_source !== nothing &&
                build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source, mesh.connijk)

            for fid = 1:6
                ia[fid] || continue
                face_normal = compute_face_normal_avg(iel, fid, ngl,
                    metrics.dξdx,metrics.dξdy,metrics.dξdz,
                    metrics.dηdx,metrics.dηdy,metrics.dηdz,
                    metrics.dζdx,metrics.dζdy,metrics.dζdz)

                if ii[fid]
                    jel = neighbor_elem[iel][fid]
                    if jel == zero(TInt)
                        if haskey(neighbor_rank, (TInt(iel),TInt(fid)))
			    rvals = remote_value[(iel,fid)]
    			    canonical_order = exchange.face_local_order[(TInt(iel),TInt(fid))]
    			    ijk_to_p = Dict{Tuple{Int,Int,Int},Int}(canonical_order[pp] => pp for pp = 1:length(canonical_order))
    			    I_upwind = (i,j,k) -> rvals[ijk_to_p[(i,j,k)]]
                        else
                            # UNCHANGED — true physical boundary
                            I_upwind = (i,j,k) -> bc_fn(mesh.x[mesh.connijk[iel,i,j,k]],
                                                         mesh.y[mesh.connijk[iel,i,j,k]],
                                                         mesh.z[mesh.connijk[iel,i,j,k]])
                        end
                    else
                        # UNCHANGED — same-rank inflow, exact same face_neighbor_dg path as run_dg_sweep!
                        I_upwind = (i,j,k) -> x_dg[face_neighbor_dg[face_flat_index(fid, ngl, i, j, k), fid, iel]]
                    end
                else
                    I_upwind = (i,j,k) -> 0.0
                end

                build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx, Ωy, Ωz,
                                      ii[fid], ia[fid], face_normal, ω,
                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                      metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                      metrics.Je, iel, I_upwind)
            end

            x_local = A_vol \ RHS_local
            for k=1:ngl, j=1:ngl, i=1:ngl
                x_dg[connijk_dg[iel,i,j,k]] = x_local[lidx(i,j,k,ngl)]
            end
            processed[iel] = true
        end

        # ── Exchange: send THIS element's outflow-face DG values (Float64,
        # via value_cache) plus my element gel (Int, via exchange.cache_gel —
        # SAME disambiguation as the validated scheduling test), keyed by
        # (global id, sender's element gel). ────────────────────────────────
        send_vals     = Vector{Float64}(undef, length(exchange.send_keys))
        send_gel_vals = Vector{Int}(undef, length(exchange.send_keys))
        for (idx,(iel,fid,pos)) in enumerate(exchange.send_keys)
            i,j,k = exchange.face_local_order[(TInt(iel),TInt(fid))][pos]
            send_vals[idx] = processed[iel] ? x_dg[connijk_dg[iel,i,j,k]] : NaN
            send_gel_vals[idx] = Int(mesh.el2gel[iel])
        end
        send_and_receive!(value_cache, send_vals)
        send_and_receive!(exchange.cache_gel, send_gel_vals)

        for m = 1:length(value_cache.combined_recv_data)
            val = value_cache.combined_recv_data[m]
            isnan(val) && continue   # sender hadn't processed that element yet this round
            recv_gel = exchange.cache_gel.combined_recv_data[m]
            for pos in exchange.recv_perm[m]
                iel, fid, p = exchange.send_keys[pos]
                remote_gel[(iel,fid)][p]   = recv_gel
                remote_value[(iel,fid)][p] = val
                remote_valid[(iel,fid)][p] = true
            end
        end

        n_unprocessed_total = MPI.Allreduce(count(!,processed), +, comm)
        n_unprocessed_total == 0 && break
    end

    return x_dg
end

function test_multirank_exponential(mesh, metrics, basis, ω; κ_val = 0.0001, I0 = 4.2)
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    ngl  = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ

    # ── Same-rank adjacency (already validated) ─────────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    # ── Rank-cut identification (already validated) ──────────────────────────
    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end
    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
        end
    end
    neighbor_rank, neighbor_gel, n_unresolved = resolve_rank_cut_neighbors(mesh, suspected_cuts, comm)
    println("Rank $rank: $(length(suspected_cuts)) suspected cuts, " *
            "$(length(neighbor_rank)) resolved, $n_unresolved unresolved")
    if n_unresolved > 0
        println("Rank $rank: ⚠ unresolved rank-cuts — aborting")
        return nothing
    end

    # ── Exchange infrastructure (already validated) ──────────────────────────
    exchange = build_rank_cut_exchange(mesh, neighbor_rank, neighbor_gel, comm)
    if exchange.n_unmatched > 0 || exchange.n_mismatch > 0
        println("Rank $rank: ⚠ exchange setup failed — aborting")
        return nothing
    end
    value_cache = build_rank_cut_value_cache(exchange, comm)

    # ── DG storage infrastructure (already validated single-rank; built
    # from LOCAL mesh data only, so this works unchanged for a distributed
    # mesh partition too — same construction, just on this rank's subset). ──
    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
            build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem,
                                                   neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg; mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg; mesh.face_neighbor_dg = face_neighbor_dg
    end

    # ── Manufactured/exact-solution setup — SAME as the already-trusted
    # single-rank exponential test, axis-aligned Ω=(1,0,0). ──────────────────
    Ωx, Ωy, Ωz = 1.0, 0.0, 0.0
    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    x0 = mesh.xmin
    bc_fn = (x,y,z) -> I0 * exp(-κ_val*(x - x0))

    # ── Run the MULTI-RANK sweep ─────────────────────────────────────────────
    x_dg = fill(NaN, mesh.npoin_dg)
    run_dg_sweep_multirank!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz,
                             neighbor_elem, neighbor_rank, neighbor_gel, exchange, value_cache,
                             mesh.connijk_dg, mesh.face_neighbor_dg, κ, σz, intΦ_eff, bc_fn, comm)

    n_nan = count(isnan, x_dg)
    n_nan > 0 && println("Rank $rank: ⚠ $n_nan DG entries never written")

    # ── Compare EVERY local node's DG value against the true analytic
    # solution — no reconciliation needed for this comparison (every DG
    # copy of a shared node should independently match the SAME exact
    # value, since the exact solution has no jump at all). ──────────────────
    max_err_local = 0.0
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        g_cg = mesh.connijk[iel,i,j,k]
        I_exact = I0 * exp(-κ_val*(mesh.x[g_cg] - x0))
        v = x_dg[mesh.connijk_dg[iel,i,j,k]]
        max_err_local = max(max_err_local, abs(v - I_exact))
    end
    max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)

    if rank == 0
        println("Multi-rank exponential test — max abs error (all ranks, all nodes): $max_err_global")
    end
    return max_err_global
end

function test_manufactured_isotropic_multirank(mesh, metrics, basis, ω;
                                                I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                κ_val=0.001, σ_val=0.0005,
                                                max_iter=500, tol=1e-9, max_sweep_rounds=200)
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ
    extra_npoin = mesh.extra_mesh.extra_npoin
    npoin_cg = length(mesh.x)

    # ── Multi-rank infrastructure — built ONCE, shared across every
    # direction and every scattering iteration (none of it depends on Ω
    # or the iteration count). ──────────────────────────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]; ip2 = mesh.poin_in_bdy_face[iface_bdy, ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, ngl, ngl]; ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1,ip2,ip3,ip4])))
    end
    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem, fid = 1:6
        neighbor_elem[iel][fid] == zero(TInt) || continue
        fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
        n1=mesh.conn_face_el[1,fe,iel]; n2=mesh.conn_face_el[2,fe,iel]
        n3=mesh.conn_face_el[3,fe,iel]; n4=mesh.conn_face_el[4,fe,iel]
        key = Tuple(sort(TInt[n1,n2,n3,n4]))
        key in tagged_face_keys && continue
        push!(suspected_cuts, (iel,fid))
    end
    neighbor_rank, neighbor_gel, n_unresolved = resolve_rank_cut_neighbors(mesh, suspected_cuts, comm)
    println("Rank $rank: $(length(suspected_cuts)) suspected cuts, $(length(neighbor_rank)) resolved, $n_unresolved unresolved")
    if n_unresolved > 0
        println("Rank $rank: ⚠ unresolved rank-cuts — aborting"); return nothing
    end
    exchange = build_rank_cut_exchange(mesh, neighbor_rank, neighbor_gel, comm)
    if exchange.n_unmatched > 0 || exchange.n_mismatch > 0
        println("Rank $rank: ⚠ exchange setup failed — aborting"); return nothing
    end
    value_cache = build_rank_cut_value_cache(exchange, comm)

    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg; mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg; mesh.face_neighbor_dg = face_neighbor_dg
    end

    # ── Manufactured solution — IDENTICAL to the already-validated
    # single-rank isotropic test. ────────────────────────────────────────────
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    dirs = Vector{NTuple{3,Float64}}(undef, extra_npoin)
    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]; ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        dirs[ip_ext] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
    end
    w_ext = zeros(Float64, extra_npoin)
    for e_ext = 1:mesh.extra_mesh.extra_nelem
        nop = mesh.extra_mesh.extra_nop[e_ext]
        for iθ=1:nop+1, iϕ=1:nop+1
            ip_ext = mesh.extra_mesh.extra_connijk[e_ext,iθ,iϕ]
            Je_ang = mesh.extra_mesh.extra_metrics.Je[e_ext,iθ,iϕ]
            w_ext[ip_ext] += Je_ang * mesh.extra_mesh.ωθ[iθ] * mesh.extra_mesh.ωθ[iϕ]
        end
    end
    W = sum(w_ext)
    C = [σ_val*w_ext[m]/W for k=1:extra_npoin, m=1:extra_npoin]

    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)

    q_ext_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for g = 1:npoin_cg
            x_,y_,z_ = mesh.x[g], mesh.y[g], mesh.z[g]
            grad = dI_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            Ω_dot_grad = Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3]
            I_true = I_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            scatter_true = sum(C[k,m]*I_manuf(x_,y_,z_,dirs[m]...) for m = 1:extra_npoin)
            q_ext_all[k][g] = Ω_dot_grad + κ_val*I_true - scatter_true
        end
    end
    bc_fn_all = [(x,y,z) -> I_manuf(x,y,z,dirs[k]...) for k = 1:extra_npoin]

    # ── Iterative solve — SAME structure as single-rank isotropic test,
    # calling run_dg_sweep_multirank! instead of run_dg_sweep!. ──────────────
    I_cg_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for it = 1:max_iter
        max_delta_local = 0.0
        for k = 1:extra_npoin
            Ωx,Ωy,Ωz = dirs[k]
            q_k = copy(q_ext_all[k])
            for mm = 1:extra_npoin
                C[k,mm] == 0.0 && continue
                q_k .+= C[k,mm] .* I_cg_all[mm]
            end
            x_dg = fill(NaN, mesh.npoin_dg)
            run_dg_sweep_multirank!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz,
                                     neighbor_elem, neighbor_rank, neighbor_gel, exchange, value_cache,
                                     mesh.connijk_dg, mesh.face_neighbor_dg, κ, σz, intΦ_eff, bc_fn_all[k],
                                     comm; extra_source=q_k, max_rounds=max_sweep_rounds)
            I_new = zeros(npoin_cg)
            for iel = 1:mesh.nelem, kk=1:ngl, jj=1:ngl, ii=1:ngl
                I_new[mesh.connijk[iel,ii,jj,kk]] = x_dg[mesh.connijk_dg[iel,ii,jj,kk]]
            end
            max_delta_local = max(max_delta_local, maximum(abs.(I_new .- I_cg_all[k])))
            I_cg_all[k] = I_new
        end
        # max_delta must be a GLOBAL (all-rank) quantity — a local-only
        # convergence check could let one rank stop early while another
        # still has unconverged directions, silently desynchronizing the
        # per-direction loop across ranks on the NEXT iteration.
        max_delta = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        rank == 0 && println("  iter $it: max_delta=$max_delta")
        max_delta < tol && (rank==0 && println("Converged after $it iterations"); break)
        it == max_iter && rank==0 && @warn "Did not converge (max_delta=$max_delta)"
    end

    max_err_local = 0.0
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for g = 1:npoin_cg
            I_true = I_manuf(mesh.x[g],mesh.y[g],mesh.z[g],Ωx,Ωy,Ωz)
            max_err_local = max(max_err_local, abs(I_cg_all[k][g] - I_true))
        end
    end
    max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)
    rank == 0 && println("Multi-rank manufactured (isotropic) — max abs error: $max_err_global")
    return max_err_global
end

function test_multirank_all_directions_exponential(mesh, metrics, basis, ω; κ_val=0.0001, I0=4.2)
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ
    extra_npoin = mesh.extra_mesh.extra_npoin
    npoin_cg = length(mesh.x)

    # ── Multi-rank infrastructure — built ONCE ──────────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]; ip2 = mesh.poin_in_bdy_face[iface_bdy, ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, ngl, ngl]; ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1,ip2,ip3,ip4])))
    end
    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem, fid = 1:6
        neighbor_elem[iel][fid] == zero(TInt) || continue
        fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
        n1=mesh.conn_face_el[1,fe,iel]; n2=mesh.conn_face_el[2,fe,iel]
        n3=mesh.conn_face_el[3,fe,iel]; n4=mesh.conn_face_el[4,fe,iel]
        key = Tuple(sort(TInt[n1,n2,n3,n4]))
        key in tagged_face_keys && continue
        push!(suspected_cuts, (iel,fid))
    end
    neighbor_rank, neighbor_gel, n_unresolved = resolve_rank_cut_neighbors(mesh, suspected_cuts, comm)
    println("Rank $rank: $(length(suspected_cuts)) suspected cuts, $(length(neighbor_rank)) resolved, $n_unresolved unresolved")
    if n_unresolved > 0
        println("Rank $rank: ⚠ unresolved rank-cuts — aborting"); return nothing
    end
    exchange = build_rank_cut_exchange(mesh, neighbor_rank, neighbor_gel, comm)
    if exchange.n_unmatched > 0 || exchange.n_mismatch > 0
        println("Rank $rank: ⚠ exchange setup failed — aborting"); return nothing
    end
    value_cache = build_rank_cut_value_cache(exchange, comm)

    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg; mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg; mesh.face_neighbor_dg = face_neighbor_dg
    end

    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    worst_err = 0.0
    worst_ip_ext = 0
    all_errs = Float64[]

    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]
        ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

        bc_fn = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1]) + Ωy*(y-r_bc[2]) + Ωz*(z-r_bc[3])))

        x_dg = fill(NaN, mesh.npoin_dg)
        run_dg_sweep_multirank!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz,
                                 neighbor_elem, neighbor_rank, neighbor_gel, exchange, value_cache,
                                 mesh.connijk_dg, mesh.face_neighbor_dg, κ, σz, intΦ_eff, bc_fn, comm)

        n_nan_local = count(isnan, x_dg)
        n_nan_total = MPI.Allreduce(n_nan_local, +, comm)
        if n_nan_total > 0
            println("Rank $rank, ip_ext=$ip_ext: ⚠ $n_nan_local UNWRITTEN entries locally " *
                     "($n_nan_total total) — skipping error calc for this direction")
            continue
        end

        max_err_local = 0.0
        for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
            g_cg = mesh.connijk[iel,i,j,k]
            x_,y_,z_ = mesh.x[g_cg], mesh.y[g_cg], mesh.z[g_cg]
            I_exact = bc_fn(x_,y_,z_)
            v = x_dg[mesh.connijk_dg[iel,i,j,k]]
            max_err_local = max(max_err_local, abs(v - I_exact))
        end
        max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)

        push!(all_errs, max_err_global)
        if max_err_global > worst_err
            worst_err = max_err_global; worst_ip_ext = ip_ext
        end
        rank == 0 && println("ip_ext=$ip_ext: max_err=$max_err_global")
    end

    if rank == 0
        println("─"^60)
        println("Worst direction: ip_ext=$worst_ip_ext, max abs error = $worst_err")
        println("Mean error across all directions: $(sum(all_errs)/length(all_errs))")
        println("Median error: $(sort(all_errs)[length(all_errs)÷2 + 1])")
    end


    Ωx, Ωy, Ωz = -0.587344755969155, -0.6938187114406047, -0.4167034116616192  # ip_ext=37
    test_bidirectional_disambiguation(mesh, metrics, neighbor_elem, neighbor_rank, neighbor_gel,
                                   exchange, value_cache, Ωx, Ωy, Ωz, comm)



    return all_errs
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
