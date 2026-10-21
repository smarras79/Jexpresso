const USE_INCREMENTAL_PACK = Ref(true)

struct GidTable <: AbstractDict{Tuple{TInt,Int,Int},Float64}
    vals::Vector{Float64}                          # aliases val_cache.combined_recv_data; NaN = not valid
    slot::Dict{Tuple{TInt,Int,Int},Int32}          # static: (remote_gel, d, gid) -> slot
    gel::Vector{Float64}; d::Vector{Float64}; gid::Vector{Float64}     # aliases of the static key caches
end
@inline gid_key_at(t::GidTable, m::Int) = (TInt(round(Int, t.gel[m])), round(Int, t.d[m]), round(Int, t.gid[m]))
@inline function Base.get(t::GidTable, key, default)
    s = get(t.slot, key, Int32(0))
    s == 0 && return default
    v = @inbounds t.vals[s]
    return isnan(v) ? default : v
end
@inline function Base.haskey(t::GidTable, key)
    s = get(t.slot, key, Int32(0))
    return s != 0 && !isnan(@inbounds t.vals[s])
end
@inline function Base.getindex(t::GidTable, key)
    s = get(t.slot, key, Int32(0))
    (s != 0 && !isnan(@inbounds t.vals[s])) || throw(KeyError(key))
    return @inbounds t.vals[s]
end
Base.length(t::GidTable) = count(!isnan, t.vals)                       # O(n), not for hot paths
function Base.iterate(t::GidTable, s::Int=1)                           # may yield a duplicate key twice
    n = length(t.vals)
    while s <= n
        v = @inbounds t.vals[s]
        isnan(v) || return (gid_key_at(t, s) => v, s + 1)
        s += 1
    end
    return nothing
end

function make_gid_table(exch)
    vc = exch.val_cache; n = length(vc.combined_recv_data)
    gelv = exch.gel_cache.combined_recv_data; dv = exch.d_cache.combined_recv_data; gidv = exch.gid_cache.combined_recv_data
    slot = Dict{Tuple{TInt,Int,Int},Int32}(); sizehint!(slot, n)
    for m in 1:n
        k = (TInt(round(Int, gelv[m])), round(Int, dv[m]), round(Int, gidv[m]))
        haskey(slot, k) || (slot[k] = Int32(m))        # duplicates of one key carry identical values
    end
    return GidTable(vc.combined_recv_data, slot, gelv, dv, gidv)
end

function augment_gid_exchange(exch, ngl::Int)
    n = exch.n_flat
    index = Dict{Tuple{Int,Int},Int}()
    pair_iel = Int[]; pair_d = Int[]; counts = Int[]
    pid = Vector{Int}(undef, n)
    for m in 1:n
        key = (exch.flat_owning_iel[m], exch.flat_d[m])
        p = get!(index, key) do
            push!(pair_iel, key[1]); push!(pair_d, key[2]); push!(counts, 0)
            length(pair_iel)
        end
        pid[m] = p; counts[p] += 1
    end
    np = length(pair_iel)
    ptr = Vector{Int}(undef, np + 1); ptr[1] = 1
    for p in 1:np; ptr[p+1] = ptr[p] + counts[p]; end
    fillpos = ptr[1:np]
    entry_slot = Vector{Int32}(undef, n); entry_lidx = Vector{Int32}(undef, n)
    for m in 1:n
        p = pid[m]; e = fillpos[p]; fillpos[p] += 1
        entry_slot[e] = Int32(m)
        entry_lidx[e] = Int32(lidx(exch.flat_i[m], exch.flat_j[m], exch.flat_k[m], ngl))
    end
    return merge(exch, (pair_iel=pair_iel, pair_d=pair_d, pair_ptr=ptr, pair_index=index, entry_slot=entry_slot, entry_lidx=entry_lidx,
                        n_pairs=np, send_buf=Vector{Float64}(undef, n), table=make_gid_table(exch)))
end

function exchange_ncf_gid_table!(exch, x_dg::Dict{Tuple{Int,Int},Vector{Float64}}, ngl::Int)
    send = exch.send_buf
    fill!(send, NaN)
    @inbounds for p in 1:exch.n_pairs
        v = get(x_dg, (exch.pair_iel[p], exch.pair_d[p]), nothing)
        v === nothing && continue
        for e in exch.pair_ptr[p]:(exch.pair_ptr[p+1]-1)
            send[exch.entry_slot[e]] = v[exch.entry_lidx[e]]
        end
    end
    send_and_receive!(exch.val_cache, send)
    return exch.table
end

# Profiled variant: same behaviour, with the pack / MPI / unpack timers and fresh-entry statistics.
#=function exchange_ncf_gid_values_profiled!(exch, x_dg, ngl, prof::SweepProfile, prev_valid::BitVector)
    t_all = time_ns()
    t0 = ptic()
    send = exch.send_buf
    fill!(send, NaN)
    @inbounds for p in 1:exch.n_pairs
        v = get(x_dg, (exch.pair_iel[p], exch.pair_d[p]), nothing)
        v === nothing && continue
        for e in exch.pair_ptr[p]:(exch.pair_ptr[p+1]-1)
            send[exch.entry_slot[e]] = v[exch.entry_lidx[e]]
        end
    end
    n_valid = 0; n_fresh = 0
    @inbounds for m in 1:exch.n_flat
        ok = !isnan(send[m]); n_valid += ok; n_fresh += (ok && !prev_valid[m]); prev_valid[m] = ok
    end
    ptoc!(prof, :ex_pack, t0)
    t0 = ptic(); send_and_receive!(exch.val_cache, send); ptoc!(prof, :ex_mpi, t0)
    t0 = ptic(); tbl = exch.table; ptoc!(prof, :ex_unpack, t0)           # nothing left to do; the row stays for comparison
    return (received=tbl, sent=exch.n_flat, valid=n_valid, fresh=n_fresh, seconds=(time_ns() - t_all)*1e-9)
end=#

function exchange_ncf_gid_values_profiled!(exch, x_dg, ngl, prof::SweepProfile, prev_valid::BitVector;
                                           newly::Union{Nothing,Vector{Tuple{Int,Int}}}=nothing,
                                           valid_count::Ref{Int}=Ref(0))
    t_all = time_ns()
    t0 = ptic()
    send = exch.send_buf
    n_valid = 0; n_fresh = 0
    if newly === nothing
        fill!(send, NaN)
        @inbounds for p in 1:exch.n_pairs
            v = get(x_dg, (exch.pair_iel[p], exch.pair_d[p]), nothing)
            v === nothing && continue
            for e in exch.pair_ptr[p]:(exch.pair_ptr[p+1]-1)
                send[exch.entry_slot[e]] = v[exch.entry_lidx[e]]
            end
        end
        @inbounds for m in 1:exch.n_flat
            ok = !isnan(send[m]); n_valid += ok; n_fresh += (ok && !prev_valid[m]); prev_valid[m] = ok
        end
    else
        for (iel, d) in newly
            p = get(exch.pair_index, (iel, d), 0)
            p == 0 && continue
            v = x_dg[(iel, d)]
            for e in exch.pair_ptr[p]:(exch.pair_ptr[p+1]-1)
                send[exch.entry_slot[e]] = v[exch.entry_lidx[e]]
                n_fresh += 1
            end
        end
        valid_count[] += n_fresh; n_valid = valid_count[]
    end
    ptoc!(prof, :ex_pack, t0)
    t0 = ptic(); send_and_receive!(exch.val_cache, send); ptoc!(prof, :ex_mpi, t0)
    t0 = ptic(); tbl = exch.table; ptoc!(prof, :ex_unpack, t0)
    return (received=tbl, sent=exch.n_flat, valid=n_valid, fresh=n_fresh, seconds=(time_ns() - t_all)*1e-9)
end

function verify_gid_table(setup, I_final, comm)
    ngl = setup.mesh.ngl; res = Float64[]
    for exch in (setup.gid_exch, setup.periodic_gid_exch)
        old = exchange_ncf_gid_values!(exch, I_final, ngl)               # old Dict path, collective
        tbl = exchange_ncf_gid_table!(exch, I_final, ngl)                # new path, collective
        n_bad = 0; n_extra = 0
        for (k, v) in old
            reinterpret(UInt64, get(tbl, k, NaN)) == reinterpret(UInt64, v) || (n_bad += 1)
        end
        for m in eachindex(tbl.vals)
            isnan(tbl.vals[m]) && continue
            haskey(old, gid_key_at(tbl, m)) || (n_extra += 1)
        end
        append!(res, Float64[length(old), n_bad, n_extra])
    end
    t = MPI.Allreduce(res, +, comm)
    return MPI.Comm_rank(comm) == 0 ?
        String[@sprintf("gid-table verification: main %d entries, %d mismatches, %d extra | periodic %d entries, %d mismatches, %d extra (mismatches and extra expected 0)",
                        Int(t[1]), Int(t[2]), Int(t[3]), Int(t[4]), Int(t[5]), Int(t[6]))] : String[]
end

function exchange_ncf_gid_table_incremental!(exch, x_dg, newly::Vector{Tuple{Int,Int}}, ngl::Int)
    send = exch.send_buf
    for (iel, d) in newly
        p = get(exch.pair_index, (iel, d), 0)
        p == 0 && continue
        v = x_dg[(iel, d)]
        for e in exch.pair_ptr[p]:(exch.pair_ptr[p+1]-1)
            send[exch.entry_slot[e]] = v[exch.entry_lidx[e]]
        end
    end
    send_and_receive!(exch.val_cache, send)
    return exch.table
end
