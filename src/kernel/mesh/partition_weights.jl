# Element weights for p4est repartitioning by cost (RT sweeps + dynamics).
#
# Flow (rt_rebalance_strategy!, or an AMR step when weights are available):
#   PART_W_OLD    per element of the mesh the weights were computed on (PART_W_MESH), in its element order;
#                 set after each RT call by rt_rebalance_check! (RT_sweep.jl).
#   mod_mesh_mesh_driver (AMR re-adapt path) maps them through the adapt step's new-to-old element map
#                 (children inherit the parent's weight), decides whether to redistribute (weighted imbalance above
#                 PART_THRESHOLD, or PART_FORCE) and hands them to GridapP4est.redistribute as cell weights.
#   The p4est partition itself is p8est_partition with the weight callback (families kept), or with
#   JEXPRESSO_P4EST_KEEP_FAMILIES=2 weighted cuts between whole trees (Jexpresso.jl).
const PART_W_OLD     = Ref{Any}(nothing)   # Vector{Float64} (µs per element per RT interval) or nothing
const PART_W_MESH    = Ref{Any}(nothing)   # the mesh PART_W_OLD belongs to
const PART_W_CELLS   = Ref{Any}(nothing)   # weights mapped to the mesh being rebuilt (element order)
const PART_FORCE     = Ref(false)          # redistribute in the next rebuild regardless of the imbalance
const PART_THRESHOLD = Ref(1.15)           # weighted max/mean above which a rebuild redistributes
const PART_WEIGHTED_LAST = Ref(false)      # the last rebuild redistributed with weights
const PART_LAST_GLUE = Ref{Any}(nothing)   # redistribution glue of the last direct redistribution (for element data)

# Start of a rebuild from `omesh`: take the weights if they belong to it
function part_weights_begin!(omesh)
    PART_WEIGHTED_LAST[] = false; PART_LAST_GLUE[] = nothing
    w = PART_W_OLD[]
    ok = w !== nothing && PART_W_MESH[] === omesh && length(w) == omesh.nelem
    PART_W_CELLS[] = ok ? copy(w) : nothing
    PART_W_OLD[] = nothing; PART_W_MESH[] = nothing            # do not keep the old mesh alive
    (w !== nothing && MPI.Comm_rank(get_mpi_comm()) == 0) &&
        println(" # p4est weighted repartition: cost weights ", ok ? "taken for the rebuild" : "DROPPED (they belong to another mesh)")
    return nothing
end

# After an adapt step with new-to-old element map n2o: weights follow the elements (children inherit the parent's).
# n2o covers the model's local cells, owned first, then ghosts; only the first nelem_new (owned) entries are used,
# and an entry pointing outside the old owned range (a ghost) gets the mean weight.
function part_weights_map!(n2o, nelem_new::Integer)
    w = PART_W_CELLS[]
    w === nothing && return nothing
    if !(n2o isa AbstractVector) || length(n2o) < nelem_new
        MPI.Comm_rank(get_mpi_comm()) == 0 && println(" # p4est weighted repartition: weights DROPPED (new-to-old map $(typeof(n2o)), length $(n2o isa AbstractVector ? length(n2o) : -1) < $nelem_new)")
        PART_W_CELLS[] = nothing; return nothing
    end
    wm = isempty(w) ? 1.0 : sum(w) / length(w); nout = 0
    wn = Vector{Float64}(undef, nelem_new)
    for k = 1:nelem_new
        o = Int(n2o[k])
        if 1 <= o <= length(w); wn[k] = w[o]; else; wn[k] = wm; nout += 1; end
    end
    no = MPI.Allreduce(nout, MPI.SUM, get_mpi_comm())
    (no > 0 && MPI.Comm_rank(get_mpi_comm()) == 0) && println(" # p4est weighted repartition: $no new elements map outside the old owned range (mean weight)")
    PART_W_CELLS[] = wn
    return nothing
end

# Decide whether the rebuilt mesh (before redistribution) must be redistributed; collective.
# With weights the decision is the weighted imbalance (element counts alone do not trigger it).
function part_weights_decide!(mesh_tmp, comm)
    w = PART_W_CELLS[]
    if w === nothing                     # no weights pending; a weighted redistribution just done stays as it is
        PART_WEIGHTED_LAST[] && (mesh_tmp.lneed_redistribute = false)
        return mesh_tmp.lneed_redistribute
    end
    length(w) == mesh_tmp.nelem || (PART_W_CELLS[] = nothing; return mesh_tmp.lneed_redistribute)
    s = MPI.Allgather(sum(w), comm); np = length(s)
    imb = maximum(s) / max(sum(s) / np, 1e-300)
    need = PART_FORCE[] || imb > PART_THRESHOLD[]
    MPI.Comm_rank(comm) == 0 && println(@sprintf(" # p4est weighted repartition: predicted cost max/mean %.3f (threshold %.2f%s) -> %s",
                                                 imb, PART_THRESHOLD[], PART_FORCE[] ? ", forced" : "",
                                                 need ? "redistributing with weights" : "kept"))
    mesh_tmp.lneed_redistribute = need
    need || (PART_W_CELLS[] = nothing)
    return need
end

# Cell weights for GridapP4est.redistribute: one Cint per local cell of the model (owned + ghost; only owned are read).
# Jexpresso's element iel is the model's owned cell own_to_local(indices)[iel]. Consumes PART_W_CELLS.
function part_cell_weights(model)
    w = PART_W_CELLS[]
    PART_W_CELLS[] = nothing; PART_FORCE[] = false
    PART_WEIGHTED_LAST[] = w !== nothing
    w === nothing && return nothing
    comm = get_mpi_comm()
    wmax = MPI.Allreduce(maximum(w; init=0.0), MPI.MAX, comm)
    scale = wmax > 0 ? 1.0e6 / wmax : 1.0                       # integer weights up to 1e6 (p4est sums in int64)
    return map(local_views(model.dmodel.models), partition(get_cell_gids(model.dmodel))) do lmodel, indices
        o2l = own_to_local(indices)
        length(o2l) == length(w) || error("part_cell_weights: $(length(w)) weights for $(length(o2l)) owned cells")
        wc = ones(Cint, num_cells(lmodel))
        for (iel, lc) in enumerate(o2l); wc[lc] = Cint(max(1, round(Int, w[iel] * scale))); end
        wc
    end
end
