using Gridap
using GridapDistributed
using PartitionedArrays
using Base.Threads
using Printf

# ===========================================================================
# USER SETTING
#
# SEND_COORDS: controls what Julia sends to Alya at every timestep.
#
#   false  →  send velocity components [u, v, ...]   (normal coupled run)
#   true   →  send coordinates         [x, y, ...]   (verification / debug)
#
# Both options produce a buffer of identical size [npoin × ndime].
# Alya receives the same buffer either way and knows nothing about this flag.
# ===========================================================================
const SEND_COORDS = true
#const SEND_COORDS = false

# ===========================================================================
# GLOBAL MPI COMMUNICATOR REFS
# ===========================================================================
const JEXPRESSO_MPI_COMM       = Ref{Union{MPI.Comm,Nothing}}(nothing)
const JEXPRESSO_MPI_COMM_WORLD = Ref{Union{MPI.Comm,Nothing}}(nothing)
const JEXPRESSO_COUPLING_DATA  = Ref{Union{Dict{Symbol,Any},Nothing}}(nothing)

# Pre-fetched JLD2 cache data — populated before with_mpi (see je_prefetch_caches!)
# to move JLD2 JIT cost and disk I/O out of the Alya-blocking critical path.
const JEXPRESSO_PREFETCHED_MESH_CACHE = Ref{Union{Nothing, Dict}}(nothing)
const JEXPRESSO_PREFETCHED_SEM_CACHE  = Ref{Union{Nothing, Tuple}}(nothing)

function set_mpi_comm(comm::MPI.Comm)
    JEXPRESSO_MPI_COMM[] = comm
end
function get_mpi_comm()
    return JEXPRESSO_MPI_COMM[] === nothing ? MPI.COMM_WORLD : JEXPRESSO_MPI_COMM[]
end

function set_mpi_comm_world(comm::MPI.Comm)
    JEXPRESSO_MPI_COMM_WORLD[] = comm
end
function get_mpi_comm_world()
    return JEXPRESSO_MPI_COMM_WORLD[] === nothing ? MPI.COMM_WORLD : JEXPRESSO_MPI_COMM_WORLD[]
end

function set_coupling_data(data::Dict{Symbol,Any})
    JEXPRESSO_COUPLING_DATA[] = data
end
function get_coupling_data()
    return JEXPRESSO_COUPLING_DATA[]
end

# ===========================================================================
# DATA STRUCTURES
# ===========================================================================

struct ElemBins
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
    bins::Vector{Vector{Int}}
end

mutable struct CouplingData
    npoin_recv::Vector{Int32}
    npoin_send::Vector{Int32}
    recv_from_ranks::Vector{Int32}
    send_to_ranks::Vector{Int32}

    comm_world::MPI.Comm
    lrank::Int
    neqs::Int
    ndime::Int

    # --------------------------------------------------------------------------
    # send_coords: user-settable flag in setup_coupling_and_mesh.
    #   false (default) → send ndime velocity components per point each step.
    #   true            → send ndime coordinate values   per point each step.
    # Both modes produce a buffer of size [npoin × ndime], so Alya's receive
    # path is identical either way. Only the physical meaning differs.
    # --------------------------------------------------------------------------
    send_coords::Bool

    send_bufs::Union{Nothing, Vector{Vector{Float64}}}

    alya_local_coords::Union{Nothing, Matrix{Float64}}
    alya_local_ids::Union{Nothing, Vector{Int32}}
    alya_owner_ranks::Union{Nothing, Vector{Int32}}

    elem_bboxes::Union{Nothing, Vector{NTuple{4,Float64}}}
    interp_bins::Union{Nothing, ElemBins}
    elem_conn::Union{Nothing, Matrix{Int}}
    elem_x::Union{Nothing, Matrix{Float64}}
    elem_y::Union{Nothing, Matrix{Float64}}
    ξ_nodes_ref::Union{Nothing, Vector{Float64}}
    ω_bary::Union{Nothing, Vector{Float64}}

    qout::Union{Nothing, Matrix{Float64}}
    u_interp::Union{Nothing, Matrix{Float64}}

    ψξ_scratch::Union{Nothing, Vector{Float64}}
    ψη_scratch::Union{Nothing, Vector{Float64}}
    dψξ_scratch::Union{Nothing, Vector{Float64}}
    dψη_scratch::Union{Nothing, Vector{Float64}}
    α_scratch::Union{Nothing, Vector{Float64}}
    x_e_scratch::Union{Nothing, Vector{Float64}}
    y_e_scratch::Union{Nothing, Vector{Float64}}

    function CouplingData(; npoin_recv, npoin_send, recv_from_ranks, send_to_ranks,
                          comm_world, lrank, neqs, ndime, send_coords=false)
        new(npoin_recv, npoin_send, recv_from_ranks, send_to_ranks,
            comm_world, lrank, neqs, ndime, send_coords,
            nothing,
            nothing, nothing, nothing,
            nothing, nothing, nothing, nothing, nothing, nothing, nothing,
            nothing, nothing,
            nothing, nothing, nothing, nothing, nothing,
            nothing, nothing)
    end
end

struct Alya_Point_Ownership
    alya_coords::Matrix{Float64}
    owner_jrank::Vector{Int32}
    owner_wrank::Vector{Int32}
    alya_point_id::Vector{Int32}
    ndime::Int32
end

# ===========================================================================
# CACHE PRE-FETCH (called before with_mpi to move JLD2 JIT + I/O off the
# Alya-blocking critical path)
# ===========================================================================
#
# Alya blocks at MPI_Barrier+Alltoall waiting for Julia to answer the matching
# call inside setup_coupling_and_mesh.  Between je_receive_alya_data returning
# and that Alltoall, Julia must:
#   1. JIT-compile driver / setup_coupling_and_mesh / sem_setup (unavoidable)
#   2. Load the mesh-topology JLD2 cache       ← moved here
#   3. Load the SEM-preprocess JLD2 cache      ← moved here
#   4. Run extract_local_alya_coordinates      (geometry, not movable)
#
# By pre-loading (2) and (3) as plain top-level calls right after
# je_receive_alya_data, we:
#   • JIT-compile the JLD2 machinery once, before the Alya wait begins
#   • put the cache data in globals so _try_load_mesh_cache! and
#     _try_load_sem_cache just copy from memory rather than re-reading disk
#   • print timing diagnostics so any remaining delay is visible
#
# _mesh_cache_path / _preprocess_cache_path are defined later in mesh.jl /
# sem_setup.jl but are looked up at call time (same module scope), so the
# forward references are fine.
# ===========================================================================
function je_prefetch_caches!(inputs::Dict, nparts::Int)
    rank = MPI.Comm_rank(get_mpi_comm())

    # ── 1. Mesh topology cache ──────────────────────────────────────────────
    if JEXPRESSO_PREFETCHED_MESH_CACHE[] === nothing
        mesh_path = _mesh_cache_path(inputs, nparts)
        if isfile(mesh_path)
            rank == 0 && (print("[prefetch] mesh cache … "); flush(stdout))
            t0 = time_ns()
            try
                raw = JLD2.load(mesh_path)
                if haskey(raw, "mesh_fields")
                    JEXPRESSO_PREFETCHED_MESH_CACHE[] = raw
                    rank == 0 && @printf("%.2f s  (%s)\n", (time_ns()-t0)/1e9, mesh_path)
                else
                    rank == 0 && println("old-format cache ignored (will rebuild)")
                end
            catch e
                rank == 0 && @warn "[prefetch] mesh cache load failed" exception=e
            end
            flush(stdout)
        end
    end

    # ── 2. SEM preprocess cache ─────────────────────────────────────────────
    if JEXPRESSO_PREFETCHED_SEM_CACHE[] === nothing
        Nξ       = inputs[:nop]
        Qξ       = get(inputs, :lexact_integration, false) ? Nξ + 1 : Nξ
        sem_path = _preprocess_cache_path(inputs, Nξ, Qξ, nparts)
        if isfile(sem_path)
            rank == 0 && (print("[prefetch] SEM cache  … "); flush(stdout))
            t0 = time_ns()
            try
                d = JLD2.load(sem_path)
                if haskey(d, "metrics") && haskey(d, "matrix")
                    JEXPRESSO_PREFETCHED_SEM_CACHE[] = (d["metrics"], d["matrix"])
                    rank == 0 && @printf("%.2f s  (%s)\n", (time_ns()-t0)/1e9, sem_path)
                end
            catch e
                rank == 0 && @warn "[prefetch] SEM cache load failed" exception=e
            end
            flush(stdout)
        end
    end
end

# ===========================================================================
# MPI HANDSHAKE
# ===========================================================================

function je_perform_coupling_handshake(world, nparts)
    wsize = MPI.Comm_size(world)
    wrank = MPI.Comm_rank(world)
    wsize <= nparts && return false
    local_chars = Vector{UInt8}(rpad("JEXPRESSO", 128, ' '))
    MPI.Gather!(local_chars, nothing, 0, world)
    if wrank == nparts
        println("[Driver] Handshake complete - Jexpresso ready.")
        flush(stdout)
    end
    return true
end

function je_receive_alya_data(world, nparts)
    wsize = MPI.Comm_size(world)
    wsize <= nparts && (@warn "je_receive_alya_data called but not in coupled mode"; return)

    # Idempotency guard: if coupling data was already received (e.g. from run.jl
    # before driver() was called), skip the MPI operations entirely.
    JEXPRESSO_COUPLING_DATA[] !== nothing && return

    ndime_buf = Vector{Int32}(undef, 1)
    MPI.Bcast!(ndime_buf, 0, world)
    ndime = Int(ndime_buf[1])

    rem_min = Vector{Float64}(undef, 3)
    rem_max = Vector{Float64}(undef, 3)
    rem_nx  = Vector{Int32}(undef, 3)
    for idime in 1:3
        MPI.Bcast!(@view(rem_min[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_max[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_nx[idime:idime]),  0, world)
    end

    nranks_alya  = wsize - nparts
    alya2world_l = zeros(Int32, nranks_alya)
    alya2world   = MPI.Allreduce(alya2world_l, MPI.SUM, world)

    set_coupling_data(Dict{Symbol,Any}(
        :ndime      => ndime,
        :rem_min    => rem_min,
        :rem_max    => rem_max,
        :rem_nx     => rem_nx,
        :alya2world => alya2world,
    ))

    lrank = MPI.Comm_rank(get_mpi_comm())
    if lrank == 0
        println("[je_receive_alya_data] ndime=$ndime  min=$rem_min  max=$rem_max  nx=$rem_nx")
        flush(stdout)
    end
end

function je_send_node_list(alya_local_ids::Vector{Int32},
                           alya_owner_ranks::Vector{Int32},
                           send_to_ranks::Vector{Int32},
                           world::MPI.Comm)
    lrank = MPI.Comm_rank(get_mpi_comm())
    wrank = MPI.Comm_rank(world)
    send_requests = MPI.Request[]
    for dest_rank in send_to_ranks
        mask    = alya_owner_ranks .== dest_rank
        gid_buf = Int32.(alya_local_ids[mask])
        push!(send_requests, MPI.Isend(gid_buf, dest_rank, 0, world))
        println("[je_send_node_list] lrank=$lrank (wrank=$wrank) → Alya world rank $dest_rank: ",
                "$(length(gid_buf)) node IDs")
    end
    isempty(send_requests) || MPI.Waitall(send_requests)
    flush(stdout)
end

# ===========================================================================
# GEOMETRY / INTERPOLATION HELPERS
# ===========================================================================

function _make_conn_accessor(mesh)
    connijk = getfield(mesh, :connijk)
    nsd = mesh.nsd
    return nsd <= 2 ? (e -> vec(@view connijk[e, :, :, 1])) :
                      (e -> vec(@view connijk[e, :, :, :]))
end

_num_elems(mesh) = mesh.nelem

function barycentric_weights(nodes::AbstractVector{<:Real})
    n = length(nodes)
    w = ones(Float64, n)
    for j in 1:n
        xj = nodes[j]
        denom = 1.0
        @inbounds for k in 1:n
            k != j && (denom *= (xj - nodes[k]))
        end
        w[j] = 1.0 / denom
    end
    return w
end

function _build_elem_bins(elem_bboxes::Vector{NTuple{4,Float64}}; bins_per_dim::Int=64)
    ne = length(elem_bboxes)
    xmin = +Inf; xmax = -Inf; ymin = +Inf; ymax = -Inf
    for (x0,x1,y0,y1) in elem_bboxes
        xmin = min(xmin, x0); xmax = max(xmax, x1)
        ymin = min(ymin, y0); ymax = max(ymax, y1)
    end
    nx = max(1, bins_per_dim); ny = max(1, bins_per_dim)
    dx = (xmax > xmin) ? (xmax - xmin) / nx : 1.0
    dy = (ymax > ymin) ? (ymax - ymin) / ny : 1.0
    bin_lists = [Int[] for _ in 1:(nx*ny)]
    for e in 1:ne
        (x0,x1,y0,y1) = elem_bboxes[e]
        ix0 = clamp(Int(floor((x0 - xmin)/dx)), 0, nx-1)
        ix1 = clamp(Int(floor((x1 - xmin)/dx)), 0, nx-1)
        iy0 = clamp(Int(floor((y0 - ymin)/dy)), 0, ny-1)
        iy1 = clamp(Int(floor((y1 - ymin)/dy)), 0, ny-1)
        for iy in iy0:iy1, ix in ix0:ix1
            push!(bin_lists[iy*nx + ix + 1], e)
        end
    end
    fallback = collect(1:ne)
    for i in eachindex(bin_lists)
        isempty(bin_lists[i]) && (bin_lists[i] = fallback)
    end
    return ElemBins(xmin, xmax, ymin, ymax, nx, ny, dx, dy, bin_lists)
end

@inline function _bin_candidates(bins::ElemBins, x::Float64, y::Float64)
    ix = clamp(Int(floor((x - bins.xmin) / bins.dx)), 0, bins.nx - 1)
    iy = clamp(Int(floor((y - bins.ymin) / bins.dy)), 0, bins.ny - 1)
    return bins.bins[iy * bins.nx + ix + 1]
end

function evaluate_lagrange_1d!(ψ::Vector{Float64}, ξ::Float64,
                                ξ_nodes::Vector{Float64}, ω::Vector{Float64})
    n = length(ξ_nodes)
    @inbounds for i in 1:n
        if abs(ξ - ξ_nodes[i]) < 1e-14
            fill!(ψ, 0.0); ψ[i] = 1.0; return
        end
    end
    sum_val = 0.0
    @inbounds for i in 1:n
        ψ[i] = ω[i] / (ξ - ξ_nodes[i])
        sum_val += ψ[i]
    end
    inv_s = 1.0 / sum_val
    @inbounds for i in 1:n; ψ[i] *= inv_s; end
end

function evaluate_lagrange_1d_derivative!(dψ::Vector{Float64}, ξ::Float64,
                                          ξ_nodes::Vector{Float64}, ω::Vector{Float64},
                                          α::Vector{Float64})
    n = length(ξ_nodes)
    @inbounds for j in 1:n
        if abs(ξ - ξ_nodes[j]) < 1e-14
            fill!(dψ, 0.0)
            for k in 1:n
                k == j && continue
                dψ[k]  =  ω[k] / (ω[j] * (ξ_nodes[j] - ξ_nodes[k]))
                dψ[j] -= 1.0 / (ξ_nodes[j] - ξ_nodes[k])
            end
            return
        end
    end
    sum_α = 0.0
    @inbounds for i in 1:n
        α[i] = ω[i] / (ξ - ξ_nodes[i])
        sum_α += α[i]
    end
    sum_dα = 0.0
    @inbounds for i in 1:n; sum_dα -= α[i] / (ξ - ξ_nodes[i]); end
    inv_sum2 = 1.0 / (sum_α * sum_α)
    @inbounds for i in 1:n
        dψ[i] = (-α[i]/(ξ - ξ_nodes[i]) * sum_α - α[i] * sum_dα) * inv_sum2
    end
end

@inline function physical_to_reference(px::Float64, py::Float64,
                                        x_elem::AbstractVector, y_elem::AbstractVector,
                                        ξ_nodes::Vector{Float64}, ω::Vector{Float64},
                                        ngl::Int,
                                        ψξ::Vector{Float64}, ψη::Vector{Float64},
                                        dψξ::Vector{Float64}, dψη::Vector{Float64},
                                        α::Vector{Float64})
    ξ, η = 0.0, 0.0
    for _ in 1:20
        evaluate_lagrange_1d!(ψξ, ξ, ξ_nodes, ω)
        evaluate_lagrange_1d!(ψη, η, ξ_nodes, ω)
        evaluate_lagrange_1d_derivative!(dψξ, ξ, ξ_nodes, ω, α)
        evaluate_lagrange_1d_derivative!(dψη, η, ξ_nodes, ω, α)
        x_c = y_c = dxdξ = dxdη = dydξ = dydη = 0.0
        idx = 1
        @inbounds for j in 1:ngl, i in 1:ngl
            ψv  = ψξ[i]*ψη[j]; dξv = dψξ[i]*ψη[j]; dηv = ψξ[i]*dψη[j]
            x_c  += ψv *x_elem[idx]; y_c  += ψv *y_elem[idx]
            dxdξ += dξv*x_elem[idx]; dxdη += dηv*x_elem[idx]
            dydξ += dξv*y_elem[idx]; dydη += dηv*y_elem[idx]
            idx += 1
        end
        rx = px - x_c; ry = py - y_c
        sqrt(rx*rx + ry*ry) < 1e-12 && return ξ, η, true
        det_J = dxdξ*dydη - dxdη*dydξ
        abs(det_J) < 1e-15 && return ξ, η, false
        inv_d = 1.0 / det_J
        ξ += inv_d*( dydη*rx - dxdη*ry)
        η += inv_d*(-dydξ*rx + dxdξ*ry)
        (abs(ξ) > 10.0 || abs(η) > 10.0) && return ξ, η, false
    end
    return ξ, η, false
end

function interpolate_solution_to_alya_coords!(u_interp::Matrix{Float64},
                                              alya_coords::Matrix{Float64},
                                              u_mat::Matrix{Float64},
                                              ξ_nodes::Vector{Float64},
                                              ω::Vector{Float64},
                                              neqs::Int,
                                              elem_bboxes::Vector{NTuple{4,Float64}},
                                              bins::ElemBins,
                                              elem_conn::Matrix{Int},
                                              elem_x::Matrix{Float64},
                                              elem_y::Matrix{Float64},
                                              ψξ::Vector{Float64}, ψη::Vector{Float64},
                                              dψξ::Vector{Float64}, dψη::Vector{Float64},
                                              α::Vector{Float64},
                                              x_e::Vector{Float64}, y_e::Vector{Float64},
                                              mesh_x::Vector{Float64},
                                              mesh_y::Vector{Float64})
    n_points = size(alya_coords, 1)
    ngl      = length(ξ_nodes)
    ngl2     = ngl * ngl
    npoin    = size(u_mat, 1)

    @inbounds for ipt in 1:n_points
        px, py = alya_coords[ipt, 1], alya_coords[ipt, 2]
        candidates = _bin_candidates(bins, px, py)
        found = false
        for e in candidates
            bb = elem_bboxes[e]
            (px < bb[1]-1e-10 || px > bb[2]+1e-10 ||
             py < bb[3]-1e-10 || py > bb[4]+1e-10) && continue
            @inbounds for k in 1:ngl2
                x_e[k] = elem_x[e, k]; y_e[k] = elem_y[e, k]
            end
            ξ_ref, η_ref, converged = physical_to_reference(
                px, py, x_e, y_e, ξ_nodes, ω, ngl, ψξ, ψη, dψξ, dψη, α)
            (!converged || abs(ξ_ref) > 1.0+1e-10 || abs(η_ref) > 1.0+1e-10) && continue
            evaluate_lagrange_1d!(ψξ, ξ_ref, ξ_nodes, ω)
            evaluate_lagrange_1d!(ψη, η_ref, ξ_nodes, ω)
            for q in 1:neqs
                val = 0.0; idx = 1
                for j in 1:ngl, i in 1:ngl
                    val += ψξ[i] * ψη[j] * u_mat[elem_conn[e, idx], q]
                    idx += 1
                end
                u_interp[ipt, q] = val
            end
            found = true; break
        end
        if !found
            nearest = 1
            min_d2  = (mesh_x[1]-px)^2 + (mesh_y[1]-py)^2
            @inbounds for ip in 2:npoin
                d2 = (mesh_x[ip]-px)^2 + (mesh_y[ip]-py)^2
                if d2 < min_d2; min_d2 = d2; nearest = ip; end
            end
            @inbounds for q in 1:neqs; u_interp[ipt, q] = u_mat[nearest, q]; end
        end
    end
end

# ===========================================================================
# ALYA COORDINATE EXTRACTION
# ===========================================================================

function extract_local_alya_coordinates(mesh, coupling_data, local_comm, world_comm;
                                        block_size::NTuple{3,Int}=(64,64,64),
                                        use_cropping::Bool=true,
                                        ξ_nodes::Union{Nothing,Vector{Float64}}=nothing)
    ndime      = coupling_data[:ndime]
    alya2world = coupling_data[:alya2world]
    @assert ndime == 2 || ndime == 3 "Only ndime==2 or ndime==3 supported"
    @assert length(alya2world) > 0   "alya2world must be non-empty"

    rem_min = Float64.(coupling_data[:rem_min])
    rem_max = Float64.(coupling_data[:rem_max])
    rem_nx  = Int.(coupling_data[:rem_nx])

    rem_dx = zeros(Float64, 3)
    for d in 1:ndime
        rem_dx[d] = rem_nx[d] > 1 ? (rem_max[d]-rem_min[d]) / (rem_nx[d]-1) : 0.0
    end

    nmax        = rem_nx[1] * rem_nx[2] * rem_nx[3]
    nranks_alya = length(alya2world)
    alya_worker_indices = [k for k in 1:nranks_alya if alya2world[k] != Int32(0)]
    nworkers_alya = length(alya_worker_indices)

    # Guard: if there are no Alya worker ranks (e.g. asize==1 so only the
    # master rank exists, which has world rank 0 and is filtered out above),
    # return empty arrays — no points can be assigned to Alya workers.
    if nworkers_alya == 0
        @warn "extract_local_alya_coordinates: no Alya worker ranks found " *
              "(alya2world=$(alya2world)). " *
              "This happens when ALYA_PROCS==1 (only the master rank exists). " *
              "Returning empty coordinate arrays."
        empty_coords  = zeros(Float64, 0, ndime)
        empty_ids     = Int32[]
        empty_owners  = Int32[]
        return empty_coords, empty_ids, empty_owners
    end

    r_w  = mod(nmax, nworkers_alya)
    np_w = div(nmax, nworkers_alya)

    xmin_local = minimum(mesh.x); xmax_local = maximum(mesh.x)
    ymin_local = minimum(mesh.y); ymax_local = maximum(mesh.y)
    zmin_local = ndime == 3 ? minimum(mesh.z) : 0.0
    zmax_local = ndime == 3 ? maximum(mesh.z) : 0.0
    tol = 1e-10

    lrank = MPI.Comm_rank(local_comm)
    wrank = MPI.Comm_rank(world_comm)

    use_elem_containment = (ξ_nodes !== nothing && ndime == 2)
    if use_elem_containment
        ω_bary    = barycentric_weights(ξ_nodes)
        nelem_loc = _num_elems(mesh)
        ngl_loc   = mesh.ngl
        ngl2_loc  = ngl_loc * ngl_loc
        gc_loc    = _make_conn_accessor(mesh)

        local_elem_bboxes = Vector{NTuple{4,Float64}}(undef, nelem_loc)
        @inbounds for e in 1:nelem_loc
            ns = gc_loc(e)
            local_elem_bboxes[e] = (minimum(mesh.x[ns]), maximum(mesh.x[ns]),
                                    minimum(mesh.y[ns]), maximum(mesh.y[ns]))
        end
        local_elem_bins = _build_elem_bins(local_elem_bboxes; bins_per_dim=64)

        local_ex = Matrix{Float64}(undef, nelem_loc, ngl2_loc)
        local_ey = Matrix{Float64}(undef, nelem_loc, ngl2_loc)
        @inbounds for e in 1:nelem_loc
            ns = gc_loc(e)
            for k in 1:ngl2_loc
                local_ex[e, k] = mesh.x[ns[k]]
                local_ey[e, k] = mesh.y[ns[k]]
            end
        end

        ψξ_s  = Vector{Float64}(undef, ngl_loc)
        ψη_s  = Vector{Float64}(undef, ngl_loc)
        dψξ_s = Vector{Float64}(undef, ngl_loc)
        dψη_s = Vector{Float64}(undef, ngl_loc)
        α_s   = Vector{Float64}(undef, ngl_loc)

        @inline function in_any_local_elem(px::Float64, py::Float64)
            for e in _bin_candidates(local_elem_bins, px, py)
                bb = local_elem_bboxes[e]
                (px < bb[1]-1e-10 || px > bb[2]+1e-10 ||
                 py < bb[3]-1e-10 || py > bb[4]+1e-10) && continue
                ξr, ηr, conv = physical_to_reference(
                    px, py, @view(local_ex[e,:]), @view(local_ey[e,:]),
                    ξ_nodes, ω_bary, ngl_loc, ψξ_s, ψη_s, dψξ_s, dψη_s, α_s)
                (conv && abs(ξr) <= 1.0+1e-10 && abs(ηr) <= 1.0+1e-10) && return true
            end
            return false
        end
    end

    @inline idx_to_xyz(i1,i2,i3) = (rem_min[1]+i1*rem_dx[1],
                                     rem_min[2]+i2*rem_dx[2],
                                     rem_min[3]+i3*rem_dx[3])
    @inline idx_to_ipoin(i1,i2,i3) = i1 + rem_nx[1]*(i2 + rem_nx[2]*i3) + 1

    @inline function inside_local(x,y,z)
        ndime == 2 ? (x >= xmin_local-tol && x <= xmax_local+tol &&
                      y >= ymin_local-tol && y <= ymax_local+tol) :
                     (x >= xmin_local-tol && x <= xmax_local+tol &&
                      y >= ymin_local-tol && y <= ymax_local+tol &&
                      z >= zmin_local-tol && z <= zmax_local+tol)
    end

    @inline function owner_alya_rank(ipoin::Int)
        iworker = ipoin <= r_w*(np_w+1) ? div(ipoin-1, np_w+1) :
                                          r_w + div(ipoin - r_w*(np_w+1) - 1, np_w)
        return alya_worker_indices[iworker + 1]
    end

    @inline function crop_1d(minA, dx, nx, minB, maxB)
        nx <= 1 && return 0, 0
        ilo = max(0,    Int(clamp(floor((minB-minA)/dx + 1e-12), 0, nx-1)) - 1)
        ihi = min(nx-1, Int(clamp( ceil((maxB-minA)/dx - 1e-12), 0, nx-1)) + 1)
        return ilo, ihi
    end

    i1_lo, i1_hi = 0, rem_nx[1]-1
    i2_lo, i2_hi = 0, rem_nx[2]-1
    i3_lo, i3_hi = 0, ndime == 3 ? rem_nx[3]-1 : 0

    if use_cropping
        rem_nx[1] > 1 && ((i1_lo, i1_hi) = crop_1d(rem_min[1], rem_dx[1], rem_nx[1], xmin_local-tol, xmax_local+tol))
        rem_nx[2] > 1 && ((i2_lo, i2_hi) = crop_1d(rem_min[2], rem_dx[2], rem_nx[2], ymin_local-tol, ymax_local+tol))
        ndime == 3 && rem_nx[3] > 1 &&
            ((i3_lo, i3_hi) = crop_1d(rem_min[3], rem_dx[3], rem_nx[3], zmin_local-tol, zmax_local+tol))
        if i1_lo > i1_hi || i2_lo > i2_hi || i3_lo > i3_hi
            println("[extract] (lrank=$lrank) no overlap with local domain.")
            flush(stdout)
            return zeros(Float64, 0, ndime), Int32[], Int32[]
        end
    end

    Bx, By, Bz = block_size
    Bx = max(1,Bx); By = max(1,By); Bz = max(1, ndime==3 ? Bz : 1)
    est = (i1_hi-i1_lo+1)*(i2_hi-i2_lo+1)*(i3_hi-i3_lo+1)
    X = Float64[]; Y = Float64[]; Z = Float64[]
    ids = Int32[]; owners = Int32[]
    sizehint!(X, est); sizehint!(Y, est)
    ndime == 3 && sizehint!(Z, est)
    sizehint!(ids, est); sizehint!(owners, est)

    @inline bin_bbox(i1s,i1e,i2s,i2e,i3s,i3e) =
        (rem_min[1]+i1s*rem_dx[1], rem_min[1]+i1e*rem_dx[1],
         rem_min[2]+i2s*rem_dx[2], rem_min[2]+i2e*rem_dx[2],
         rem_min[3]+i3s*rem_dx[3], rem_min[3]+i3e*rem_dx[3])

    @inline function outside_local(xmn,xmx,ymn,ymx,zmn,zmx)
        ndime==2 ? (xmx<xmin_local-tol || xmn>xmax_local+tol ||
                    ymx<ymin_local-tol || ymn>ymax_local+tol) :
                   (xmx<xmin_local-tol || xmn>xmax_local+tol ||
                    ymx<ymin_local-tol || ymn>ymax_local+tol ||
                    zmx<zmin_local-tol || zmn>zmax_local+tol)
    end

    @inline function inside_box(xmn,xmx,ymn,ymx,zmn,zmx)
        ndime==2 ? (xmn>=xmin_local-tol && xmx<=xmax_local+tol &&
                    ymn>=ymin_local-tol && ymx<=ymax_local+tol) :
                   (xmn>=xmin_local-tol && xmx<=xmax_local+tol &&
                    ymn>=ymin_local-tol && ymx<=ymax_local+tol &&
                    zmn>=zmin_local-tol && zmx<=zmax_local+tol)
    end

    for i3s in i3_lo:Bz:i3_hi
        i3e = min(i3s+Bz-1, i3_hi)
        for i2s in i2_lo:By:i2_hi
            i2e = min(i2s+By-1, i2_hi)
            for i1s in i1_lo:Bx:i1_hi
                i1e = min(i1s+Bx-1, i1_hi)
                bxmn,bxmx,bymn,bymx,bzmn,bzmx = bin_bbox(i1s,i1e,i2s,i2e,i3s,i3e)
                outside_local(bxmn,bxmx,bymn,bymx,bzmn,bzmx) && continue
                if !use_elem_containment && inside_box(bxmn,bxmx,bymn,bymx,bzmn,bzmx)
                    @inbounds for i3 in i3s:i3e, i2 in i2s:i2e, i1 in i1s:i1e
                        x,y,z = idx_to_xyz(i1,i2,i3)
                        ipoin = idx_to_ipoin(i1,i2,i3)
                        push!(X,x); push!(Y,y); ndime==3 && push!(Z,z)
                        push!(ids, Int32(ipoin))
                        push!(owners, Int32(alya2world[owner_alya_rank(ipoin)]))
                    end
                else
                    @inbounds for i3 in i3s:i3e, i2 in i2s:i2e, i1 in i1s:i1e
                        x,y,z = idx_to_xyz(i1,i2,i3)
                        accept = use_elem_containment ? in_any_local_elem(x, y) :
                                                        inside_local(x,y,z)
                        accept || continue
                        ipoin = idx_to_ipoin(i1,i2,i3)
                        push!(X,x); push!(Y,y); ndime==3 && push!(Z,z)
                        push!(ids, Int32(ipoin))
                        push!(owners, Int32(alya2world[owner_alya_rank(ipoin)]))
                    end
                end
            end
        end
    end

    n_local = length(ids)

    local_claim = fill(Int32(-1), nmax)
    @inbounds for i in 1:n_local; local_claim[Int(ids[i])] = Int32(lrank); end
    global_claim = MPI.Allreduce(local_claim, MPI.MAX, local_comm)

    if n_local > 0
        keep = [global_claim[Int(ids[i])] == Int32(lrank) for i in 1:n_local]
        if !all(keep)
            X = X[keep]; Y = Y[keep]
            ndime == 3 && (Z = Z[keep])
            ids = ids[keep]; owners = owners[keep]
            n_local = length(ids)
            println("[extract] (lrank=$lrank) after phase-1 dedup: $n_local points kept.")
            flush(stdout)
        end
    end

    rescue_claim = fill(Int32(typemax(Int32)), nmax)
    for gid in 1:nmax
        global_claim[gid] != -1 && continue
        i1 = (gid-1) % rem_nx[1]
        i2 = div(gid-1, rem_nx[1]) % rem_nx[2]
        i3 = div(gid-1, rem_nx[1]*rem_nx[2])
        x,y,z = idx_to_xyz(i1,i2,i3)
        inside_local(x,y,z) || continue
        rescue_claim[gid] = Int32(lrank)
    end
    global_rescue = MPI.Allreduce(rescue_claim, MPI.MIN, local_comm)

    n_rescued = 0
    for gid in 1:nmax
        global_claim[gid] != -1             && continue
        global_rescue[gid] == typemax(Int32) && continue
        global_rescue[gid] != Int32(lrank)   && continue
        i1 = (gid-1) % rem_nx[1]
        i2 = div(gid-1, rem_nx[1]) % rem_nx[2]
        i3 = div(gid-1, rem_nx[1]*rem_nx[2])
        x,y,z = idx_to_xyz(i1,i2,i3)
        push!(X,x); push!(Y,y); ndime==3 && push!(Z,z)
        push!(ids, Int32(gid))
        push!(owners, Int32(alya2world[owner_alya_rank(gid)]))
        n_rescued += 1
    end
    n_local = length(ids)
    n_rescued > 0 && println("[extract] (lrank=$lrank) rescued $n_rescued unclaimed points.")

    alya_local_coords = zeros(Float64, n_local, ndime)
    @inbounds for i in 1:n_local
        alya_local_coords[i,1] = X[i]
        alya_local_coords[i,2] = Y[i]
        ndime == 3 && (alya_local_coords[i,3] = Z[i])
    end

    if n_local > 1
        perm  = sortperm(collect(zip(owners, ids)))
        alya_local_coords = alya_local_coords[perm, :]
        ids    = ids[perm]
        owners = owners[perm]
    end

    println("[extract] (lrank=$lrank, wrank=$wrank) collected $n_local Alya points.")
    flush(stdout)
    return alya_local_coords, ids, owners
end

# ===========================================================================
# DIAGNOSTICS
# ===========================================================================

function build_alya_point_ownership_map(mesh, coupling_data, local_comm, world_comm)
    ndime      = coupling_data[:ndime]
    rem_min    = coupling_data[:rem_min]
    rem_max    = coupling_data[:rem_max]
    rem_nx     = coupling_data[:rem_nx]
    alya2world = coupling_data[:alya2world]
    wrank = MPI.Comm_rank(world_comm)
    lrank = MPI.Comm_rank(local_comm)
    rem_dx = zeros(Float64, 3)
    for idim in 1:ndime
        rem_dx[idim] = rem_nx[idim] > 1 ?
            (rem_max[idim]-rem_min[idim])/(rem_nx[idim]-1) : 0.0
    end
    nmax   = rem_nx[1]*rem_nx[2]*rem_nx[3]
    xmin_l = minimum(mesh.x); xmax_l = maximum(mesh.x)
    ymin_l = minimum(mesh.y); ymax_l = maximum(mesh.y)
    zmin_l = ndime==3 ? minimum(mesh.z) : 0.0
    zmax_l = ndime==3 ? maximum(mesh.z) : 0.0
    tol    = 1e-10
    local_coords      = zeros(Float64, nmax, ndime)
    local_owner_jrank = fill(Int32(-1), nmax)
    local_owner_wrank = fill(Int32(-1), nmax)
    for ipoin in 1:nmax
        i0  = ipoin - 1
        ri3 = div(i0, rem_nx[1]*rem_nx[2])
        ri2 = div(i0 - ri3*rem_nx[1]*rem_nx[2], rem_nx[1])
        ri1 = mod(i0 - ri3*rem_nx[1]*rem_nx[2] - ri2*rem_nx[1], rem_nx[1])
        x = [rem_min[1]+ri1*rem_dx[1], rem_min[2]+ri2*rem_dx[2], rem_min[3]+ri3*rem_dx[3]]
        for d in 1:ndime; local_coords[ipoin,d] = x[d]; end
        in_r = ndime==2 ? (x[1]>=xmin_l-tol && x[1]<=xmax_l+tol &&
                           x[2]>=ymin_l-tol && x[2]<=ymax_l+tol) :
                          (x[1]>=xmin_l-tol && x[1]<=xmax_l+tol &&
                           x[2]>=ymin_l-tol && x[2]<=ymax_l+tol &&
                           x[3]>=zmin_l-tol && x[3]<=zmax_l+tol)
        if in_r
            local_owner_jrank[ipoin] = lrank
            local_owner_wrank[ipoin] = wrank
        end
    end
    global_owner_jrank = MPI.Allreduce(local_owner_jrank, MPI.MAX, local_comm)
    global_owner_wrank = MPI.Allreduce(local_owner_wrank, MPI.MAX, local_comm)
    return Alya_Point_Ownership(local_coords, global_owner_jrank, global_owner_wrank,
                                collect(Int32(1):Int32(nmax)), Int32(ndime))
end

function print_alya_point_ownership(ownership, coupling_data; max_points_to_print=20, only_owned=false)
    MPI.Comm_rank(get_mpi_comm()) != 0 && return
    npoints = length(ownership.alya_point_id)
    ndime   = ownership.ndime
    n_owned = count(x -> x >= 0, ownership.owner_jrank)
    println("\n" * "="^70); println("Alya Point Ownership Map"); println("="^70)
    println("Total Alya points: $npoints  (owned: $n_owned, unowned: $(npoints-n_owned))")
    println("-"^70)
    n_to_print = max_points_to_print < 0 ? npoints : max_points_to_print
    printed = 0
    for i in 1:npoints
        only_owned && ownership.owner_jrank[i] < 0 && continue
        if printed >= n_to_print
            remaining = only_owned ? count(x->x>=0, ownership.owner_jrank[i+1:end]) : npoints-i
            remaining > 0 && println("... ($remaining more not shown)")
            break
        end
        coord_str = ndime==2 ?
            @sprintf("(%10.2f, %10.2f)", ownership.alya_coords[i,1], ownership.alya_coords[i,2]) :
            @sprintf("(%10.2f, %10.2f, %10.2f)", ownership.alya_coords[i,1],
                     ownership.alya_coords[i,2], ownership.alya_coords[i,3])
        jr = ownership.owner_jrank[i]; wr = ownership.owner_wrank[i]
        owner_str = jr >= 0 ? @sprintf("Jexpresso rank %d (world rank %d)", jr, wr) : "NOT OWNED"
        println(@sprintf("%8d | %s | %s", ownership.alya_point_id[i], coord_str, owner_str))
        printed += 1
    end
    println("="^70); println(); flush(stdout)
end

function verify_alya_point_distribution(ownership, coupling_data)
    MPI.Comm_rank(get_mpi_comm()) != 0 && return
    lsize   = MPI.Comm_size(get_mpi_comm())
    npoints = length(ownership.alya_point_id)
    println("\n" * "="^70); println("Alya Point Distribution Statistics"); println("="^70)
    pts_per_rank = zeros(Int, lsize)
    for i in 1:npoints
        jr = ownership.owner_jrank[i]
        jr >= 0 && (pts_per_rank[jr+1] += 1)
    end
    for r in 0:lsize-1
        println(@sprintf("  Jexpresso rank %d: %d points (%.1f%%)",
                r, pts_per_rank[r+1], 100.0*pts_per_rank[r+1]/npoints))
    end
    println("-"^70)
    n_unowned = count(x->x<0, ownership.owner_jrank)
    n_unowned > 0 ? println("WARNING: $n_unowned Alya points are not owned by any Jexpresso rank!") :
                    println("✓ All Alya points are owned by Jexpresso ranks")
    println("="^70); println(); flush(stdout)
end

# ===========================================================================
# COUPLING SETUP ORCHESTRATOR
# ===========================================================================
#    setup_coupling_and_mesh(world, lsize, inputs, nranks, distribute, rank,
#                            OUTPUT_DIR, TFloat; send_coords=false)
#
# Set up the coupled simulation.
#
# Keyword argument
# - `send_coords=false` *(default)*: Julia sends **velocity** components
#  `[u, v, ...]` (columns `2:neqs-1` of the solution) to Alya each step.
# - `send_coords=true`: Julia sends the **coordinates** `[x, y, ...]`
#  of the Alya grid points it computed (useful for verifying point location).
#
# Both modes produce a buffer of shape `[npoin × ndime]` with the same
#field-fastest layout, so Alya's receive path is identical either way.
# ===========================================================================
function setup_coupling_and_mesh(world, lsize, inputs, nranks, distribute, rank,
                                 OUTPUT_DIR, TFloat;
                                 send_coords::Bool = SEND_COORDS)
    
    if rank == 1
        println("[setup_coupling] JIT ... — starting mesh setup")
    end
    je_receive_alya_data(world, lsize)
    coupling_data = get_coupling_data()
    local_comm    = get_mpi_comm()
    lrank         = MPI.Comm_rank(local_comm)
    wsize         = MPI.Comm_size(world)

    # Sign-of-life print: appears as soon as setup_coupling_and_mesh has been
    # JIT-compiled.  This breaks up the otherwise-silent period between
    # je_receive_alya_data and the "Read gmsh grid" message from sem_setup.
    if lrank == 0
        println("[setup_coupling] JIT done — starting mesh setup")
        flush(stdout)
    end

    t_sem = @elapsed begin
        sem, partitioned_model = sem_setup(inputs, nranks, distribute, rank)
        if inputs[:backend] != CPU()
            convert_mesh_arrays!(sem.mesh.SD, sem.mesh, inputs[:backend], inputs)
        end
    end
    lrank == 0 && (@printf("[setup_coupling] sem_setup done in %.2f s\n", t_sem); flush(stdout))

    t_geom = @elapsed alya_local_coords, alya_local_ids, alya_owner_ranks =
        extract_local_alya_coordinates(sem.mesh, coupling_data, local_comm, world;
                                       block_size=(64,64,64), use_cropping=true,
                                       ξ_nodes=Vector{Float64}(sem.ξ))
    lrank == 0 && (@printf("[setup_coupling] extract_alya_coords done in %.2f s\n", t_geom); flush(stdout))

    npoin_send = zeros(Int32, wsize)
    for wr in alya_owner_ranks; npoin_send[wr+1] += 1; end
    send_to_ranks = Int32[i-1 for i in 1:wsize if npoin_send[i] > 0]

    lrank == 0 && (println("[setup_coupling] Alltoall to exchange point counts…"); flush(stdout))

    recv_counts_from_alya = zeros(Int32, wsize)
    MPI.Barrier(world)
    MPI.Alltoall!(Vector{Int32}(npoin_send), recv_counts_from_alya, 1, world)

    je_send_node_list(Int32.(alya_local_ids), alya_owner_ranks, send_to_ranks, world)

    verify_coupling_communication_pattern(npoin_send, alya_owner_ranks, alya_local_ids,
                                          coupling_data, local_comm, world)

    if lrank == 0
        println("[setup_coupling] Send pattern: ",
                [(i-1, npoin_send[i]) for i in 1:wsize if npoin_send[i]>0])
        println("[setup_coupling] send_coords=$send_coords  ",
                send_coords ? "→ sending COORDINATES each step" :
                              "→ sending VELOCITY each step")
        flush(stdout)
    end

    ownership = build_alya_point_ownership_map(sem.mesh, coupling_data, local_comm, world)
    print_alya_point_ownership(ownership, coupling_data; max_points_to_print=20, only_owned=true)
    verify_alya_point_distribution(ownership, coupling_data)

    qp    = initialize(sem.mesh.SD, 0, sem.mesh, inputs, OUTPUT_DIR, TFloat)
    ndime = coupling_data[:ndime]

    # In velocity mode nfields = neqs-2 must equal ndime (Alya's assumption).
    # In coords mode nfields is always ndime, so no additional check is needed.
    if !send_coords
        @assert qp.neqs - 2 == ndime """
            nfields mismatch: Julia sends $(qp.neqs-2) velocity fields (neqs-2)
            but Alya expects ndime=$ndime.
        """
    end

    coupling = CouplingData(
        npoin_recv      = zeros(Int32, wsize),
        npoin_send      = npoin_send,
        recv_from_ranks = Int32[],
        send_to_ranks   = send_to_ranks,
        comm_world      = world,
        lrank           = lrank,
        neqs            = qp.neqs,
        ndime           = ndime,
        send_coords     = send_coords,
    )

    coupling.send_bufs = [zeros(Float64, npoin_send[i] * ndime) for i in 1:wsize]

    coupling.alya_local_coords = alya_local_coords
    coupling.alya_local_ids    = alya_local_ids
    coupling.alya_owner_ranks  = alya_owner_ranks

    let mesh  = sem.mesh,
        nelem = _num_elems(sem.mesh),
        ngl   = sem.mesh.ngl,
        ngl2  = sem.mesh.ngl * sem.mesh.ngl,
        gc    = _make_conn_accessor(sem.mesh)

        bb = Vector{NTuple{4,Float64}}(undef, nelem)
        @inbounds for e in 1:nelem
            ns = gc(e)
            bb[e] = (minimum(mesh.x[ns]), maximum(mesh.x[ns]),
                     minimum(mesh.y[ns]), maximum(mesh.y[ns]))
        end
        coupling.elem_bboxes = bb
        coupling.interp_bins = _build_elem_bins(bb; bins_per_dim=64)

        conn_mat = Matrix{Int}(undef, nelem, ngl2)
        ex_mat   = Matrix{Float64}(undef, nelem, ngl2)
        ey_mat   = Matrix{Float64}(undef, nelem, ngl2)
        @inbounds for e in 1:nelem
            ns = gc(e)
            for k in 1:ngl2
                conn_mat[e,k] = ns[k]
                ex_mat[e,k]   = mesh.x[ns[k]]
                ey_mat[e,k]   = mesh.y[ns[k]]
            end
        end
        coupling.elem_conn = conn_mat
        coupling.elem_x    = ex_mat
        coupling.elem_y    = ey_mat

        ξ_nodes_v        = Vector{Float64}(sem.ξ)
        coupling.ξ_nodes_ref = ξ_nodes_v
        coupling.ω_bary      = barycentric_weights(ξ_nodes_v)

        coupling.ψξ_scratch  = Vector{Float64}(undef, ngl)
        coupling.ψη_scratch  = Vector{Float64}(undef, ngl)
        coupling.dψξ_scratch = Vector{Float64}(undef, ngl)
        coupling.dψη_scratch = Vector{Float64}(undef, ngl)
        coupling.α_scratch   = Vector{Float64}(undef, ngl)
        coupling.x_e_scratch = Vector{Float64}(undef, ngl2)
        coupling.y_e_scratch = Vector{Float64}(undef, ngl2)

        coupling.qout     = zeros(Float64, mesh.npoin,             coupling.neqs)
        coupling.u_interp = zeros(Float64, length(alya_local_ids), coupling.neqs)
    end

    return coupling, sem, partitioned_model, qp
end

# ===========================================================================
# TIME-LOOP DATA EXCHANGE
# ===========================================================================

# Pack Alya coordinates [x, y, ...] into send buffers — coord-fastest layout,
# identical to the field-fastest velocity layout Alya already expects.
function pack_coord_data!(cpg::CouplingData,
                          alya_coords::Matrix{Float64},
                          alya_owner_ranks::Vector{Int32})
    n_local = size(alya_coords, 1)
    ndime   = cpg.ndime
    wsize   = length(cpg.npoin_send)
    for i in 1:wsize; fill!(cpg.send_bufs[i], 0.0); end
    send_offsets = zeros(Int, wsize)
    @inbounds for i in 1:n_local
        bidx = alya_owner_ranks[i] + 1
        cpg.npoin_send[bidx] > 0 || continue
        off = send_offsets[bidx]
        for d in 1:ndime
            cpg.send_bufs[bidx][off+d] = alya_coords[i, d]
        end
        send_offsets[bidx] += ndime
    end
end

# Pack interpolated velocity components [u, v, ...] into send buffers.
function pack_velocity_data!(cpg::CouplingData,
                              interp_values::Matrix{Float64},
                              alya_owner_ranks::Vector{Int32})
    n_local = size(interp_values, 1)
    nfields = size(interp_values, 2)
    wsize   = length(cpg.npoin_send)
    for i in 1:wsize; fill!(cpg.send_bufs[i], 0.0); end
    send_offsets = zeros(Int, wsize)
    @inbounds for i in 1:n_local
        bidx = alya_owner_ranks[i] + 1
        cpg.npoin_send[bidx] > 0 || continue
        off = send_offsets[bidx]
        for q in 1:nfields
            cpg.send_bufs[bidx][off+q] = interp_values[i, q]
        end
        send_offsets[bidx] += nfields
    end
end

function coupling_exchange_data!(cpg::CouplingData)
    send_requests = MPI.Request[]
    for dest_rank in cpg.send_to_ranks
        cpg.npoin_send[dest_rank+1] > 0 || continue
        push!(send_requests, MPI.Isend(cpg.send_bufs[dest_rank+1], dest_rank, 0, cpg.comm_world))
    end
    isempty(send_requests) || MPI.Waitall(send_requests)
end

# Full per-step exchange.  Branches on cpg.send_coords (set once at startup):
#   false -> interpolate solution, pack velocity columns 2:neqs-1, send.
#   true  -> pack Alya coordinates directly, send. No interpolation needed.
function je_perform_coupling_exchange(u, u_mat, t, cpg::CouplingData,
                                      qout::Matrix{Float64},
                                      u_interp::Matrix{Float64},
                                      ξ_nodes::Vector{Float64},
                                      ω::Vector{Float64},
                                      e_conn::Matrix{Int},
                                      elem_x::Matrix{Float64},
                                      elem_y::Matrix{Float64},
                                      ψξ::Vector{Float64},  ψη::Vector{Float64},
                                      dψξ::Vector{Float64}, dψη::Vector{Float64},
                                      α::Vector{Float64},
                                      x_e::Vector{Float64}, y_e::Vector{Float64},
                                      alya_coords::Matrix{Float64},
                                      owner_ranks::Vector{Int32},
                                      mesh_x::Vector{Float64},
                                      mesh_y::Vector{Float64},
                                      inputs, neqs::Int,
                                      elem_bboxes::Vector{NTuple{4,Float64}},
                                      bins::ElemBins)
    if cpg.send_coords
        pack_coord_data!(cpg, alya_coords, owner_ranks)
    else
        npoin = size(qout, 1)
        u2uaux!(u_mat, u, neqs, npoin)
        call_user_uout(qout, u_mat, u_mat, 0, inputs[:SOL_VARS_TYPE], npoin, neqs, neqs)
        interpolate_solution_to_alya_coords!(
            u_interp, alya_coords, qout,
            ξ_nodes, ω, neqs,
            elem_bboxes, bins,
            e_conn, elem_x, elem_y,
            ψξ, ψη, dψξ, dψη, α, x_e, y_e,
            mesh_x, mesh_y)
        pack_velocity_data!(cpg, u_interp[:, 2:neqs-1], owner_ranks)
    end
    coupling_exchange_data!(cpg)
end

# ===========================================================================
# SETUP-TIME COMMUNICATION VERIFICATION
# ===========================================================================

function verify_coupling_communication_pattern(npoin_send, alya_owner_ranks, alya_local_ids,
                                               coupling_data, local_comm, world_comm)
    lrank       = MPI.Comm_rank(local_comm)
    wrank       = MPI.Comm_rank(world_comm)
    wsize       = MPI.Comm_size(world_comm)
    nranks_alya = length(coupling_data[:alya2world])
    total_local = length(alya_owner_ranks)

    println("="^80)
    println("[VERIFY] Julia lrank=$lrank, wrank=$wrank")
    println("="^80)

    println("\n[CHECK 1] Verifying npoin_send against alya_owner_ranks...")
    local_counts = zeros(Int32, wsize)
    for wr in alya_owner_ranks; local_counts[wr+1] += 1; end
    mismatch = false
    for i in 1:wsize
        if local_counts[i] != npoin_send[i]
            println("  MISMATCH at world rank $(i-1): extracted=$(local_counts[i]) vs npoin_send=$(npoin_send[i])")
            mismatch = true
        end
    end
    mismatch ? (@warn "npoin_send does NOT match alya_owner_ranks!") :
               println("  ✓ npoin_send matches alya_owner_ranks distribution")

    println("\n[CHECK 2] Global point conservation...")
    rem_nx       = coupling_data[:rem_nx]
    total_alya   = rem_nx[1]*rem_nx[2]*rem_nx[3]
    global_owned = MPI.Allreduce(total_local, MPI.SUM, local_comm)
    println("  Total Alya grid points: $total_alya")
    println("  Owned by all Julia ranks: $global_owned")
    global_owned == total_alya ? println("  ✓ All Alya points accounted for") :
                                 (@warn "Point conservation failed! Expected $total_alya, got $global_owned")

    println("\n[CHECK 3] Send pattern (npoin_send):")
    for i in 1:wsize
        npoin_send[i] == 0 && continue
        wr  = i-1
        tag = wr < nranks_alya ? "Alya rank $wr" : "Julia rank $(wr-nranks_alya)"
        println("    $(npoin_send[i]) pts TO world rank $wr ($tag)")
    end

    println("\n[CHECK 4] Sample point ownership (first $(min(10,length(alya_local_ids)))):")
    for i in 1:min(10, length(alya_local_ids))
        wr  = alya_owner_ranks[i]
        tag = wr < nranks_alya ? "Alya rank $wr" : "Julia rank $(wr-nranks_alya)"
        println("    ID $(alya_local_ids[i]) → world rank $wr ($tag)")
    end

    println("\n[CHECK 5] Julia-wide send totals...")
    send_to_alya = Int32[npoin_send[i] for i in 1:nranks_alya]
    g_send       = MPI.Allreduce(send_to_alya, MPI.SUM, local_comm)
    if lrank == 0
        for i in 1:nranks_alya
            g_send[i] > 0 && println("  Alya rank $(i-1): total Julia→Alya = $(g_send[i]) pts")
        end
    end

    println("\n" * "="^80)
    println("[VERIFY] Complete for Julia lrank=$lrank")
    println("="^80)
    return true
end

# ===========================================================================
# DIFFERENTIALEQUATIONS.JL CALLBACK
# ===========================================================================

function setup_coupling_callback(is_coupled, params, inputs)
    is_coupled === false && return nothing

    cpg = params.coupling
    @assert cpg !== nothing "params.coupling must be set during setup."

    t0   = params.tspan[1]
    tol0 = get(inputs, :couple_time_tol, 1e-12)

    @inline coupling_condition(u_state, t, integrator) = t > t0 + tol0

    neqs = params.neqs
    mesh = params.mesh

    _qout        = cpg.qout::Matrix{Float64}
    _u_interp    = cpg.u_interp::Matrix{Float64}
    _ξ_nodes     = cpg.ξ_nodes_ref::Vector{Float64}
    _ω           = cpg.ω_bary::Vector{Float64}
    _e_conn      = cpg.elem_conn::Matrix{Int}
    _elem_x      = cpg.elem_x::Matrix{Float64}
    _elem_y      = cpg.elem_y::Matrix{Float64}
    _ψξ          = cpg.ψξ_scratch::Vector{Float64}
    _ψη          = cpg.ψη_scratch::Vector{Float64}
    _dψξ         = cpg.dψξ_scratch::Vector{Float64}
    _dψη         = cpg.dψη_scratch::Vector{Float64}
    _α           = cpg.α_scratch::Vector{Float64}
    _x_e         = cpg.x_e_scratch::Vector{Float64}
    _y_e         = cpg.y_e_scratch::Vector{Float64}
    _alya_coords = cpg.alya_local_coords::Matrix{Float64}
    _owner_ranks = cpg.alya_owner_ranks::Vector{Int32}
    _elem_bboxes = cpg.elem_bboxes::Vector{NTuple{4,Float64}}
    _bins        = cpg.interp_bins::ElemBins
    _mesh_x      = mesh.x::Vector{Float64}
    _mesh_y      = mesh.y::Vector{Float64}

    function do_coupling_exchange!(integrator)
        je_perform_coupling_exchange(
            integrator.u, integrator.p.uaux, integrator.t, cpg,
            _qout, _u_interp, _ξ_nodes, _ω, _e_conn, _elem_x, _elem_y,
            _ψξ, _ψη, _dψξ, _dψη, _α, _x_e, _y_e,
            _alya_coords, _owner_ranks, _mesh_x, _mesh_y,
            inputs, neqs, _elem_bboxes, _bins)
    end

    return DiscreteCallback(coupling_condition, do_coupling_exchange!)
end
