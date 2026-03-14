using Gridap
using GridapDistributed
using PartitionedArrays
using Base.Threads

# ---------------------------------------------------------------------------
# Global MPI communicator refs
# Default is COMM_WORLD; in coupled mode these are set by the driver.
# ---------------------------------------------------------------------------
const JEXPRESSO_MPI_COMM       = Ref{Union{MPI.Comm,Nothing}}(nothing)
const JEXPRESSO_MPI_COMM_WORLD = Ref{Union{MPI.Comm,Nothing}}(nothing)
const JEXPRESSO_COUPLING_DATA  = Ref{Union{Dict{Symbol,Any},Nothing}}(nothing)

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

# ---------------------------------------------------------------------------
# CouplingData — all communication metadata for one Julia rank
# ---------------------------------------------------------------------------
mutable struct CouplingData
    npoin_recv::Vector{Int32}
    npoin_send::Vector{Int32}
    recv_from_ranks::Vector{Int32}
    send_to_ranks::Vector{Int32}

    comm_world::MPI.Comm
    lrank::Int
    neqs::Int
    ndime::Int

    send_bufs::Union{Nothing, Vector{Vector{Float64}}}
    recv_bufs::Union{Nothing, Vector{Vector{Float64}}}
    send_coord_bufs::Union{Nothing, Vector{Vector{Float64}}}

    alya_local_coords::Union{Nothing, Matrix{Float64}}
    alya_local_ids::Union{Nothing, Vector{Int32}}
    alya_owner_ranks::Union{Nothing, Vector{Int32}}

    # Precomputed interpolation data (mesh geometry never changes)
    elem_bboxes::Union{Nothing, Vector{NTuple{4,Float64}}}
    interp_bins::Union{Nothing, Dict{Symbol,Any}}

    function CouplingData(; npoin_recv, npoin_send, recv_from_ranks, send_to_ranks,
                          comm_world, lrank, neqs, ndime)
        new(npoin_recv, npoin_send, recv_from_ranks, send_to_ranks,
            comm_world, lrank, neqs, ndime,
            nothing, nothing, nothing,
            nothing, nothing, nothing,
            nothing, nothing)
    end
end

# ---------------------------------------------------------------------------
# Alya_Point_Ownership — which Jexpresso rank owns each Alya grid point
# ---------------------------------------------------------------------------
struct Alya_Point_Ownership
    alya_coords::Matrix{Float64}
    owner_jrank::Vector{Int32}
    owner_wrank::Vector{Int32}
    alya_point_id::Vector{Int32}
    ndime::Int32
end

# ===========================================================================
# MPI HANDSHAKE  (Julia ← Alya)
# ===========================================================================

# Initial identity exchange with Alya. Returns true when running in coupled mode
# (world size > nparts), false in standalone mode.
function je_perform_coupling_handshake(world, nparts)
    wsize = MPI.Comm_size(world)
    wrank = MPI.Comm_rank(world)
    if wsize <= nparts
        return false
    end
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

    if wsize <= nparts
        @warn "je_receive_alya_data called but not in coupled mode"
        return
    end

    # 1. ndime
    ndime_buf = Vector{Int32}(undef, 1)
    MPI.Bcast!(ndime_buf, 0, world)
    ndime = Int(ndime_buf[1])

    # 2. rem_min / rem_max / rem_nx (element-by-element to avoid SubArray issues)
    rem_min = Vector{Float64}(undef, 3)
    rem_max = Vector{Float64}(undef, 3)
    rem_nx  = Vector{Int32}(undef, 3)
    for idime in 1:3
        MPI.Bcast!(@view(rem_min[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_max[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_nx[idime:idime]),  0, world)
    end

    # 3. Alya→world rank map via Allreduce
    nranks_alya  = wsize - nparts
    alya2world_l = zeros(Int32, nranks_alya)
    alya2world   = MPI.Allreduce(alya2world_l, MPI.SUM, world)

    set_coupling_data(Dict{Symbol,Any}(
        :ndime      => ndime,
        :neqs       => 0,       # filled later by neqs Allreduce
        :nsteps     => 0,
        :rem_min    => rem_min,
        :rem_max    => rem_max,
        :rem_nx     => rem_nx,
        :alya2world => alya2world,
    ))

    lcomm = get_mpi_comm()
    lrank = MPI.Comm_rank(lcomm)
    if lrank == 0
        println("[je_receive_alya_data] ndime=$ndime")
        println("  min=$rem_min, max=$rem_max, nx=$rem_nx")
        flush(stdout)
    end
end

# Send the Alya grid-point IDs for each partner: only those whose owner
# equals dest_rank, preserving the same order used for data packing.
function je_send_node_list(alya_local_ids::Vector{Int64},
                           alya_owner_ranks::Vector{Int32},
                           send_to_ranks::Vector{Int32},
                           world::MPI.Comm)
    lcomm = get_mpi_comm()
    lrank = MPI.Comm_rank(lcomm)
    wrank = MPI.Comm_rank(world)

    send_requests = MPI.Request[]
    for dest_rank in send_to_ranks
        mask    = alya_owner_ranks .== dest_rank
        gid_buf = Int64.(alya_local_ids[mask])
        push!(send_requests, MPI.Isend(gid_buf, dest_rank, 0, world))
        println("[je_send_node_list] Jexpresso lrank=$lrank (wrank=$wrank) → Alya world rank $dest_rank: ",
                "$(length(gid_buf)) node IDs")
        for (k, gid) in enumerate(gid_buf)
            println("  [$k] gid = $gid")
        end
    end
    isempty(send_requests) || MPI.Waitall(send_requests)
    flush(stdout)
end

# ===========================================================================
# GEOMETRY / INTERPOLATION HELPERS
# ===========================================================================

# Connectivity accessor using connijk (tensor-product ordered, i fastest)
function _make_conn_accessor(mesh)
    connijk = getfield(mesh, :connijk)
    nsd = mesh.nsd
    if nsd <= 2
        return e -> vec(@view connijk[e, :, :, 1])
    else
        return e -> vec(@view connijk[e, :, :, :])
    end
end

_num_elems(mesh) = mesh.nelem

# Barycentric weights:  w_j = 1 / ∏_{k≠j} (x_j - x_k)
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

# Uniform binning over element bounding boxes
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

    bins = [Int[] for _ in 1:(nx*ny)]
    for e in 1:ne
        (x0,x1,y0,y1) = elem_bboxes[e]
        ix0 = clamp(Int(floor((x0 - xmin)/dx)), 0, nx-1)
        ix1 = clamp(Int(floor((x1 - xmin)/dx)), 0, nx-1)
        iy0 = clamp(Int(floor((y0 - ymin)/dy)), 0, ny-1)
        iy1 = clamp(Int(floor((y1 - ymin)/dy)), 0, ny-1)
        for iy in iy0:iy1, ix in ix0:ix1
            push!(bins[iy*nx + ix + 1], e)
        end
    end
    return Dict(:xmin=>xmin, :xmax=>xmax, :ymin=>ymin, :ymax=>ymax,
                :nx=>nx, :ny=>ny, :dx=>dx, :dy=>dy, :bins=>bins)
end

function _bin_candidates(bins, x, y, elem_bboxes)
    nx = bins[:nx]; dx = bins[:dx]; xmin = bins[:xmin]
    ny = bins[:ny]; dy = bins[:dy]; ymin = bins[:ymin]
    ix = clamp(Int(floor((x - xmin)/dx)), 0, nx-1)
    iy = clamp(Int(floor((y - ymin)/dy)), 0, ny-1)
    cand = bins[:bins][iy*nx + ix + 1]
    return isempty(cand) ? (1:length(elem_bboxes)) : cand
end

# 1D Lagrange basis at ξ using barycentric formula
function evaluate_lagrange_1d(ξ::Float64, ξ_nodes::Vector{Float64}, ω::Vector{Float64})
    n = length(ξ_nodes)
    ψ = zeros(Float64, n)
    @inbounds for i in 1:n
        if abs(ξ - ξ_nodes[i]) < 1e-14
            ψ[i] = 1.0
            return ψ
        end
    end
    sum_val = 0.0
    @inbounds for i in 1:n
        ψ[i] = ω[i] / (ξ - ξ_nodes[i])
        sum_val += ψ[i]
    end
    inv_s = 1.0 / sum_val
    @inbounds for i in 1:n; ψ[i] *= inv_s; end
    return ψ
end

# 1D Lagrange derivative at ξ using barycentric formula
function evaluate_lagrange_1d_derivative(ξ::Float64, ξ_nodes::Vector{Float64},
                                         ω::Vector{Float64})
    n = length(ξ_nodes)
    dψ = zeros(Float64, n)
    @inbounds for j in 1:n
        if abs(ξ - ξ_nodes[j]) < 1e-14
            for k in 1:n
                k == j && continue
                dψ[k]  =  ω[k] / (ω[j] * (ξ_nodes[j] - ξ_nodes[k]))
                dψ[j] -= 1.0 / (ξ_nodes[j] - ξ_nodes[k])
            end
            return dψ
        end
    end
    α = zeros(Float64, n)
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
    return dψ
end

# Map physical (px, py) → reference (ξ, η) via Newton iteration
function physical_to_reference(px::Float64, py::Float64,
                                x_elem::AbstractVector, y_elem::AbstractVector,
                                ξ_nodes::Vector{Float64}, ω::Vector{Float64}, ngl::Int)
    ξ, η = 0.0, 0.0
    for _ in 1:20
        ψξ  = evaluate_lagrange_1d(ξ, ξ_nodes, ω)
        ψη  = evaluate_lagrange_1d(η, ξ_nodes, ω)
        dψξ = evaluate_lagrange_1d_derivative(ξ, ξ_nodes, ω)
        dψη = evaluate_lagrange_1d_derivative(η, ξ_nodes, ω)

        x_c = y_c = dxdξ = dxdη = dydξ = dydη = 0.0
        idx = 1
        @inbounds for j in 1:ngl, i in 1:ngl
            ψv = ψξ[i]*ψη[j]; dξv = dψξ[i]*ψη[j]; dηv = ψξ[i]*dψη[j]
            x_c += ψv*x_elem[idx];  y_c  += ψv*y_elem[idx]
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

# Interpolate SEM solution to arbitrary coordinates.
# Pass precomp_bboxes and precomp_bins (precomputed once at setup) to avoid
# rebuilding them on every timestep call.
function interpolate_solution_to_alya_coords(alya_coords::Matrix{Float64}, mesh,
                                             u_mat::AbstractMatrix, basis, ξ, neqs, inputs;
                                             use_bins::Bool=true, bins_per_dim::Int=64,
                                             precomp_bboxes=nothing, precomp_bins=nothing)
    n_points = size(alya_coords, 1)
    u_interp = zeros(Float64, n_points, neqs)

    get_conn = _make_conn_accessor(mesh)
    nelem = _num_elems(mesh)
    ngl = mesh.ngl
    ξ_nodes = Vector{Float64}(ξ)
    ω = barycentric_weights(ξ_nodes)

    elem_bboxes = if precomp_bboxes !== nothing
        precomp_bboxes
    else
        bb = Vector{NTuple{4,Float64}}(undef, nelem)
        @inbounds for e in 1:nelem
            ns = get_conn(e)
            bb[e] = (minimum(mesh.x[ns]), maximum(mesh.x[ns]),
                     minimum(mesh.y[ns]), maximum(mesh.y[ns]))
        end
        bb
    end

    bins = if precomp_bins !== nothing
        precomp_bins
    elseif use_bins
        _build_elem_bins(elem_bboxes; bins_per_dim=bins_per_dim)
    else
        nothing
    end

    @inbounds for ipt in 1:n_points
        px, py = alya_coords[ipt,1], alya_coords[ipt,2]
        candidates = bins !== nothing ? _bin_candidates(bins, px, py, elem_bboxes) : (1:nelem)
        found = false
        for e in candidates
            ns = get_conn(e)
            x_e = mesh.x[ns]; y_e = mesh.y[ns]
            (px < minimum(x_e)-1e-10 || px > maximum(x_e)+1e-10 ||
             py < minimum(y_e)-1e-10 || py > maximum(y_e)+1e-10) && continue
            ξr, ηr, conv = physical_to_reference(px, py, x_e, y_e, ξ_nodes, ω, ngl)
            (conv && abs(ξr) <= 1+1e-10 && abs(ηr) <= 1+1e-10) || continue
            ψξ = evaluate_lagrange_1d(ξr, ξ_nodes, ω)
            ψη = evaluate_lagrange_1d(ηr, ξ_nodes, ω)
            for q in 1:neqs
                val = 0.0; idx = 1
                for j in 1:ngl, i in 1:ngl
                    val += ψξ[i]*ψη[j]*u_mat[ns[idx], q]; idx += 1
                end
                u_interp[ipt,q] = val
            end
            found = true; break
        end
        if !found   # nearest-neighbour fallback
            nearest = argmin((mesh.x .- px).^2 .+ (mesh.y .- py).^2)
            for q in 1:neqs; u_interp[ipt,q] = u_mat[nearest,q]; end
        end
    end
    return u_interp
end

# ===========================================================================
# ALYA COORDINATE EXTRACTION
# ===========================================================================

function extract_local_alya_coordinates(mesh, coupling_data, local_comm, world_comm;
                                        block_size::NTuple{3,Int}=(64,64,64),
                                        use_cropping::Bool=true)
    ndime      = coupling_data[:ndime]
    rem_min_f  = coupling_data[:rem_min]
    rem_max_f  = coupling_data[:rem_max]
    rem_nx_i   = coupling_data[:rem_nx]
    alya2world = coupling_data[:alya2world]
    @assert ndime == 2 || ndime == 3 "Only ndime==2 or ndime==3 supported"
    @assert length(alya2world) > 0 "alya2world must be non-empty"

    rem_min = Float64.(rem_min_f); rem_max = Float64.(rem_max_f)
    rem_nx  = Int.(rem_nx_i)

    rem_dx = zeros(Float64, 3)
    for d in 1:ndime
        rem_dx[d] = rem_nx[d] > 1 ? (rem_max[d] - rem_min[d]) / (rem_nx[d] - 1) : 0.0
    end

    nmax = rem_nx[1] * rem_nx[2] * rem_nx[3]
    nranks_alya = length(alya2world)
    alya_driving_world_rank = Int32(0)
    alya_worker_indices = [k for k in 1:nranks_alya if alya2world[k] != alya_driving_world_rank]
    nworkers_alya = length(alya_worker_indices)
    r_w  = mod(nmax, nworkers_alya)
    np_w = div(nmax, nworkers_alya)

    xmin_local = minimum(mesh.x); xmax_local = maximum(mesh.x)
    ymin_local = minimum(mesh.y); ymax_local = maximum(mesh.y)
    zmin_local = ndime == 3 ? minimum(mesh.z) : 0.0
    zmax_local = ndime == 3 ? maximum(mesh.z) : 0.0
    tol = 1e-10

    lrank = MPI.Comm_rank(local_comm)
    wrank = MPI.Comm_rank(world_comm)

    @inline idx_to_xyz(i1,i2,i3) = (rem_min[1]+i1*rem_dx[1],
                                    rem_min[2]+i2*rem_dx[2],
                                    rem_min[3]+i3*rem_dx[3])
    @inline idx_to_ipoin(i1,i2,i3) = i1 + rem_nx[1]*(i2 + rem_nx[2]*i3) + 1  # 1-based

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
        ilo = max(0,   Int(clamp(floor((minB - minA)/dx + 1e-12), 0, nx-1)) - 1)
        ihi = min(nx-1, Int(clamp( ceil((maxB - minA)/dx - 1e-12), 0, nx-1)) + 1)
        return ilo, ihi
    end

    i1_lo, i1_hi = 0, rem_nx[1]-1
    i2_lo, i2_hi = 0, rem_nx[2]-1
    i3_lo, i3_hi = 0, ndime == 3 ? rem_nx[3]-1 : 0

    if use_cropping
        rem_nx[1] > 1 && ((i1_lo, i1_hi) = crop_1d(rem_min[1], rem_dx[1], rem_nx[1], xmin_local-tol, xmax_local+tol))
        rem_nx[2] > 1 && ((i2_lo, i2_hi) = crop_1d(rem_min[2], rem_dx[2], rem_nx[2], ymin_local-tol, ymax_local+tol))
        ndime == 3 && rem_nx[3] > 1 && ((i3_lo, i3_hi) = crop_1d(rem_min[3], rem_dx[3], rem_nx[3], zmin_local-tol, zmax_local+tol))
        if i1_lo > i1_hi || i2_lo > i2_hi || i3_lo > i3_hi
            println("[extract_local_alya_coordinates] (lrank=$lrank) no overlap with local domain.")
            flush(stdout)
            return zeros(Float64, 0, ndime), Int32[], Int32[]
        end
    end

    Bx, By, Bz = block_size
    Bx = max(1,Bx); By = max(1,By); Bz = max(1, ndime==3 ? Bz : 1)
    est = (i1_hi-i1_lo+1)*(i2_hi-i2_lo+1)*(i3_hi-i3_lo+1)
    X = Float64[]; Y = Float64[]; Z = ndime==3 ? Float64[] : Float64[]
    ids = Int32[]; owners = Int32[]
    sizehint!(X, est); sizehint!(Y, est)
    ndime==3 && sizehint!(Z, est)
    sizehint!(ids, est); sizehint!(owners, est)

    @inline function bin_bbox(i1s,i1e,i2s,i2e,i3s,i3e)
        (rem_min[1]+i1s*rem_dx[1], rem_min[1]+i1e*rem_dx[1],
         rem_min[2]+i2s*rem_dx[2], rem_min[2]+i2e*rem_dx[2],
         rem_min[3]+i3s*rem_dx[3], rem_min[3]+i3e*rem_dx[3])
    end
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
                if inside_box(bxmn,bxmx,bymn,bymx,bzmn,bzmx)
                    @inbounds for i3 in i3s:i3e, i2 in i2s:i2e, i1 in i1s:i1e
                        x,y,z = idx_to_xyz(i1,i2,i3)
                        ipoin = idx_to_ipoin(i1,i2,i3)
                        aw = alya2world[owner_alya_rank(ipoin)]
                        push!(X,x); push!(Y,y); ndime==3 && push!(Z,z)
                        push!(ids, Int32(ipoin)); push!(owners, Int32(aw))
                    end
                else
                    @inbounds for i3 in i3s:i3e, i2 in i2s:i2e, i1 in i1s:i1e
                        x,y,z = idx_to_xyz(i1,i2,i3)
                        if inside_local(x,y,z)
                            ipoin = idx_to_ipoin(i1,i2,i3)
                            aw = alya2world[owner_alya_rank(ipoin)]
                            push!(X,x); push!(Y,y); ndime==3 && push!(Z,z)
                            push!(ids, Int32(ipoin)); push!(owners, Int32(aw))
                        end
                    end
                end
            end
        end
    end

    n_local = length(ids)
    alya_local_coords = zeros(Float64, n_local, ndime)
    if ndime == 2
        @inbounds for i in 1:n_local
            alya_local_coords[i,1] = X[i]; alya_local_coords[i,2] = Y[i]
        end
    else
        @inbounds for i in 1:n_local
            alya_local_coords[i,1] = X[i]; alya_local_coords[i,2] = Y[i]
            alya_local_coords[i,3] = Z[i]
        end
    end

    # Sort by (owner_rank, global_id) so data arrives at Alya in ascending index order
    if n_local > 1
        perm = sortperm(collect(zip(owners, ids)))
        alya_local_coords = alya_local_coords[perm, :]
        ids    = ids[perm]
        owners = owners[perm]
    end

    println("[extract_local_alya_coordinates] (lrank=$lrank, wrank=$wrank) ",
            "collected $n_local Alya points (cropping=$(use_cropping)).")
    flush(stdout)
    return alya_local_coords, ids, owners
end

# ===========================================================================
# BUFFER ALLOCATION
# ===========================================================================

function allocate_coupling_buffers(npoin_recv, npoin_send, neqs, ndime)
    wsize = length(npoin_recv)
    send_bufs       = [zeros(Float64, npoin_send[i] * neqs)  for i in 1:wsize]
    recv_bufs       = [zeros(Float64, npoin_recv[i] * neqs)  for i in 1:wsize]
    send_coord_bufs = [zeros(Float64, npoin_send[i] * ndime) for i in 1:wsize]
    return send_bufs, recv_bufs, send_coord_bufs
end

# ===========================================================================
# OWNERSHIP MAP  (diagnostics)
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

    nmax        = rem_nx[1]*rem_nx[2]*rem_nx[3]
    nranks_alya = length(alya2world)
    nworkers    = nranks_alya - 1          # rank 0 is the master
    r     = mod(nmax, nworkers)
    npoin = div(nmax, nworkers)

    xmin_l = minimum(mesh.x); xmax_l = maximum(mesh.x)
    ymin_l = minimum(mesh.y); ymax_l = maximum(mesh.y)
    zmin_l = ndime==3 ? minimum(mesh.z) : 0.0
    zmax_l = ndime==3 ? maximum(mesh.z) : 0.0
    tol = 1e-10

    local_coords      = zeros(Float64, nmax, ndime)
    local_owner_jrank = fill(Int32(-1), nmax)
    local_owner_wrank = fill(Int32(-1), nmax)

    for ipoin in 1:nmax
        i0 = ipoin - 1
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

    return Alya_Point_Ownership(local_coords,
                                global_owner_jrank, global_owner_wrank,
                                collect(Int32(1):Int32(nmax)), Int32(ndime))
end

function print_alya_point_ownership(ownership, coupling_data;
                                    max_points_to_print=20, only_owned=false)
    lrank = MPI.Comm_rank(get_mpi_comm())
    lrank != 0 && return

    npoints = length(ownership.alya_point_id)
    ndime   = ownership.ndime
    n_owned   = count(x -> x >= 0, ownership.owner_jrank)
    n_unowned = npoints - n_owned

    println("\n" * "="^70)
    println("Alya Point Ownership Map")
    println("="^70)
    println("Total Alya points: $npoints  (owned: $n_owned, unowned: $n_unowned)")
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
            @sprintf("(%10.2f, %10.2f)",   ownership.alya_coords[i,1], ownership.alya_coords[i,2]) :
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
    lrank = MPI.Comm_rank(get_mpi_comm())
    lrank != 0 && return

    lsize   = MPI.Comm_size(get_mpi_comm())
    npoints = length(ownership.alya_point_id)

    println("\n" * "="^70)
    println("Alya Point Distribution Statistics")
    println("="^70)

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
    if n_unowned > 0
        println("WARNING: $n_unowned Alya points are not owned by any Jexpresso rank!")
    else
        println("✓ All Alya points are owned by Jexpresso ranks")
    end
    println("="^70); println(); flush(stdout)
end

# ===========================================================================
# COUPLING SETUP ORCHESTRATOR
# ===========================================================================

function setup_coupling_and_mesh(world, lsize, inputs, nranks, distribute, rank, OUTPUT_DIR, TFloat)

    je_receive_alya_data(world, lsize)
    coupling_data = get_coupling_data()
    local_comm = get_mpi_comm()
    lrank = MPI.Comm_rank(local_comm)

    sem, partitioned_model = sem_setup(inputs, nranks, distribute, rank)
    if inputs[:backend] != CPU()
        convert_mesh_arrays!(sem.mesh.SD, sem.mesh, inputs[:backend], inputs)
    end

    alya_local_coords, alya_local_ids, alya_owner_ranks =
        extract_local_alya_coordinates(sem.mesh, coupling_data, local_comm, world;
                                       block_size=(64,64,64), use_cropping=true)

    wsize = MPI.Comm_size(world)
    npoin_recv = zeros(Int32, wsize)
    for owner_wrank in alya_owner_ranks
        npoin_recv[owner_wrank + 1] += 1
    end

    recv_from_ranks = Int32[i-1 for i in 1:wsize if npoin_recv[i] > 0]

    if lrank == 0
        println("[setup_coupling] Alltoall to exchange point counts…"); flush(stdout)
    end

    send_counts_to_alya  = Vector{Int32}(npoin_recv)
    recv_counts_from_alya = zeros(Int32, wsize)
    MPI.Barrier(world)
    MPI.Alltoall!(send_counts_to_alya, recv_counts_from_alya, 1, world)

    npoin_send    = copy(npoin_recv)
    send_to_ranks = copy(recv_from_ranks)

    je_send_node_list(Int64.(alya_local_ids), alya_owner_ranks, send_to_ranks, world)

    verify_coupling_communication_pattern(npoin_recv, npoin_send,
                                          alya_owner_ranks, alya_local_ids,
                                          coupling_data, local_comm, world)

    if lrank == 0
        println("[setup_coupling] Communication pattern:")
        println("  npoin_recv: ", [(i-1, npoin_recv[i]) for i in 1:wsize if npoin_recv[i]>0])
        println("  npoin_send: ", [(i-1, npoin_send[i]) for i in 1:wsize if npoin_send[i]>0])
        flush(stdout)
    end

    ownership = build_alya_point_ownership_map(sem.mesh, coupling_data, local_comm, world)
    print_alya_point_ownership(ownership, coupling_data; max_points_to_print=20, only_owned=true)
    verify_alya_point_distribution(ownership, coupling_data)

    qp = initialize(sem.mesh.SD, 0, sem.mesh, inputs, OUTPUT_DIR, TFloat)

    # Broadcast neqs to Alya via Allreduce(MAX): Alya contributes 0, Julia contributes qp.neqs
    neqs_buf = Ref(Int32(qp.neqs))
    MPI.Allreduce!(neqs_buf, MPI.MAX, world)

    ndime = coupling_data[:ndime]
    coupling = CouplingData(
        npoin_recv    = npoin_recv,
        npoin_send    = npoin_send,
        recv_from_ranks = recv_from_ranks,
        send_to_ranks   = send_to_ranks,
        comm_world    = world,
        lrank         = lrank,
        neqs          = qp.neqs,
        ndime         = ndime,
    )

    coupling.send_bufs, coupling.recv_bufs, coupling.send_coord_bufs =
        allocate_coupling_buffers(npoin_recv, npoin_send, coupling.neqs, ndime)

    coupling.alya_local_coords = alya_local_coords
    coupling.alya_local_ids    = alya_local_ids
    coupling.alya_owner_ranks  = alya_owner_ranks

    # Precompute element bounding boxes and spatial bins once — mesh never changes.
    let get_conn = _make_conn_accessor(sem.mesh),
        nelem    = _num_elems(sem.mesh)
        bb = Vector{NTuple{4,Float64}}(undef, nelem)
        @inbounds for e in 1:nelem
            ns = get_conn(e)
            bb[e] = (minimum(sem.mesh.x[ns]), maximum(sem.mesh.x[ns]),
                     minimum(sem.mesh.y[ns]), maximum(sem.mesh.y[ns]))
        end
        coupling.elem_bboxes = bb
        coupling.interp_bins = _build_elem_bins(bb; bins_per_dim=64)
    end

    return coupling, sem, partitioned_model, qp
end

# ===========================================================================
# TIME-LOOP DATA EXCHANGE
# ===========================================================================

function pack_interpolated_data!(cpg::CouplingData, interp_values::Matrix{Float64},
                                 alya_owner_ranks::Vector{Int32},
                                 alya_local_coords::Matrix{Float64})
    n_local = size(interp_values, 1)
    neqs    = cpg.neqs
    ndime   = cpg.ndime
    wsize   = length(cpg.npoin_send)

    for i in 1:wsize
        fill!(cpg.send_bufs[i], 0.0)
        fill!(cpg.send_coord_bufs[i], 0.0)
    end

    send_offsets  = zeros(Int, wsize)
    coord_offsets = zeros(Int, wsize)

    @inbounds for i in 1:n_local
        bidx = alya_owner_ranks[i] + 1
        cpg.npoin_send[bidx] > 0 || continue

        off = send_offsets[bidx]
        for q in 1:neqs; cpg.send_bufs[bidx][off+q] = interp_values[i,q]; end
        send_offsets[bidx] += neqs

        coff = coord_offsets[bidx]
        for d in 1:ndime; cpg.send_coord_bufs[bidx][coff+d] = alya_local_coords[i,d]; end
        coord_offsets[bidx] += ndime
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

function unpack_received_data!(cpg::CouplingData, u::AbstractVector, mesh,
                               alya_local_coords::Matrix{Float64},
                               alya_local_ids::Vector{Int32})
    # Placeholder: implement physics-specific coupling (forcing, BCs, etc.)
end

function je_perform_coupling_exchange(u, u_mat, t, cpg::CouplingData,
                                      mesh, basis, inputs, ξ, neqs)
    npoin = mesh.npoin
    qout  = zeros(Float64, npoin, neqs)
    u2uaux!(u_mat, u, neqs, npoin)
    call_user_uout(qout, u_mat, u_mat, 0, inputs[:SOL_VARS_TYPE], npoin, neqs, neqs)

    u_interp = interpolate_solution_to_alya_coords(
        cpg.alya_local_coords, mesh, qout, basis, ξ, neqs, inputs;
        use_bins=true, bins_per_dim=64,
        precomp_bboxes=cpg.elem_bboxes, precomp_bins=cpg.interp_bins)

    pack_interpolated_data!(cpg, u_interp, cpg.alya_owner_ranks, cpg.alya_local_coords)
    coupling_exchange_data!(cpg)
    unpack_received_data!(cpg, u, mesh, cpg.alya_local_coords, cpg.alya_local_ids)
end

# ===========================================================================
# VERIFICATION  (setup-time diagnostic — no Alya participation needed)
# ===========================================================================

function verify_coupling_communication_pattern(npoin_recv, npoin_send,
                                               alya_owner_ranks, alya_local_ids,
                                               coupling_data, local_comm, world_comm)
    lrank = MPI.Comm_rank(local_comm)
    wrank = MPI.Comm_rank(world_comm)
    wsize = MPI.Comm_size(world_comm)
    lsize = MPI.Comm_size(local_comm)
    nranks_alya = length(coupling_data[:alya2world])

    println("="^80)
    println("[VERIFY] Julia lrank=$lrank, wrank=$wrank")
    println("="^80)

    # CHECK 1: npoin_recv matches alya_owner_ranks distribution
    println("\n[CHECK 1] Verifying npoin_recv against alya_owner_ranks...")
    local_counts = zeros(Int32, wsize)
    for wr in alya_owner_ranks; local_counts[wr+1] += 1; end
    mismatch = false
    for i in 1:wsize
        if local_counts[i] != npoin_recv[i]
            println("  MISMATCH at world rank $(i-1): extracted=$(local_counts[i]) vs npoin_recv=$(npoin_recv[i])")
            mismatch = true
        end
    end
    mismatch ? (@warn "npoin_recv does NOT match alya_owner_ranks!") :
               println("  ✓ npoin_recv matches alya_owner_ranks distribution")

    # CHECK 2: local totals
    println("\n[CHECK 2] Verifying total points to receive...")
    total_local = length(alya_owner_ranks)
    total_recv  = sum(npoin_recv)
    println("  Local Alya points extracted: $total_local")
    println("  Sum of npoin_recv: $total_recv")
    total_local == total_recv ? println("  ✓ Total counts match") :
                                (@warn "Total counts do NOT match!")

    # CHECK 3: global conservation (Julia ranks only)
    println("\n[CHECK 3] Checking global point conservation...")
    rem_nx = coupling_data[:rem_nx]
    total_alya = rem_nx[1]*rem_nx[2]*rem_nx[3]
    global_owned = MPI.Allreduce(total_local, MPI.SUM, local_comm)
    println("  Total Alya grid points: $total_alya")
    println("  Total points owned by all Julia ranks: $global_owned")
    global_owned == total_alya ? println("  ✓ All Alya points accounted for") :
                                 (@warn "Point conservation failed! Expected $total_alya, got $global_owned")

    # CHECK 4: send/recv breakdown
    println("\n[CHECK 4] Communication pattern:")
    println("  RECEIVE (npoin_recv):")
    for i in 1:wsize
        npoin_recv[i] == 0 && continue
        wr = i-1
        tag = wr < nranks_alya ? "Alya rank $wr" : "Julia rank $(wr-nranks_alya)"
        println("    $(npoin_recv[i]) pts FROM world rank $wr ($tag)")
    end
    println("  SEND (npoin_send):")
    for i in 1:wsize
        npoin_send[i] == 0 && continue
        wr = i-1
        tag = wr < nranks_alya ? "Alya rank $wr" : "Julia rank $(wr-nranks_alya)"
        println("    $(npoin_send[i]) pts TO world rank $wr ($tag)")
    end

    # CHECK 5: symmetry
    println("\n[CHECK 5] Communication partners:")
    println("  Recv from: $([i-1 for i in 1:wsize if npoin_recv[i]>0])")
    println("  Send to:   $([i-1 for i in 1:wsize if npoin_send[i]>0])")

    # CHECK 6: sample IDs
    println("\n[CHECK 6] Sample Alya point ownership (first $(min(10,length(alya_local_ids)))):")
    for i in 1:min(10, length(alya_local_ids))
        wr = alya_owner_ranks[i]
        tag = wr < nranks_alya ? "Alya rank $wr" : "Julia rank $(wr-nranks_alya)"
        println("    ID $(alya_local_ids[i]) → world rank $wr ($tag)")
    end

    # CHECK 7: Julia-wide send/recv totals (Julia ranks only)
    println("\n[CHECK 7] Julia-wide totals (via local_comm Allreduce)...")
    send_to_alya  = Int32[npoin_send[i]  for i in 1:nranks_alya]
    recv_fr_alya  = Int32[npoin_recv[i]  for i in 1:nranks_alya]
    g_send = MPI.Allreduce(send_to_alya,  MPI.SUM, local_comm)
    g_recv = MPI.Allreduce(recv_fr_alya,  MPI.SUM, local_comm)
    if lrank == 0
        for i in 1:nranks_alya
            (g_send[i]>0 || g_recv[i]>0) &&
                println("  Alya rank $(i-1): Julia sends $(g_send[i]), receives $(g_recv[i])")
        end
    end

    println("\n" * "="^80)
    println("[VERIFY] Complete for Julia lrank=$lrank")
    println("="^80)
    return true
end
