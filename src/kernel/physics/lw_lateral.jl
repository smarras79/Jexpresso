# Longwave lateral inflow from plane-parallel solutions of the LOCAL column at each lateral-boundary point.
#
# For a point (x,y,z) on a lateral boundary (x = xmin/xmax or y = ymin/ymax) and a direction with μ = cos θ
# (θ polar from +z), the inflow radiance is the 1D plane-parallel radiance of the column through (x,y):
#   μ < 0 (downward): integrate from the top, I_top = σ T_space⁴ / π
#   μ > 0 (upward):   integrate from the surface, I_sfc = ε σ T_sfc⁴ / π  (no reflection, as in the 3D bottom BC)
#   |μ| → 0:          I → B(z)
# Extinction is κ_abs + σ_scat and the source is (κ_abs + σ_scat) B, i.e. the scattered part assumes J ≈ B
# (LW scattering is weak; in optically thick cloud J → B). Within each layer between column levels, B is linear in
# optical depth and κ is the layer mean; the layer update is exact for that model and stable for any optical depth.
#
# Uses:
#   - real walls (non-periodic meshes): the lateral branch of the LW boundary condition (:RT_lw_lateral_bc => :plane_parallel)
#   - periodic meshes: lw_lateral_prev_override builds, for every (element, direction), the values that replace the
#     lagged periodic inflow, which turns the periodic faces into walls with this inflow (:RT_lw_periodic_mode => :walls)
#
# The columns are gathered once per RT call (collective): every rank packs its nodes on the lateral planes and
# all ranks receive all of them (Allgatherv). Size ~ 4 · n_horizontal · n_vertical nodes · 5 doubles.

const LW_σ_SB = 5.670374419e-8

struct LWColumn
    z::Vector{Float64}          # ascending
    κ::Vector{Float64}          # extinction at the levels
    B::Vector{Float64}          # Planck radiance σT⁴/π at the levels
    I_sfc::Float64
    I_top::Float64
end

mutable struct LWLateralColumns
    cols::Dict{Tuple{Int,Int},LWColumn}                         # key: (x, y) rounded to mm
    cache::Dict{Tuple{Tuple{Int,Int},Float64},Vector{Float64}}  # (column, μ) -> radiance at the column levels
    xmin::Float64; xmax::Float64; ymin::Float64; ymax::Float64
    tol::Float64
end

lw_mm(v) = round(Int, v * 1000)

lw_on_lateral(L::LWLateralColumns, x, y) =
    abs(x - L.xmin) < L.tol || abs(x - L.xmax) < L.tol || abs(y - L.ymin) < L.tol || abs(y - L.ymax) < L.tol

"""
    build_lw_lateral_columns(mesh, κ_ext, T, comm; ε_surface, T_space) -> LWLateralColumns

Collective. κ_ext and T are nodal (length npoin). Gathers every column that lies on a lateral plane.
"""
function build_lw_lateral_columns(mesh, κ_ext::AbstractVector, T::AbstractVector, comm;
                                  ε_surface::Float64=0.97, T_space::Float64=0.0)
    npoin = length(κ_ext)
    x = mesh.x; y = mesh.y; z = mesh.z
    xmin = MPI.Allreduce(npoin > 0 ? minimum(view(x, 1:npoin)) :  Inf, MPI.MIN, comm)
    xmax = MPI.Allreduce(npoin > 0 ? maximum(view(x, 1:npoin)) : -Inf, MPI.MAX, comm)
    ymin = MPI.Allreduce(npoin > 0 ? minimum(view(y, 1:npoin)) :  Inf, MPI.MIN, comm)
    ymax = MPI.Allreduce(npoin > 0 ? maximum(view(y, 1:npoin)) : -Inf, MPI.MAX, comm)
    tol = 1e-6 * max(xmax - xmin, ymax - ymin, 1.0)
    L = LWLateralColumns(Dict{Tuple{Int,Int},LWColumn}(), Dict{Tuple{Tuple{Int,Int},Float64},Vector{Float64}}(),
                         xmin, xmax, ymin, ymax, tol)

    buf = Float64[]
    for ip = 1:npoin
        lw_on_lateral(L, x[ip], y[ip]) || continue
        push!(buf, x[ip], y[ip], z[ip], κ_ext[ip], lw_B(T[ip]))
    end
    counts = MPI.Allgather(Int32(length(buf)), comm)
    allbuf = Vector{Float64}(undef, sum(counts))
    MPI.Allgatherv!(buf, MPI.VBuffer(allbuf, counts), comm)

    # group by column, merge duplicate levels (shared nodes on several ranks) by averaging
    acc = Dict{Tuple{Int,Int},Dict{Int,NTuple{4,Float64}}}()    # col -> zmm -> (z, Σκ, ΣB, n)
    for r = 1:5:length(allbuf)
        key = (lw_mm(allbuf[r]), lw_mm(allbuf[r+1]))
        zm = lw_mm(allbuf[r+2])
        d = get!(acc, key, Dict{Int,NTuple{4,Float64}}())
        a = get(d, zm, (allbuf[r+2], 0.0, 0.0, 0.0))
        d[zm] = (a[1], a[2] + allbuf[r+3], a[3] + allbuf[r+4], a[4] + 1.0)
    end
    I_top = lw_B(T_space)
    for (key, d) in acc
        zs = sort(collect(keys(d)))
        zv = [d[k][1] for k in zs]; κv = [d[k][2]/d[k][4] for k in zs]; Bv = [d[k][3]/d[k][4] for k in zs]
        L.cols[key] = LWColumn(zv, κv, Bv, ε_surface * Bv[1], I_top)
    end
    return L
end

# Linear-in-optical-depth source across a layer of slant optical thickness t: contribution at the exit.
@inline function lw_layer(B_in::Float64, B_out::Float64, t::Float64)
    if t < 1e-6
        return B_in * t + (B_out - B_in) * t / 2
    end
    e = exp(-t); a = 1 - e
    return B_in * a + (B_out - B_in) * (1 - a / t)
end

# top_fn = false: start from c.I_top even with :RT_lw_top_bc => :column (used for the column above the lid itself)
function lw_column_radiance(c::LWColumn, μ::Float64; top_fn::Bool=true)
    n = length(c.z); I = zeros(n)
    if abs(μ) < 1e-6
        I .= c.B; return I
    end
    if μ < 0
        ftop = top_fn ? LW_TOP_FN[] : nothing    # :RT_lw_top_bc => :column (lw_top_column.jl); nothing = vacuum (I_top)
        I[n] = ftop === nothing ? c.I_top : ftop(abs(μ))
        for k = n-1:-1:1
            t = 0.5 * (c.κ[k] + c.κ[k+1]) * (c.z[k+1] - c.z[k]) / abs(μ)
            I[k] = I[k+1] * exp(-t) + lw_layer(c.B[k+1], c.B[k], t)
        end
    else
        I[1] = c.I_sfc
        for k = 2:n
            t = 0.5 * (c.κ[k-1] + c.κ[k]) * (c.z[k] - c.z[k-1]) / μ
            I[k] = I[k-1] * exp(-t) + lw_layer(c.B[k-1], c.B[k], t)
        end
    end
    return I
end

"""
    lw_lateral_radiance(L, x, y, z, μ) -> Float64

Plane-parallel radiance of the column through (x, y), at height z, for direction cosine μ.
"""
function lw_lateral_radiance(L::LWLateralColumns, x, y, z, μ::Float64)
    key = (lw_mm(x), lw_mm(y))
    c = get(L.cols, key, nothing)
    c === nothing && error("LW lateral column not found at (x,y) = ($x, $y)")
    I = get!(() -> lw_column_radiance(c, μ), L.cache, (key, μ))
    k = searchsortedfirst(c.z, z - L.tol)
    if k <= length(c.z) && abs(c.z[k] - z) < max(L.tol, 1e-3)
        return I[k]
    end
    # not a column level (should not happen for nodes of the mesh): linear interpolation
    k = clamp(searchsortedlast(c.z, z), 1, length(c.z) - 1)
    w = (z - c.z[k]) / (c.z[k+1] - c.z[k])
    return (1 - w) * I[k] + w * I[k+1]
end

"""
    lw_lateral_prev_override(mesh, L) -> Dict{(iel,d) => (lidx::Vector{Int}, vals::Vector{Float64})}

For every local element with nodes on a lateral plane and every direction of its angular mesh, the plane-parallel
radiance at those nodes. Replacing these entries in the lagged iterate that the sweep reads across periodic faces
turns the periodic faces into walls with plane-parallel inflow.
"""
function lw_lateral_prev_override(mesh, L::LWLateralColumns)
    ngl = mesh.ngl
    out = Dict{Tuple{Int,Int},Tuple{Vector{Int},Vector{Float64}}}()
    for iel = 1:mesh.nelem
        lids = Int[]; pts = NTuple{3,Float64}[]
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel,i,j,k]
            lw_on_lateral(L, mesh.x[ip], mesh.y[ip]) || continue
            push!(lids, lidx(i,j,k,ngl)); push!(pts, (mesh.x[ip], mesh.y[ip], mesh.z[ip]))
        end
        isempty(lids) && continue
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            μ = cos(em.extra_coords[1,d])
            out[(iel,d)] = (lids, [lw_lateral_radiance(L, p[1], p[2], p[3], μ) for p in pts])
        end
    end
    return out
end

"""
    lw_column_mean_init(mesh, κ_ext, T, comm; ε_surface, T_space) -> Dict{(iel,d) => Vector (ngl^3, lidx order)}

Collective. Initial guess for a cold LW solve (:RT_lw_cold_init => :column_mean): the plane-parallel radiance of the
HORIZONTAL-MEAN column (κ_ext and B averaged over all nodes at each height, node-count weighted), at every node and
every direction. Surface emission ε·B̄(z_min), top σT_space⁴/π, as for the lateral columns.
"""
function lw_column_mean_init(mesh, κ_ext::AbstractVector, T::AbstractVector, comm;
                             ε_surface::Float64=0.97, T_space::Float64=0.0)
    npoin = length(κ_ext); ngl = mesh.ngl
    acc = Dict{Int,NTuple{4,Float64}}()                       # zmm -> (z, Σκ, ΣB, n)
    for ip = 1:npoin
        zm = lw_mm(mesh.z[ip])
        a = get(acc, zm, (mesh.z[ip], 0.0, 0.0, 0.0))
        acc[zm] = (a[1], a[2] + κ_ext[ip], a[3] + lw_B(T[ip]), a[4] + 1.0)
    end
    buf = Float64[]
    for (zm, a) in acc; push!(buf, a[1], a[2], a[3], a[4]); end
    counts = MPI.Allgather(Int32(length(buf)), comm)
    allbuf = Vector{Float64}(undef, sum(counts))
    MPI.Allgatherv!(buf, MPI.VBuffer(allbuf, counts), comm)
    g = Dict{Int,NTuple{4,Float64}}()
    for r = 1:4:length(allbuf)
        zm = lw_mm(allbuf[r]); a = get(g, zm, (allbuf[r], 0.0, 0.0, 0.0))
        g[zm] = (a[1], a[2] + allbuf[r+1], a[3] + allbuf[r+2], a[4] + allbuf[r+3])
    end
    zs = sort(collect(keys(g)))
    zv = [g[k][1] for k in zs]; κv = [g[k][2]/g[k][4] for k in zs]; Bv = [g[k][3]/g[k][4] for k in zs]
    col = LWColumn(zv, κv, Bv, ε_surface * Bv[1], lw_B(T_space))
    lev = Dict(zs[i] => i for i in eachindex(zs))
    cache = Dict{Float64,Vector{Float64}}()
    out = Dict{Tuple{Int,Int},Vector{Float64}}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        kz = Vector{Int}(undef, ngl^3)
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            kz[lidx(i,j,k,ngl)] = lev[lw_mm(mesh.z[mesh.connijk[iel,i,j,k]])]
        end
        for d = 1:em.extra_npoin
            I = get!(() -> lw_column_radiance(col, cos(em.extra_coords[1,d])), cache, cos(em.extra_coords[1,d]))
            out[(iel,d)] = I[kz]
        end
    end
    return out
end

"""
    lw_column_local_init(mesh, κ_ext, T, comm; ε_surface, T_space) -> Dict{(iel,d) => Vector (ngl^3, lidx order)}

Collective. Initial guess for a cold LW solve (:RT_lw_cold_init => :column_local): at every node, the plane-parallel
radiance of the node's OWN (x,y) column (independent column approximation), same model as the lateral columns.
Columns are split across ranks, so each column is assembled on an owner rank (chosen from its (x,y) key) and the full
profile is sent back to every rank holding part of it: two Alltoallv, memory ~ the rank's own columns.
"""
function lw_column_local_init(mesh, κ_ext::AbstractVector, T::AbstractVector, comm;
                              ε_surface::Float64=0.97, T_space::Float64=0.0)
    np = MPI.Comm_size(comm); npoin = length(κ_ext); ngl = mesh.ngl
    owner(k) = mod(k[1] * 1_000_003 + k[2], np) + 1
    # 1) nodes -> column owners: (x_mm, y_mm, z, κ, B)
    sbufs = [Float64[] for _ = 1:np]
    for ip = 1:npoin
        key = (lw_mm(mesh.x[ip]), lw_mm(mesh.y[ip]))
        push!(sbufs[owner(key)], key[1], key[2], mesh.z[ip], κ_ext[ip], lw_B(T[ip]))
    end
    recv = lw_alltoallv(sbufs, comm)
    # 2) owners assemble columns; remember which ranks asked for which column
    acc = Dict{Tuple{Int,Int},Dict{Int,NTuple{4,Float64}}}(); askers = Dict{Tuple{Int,Int},Set{Int}}()
    for (src, buf) in enumerate(recv), r = 1:5:length(buf)
        key = (round(Int, buf[r]), round(Int, buf[r+1])); zm = lw_mm(buf[r+2])
        d = get!(acc, key, Dict{Int,NTuple{4,Float64}}())
        a = get(d, zm, (buf[r+2], 0.0, 0.0, 0.0))
        d[zm] = (a[1], a[2] + buf[r+3], a[3] + buf[r+4], a[4] + 1.0)
        push!(get!(askers, key, Set{Int}()), src)
    end
    sbufs = [Float64[] for _ = 1:np]
    for (key, d) in acc
        zs = sort(collect(keys(d))); n = length(zs)
        prof = vcat(Float64[key[1], key[2], n], [d[k][1] for k in zs], [d[k][2]/d[k][4] for k in zs],
                    [d[k][3]/d[k][4] for k in zs])
        for dst in askers[key]; append!(sbufs[dst], prof); end
    end
    recv = lw_alltoallv(sbufs, comm)
    # 3) local columns -> radiance at every node and direction
    I_top = lw_B(T_space)
    cols = Dict{Tuple{Int,Int},LWColumn}()
    for buf in recv
        r = 1
        while r <= length(buf)
            key = (round(Int, buf[r]), round(Int, buf[r+1])); n = round(Int, buf[r+2]); r += 3
            zv = buf[r:r+n-1]; κv = buf[r+n:r+2n-1]; Bv = buf[r+2n:r+3n-1]; r += 3n
            cols[key] = LWColumn(zv, κv, Bv, ε_surface * Bv[1], I_top)
        end
    end
    cache = Dict{Tuple{Tuple{Int,Int},Float64},Vector{Float64}}()
    out = Dict{Tuple{Int,Int},Vector{Float64}}()
    keyk = Vector{Tuple{Int,Int}}(undef, ngl^3); lev = Vector{Int}(undef, ngl^3)
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel,i,j,k]; q = lidx(i,j,k,ngl)
            keyk[q] = (lw_mm(mesh.x[ip]), lw_mm(mesh.y[ip]))
            c = cols[keyk[q]]
            lev[q] = searchsortedfirst(c.z, mesh.z[ip] - 1e-4)
        end
        for d = 1:em.extra_npoin
            μ = cos(em.extra_coords[1,d]); v = Vector{Float64}(undef, ngl^3)
            for q = 1:ngl^3
                I = get!(() -> lw_column_radiance(cols[keyk[q]], μ), cache, (keyk[q], μ))
                v[q] = I[lev[q]]
            end
            out[(iel,d)] = v
        end
    end
    return out
end

# Variable-size all-to-all of Float64 buffers (one per destination rank, 1-based); returns one buffer per source rank.
function lw_alltoallv(sbufs::Vector{Vector{Float64}}, comm)
    np = MPI.Comm_size(comm)
    scounts = Int32[length(b) for b in sbufs]
    rcounts = MPI.Alltoall(MPI.UBuffer(scounts, 1), comm)
    sall = reduce(vcat, sbufs; init=Float64[]); rall = Vector{Float64}(undef, sum(rcounts))
    MPI.Alltoallv!(MPI.VBuffer(sall, scounts), MPI.VBuffer(rall, rcounts), comm)
    out = Vector{Vector{Float64}}(undef, np); o = 0
    for s = 1:np; out[s] = rall[o+1:o+rcounts[s]]; o += rcounts[s]; end
    return out
end

"""
    lateral_prev_override(mesh, f, comm) -> Dict{(iel,d) => (lidx, vals)}

Collective (domain extents). Like lw_lateral_prev_override, for any inflow radiance f(x, y, z, θ, ϕ): turns periodic
faces into walls with that inflow (RadiativeCase.prev_override). Used for SW walls mode (:RT_sw_periodic_mode).
"""
function lateral_prev_override(mesh, f::Function, comm)
    np_ = mesh.npoin
    ext(v, op, d) = MPI.Allreduce(np_ > 0 ? op(view(v, 1:np_)) : d, op === minimum ? MPI.MIN : MPI.MAX, comm)
    xmin, xmax = ext(mesh.x, minimum, Inf), ext(mesh.x, maximum, -Inf)
    ymin, ymax = ext(mesh.y, minimum, Inf), ext(mesh.y, maximum, -Inf)
    tol = 1e-6 * max(xmax - xmin, ymax - ymin, 1.0)
    onlat(x, y) = abs(x - xmin) < tol || abs(x - xmax) < tol || abs(y - ymin) < tol || abs(y - ymax) < tol
    ngl = mesh.ngl
    out = Dict{Tuple{Int,Int},Tuple{Vector{Int},Vector{Float64}}}()
    for iel = 1:mesh.nelem
        lids = Int[]; pts = NTuple{3,Float64}[]
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel,i,j,k]
            onlat(mesh.x[ip], mesh.y[ip]) || continue
            push!(lids, lidx(i,j,k,ngl)); push!(pts, (mesh.x[ip], mesh.y[ip], mesh.z[ip]))
        end
        isempty(lids) && continue
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            out[(iel,d)] = (lids, [f(p[1], p[2], p[3], θ, ϕ) for p in pts])
        end
    end
    return out
end
