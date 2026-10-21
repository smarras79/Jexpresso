const USE_UPWIND_ORDER = Ref(true)

struct ElemFactor
    F::Union{LinearAlgebra.LU,Nothing}       # LU factors, or nothing when the matrix is triangular
    A::Union{Matrix{Float64},Nothing}        # kept only in the triangular case (A \ x, as the original)
    faces::Vector{Int}                       # inflow+active faces, increasing fid
    weights::Vector{Vector{Float64}}         # flux weight per face node, face_ip_range order
end

mutable struct ElemCache
    faceinfo::Vector{Vector{UInt16}}
    factors::Dict{Tuple{Int,Int},ElemFactor}
    budget_bytes::Float64
    bytes::Float64
    bytes_per_pair::Float64
    t_build::Float64
    ws_A::Matrix{Float64}
    ws_faces::Vector{Any}
    ws_br::Vector{Symbol}
    ord_iel::Vector{Int32}          # NEW: scan order of the (iel,d) pairs
    ord_d::Vector{Int32}
end

@inline mask_inflow(m::UInt16, fid::Int) = ((m >> (fid-1)) & UInt16(1)) == UInt16(1)
@inline mask_active(m::UInt16, fid::Int) = ((m >> (5+fid)) & UInt16(1)) == UInt16(1)

function make_elem_cache(mesh, metrics, n_dirs::Vector{Int}, ngl::Int; budget_mb::Float64=0.0)
    n3 = ngl^3
    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
    faceinfo = Vector{Vector{UInt16}}(undef, mesh.nelem)
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        masks = Vector{UInt16}(undef, n_dirs[iel])
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz)
            m = UInt16(0)
            for f = 1:6
                is_active[f] && (m |= UInt16(1) << (5+f))
                (is_inflow[f] && is_active[f]) && (m |= UInt16(1) << (f-1))
            end
            masks[d] = m
        end
        faceinfo[iel] = masks
    end

    # ---- scan order: natural (iel outer, d inner) or upstream-first by centroid·Ω ----
    npairs = sum(n_dirs)
    ord_iel = Vector{Int32}(undef, npairs); ord_d = Vector{Int32}(undef, npairs)
    upkeys = zeros(Float64, npairs)
    p = 0
    for iel = 1:mesh.nelem
        cx = 0.0; cy = 0.0; cz = 0.0
        if USE_UPWIND_ORDER[]
            for k=1:ngl, j=1:ngl, i=1:ngl
                ip = mesh.connijk[iel,i,j,k]; cx += mesh.x[ip]; cy += mesh.y[ip]; cz += mesh.z[ip]
            end
            cx /= n3; cy /= n3; cz /= n3
        end
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            p += 1
            ord_iel[p] = Int32(iel); ord_d[p] = Int32(d)
            if USE_UPWIND_ORDER[]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                upkeys[p] = cx*sin(θ)*cos(ϕ) + cy*sin(θ)*sin(ϕ) + cz*cos(θ)
            end
        end
    end
    if USE_UPWIND_ORDER[]
        perm = sortperm(upkeys; alg=MergeSort)                  # stable: ties keep the natural order
        ord_iel = ord_iel[perm]; ord_d = ord_d[perm]
    end

    bpp = 8.0*(n3*n3 + n3 + 6*ngl*ngl) + 200.0
    return ElemCache(faceinfo, Dict{Tuple{Int,Int},ElemFactor}(), budget_mb*1e6, 0.0, bpp, 0.0,
                     zeros(n3, n3), Vector{Any}(nothing, 6), fill(:none, 6), ord_iel, ord_d)
end

function elem_factor_build(mesh, metrics, ψ, dψ, ω, κ, iel::Int, d::Int, ngl::Int, m::UInt16)
    n3 = ngl^3; em = mesh.extra_mesh[iel]
    θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
    Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
    A = zeros(n3, n3); rhs = zeros(n3)                     # rhs only absorbs the zero-upwind flux terms
    build_dg_volume_matrix!(A, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
        metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz,
        metrics.Je, κ, mesh.connijk)
    faces = Int[]; weights = Vector{Vector{Float64}}()
    for fid = 1:6
        mask_active(m, fid) || continue
        inflow = mask_inflow(m, fid)
        fn = compute_face_normal_avg(iel, fid, ngl,
            metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz)
        w = zeros(ngl*ngl)
        build_dg_flux_terms!(A, rhs, fid, ngl, Ωx,Ωy,Ωz, inflow, true, fn, ω,
            metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz,
            metrics.Je, iel, (i,j,k)->0.0; weight_out=w)
        if inflow; push!(faces, fid); push!(weights, w); end
    end
    (istril(A) || istriu(A)) && return ElemFactor(nothing, A, faces, weights)
    return ElemFactor(lu!(A), nothing, faces, weights)
end

# Solves one ready pair. Upwind closures are in cache.ws_faces. Returns the new local solution vector.

function solve_pair!(cache::ElemCache, mesh, metrics, ψ, dψ, ω, κ, iel::Int, d::Int, ngl::Int, m::UInt16,
                     extra_source_all, prof)
    n3 = ngl^3; em = mesh.extra_mesh[iel]
    θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
    Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
    x = zeros(n3)
    t0 = gtic(prof)
    if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
        build_dg_source_term!(x, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
    end
    gtoc!(prof, :source_term, t0)

    ef = get(cache.factors, (iel,d), nothing)
    if ef === nothing && cache.bytes + cache.bytes_per_pair <= cache.budget_bytes
        t0 = time_ns()
        ef = elem_factor_build(mesh, metrics, ψ, dψ, ω, κ, iel, d, ngl, m)
        cache.factors[(iel,d)] = ef; cache.bytes += cache.bytes_per_pair
        cache.t_build += (time_ns() - t0) * 1e-9
        gtoc!(prof, :lu_build, t0)
    end

    if ef !== nothing                                      # cached path
        t0 = gtic(prof)
        for (q, fid) in enumerate(ef.faces)
            clo = cache.ws_faces[fid]; w = ef.weights[q]
            tf = gtic(prof)
            irange, jrange, krange = face_ip_range(fid, ngl)
            idx = 0
            for k in krange, j in jrange, i in irange
                idx += 1
                x[lidx(i,j,k,ngl)] += w[idx] * clo(i,j,k)
            end
            prof === nothing || padd!(prof, BRANCH_EVAL[cache.ws_br[fid]], (time_ns() - tf) * 1e-9)
        end
        gtoc!(prof, :flux_terms, t0)
        t0 = gtic(prof)
        if ef.F === nothing
            x = ef.A \ x
        else
            ldiv!(ef.F, x)
        end
        gtoc!(prof, :solve, t0)
        return x
    end

    A = cache.ws_A                                         # workspace path (budget exhausted or 0)
    t0 = gtic(prof)
    build_dg_volume_matrix!(A, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
        metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz,
        metrics.Je, κ, mesh.connijk)
    gtoc!(prof, :vol_matrix, t0)
    t0 = gtic(prof)
    for fid = 1:6
        mask_active(m, fid) || continue
        tf = gtic(prof)
        inflow = mask_inflow(m, fid)
        fn = compute_face_normal_avg(iel, fid, ngl,
            metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz)
        I_up = inflow ? cache.ws_faces[fid] : ((i,j,k)->0.0)
        build_dg_flux_terms!(A, x, fid, ngl, Ωx,Ωy,Ωz, inflow, true, fn, ω,
            metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz,
            metrics.Je, iel, I_up)
        prof === nothing || padd!(prof, BRANCH_EVAL[cache.ws_br[fid]], (time_ns() - tf) * 1e-9)
    end
    gtoc!(prof, :flux_terms, t0)
    t0 = gtic(prof)
    if istril(A) || istriu(A)
        x = A \ x
    else
        F = lu!(A); ldiv!(F, x)
    end
    gtoc!(prof, :solve, t0)
    return x
end

append!(PROF_SECTIONS, [:lu_build])

function elem_cache_report(cache::ElemCache, n_pairs::Int, comm)          # collective
    np = MPI.Comm_size(comm)
    v = Float64[length(cache.factors), n_pairs, cache.bytes, cache.t_build]
    s = MPI.Allreduce(v, +, comm); mx = MPI.Allreduce(v, MPI.MAX, comm)
    MPI.Comm_rank(comm) == 0 || return String[]
    pct = round(100*s[1]/max(s[2],1); digits=1)
    return String["elem cache: budget $(round(cache.budget_bytes/1e6; digits=0)) MB/rank, cached $(Int(s[1])) of $(Int(s[2])) pairs ($pct%), " *
        "memory mean $(round(s[3]/np/1e9; digits=2)) GB max $(round(mx[3]/1e9; digits=2)) GB per rank, " *
        "factor build time mean $(round(s[4]/np; digits=2)) s max $(round(mx[4]; digits=2)) s"]
end
