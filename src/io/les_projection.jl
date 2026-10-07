#=============================================================================
 les_projection.jl -- LES statistics on a user-defined structured plane.

 WHY
 ---
 horizontal_mean! and compute_xz_cross_section! average NODE values with equal
 weights and report one value per LGL node. LGL nodes cluster at element
 edges, so the output grid is non-uniform, depends on the mesh (and on AMR),
 and anything that varies with a node's position inside its element shows up
 as a "griddy" element-periodic pattern. Here the solution polynomial is
 instead EVALUATED at the points of a uniform structured grid, the statistics
 are formed there, and the average along the normal direction is taken over
 uniformly spaced samples. The output is mesh-independent and can be put
 directly on the protocol grid.

 INPUT (user_inputs.jl)
 ----------------------
     :les_projection => [
         (plane = "xz", npts = (513, 300), range1 = (0.0, 10240.0), range2 = (5.0, 2995.0)),
         (plane = "xy", npts = (512, 512), at = 100.0),
     ]

   plane   "xz", "xy" or "yz". The two letters are the plane axes, in order.
   npts    (n1, n2) points along the two plane axes, endpoints included.
   range1  (a, b) extent along the first axis;  default: the whole domain.
   range2  (a, b) extent along the second axis; default: the whole domain.
   at      position along the normal axis. Given: a SLICE, time average only.
           Omitted: the plane is AVERAGED along the normal axis.
   navg    number of samples along the normal axis when averaging; default
           the domain length over the first-axis spacing. The samples are cell
           centred, L*(k-1/2)/navg, so none sits on a periodic boundary.

 The quantities are exactly those of :lesprofile_vars / :lesstress_vars: each
 target point gets the conserved variables, the reference state and the SGS
 cache interpolated to it, and the deck's own user_les_profiles! and
 user_les_stress! are called on those values.

 OUTPUT (at the end of the statistics window, rank 0)
 ------
   les_proj_<k>_<plane>_tavg.vtr        the time-averaged plane, all variables
   les_proj_<k>_<plane>_profile_tavg.dat  for xz and yz: the plane further
                                        averaged over its first (horizontal) axis,
                                        i.e. the z profile, same columns as
                                        les_statistics_tavg.dat + les_stress_tavg.dat

 COST
 ----
 Setup locates every target point once (element bounding box, then Newton on
 the element map, so warped and AMR meshes work). Each statistics call is one
 tensor-product interpolation per point and NO communication; the sums are
 reduced to rank 0 once, in chunks, at the end.
=============================================================================#

mutable struct LESProjPlane
    name     ::String                  # "les_proj_1_xz"
    axes     ::NTuple{3,Int}           # (axis1, axis2, normal) as 1=x, 2=y, 3=z
    coord1   ::Vector{Float64}
    coord2   ::Vector{Float64}
    normal_at::Float64                 # NaN when averaging
    navg     ::Int                     # samples per (i1,i2) along the normal (1 for a slice)
    # located points owned by this rank
    elem     ::Vector{Int32}
    slot     ::Vector{Int32}           # local cell slot of each point
    Lx       ::Matrix{Float64}         # ngl × npts, Lagrange weights in ξ, η, ζ
    Ly       ::Matrix{Float64}
    Lz       ::Matrix{Float64}
    cells    ::Vector{Int32}           # global linear cell index of each local slot (sorted)
    count    ::Vector{Int32}           # points per local slot
    sum_mean ::Matrix{Float64}         # nslot × nprofiles
    sum_str  ::Matrix{Float64}         # nslot × nstress
    nsamples ::Int
end

const LES_PROJ = Ref{Any}(nothing)

# ---- 1D Lagrange basis on the LGL nodes at an arbitrary point ----
@inline function _lagrange!(L, ξn, x)
    n = length(ξn)
    @inbounds for j in 1:n
        v = 1.0
        for k in 1:n
            k == j && continue
            v *= (x - ξn[k]) / (ξn[j] - ξn[k])
        end
        L[j] = v
    end
    return L
end

@inline function _lagrange_deriv!(dL, ξn, x)
    n = length(ξn)
    @inbounds for j in 1:n
        s = 0.0
        for m in 1:n
            m == j && continue
            p = 1.0 / (ξn[j] - ξn[m])
            for k in 1:n
                (k == j || k == m) && continue
                p *= (x - ξn[k]) / (ξn[j] - ξn[k])
            end
            s += p
        end
        dL[j] = s
    end
    return dL
end

# Newton inverse of the element map x(ξ) = Σ La(ξ)Lb(η)Lc(ζ) X_abc.
# Returns (ok, ξ, η, ζ). X is ngl×ngl×ngl×3.
function _invert_element(X, p, ξn, ξ0, bufs)
    Lx, Ly, Lz, dLx, dLy, dLz = bufs
    ngl = length(ξn)
    ξ = ξ0[1]; η = ξ0[2]; ζ = ξ0[3]
    for it in 1:30
        _lagrange!(Lx, ξn, ξ); _lagrange!(Ly, ξn, η); _lagrange!(Lz, ξn, ζ)
        _lagrange_deriv!(dLx, ξn, ξ); _lagrange_deriv!(dLy, ξn, η); _lagrange_deriv!(dLz, ξn, ζ)
        f1 = f2 = f3 = 0.0
        J11 = J12 = J13 = J21 = J22 = J23 = J31 = J32 = J33 = 0.0
        @inbounds for c in 1:ngl, b in 1:ngl, a in 1:ngl
            w   = Lx[a]*Ly[b]*Lz[c]
            wξ  = dLx[a]*Ly[b]*Lz[c]
            wη  = Lx[a]*dLy[b]*Lz[c]
            wζ  = Lx[a]*Ly[b]*dLz[c]
            x1 = X[a,b,c,1]; x2 = X[a,b,c,2]; x3 = X[a,b,c,3]
            f1 += w*x1;  f2 += w*x2;  f3 += w*x3
            J11 += wξ*x1; J12 += wη*x1; J13 += wζ*x1
            J21 += wξ*x2; J22 += wη*x2; J23 += wζ*x2
            J31 += wξ*x3; J32 += wη*x3; J33 += wζ*x3
        end
        r1 = f1 - p[1]; r2 = f2 - p[2]; r3 = f3 - p[3]
        det = J11*(J22*J33 - J23*J32) - J12*(J21*J33 - J23*J31) + J13*(J21*J32 - J22*J31)
        abs(det) < 1e-300 && return (false, ξ, η, ζ)
        d1 = ( (J22*J33 - J23*J32)*r1 - (J12*J33 - J13*J32)*r2 + (J12*J23 - J13*J22)*r3) / det
        d2 = (-(J21*J33 - J23*J31)*r1 + (J11*J33 - J13*J31)*r2 - (J11*J23 - J13*J21)*r3) / det
        d3 = ( (J21*J32 - J22*J31)*r1 - (J11*J32 - J12*J31)*r2 + (J11*J22 - J12*J21)*r3) / det
        ξ = clamp(ξ - d1, -1.5, 1.5); η = clamp(η - d2, -1.5, 1.5); ζ = clamp(ζ - d3, -1.5, 1.5)
        abs(d1) + abs(d2) + abs(d3) < 1e-12 && return (true, ξ, η, ζ)
    end
    return (false, ξ, η, ζ)
end

_axis_index(c::Char) = c == 'x' ? 1 : c == 'y' ? 2 : c == 'z' ? 3 :
    error(":les_projection: plane axes must be x, y or z; got '$c'")

function _uniform(a, b, n)
    n == 1 && return [0.5*(a + b)]
    return collect(range(a, b; length = n))
end

"""
    build_les_projection(params) -> Vector{LESProjPlane} or nothing

Collective. Reads :les_projection, builds the target grids and locates every
target point in a local element. Called lazily on the first statistics call.
"""
function build_les_projection(params)
    specs = get(params.inputs, :les_projection, nothing)
    (specs === nothing || isempty(specs)) && return nothing
    specs isa NamedTuple && (specs = [specs])

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    mesh = params.mesh
    ngl  = Int(mesh.ngl)
    nelem = Int(mesh.nelem)
    connijk = Array(mesh.connijk)
    xyz = (Array(mesh.x), Array(mesh.y), Array(mesh.z))
    ξn  = Array(basis_structs_ξ_ω!(LGL(), ngl - 1, CPU()).ξ)

    lo = [MPI.Allreduce(minimum(xyz[d]), MPI.MIN, comm) for d in 1:3]
    hi = [MPI.Allreduce(maximum(xyz[d]), MPI.MAX, comm) for d in 1:3]
    Lper = hi .- lo

    nprof = length(params.inputs[:lesprofile_vars])
    nstr  = length(params.inputs[:lesstress_vars])

    planes = LESProjPlane[]
    X = zeros(ngl, ngl, ngl, 3)
    bufs = ntuple(_ -> zeros(ngl), 6)
    Lb = (zeros(ngl), zeros(ngl), zeros(ngl))
    tolr = 1e-8

    for (k, sp) in enumerate(specs)
        pl = lowercase(String(sp.plane))
        length(pl) == 2 || error(":les_projection: plane must be \"xz\", \"xy\" or \"yz\"; got \"$pl\"")
        a1 = _axis_index(pl[1]); a2 = _axis_index(pl[2])
        a1 == a2 && error(":les_projection: plane \"$pl\" repeats an axis")
        an = 6 - a1 - a2
        n1, n2 = Int(sp.npts[1]), Int(sp.npts[2])
        r1 = haskey(sp, :range1) ? Float64.(sp.range1) : (lo[a1], hi[a1])
        r2 = haskey(sp, :range2) ? Float64.(sp.range2) : (lo[a2], hi[a2])
        c1 = _uniform(r1[1], r1[2], n1)
        c2 = _uniform(r2[1], r2[2], n2)
        slice = haskey(sp, :at)
        at = slice ? Float64(sp.at) : NaN
        if slice
            cn = [at]
        else
            h1 = n1 > 1 ? (c1[end] - c1[1]) / (n1 - 1) : Lper[an]
            navg = haskey(sp, :navg) ? Int(sp.navg) : max(1, round(Int, Lper[an] / h1))
            cn = [lo[an] + Lper[an]*(m - 0.5)/navg for m in 1:navg]
        end
        coords = Vector{Vector{Float64}}(undef, 3)
        coords[a1] = c1; coords[a2] = c2; coords[an] = cn

        elem_v = Int32[]; cell_v = Int32[]
        Lxv = Float64[]; Lyv = Float64[]; Lzv = Float64[]

        for e in 1:nelem
            @inbounds for c in 1:ngl, b in 1:ngl, a in 1:ngl
                ip = connijk[e, a, b, c]
                X[a,b,c,1] = xyz[1][ip]; X[a,b,c,2] = xyz[2][ip]; X[a,b,c,3] = xyz[3][ip]
            end
            # periodic wrap: an element whose extent exceeds half the domain has a
            # face that the mesh stores at the other end -- unwrap it
            for d in 1:2
                emin = minimum(@view X[:,:,:,d]); emax = maximum(@view X[:,:,:,d])
                if emax - emin > 0.5*Lper[d]
                    Xd = @view X[:,:,:,d]
                    mid = 0.5*(lo[d] + hi[d])
                    Xd[Xd .< mid] .+= Lper[d]
                end
            end
            bmin = ntuple(d -> minimum(@view X[:,:,:,d]), 3)
            bmax = ntuple(d -> maximum(@view X[:,:,:,d]), 3)
            # candidate target indices inside the bounding box, per axis
            rng = ntuple(3) do d
                cd = coords[d]
                i0 = searchsortedfirst(cd, bmin[d] - 1e-9*max(1.0, abs(bmin[d])))
                i1 = searchsortedlast(cd,  bmax[d] + 1e-9*max(1.0, abs(bmax[d])))
                i0:i1
            end
            (isempty(rng[1]) || isempty(rng[2]) || isempty(rng[3])) && continue
            for i3 in rng[3], i2 in rng[2], i1 in rng[1]
                p = (coords[1][i1], coords[2][i2], coords[3][i3])
                g = ntuple(d -> 2*(p[d] - bmin[d])/max(bmax[d] - bmin[d], eps()) - 1, 3)
                ok, ξ, η, ζ = _invert_element(X, p, ξn, g, bufs)
                ok || continue
                # half-open ownership: a point on a shared face belongs to the
                # element on its high side, except on the domain's upper boundary
                own = true
                for (d, s) in ((1, ξ), (2, η), (3, ζ))
                    if s < -1 - tolr || s > 1 + tolr
                        own = false
                    elseif s > 1 - tolr && p[d] < hi[d] - 1e-9*max(1.0, abs(hi[d]))
                        own = false
                    end
                end
                own || continue
                idx = (i1, i2, i3)
                j1 = idx[a1]; j2 = idx[a2]
                push!(elem_v, Int32(e))
                push!(cell_v, Int32((j2 - 1)*n1 + j1))
                append!(Lxv, _lagrange!(Lb[1], ξn, clamp(ξ, -1.0, 1.0)))
                append!(Lyv, _lagrange!(Lb[2], ξn, clamp(η, -1.0, 1.0)))
                append!(Lzv, _lagrange!(Lb[3], ξn, clamp(ζ, -1.0, 1.0)))
            end
        end

        cells = sort(unique(cell_v))
        slotof = Dict{Int32,Int32}(c => Int32(i) for (i, c) in enumerate(cells))
        slot = Int32[slotof[c] for c in cell_v]
        count = zeros(Int32, length(cells))
        for s in slot; count[s] += 1; end

        npts = length(elem_v)
        pln = LESProjPlane("les_proj_$(k)_$(pl)", (a1, a2, an), c1, c2, at,
                           length(cn), elem_v, slot,
                           reshape(Lxv, ngl, npts), reshape(Lyv, ngl, npts), reshape(Lzv, ngl, npts),
                           cells, count,
                           zeros(length(cells), nprof), zeros(length(cells), nstr), 0)

        # every target point must be found exactly once
        ntot = MPI.Allreduce(npts, +, comm)
        nexp = n1 * n2 * length(cn)
        if rank == 0
            if ntot == nexp
                @info "$(pln.name): $n1 x $n2 plane, $(slice ? "slice at $("xyz"[an]) = $at" : "averaged along $("xyz"[an]) over $(length(cn)) samples"), all $nexp points located"
            else
                @warn "$(pln.name): located $ntot of $nexp target points -- points outside the mesh or on an unowned boundary are missing"
            end
        end
        push!(planes, pln)
    end
    return planes
end

"""
    les_projection_accumulate!(params)

Called from les_statistics after uaux and the SGS cache are filled. Local
only, no communication.
"""
function les_projection_accumulate!(params)
    if LES_PROJ[] === nothing
        LES_PROJ[] = something(build_les_projection(params), false)
    end
    LES_PROJ[] === false && return

    uaux = params.uaux
    qe   = params.qp.qe
    sgs  = params.sgs_stress
    ET   = params.SOL_VARS_TYPE
    connijk = params.mesh.connijk
    ngl  = Int(params.mesh.ngl)
    nu, nq, ns = size(uaux, 2), size(qe, 2), size(sgs, 2)
    qp  = zeros(nu); qep = zeros(nq); sp = zeros(ns)
    nprof = length(params.inputs[:lesprofile_vars])
    nstr  = length(params.inputs[:lesstress_vars])
    mbuf = zeros(nprof); sbuf = zeros(nstr)

    # LES_PROJ is a Ref{Any}: `pl` is not inferable here, so the hot loop is
    # behind a function barrier. Without it every access in the loop is dynamic
    # and one statistics call cost ~40 s at 1024 ranks instead of ~0.2 s.
    for pl in LES_PROJ[]
        t0 = time()
        _proj_accumulate_plane!(pl::LESProjPlane, uaux, qe, sgs, connijk, ngl,
                                qp, qep, sp, mbuf, sbuf, ET)
        # timing of the first few calls on rank 0, so the cost is on record
        if pl.nsamples <= 3 && MPI.Comm_rank(get_mpi_comm()) == 0
            @printf(" # %s: statistics call %d took %.3f s on rank 0 (%d points)\n",
                    pl.name, pl.nsamples, time() - t0, length(pl.elem))
            flush(stdout)
        end
    end
end

function _proj_accumulate_plane!(pl::LESProjPlane, uaux, qe, sgs, connijk, ngl::Int,
                                 qp, qep, sp, mbuf, sbuf, ET)
    nu, nq, ns = length(qp), length(qep), length(sp)
    nprof, nstr = length(mbuf), length(sbuf)
    elem = pl.elem; slot = pl.slot
    Lx = pl.Lx; Ly = pl.Ly; Lz = pl.Lz
    sum_mean = pl.sum_mean; sum_str = pl.sum_str
    @inbounds for p in eachindex(elem)
        e = elem[p]
        fill!(qp, 0.0); fill!(qep, 0.0); fill!(sp, 0.0)
        for c in 1:ngl
            wz = Lz[c, p]
            for b in 1:ngl
                wyz = Ly[b, p] * wz
                for a in 1:ngl
                    w  = Lx[a, p] * wyz
                    ip = connijk[e, a, b, c]
                    for v in 1:nu; qp[v]  += w * uaux[ip, v]; end
                    for v in 1:nq; qep[v] += w * qe[ip, v];   end
                    for v in 1:ns; sp[v]  += w * sgs[ip, v];  end
                end
            end
        end
        user_les_profiles!(mbuf, sbuf, qp, qep, sp, ET)
        s = slot[p]
        for v in 1:nprof; sum_mean[s, v] += mbuf[v]; end
        for v in 1:nstr;  sum_str[s, v]  += sbuf[v]; end
    end
    pl.nsamples += 1
    return nothing
end

"""
    les_projection_finalize!(params, t)

Collective. Reduces the sums to rank 0 in chunks (so no rank ever holds a
second full-size copy), applies user_les_stress! with the time-and-space mean,
and writes the plane VTK and, for vertical planes, the z profile.
"""
function les_projection_finalize!(params, t)
    (LES_PROJ[] === nothing || LES_PROJ[] === false) && return
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    pvars = params.inputs[:lesprofile_vars]
    svars = params.inputs[:lesstress_vars]
    nprof = length(pvars); nstr = length(svars)
    nf    = 1 + nprof + nstr                    # count, means, raw products
    outdir = params.inputs[:output_dir]

    for pl in LES_PROJ[]
        ns = MPI.Allreduce(pl.nsamples, MPI.MAX, comm)
        ns == 0 && continue
        n1, n2 = length(pl.coord1), length(pl.coord2)
        ncell  = n1 * n2
        G = rank == 0 ? zeros(ncell, nf) : zeros(0, 0)
        chunk = max(1, 1_000_000 ÷ nf)
        buf = zeros(chunk * nf)
        ptr = 1
        for c0 in 1:chunk:ncell
            c1 = min(ncell, c0 + chunk - 1); nc = c1 - c0 + 1
            fill!(buf, 0.0)
            while ptr <= length(pl.cells) && pl.cells[ptr] <= c1
                cl = pl.cells[ptr] - c0 + 1
                buf[cl] += pl.count[ptr]
                for v in 1:nprof; buf[(v)*nc + cl]        += pl.sum_mean[ptr, v]; end
                for v in 1:nstr;  buf[(nprof + v)*nc + cl] += pl.sum_str[ptr, v];  end
                ptr += 1
            end
            view_buf = @view buf[1:nc*nf]
            MPI.Reduce!(view_buf, +, 0, comm)
            if rank == 0
                G[c0:c1, :] .= reshape(view_buf, nc, nf)
            end
        end
        rank == 0 || continue

        cnt = G[:, 1]
        any(cnt .== 0) && @info "$(pl.name): $(count(==(0), cnt)) cells received no samples (outside the mesh, e.g. under terrain) -- written as NaN"
        den = max.(cnt, 1) .* ns                   # count was per sample
        means = G[:, 2:1+nprof] ./ den
        raw   = G[:, 2+nprof:end] ./ den
        stress = zeros(ncell, nstr)
        pbuf = zeros(nstr)
        for i in 1:ncell
            user_les_stress!(pbuf, @view(raw[i, :]), @view(means[i, :]))
            stress[i, :] .= pbuf
        end
        # Cells no rank sampled lie outside the mesh -- under the terrain on a
        # warped mesh. They have no value: NaN, as the official TABLES
        # interpolation writes below the surface, not a misleading 0.
        empty = cnt .== 0
        means[empty, :]  .= NaN
        stress[empty, :] .= NaN

        # ---- plane: rectilinear VTK, the normal axis a single coordinate ----
        a1, a2, an = pl.axes
        cs = Vector{Vector{Float64}}(undef, 3)
        cs[a1] = pl.coord1; cs[a2] = pl.coord2
        cs[an] = [isnan(pl.normal_at) ? 0.0 : pl.normal_at]
        dims = (length(cs[1]), length(cs[2]), length(cs[3]))
        # map (i1,i2) cell order onto (x,y,z) array order
        function toarr(col)
            A = zeros(dims)
            for j2 in 1:n2, j1 in 1:n1
                idx = [1, 1, 1]; idx[a1] = j1; idx[a2] = j2
                A[idx...] = col[(j2 - 1)*n1 + j1]
            end
            return A
        end
        vtk = vtk_grid(joinpath(outdir, "$(pl.name)_tavg"), cs[1], cs[2], cs[3])
        for v in 1:nprof; vtk[pvars[v], VTKPointData()] = toarr(means[:, v]);  end
        for v in 1:nstr;  vtk[svars[v], VTKPointData()] = toarr(stress[:, v]); end
        vtk["n_samples_per_point", VTKPointData()] = toarr(den)
        vtk_save(vtk)

        # ---- vertical plane: z profile over the first (horizontal) axis ----
        if a2 == 3
            prof_file = joinpath(outdir, "$(pl.name)_profile_tavg.dat")
            open(prof_file, "w") do io
                print(io, "# time_end=", @sprintf("%.6e", t), "  n_samples=", ns,
                      "  plane=$(pl.name) averaged over its first axis  z")
                for v in 1:nprof; print(io, "  ", pvars[v]); end
                for v in 1:nstr;  print(io, "  ", svars[v]); end
                println(io)
                mz = zeros(nprof); rz = zeros(nstr); sz = zeros(nstr)
                for j2 in 1:n2
                    fill!(mz, 0.0); fill!(rz, 0.0); dsum = 0.0
                    for j1 in 1:n1
                        i = (j2 - 1)*n1 + j1
                        empty[i] && continue           # below the terrain
                        mz .+= @view G[i, 2:1+nprof]
                        rz .+= @view G[i, 2+nprof:end]
                        dsum += den[i]
                    end
                    if dsum > 0
                        mz ./= dsum; rz ./= dsum
                        user_les_stress!(sz, rz, mz)
                    else
                        fill!(mz, NaN); fill!(sz, NaN)
                    end
                    @printf(io, "%.6e", pl.coord2[j2])
                    for v in 1:nprof; @printf(io, "  %.6e", mz[v]); end
                    for v in 1:nstr;  @printf(io, "  %.6e", sz[v]); end
                    println(io)
                end
            end
        end
        @info "$(pl.name): wrote $(pl.name)_tavg.vtr$(a2 == 3 ? " and $(pl.name)_profile_tavg.dat" : "") ($ns samples)"
    end
end
