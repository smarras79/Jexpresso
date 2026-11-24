# ── Shortwave direct beam: slant-path optical depth by polynomial ray tracing through the elements ───────────────────
#
# F_dir(x) = S₀ μ₀ exp(−τ(x)/μ₀),  τ(x) = ∫_{z}^{z_top} κ_ext(x_h + (z′−z) s, z′) dz′,  s = −tanθ₀ (cosφ₀, sinφ₀)
#
# (photons travel along Ω_sun = (sinθ cosφ₀, sinθ sinφ₀, cosθ), θ = π − acos μ₀, so the ray from a node TOWARD the sun
# moves −tanθ₀ (cosφ₀, sinφ₀) horizontally per metre of height). Every node's ray is traced element by element: in each
# element κ_ext is the element's own tensor Lagrange interpolant of its nodal values (the same polynomial the SEM
# uses), the segment between entry and exit is integrated with 7-point Gauss–Legendre (exact for the degree ≤ 3(ngl−1)
# integrand when ngl ≤ 5; segments in steep elements (nodal max/min > 4) or where a Gauss point dips below 0 are
# integrated on 64 pieces), and the ray moves on to the element containing the next point. Elements are located by
# point (a uniform x–y bucket grid, z-sorted), never through face connectivity, so conforming, non-conforming
# (hanging-node / p4est-adapted) and periodic meshes are all handled the same way. The integrand is max(κ_poly, 0):
# near a sharp cloud edge the degree-(ngl−1) interpolant can undershoot below zero between nodes.
#
# Requirements: hexahedral elements that are axis-aligned boxes (flat-bottom box domains: BOMEX, bubble, p4est leaves of
# a box); an error is raised otherwise. Every rank gathers every element's box, 1D node coordinates and nodal κ_ext
# (Allgatherv, ~(6 + 3 ngl + ngl³)·8 bytes per element), so memory scales with the global element count.
# Lateral wrap is periodic along axes tagged periodicx / periodicy, else the ray stays on the domain edge (the
# atmosphere outside is taken to be the edge column). COLLECTIVE.

struct SWBeamElements
    lo   :: Matrix{Float64}        # 3 × ne  box lower corner
    hi   :: Matrix{Float64}        # 3 × ne  box upper corner
    xn   :: Array{Float64,3}       # ngl × 3 × ne  1D node coordinates per axis
    kap  :: Array{Float64,4}       # ngl × ngl × ngl × ne  nodal κ_ext (i ↔ x, j ↔ y, k ↔ z)
    # x–y bucket grid: elements overlapping each cell, sorted by lo z
    bx0  :: Float64; by0 :: Float64; bdx :: Float64; bdy :: Float64; nbx :: Int; nby :: Int
    cell :: Vector{Vector{Int}}
    steep :: BitVector             # nodal κ_ext max/min > 4 in the element (cloud edge): always refine segments
    dom_lo :: NTuple{3,Float64}; dom_hi :: NTuple{3,Float64}
    bw :: Vector{Float64}          # barycentric weights of the reference nodes on [-1, 1] (box elements: affine map)
end

# reference-node barycentric weights 1/Π_{m≠i}(ξ_i − ξ_m), from one element's 1D nodes mapped to [-1, 1]
function sw_beam_ref_weights(xn1::AbstractVector{Float64})
    n = length(xn1); ξ = 2 .* (xn1 .- xn1[1]) ./ (xn1[end] - xn1[1]) .- 1
    return [1 / prod(ξ[i] - ξ[m] for m = 1:n if m != i) for i = 1:n]
end

# 7-point Gauss–Legendre on [-1, 1] (Golub–Welsch)
const SW_BEAM_GL = let n = 7
    β = [k / sqrt(4k^2 - 1) for k = 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), β))
    (E.values, 2 .* E.vectors[1, :] .^ 2)
end

function sw_beam_gather_elements(mesh, κ_ext::Vector{Float64}, comm)
    ngl = mesh.ngl; ne = mesh.nelem
    rec = 6 + 3ngl + ngl^3
    buf = Vector{Float64}(undef, rec * ne)
    tolrel = 1e-6
    for iel = 1:ne
        o = (iel - 1) * rec
        ip(i, j, k) = mesh.connijk[iel, i, j, k]
        xs = [mesh.x[ip(i, 1, 1)] for i = 1:ngl]; ys = [mesh.y[ip(1, j, 1)] for j = 1:ngl]; zs = [mesh.z[ip(1, 1, k)] for k = 1:ngl]
        # axis-aligned box check: x depends only on i, y only on j, z only on k
        h = max(xs[end] - xs[1], ys[end] - ys[1], zs[end] - zs[1])
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            p = ip(i, j, k)
            if abs(mesh.x[p] - xs[i]) > tolrel * h || abs(mesh.y[p] - ys[j]) > tolrel * h || abs(mesh.z[p] - zs[k]) > tolrel * h
                error("compute_direct_beam_raytrace: element $iel is not an axis-aligned box " *
                      "(only box elements are supported); use :RT_sw_direct_beam => :reference_column")
            end
        end
        (xs[end] > xs[1] && ys[end] > ys[1] && zs[end] > zs[1]) ||
            error("compute_direct_beam_raytrace: element $iel has reversed node ordering")
        buf[o+1:o+3] .= (xs[1], ys[1], zs[1]); buf[o+4:o+6] .= (xs[end], ys[end], zs[end])
        buf[o+7:o+6+ngl] .= xs; buf[o+7+ngl:o+6+2ngl] .= ys; buf[o+7+2ngl:o+6+3ngl] .= zs
        q = o + 6 + 3ngl
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            q += 1; buf[q] = κ_ext[ip(i, j, k)]
        end
    end
    cnt = MPI.Allgather(Int32(length(buf)), comm)
    all = MPI.Allgatherv(buf, cnt, comm)
    neg = length(all) ÷ rec
    lo = Matrix{Float64}(undef, 3, neg); hi = similar(lo)
    xn = Array{Float64,3}(undef, ngl, 3, neg); kap = Array{Float64,4}(undef, ngl, ngl, ngl, neg)
    for e = 1:neg
        o = (e - 1) * rec
        lo[:, e] .= all[o+1:o+3]; hi[:, e] .= all[o+4:o+6]
        for a = 1:3; xn[:, a, e] .= all[o+7+(a-1)*ngl:o+6+a*ngl]; end
        kap[:, :, :, e] .= reshape(all[o+7+3ngl:o+6+3ngl+ngl^3], ngl, ngl, ngl)
    end
    dlo = (minimum(lo[1, :]), minimum(lo[2, :]), minimum(lo[3, :]))
    dhi = (maximum(hi[1, :]), maximum(hi[2, :]), maximum(hi[3, :]))
    # bucket grid: cell size = smallest element width per axis
    bdx = minimum(hi[1, :] .- lo[1, :]); bdy = minimum(hi[2, :] .- lo[2, :])
    nbx = max(1, round(Int, (dhi[1] - dlo[1]) / bdx)); nby = max(1, round(Int, (dhi[2] - dlo[2]) / bdy))
    bdx = (dhi[1] - dlo[1]) / nbx; bdy = (dhi[2] - dlo[2]) / nby
    cell = [Int[] for _ = 1:nbx*nby]
    ϵ = 1e-9 * max(dhi[1] - dlo[1], dhi[2] - dlo[2])
    for e = 1:neg
        i0 = clamp(floor(Int, (lo[1, e] - dlo[1] + ϵ) / bdx) + 1, 1, nbx); i1 = clamp(ceil(Int, (hi[1, e] - dlo[1] - ϵ) / bdx), 1, nbx)
        j0 = clamp(floor(Int, (lo[2, e] - dlo[2] + ϵ) / bdy) + 1, 1, nby); j1 = clamp(ceil(Int, (hi[2, e] - dlo[2] - ϵ) / bdy), 1, nby)
        for j = j0:j1, i = i0:i1; push!(cell[(j-1)*nbx+i], e); end
    end
    for c in cell; sort!(c; by = e -> lo[3, e]); end
    steep = BitVector([maximum(view(kap, :, :, :, e)) > 4 * max(minimum(view(kap, :, :, :, e)), 1e-300) for e = 1:neg])
    return SWBeamElements(lo, hi, xn, kap, dlo[1], dlo[2], bdx, bdy, nbx, nby, cell, steep, dlo, dhi,
                          sw_beam_ref_weights(xn[:, 1, 1]))
end

# element containing point p (first match; points on shared faces may belong to either neighbour)
function sw_beam_locate(B::SWBeamElements, p)
    i = clamp(floor(Int, (p[1] - B.bx0) / B.bdx) + 1, 1, B.nbx); j = clamp(floor(Int, (p[2] - B.by0) / B.bdy) + 1, 1, B.nby)
    c = B.cell[(j-1)*B.nbx+i]
    tz = 1e-9 * (B.dom_hi[3] - B.dom_lo[3]); tx = 1e-9 * (B.dom_hi[1] - B.dom_lo[1]); ty = 1e-9 * (B.dom_hi[2] - B.dom_lo[2])
    # first element whose lo z exceeds p z: candidates are before it
    n = searchsortedlast(c, p[3] + tz; by = e -> e isa Int ? B.lo[3, e] : e)
    for q = n:-1:1
        e = c[q]
        B.hi[3, e] + tz < p[3] && continue
        (B.lo[1, e] - tx <= p[1] <= B.hi[1, e] + tx && B.lo[2, e] - ty <= p[2] <= B.hi[2, e] + ty) && return e
    end
    return 0
end

# Lagrange basis values at x on element e's axis a: l_i = w_i Π_{m≠i}(ξ − ξ_m) in the reference coordinate ξ, with
# the products from prefix/suffix sweeps (no divisions; nodes are an affine image of the reference nodes)
@inline function sw_beam_lagrange(xn, a, e, x, ngl, out, bw)
    @inbounds begin
        x1 = xn[1, a, e]; sc = 2 / (xn[ngl, a, e] - x1)
        p = 1.0
        for i = 1:ngl                                   # prefix products Π_{m<i}
            out[i] = p
            p *= (x - xn[i, a, e]) * sc
        end
        p = 1.0
        for i = ngl:-1:1                                # times suffix products Π_{m>i}, times w_i
            out[i] *= p * bw[i]
            p *= (x - xn[i, a, e]) * sc
        end
    end
    return out
end

function sw_beam_kappa(B::SWBeamElements, e, p, ngl, lx, ly, lz)
    sw_beam_lagrange(B.xn, 1, e, p[1], ngl, lx, B.bw); sw_beam_lagrange(B.xn, 2, e, p[2], ngl, ly, B.bw)
    sw_beam_lagrange(B.xn, 3, e, p[3], ngl, lz, B.bw)
    s = 0.0
    @inbounds for k = 1:ngl, j = 1:ngl
        t = ly[j] * lz[k]
        for i = 1:ngl; s += B.kap[i, j, k, e] * lx[i] * t; end
    end
    return s                       # raw polynomial value (may undershoot below 0 near sharp cloud edges)
end

# ∫ max(κ_poly, 0) dz over [za, zb] in element e along the ray q + s (z − z0): 7-point Gauss on `nsub` equal pieces.
# Returns (integral, clipped) where clipped = some Gauss point had κ_poly < 0.
function sw_beam_segment(B, e, q, s, z0, za, zb, nsub, ngl, lx, ly, lz)
    ξg, wg = SW_BEAM_GL
    acc = 0.0; clipped = false
    h = (zb - za) / nsub
    for m = 1:nsub
        a = za + (m - 1) * h; hz = 0.5h; cz = a + hz
        for g = 1:length(ξg)
            zg = cz + hz * ξg[g]
            v = sw_beam_kappa(B, e, (q[1] + s[1] * (zg - z0), q[2] + s[2] * (zg - z0), zg), ngl, lx, ly, lz)
            v < 0 && (clipped = true; v = 0.0)
            acc += wg[g] * v * hz
        end
    end
    return acc, clipped
end

# ── Exact segment integral ∫ max(κ_poly, 0) dz ────────────────────────────────────────────────────────────────────────
# Along a straight ray inside one element, κ_poly is a single polynomial of degree ≤ 3(ngl−1) in the ray parameter.
# Its values at the n+1 Chebyshev–Lobatto points (n = 3(ngl−1)) give its Chebyshev coefficients exactly; if
# c₀ − Σ|c_k| ≥ 0 it is non-negative on the segment and the integral is the exact Chebyshev integral. Otherwise the
# sign is scanned on 4n sub-intervals (cheap Clenshaw evaluations), roots are refined by bisection, and only the
# positive parts are integrated exactly. Replaces 64 pieces × 7 Gauss points (448 evaluations of the 3D element
# polynomial) by n+1 = 13 evaluations at ngl = 5. :RT_sw_beam_exact => false restores the 64-piece rule.
const SW_BEAM_EXACT = Ref(true)

struct SWBeamCheb
    n :: Int
    t :: Vector{Float64}            # Chebyshev–Lobatto points cos(πj/n), j = 0..n
    M :: Matrix{Float64}            # values at t → Chebyshev coefficients (p = Σ c_k T_k)
    f :: Vector{Float64}; c :: Vector{Float64}; T :: Vector{Float64}
end
function SWBeamCheb(ngl::Int)
    n = 3 * (ngl - 1)
    M = zeros(n + 1, n + 1)
    for k = 0:n, j = 0:n
        M[k+1, j+1] = (2 / n) * ((j == 0 || j == n) ? 0.5 : 1.0) * cos(π * j * k / n)
    end
    M[1, :] .*= 0.5; M[n+1, :] .*= 0.5
    return SWBeamCheb(n, [cos(π * j / n) for j = 0:n], M, zeros(n + 1), zeros(n + 1), zeros(n + 2))
end

@inline function sw_cheb_eval(C::SWBeamCheb, t)
    c = C.c; b1 = 0.0; b2 = 0.0
    @inbounds for k = C.n:-1:1
        b1, b2 = 2t * b1 - b2 + c[k+1], b1
    end
    return t * b1 - b2 + c[1]
end

# antiderivative of Σ c_k T_k at t:  ∫T₀ = t, ∫T₁ = t²/2, ∫T_k = (T_{k+1}/(k+1) − T_{k−1}/(k−1))/2
@inline function sw_cheb_prim(C::SWBeamCheb, t)
    c = C.c; T = C.T; n = C.n
    @inbounds begin
        T[1] = 1.0; T[2] = t
        for k = 2:n+1; T[k+1] = 2t * T[k] - T[k-1]; end
        s = c[1] * t + (n >= 1 ? c[2] * t * t / 2 : 0.0)
        for k = 2:n
            s += c[k+1] * (T[k+2] / (k + 1) - T[k] / (k - 1)) / 2
        end
    end
    return s
end

function sw_beam_segment_exact(B, e, q, s, z0, za, zb, ngl, lx, ly, lz, C::SWBeamCheb)
    hz = 0.5 * (zb - za); cz = 0.5 * (za + zb); n = C.n
    @inbounds for j = 0:n
        zg = cz + hz * C.t[j+1]
        C.f[j+1] = sw_beam_kappa(B, e, (q[1] + s[1] * (zg - z0), q[2] + s[2] * (zg - z0), zg), ngl, lx, ly, lz)
    end
    mul!(C.c, C.M, C.f)
    bound = C.c[1]
    @inbounds for k = 2:n+1; bound -= abs(C.c[k]); end
    bound >= 0 && return hz * (sw_cheb_prim(C, 1.0) - sw_cheb_prim(C, -1.0)), false
    m = 4n; acc = 0.0
    a = -1.0; pa = sw_cheb_eval(C, a); Fa = sw_cheb_prim(C, a)
    for i = 1:m
        b = -1.0 + 2.0 * i / m; pb = sw_cheb_eval(C, b); Fb = sw_cheb_prim(C, b)
        if pa >= 0 && pb >= 0
            acc += Fb - Fa
        elseif (pa >= 0) != (pb >= 0)
            lo = a; hi = b; plo = pa                       # bisection for the sign change
            for _ = 1:52
                mid = 0.5 * (lo + hi); pm = sw_cheb_eval(C, mid)
                ((pm >= 0) == (plo >= 0)) ? (lo = mid; plo = pm) : (hi = mid)
                hi - lo < 1e-14 && break
            end
            r = 0.5 * (lo + hi); Fr = sw_cheb_prim(C, r)
            acc += pa >= 0 ? Fr - Fa : Fb - Fr
        end
        a = b; pa = pb; Fa = Fb
    end
    return hz * max(acc, 0.0), true
end

# Ray context for τ at arbitrary points (nodes, or quadrature points for the projected beam source)
mutable struct SWBeamRay
    B :: SWBeamElements
    per :: Tuple{Bool,Bool}; L :: NTuple{2,Float64}; ztop :: Float64; δ :: Float64; s0 :: NTuple{2,Float64}
    ngl :: Int; lx :: Vector{Float64}; ly :: Vector{Float64}; lz :: Vector{Float64}
    nfail :: Int; nclip :: Int; nseg_max :: Int
    cheb :: SWBeamCheb
end

function sw_beam_ray(mesh, ext::Vector{Float64}, sw::SWParams, comm)        # COLLECTIVE
    B = sw_beam_gather_elements(mesh, ext, comm)
    pax = detect_periodic_axes(mesh)
    pf = MPI.Allreduce(Int[1 in pax, 2 in pax], MPI.MAX, comm)
    μ0 = sw.μ₀; tanθ = sqrt(max(1 - μ0^2, 0.0)) / max(μ0, 1e-12)
    ngl = mesh.ngl
    return SWBeamRay(B, (pf[1] == 1, pf[2] == 1), (B.dom_hi[1] - B.dom_lo[1], B.dom_hi[2] - B.dom_lo[2]), B.dom_hi[3],
                     1e-7 * (B.dom_hi[3] - B.dom_lo[3]), (-tanθ * cos(sw.φ₀), -tanθ * sin(sw.φ₀)),
                     ngl, zeros(ngl), zeros(ngl), zeros(ngl), 0, 0, 0, SWBeamCheb(ngl))
end

# Vertical-path-integrated extinction τ from (x, y, z) to the model top along the ray toward the sun.
function sw_beam_tau(R::SWBeamRay, x::Float64, y::Float64, z::Float64)
    B = R.B; per = R.per; L = R.L; δ = R.δ; ngl = R.ngl
    wrap(v, a) = per[a] ? B.dom_lo[a] + mod(v - B.dom_lo[a], L[a]) : clamp(v, B.dom_lo[a], B.dom_hi[a])
    q = [x, y, z]; s = [R.s0[1], R.s0[2]]
    τ = 0.0; nseg = 0
    while z < R.ztop - δ
        # stop drifting along a non-periodic axis once the ray reaches the domain edge
        for a = 1:2
            if !per[a] && ((s[a] < 0 && q[a] <= B.dom_lo[a] + δ) || (s[a] > 0 && q[a] >= B.dom_hi[a] - δ)); s[a] = 0.0; end
        end
        probe = (wrap(q[1] + s[1] * δ, 1), wrap(q[2] + s[2] * δ, 2), z + δ)
        e = sw_beam_locate(B, probe)
        if e == 0; R.nfail += 1; break; end
        # express q in the periodic image of e (q may sit on the far side of a periodic boundary)
        for a = 1:2
            if per[a]
                q[a] < B.lo[a, e] - δ && (q[a] += L[a])
                q[a] > B.hi[a, e] + δ && (q[a] -= L[a])
            end
        end
        # exit height from the element box along the ray (parametrised by z)
        zb = min(B.hi[3, e], R.ztop)
        for a = 1:2
            s[a] > 0 && (zb = min(zb, z + (B.hi[a, e] - q[a]) / s[a]))
            s[a] < 0 && (zb = min(zb, z + (B.lo[a, e] - q[a]) / s[a]))
        end
        zb = max(zb, z + δ)                     # guard: always advance
        # 7-point Gauss is exact for the polynomial. Where it can dip below 0 (clipping makes kinks) integrate on 64
        # pieces: elements with a steep nodal κ_ext (max/min > 4: cloud edges; dips can fall between Gauss points)
        # and any segment where a Gauss point was clipped
        if SW_BEAM_EXACT[]
            # exact positive-part integral (steep elements directly; others when a Gauss point was clipped)
            acc, clipped = B.steep[e] ? sw_beam_segment_exact(B, e, q, s, z, z, zb, ngl, R.lx, R.ly, R.lz, R.cheb) :
                                        sw_beam_segment(B, e, q, s, z, z, zb, 1, ngl, R.lx, R.ly, R.lz)
            if clipped && !B.steep[e]
                acc, clipped = sw_beam_segment_exact(B, e, q, s, z, z, zb, ngl, R.lx, R.ly, R.lz, R.cheb)
            end
            clipped && (R.nclip += 1)
        else
            acc, clipped = B.steep[e] ? (0.0, true) : sw_beam_segment(B, e, q, s, z, z, zb, 1, ngl, R.lx, R.ly, R.lz)
            if clipped
                acc, _ = sw_beam_segment(B, e, q, s, z, z, zb, 64, ngl, R.lx, R.ly, R.lz); R.nclip += 1
            end
        end
        τ += acc
        q[1] = wrap(q[1] + s[1] * (zb - z), 1); q[2] = wrap(q[2] + s[2] * (zb - z), 2)
        z = zb; nseg += 1
        nseg > 100_000 && error("sw_beam_tau: ray from ($x, $y) did not reach the top")
    end
    R.nseg_max = max(R.nseg_max, nseg)
    return τ
end

const SW_BEAM_GQ3 = ([-sqrt(3 / 5), 0.0, sqrt(3 / 5)], [5 / 9, 8 / 9, 5 / 9])     # 3-point Gauss on [-1, 1]

"""
    sw_beam_projected_source(mesh, R, σ_sw, F_dir, sw; nsub=4, steep_ratio=2.0) -> (C, info)

Element-local nodal values C[l, iel] (l = lidx(i,j,k)) of the scattered-beam source coefficient c = σ F_dir / μ₀, such
that the scheme's lumped (GLL-collocated) load vector  ω_i ω_j ω_k J c_l  equals the EXACT load ∫ φ_l c dV. Elements
where c varies steeply (nodal F_dir or σ max/min > steep_ratio) are integrated with composite 3-point Gauss on nsub³
sub-cells, F_dir traced at every quadrature point; the rest keep the nodal values (collocation is accurate there).
Then a Zhang–Shu positivity limiter scales the nodal values toward their (exact, ≥ 0) mean where any is negative.
Requires box elements (as the ray tracer). Not collective.
"""
function sw_beam_projected_source(mesh, R::SWBeamRay, σ_sw::Vector{Float64}, F_dir::Vector{Float64}, sw::SWParams;
                                  nsub::Int=4, steep_ratio::Float64=2.0)
    ngl = mesh.ngl; n3 = ngl^3; μ0 = sw.μ₀
    C = zeros(n3, mesh.nelem)
    ω = mesh.ngl == length(R.lx) ? gll_weights_1d(ngl) : error("ngl mismatch")
    xg, wg = SW_BEAM_GQ3
    xs = zeros(ngl); ys = zeros(ngl); zs = zeros(ngl); lx = zeros(ngl); ly = zeros(ngl); lz = zeros(ngl)
    lag!(out, nodes, x) = (for i = 1:ngl; l = 1.0; for m = 1:ngl; m == i && continue; l *= (x - nodes[m]) / (nodes[i] - nodes[m]); end; out[i] = l; end; out)
    nproj = 0; nlim = 0; moved = 0.0; total = 0.0
    load = zeros(n3)
    for iel = 1:mesh.nelem
        for i = 1:ngl
            xs[i] = mesh.x[mesh.connijk[iel, i, 1, 1]]; ys[i] = mesh.y[mesh.connijk[iel, 1, i, 1]]
            zs[i] = mesh.z[mesh.connijk[iel, 1, 1, i]]
        end
        Fmn = Inf; Fmx = 0.0; smn = Inf; smx = 0.0
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            g = mesh.connijk[iel, i, j, k]; l = lidx(i, j, k, ngl)
            C[l, iel] = μ0 > 0 ? σ_sw[g] * F_dir[g] / μ0 : 0.0
            Fmn = min(Fmn, F_dir[g]); Fmx = max(Fmx, F_dir[g]); smn = min(smn, σ_sw[g]); smx = max(smx, σ_sw[g])
        end
        steep = μ0 > 0 && smx > 0 && (Fmx > steep_ratio * max(Fmn, 1e-300) || smx > steep_ratio * max(smn, 1e-300))
        V = (xs[end] - xs[1]) * (ys[end] - ys[1]) * (zs[end] - zs[1])
        if steep
            nproj += 1; fill!(load, 0.0)
            hx = (xs[end] - xs[1]) / nsub; hy = (ys[end] - ys[1]) / nsub; hz = (zs[end] - zs[1]) / nsub
            for cz = 1:nsub, gz = 1:3
                z = zs[1] + (cz - 0.5 + 0.5xg[gz]) * hz; lag!(lz, zs, z)
                for cy = 1:nsub, gy = 1:3
                    y = ys[1] + (cy - 0.5 + 0.5xg[gy]) * hy; lag!(ly, ys, y)
                    for cx = 1:nsub, gx = 1:3
                        x = xs[1] + (cx - 0.5 + 0.5xg[gx]) * hx; lag!(lx, xs, x)
                        σp = 0.0
                        for k = 1:ngl, j = 1:ngl
                            t = ly[j] * lz[k]
                            for i = 1:ngl; σp += σ_sw[mesh.connijk[iel, i, j, k]] * lx[i] * t; end
                        end
                        σp <= 0 && continue
                        c = σp * sw.S₀_flux * exp(-sw_beam_tau(R, x, y, z) / μ0)          # σ F_dir / μ₀
                        w = c * (wg[gx] * wg[gy] * wg[gz] / 8) * hx * hy * hz
                        for k = 1:ngl, j = 1:ngl
                            t = w * ly[j] * lz[k]
                            for i = 1:ngl; load[lidx(i, j, k, ngl)] += lx[i] * t; end
                        end
                    end
                end
            end
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                l = lidx(i, j, k, ngl); C[l, iel] = load[l] / (ω[i] * ω[j] * ω[k] * V / 8)
            end
        end
        # Zhang–Shu positivity: scale toward the element mean (weights ω_i ω_j ω_k; the mean is the exact one, ≥ 0)
        sw_ = 0.0; sm = 0.0; mn = Inf; tot = 0.0
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            l = lidx(i, j, k, ngl); w = ω[i] * ω[j] * ω[k]
            sw_ += w; sm += w * C[l, iel]; mn = min(mn, C[l, iel]); tot += w * abs(C[l, iel])
        end
        total += tot * V / 8
        if mn < 0
            m = max(sm / sw_, 0.0); θ = m > 0 ? m / (m - mn) : 0.0; mv = 0.0
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                l = lidx(i, j, k, ngl); cn = m + θ * (C[l, iel] - m)
                mv += ω[i] * ω[j] * ω[k] * abs(cn - C[l, iel]); C[l, iel] = cn
            end
            nlim += 1; moved += mv * V / 8
        end
    end
    return C, (projected = nproj, limited = nlim, moved = moved, total = total)
end

gll_weights_1d(ngl) = begin
    # Gauss–Lobatto–Legendre weights on [-1, 1]: w_i = 2 / (n (n−1) P_{n−1}(x_i)²), nodes = ±1 and roots of P'_{n−1}
    n = ngl; N = n - 1
    x = vcat(-1.0, n > 2 ? eigen(SymTridiagonal(zeros(N - 1), [sqrt(k * (k + 2) / ((2k + 1) * (2k + 3))) for k = 1:N-2])).values : Float64[], 1.0)
    P(t) = (p0 = 1.0; p1 = t; for k = 2:N; p0, p1 = p1, ((2k - 1) * t * p1 - (k - 1) * p0) / k; end; N == 0 ? p0 : p1)
    [2 / (n * N * P(t)^2) for t in x]
end

"""
    compute_direct_beam_raytrace(mesh, κ_sw, σ_sw, sw, comm; project_source=false) -> (F_dir, G_dir, Q_dir, info)

Direct solar beam with the 3D slant path, traced element by element through the polynomial κ_ext = κ + σ (see the
header of sw_direct_beam.jl). Works on conforming, non-conforming (hanging-node) and periodic box-element meshes.
project_source = true also returns info.source_elem: the projected, positivity-limited scattered-beam source
(sw_beam_projected_source), and its statistics in info.source_stats.
"""
function compute_direct_beam_raytrace(mesh, κ_sw::Vector{Float64}, σ_sw::Vector{Float64}, sw::SWParams, comm;
                                      project_source::Bool=false)
    npoin = length(mesh.x)
    R = sw_beam_ray(mesh, κ_sw .+ σ_sw, sw, comm)
    μ0 = sw.μ₀
    F_dir = zeros(npoin)
    if μ0 > 0
        for ip = 1:npoin
            F_dir[ip] = sw.S₀_flux * μ0 * exp(-sw_beam_tau(R, mesh.x[ip], mesh.y[ip], mesh.z[ip]) / μ0)
        end
    end
    src = nothing; sst = nothing
    project_source && ((src, sst) = sw_beam_projected_source(mesh, R, σ_sw, F_dir, sw))
    nf = MPI.Allreduce(R.nfail, MPI.SUM, comm)
    nf == 0 || error("compute_direct_beam_raytrace: $nf rays left the element set (point location failed)")
    G_dir = μ0 > 0 ? F_dir ./ μ0 : zeros(npoin)
    Q_dir = κ_sw .* G_dir
    B = R.B
    sstats = sst === nothing ? nothing :
        (v = MPI.Allreduce(Float64[sst.projected, sst.limited, sst.moved, sst.total, mesh.nelem], MPI.SUM, comm);
         (projected_elements = Int(v[1]), limited_elements = Int(v[2]), elements = Int(v[5]),
          moved_fraction = v[3] / max(v[4], 1e-300)))
    info = (nelem_global = size(B.lo, 2), periodic = R.per, max_segments = MPI.Allreduce(R.nseg_max, MPI.MAX, comm),
            refined_segments = MPI.Allreduce(R.nclip, MPI.SUM, comm),
            gathered_MB = 8 * (length(B.kap) + length(B.xn) + 2length(B.lo)) / 1e6,
            shift_km_per_km = sqrt(max(1 - μ0^2, 0.0)) / max(μ0, 1e-12),
            source_elem = src, source_stats = sstats)
    return F_dir, G_dir, Q_dir, info
end
