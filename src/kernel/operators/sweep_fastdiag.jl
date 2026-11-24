# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  Fast-diagonalisation element solver (opt-in: SweepOptions(elem_solver = :fastdiag), :RT_elem_solver, or
#  ENV JEXPRESSO_RT_FASTDIAG = 1). Used only on the workspace path of solve_pair! (no cached LU factor).
#
#  For an AFFINE element (constant metric terms and Jacobian), the DG transport element matrix built by
#  build_dg_volume_matrix! + build_dg_flux_terms! satisfies, with W = diag(ω_i ω_j ω_k J) and s_a = Ω·∇ξ_a,
#
#       W⁻¹ A = L_ξ ⊕ L_η ⊕ L_ζ + diag(κ),      L_a = |s_a| M_±,   M_+ = Dᵀ + E_1/ω_1,   M_- = -Dᵀ + E_n/ω_n,
#
#  (Dᵀ[i,m] = dψ[m,i]; E_e the unit diagonal at the inflow end; ± = sign(s_a): the inflow face is i = 1 when s_a > 0).
#  Only M_+ and M_- are eigendecomposed (once). With κ̄ = mean κ the Kronecker sum + κ̄ I is inverted exactly by 1-D
#  transforms; variable κ is handled by right-preconditioned GMRES on the exact operator (applied in tensor form) with
#  that inverse as preconditioner. Not converged within FD_MAXIT, a non-affine element, an inactive or inconsistent
#  inflow face, or a large κ contrast: the caller falls back to the dense LU path, unchanged.
#  Self-check: the first FD_CHECK pairs of a process are also solved densely; a relative difference above 1e-9
#  disables the fast path for the rest of the run (with a warning).
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
const USE_FASTDIAG   = Ref(false)
const FD_MAXIT       = Ref(12)
const FD_TOL         = Ref(1e-12)
const FD_KAPPA_SKIP  = Ref(0.5)      # skip when (κmax - κmin) > FD_KAPPA_SKIP * (κ̄ + min_a |s_a|): GMRES would stall
const FD_CHECK       = Ref(64)       # self-check comparisons per process
const FD_DISABLED    = Ref(false)
mutable struct FDStats
    n_fd::Int; n_fallback::Int; n_skip::Int; its::Int; n_checked::Int; max_check::Float64
    n_tan::Int; its_tan::Int; n_fallback_tan::Int; n_checked_tan::Int
end
const FD_STATS = FDStats(0, 0, 0, 0, 0, 0.0, 0, 0, 0, 0)

struct FDBase
    Mp::Matrix{Float64}; Mm::Matrix{Float64}
    Vp::Matrix{ComplexF64}; Vip::Matrix{ComplexF64}; μp::Vector{ComplexF64}
    Vm::Matrix{ComplexF64}; Vim::Matrix{ComplexF64}; μm::Vector{ComplexF64}
    D::Matrix{Float64}                                  # transport only (tangential axis: s_a D, no face terms)
    Ic::Matrix{ComplexF64}; μ0::Vector{ComplexF64}      # identity transform / zero spectrum for a tangential axis
end
mutable struct FDWork
    n::Int; maxit::Int
    V::Matrix{Float64}; Z::Matrix{Float64}; H::Matrix{Float64}; g::Vector{Float64}; cs::Vector{Float64}; sn::Vector{Float64}
    w::Vector{Float64}; f::Vector{Float64}; κl::Vector{Float64}; W::Vector{Float64}; den::Vector{ComplexF64}
    c1::Vector{ComplexF64}; c2::Vector{ComplexF64}; c3::Vector{ComplexF64}; r1::Vector{Float64}; r2::Vector{Float64}
    xsrc::Vector{Float64}; ucheck::Vector{Float64}      # source-only rhs (restored on fallback); self-check solution
    # preconditioner transforms of the current pair (FDBase's, per axis mode)
    PVx::Matrix{ComplexF64}; PVy::Matrix{ComplexF64}; PVz::Matrix{ComplexF64}
    PVix::Matrix{ComplexF64}; PViy::Matrix{ComplexF64}; PViz::Matrix{ComplexF64}
    # exact operator of the current pair: axis a contributes c_a M_a (|s_a| M_± if active, s_a D if tangential)
    AMx::Matrix{Float64}; AMy::Matrix{Float64}; AMz::Matrix{Float64}; acx::Float64; acy::Float64; acz::Float64
    tan::Bool                                           # at least one tangential (inactive-face) axis
end
const FDAffine = Union{Nothing,NTuple{3,NTuple{3,Float64}}}
const FD_AFFINE_G   = FDAffine[]      # per element: contravariant gradients if affine, else nothing
const FD_AFFINE_SET = BitVector()     # per element: already classified
const FD_BASE = Ref{Any}(nothing)
const FD_WORK = Ref{Any}(nothing)
const FD_AFFINE = Ref{Any}(nothing)   # the mesh the per-element affine cache below belongs to (reset: = nothing)

function fd_base(dψ, ω, ngl::Int)
    b = FD_BASE[]
    (b !== nothing && size(b.Mp, 1) == ngl) && return b
    Dt = [Float64(dψ[m,i]) for i = 1:ngl, m = 1:ngl]
    Mp = copy(Dt);  Mp[1,1]     += 1 / ω[1]
    Mm = -copy(Dt); Mm[ngl,ngl] += 1 / ω[ngl]
    Ep = eigen(Mp); Em = eigen(Mm)
    b = FDBase(Mp, Mm, Ep.vectors, inv(Ep.vectors), Ep.values, Em.vectors, inv(Em.vectors), Em.values,
               Dt, Matrix{ComplexF64}(LinearAlgebra.I, ngl, ngl), zeros(ComplexF64, ngl))
    FD_BASE[] = b
    return b
end

function fd_work(ngl::Int)
    n = ngl^3
    w = FD_WORK[]
    (w !== nothing && w.n == n && w.maxit == FD_MAXIT[]) && return w
    m = FD_MAXIT[]; c() = zeros(ComplexF64, ngl, ngl)
    w = FDWork(n, m, zeros(n, m+1), zeros(n, m), zeros(m+1, m), zeros(m+1), zeros(m), zeros(m),
               zeros(n), zeros(n), zeros(n), zeros(n), zeros(ComplexF64, n), zeros(ComplexF64, n), zeros(ComplexF64, n),
               zeros(ComplexF64, n), zeros(n), zeros(n), zeros(n), zeros(n),
               c(), c(), c(), c(), c(), c(),
               zeros(ngl, ngl), zeros(ngl, ngl), zeros(ngl, ngl), 0.0, 0.0, 0.0, false)
    FD_WORK[] = w
    return w
end

# per element: the three contravariant gradients (∇ξ, ∇η, ∇ζ) if the element is affine (constant metric terms and Je)
function fd_affine(mesh, metrics, iel::Int, ngl::Int)
    G = FD_AFFINE_G; S = FD_AFFINE_SET
    if FD_AFFINE[] !== mesh
        resize!(G, mesh.nelem); fill!(G, nothing); resize!(S, mesh.nelem); fill!(S, false); FD_AFFINE[] = mesh
    end
    @inbounds S[iel] && return G[iel]
    S[iel] = true
    m = metrics           # St_metrics fields are untyped (Any): pass the arrays through a function barrier
    g = fd_affine_compute(m.dξdx, m.dξdy, m.dξdz, m.dηdx, m.dηdy, m.dηdz, m.dζdx, m.dζdy, m.dζdz, m.Je, iel, ngl)::FDAffine
    g === nothing && return nothing
    G[iel] = g
    return g
end

# affine test tolerance (relative to the largest metric term / to Je): gmsh node coordinates carry ~1e-9 m rounding
# (e.g. the column next to the periodic x seam of the BOMEX meshes: dζdx ≈ -8e-12 instead of 0, i.e. 4e-10 × scale)
const FD_AFFINE_TOL = 1e-9

function fd_affine_compute(dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, Je, iel::Int, ngl::Int)
    g = ((dξdx[iel,1,1,1], dξdy[iel,1,1,1], dξdz[iel,1,1,1]),
         (dηdx[iel,1,1,1], dηdy[iel,1,1,1], dηdz[iel,1,1,1]),
         (dζdx[iel,1,1,1], dζdy[iel,1,1,1], dζdz[iel,1,1,1]))
    J0 = Je[iel,1,1,1]
    scale = maximum(abs, (g[1]..., g[2]..., g[3]...))
    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl
        gg = ((dξdx[iel,i,j,k], dξdy[iel,i,j,k], dξdz[iel,i,j,k]),
              (dηdx[iel,i,j,k], dηdy[iel,i,j,k], dηdz[iel,i,j,k]),
              (dζdx[iel,i,j,k], dζdy[iel,i,j,k], dζdz[iel,i,j,k]))
        for a = 1:3, c = 1:3
            abs(gg[a][c] - g[a][c]) <= FD_AFFINE_TOL * scale || return nothing
        end
        abs(Je[iel,i,j,k] - J0) <= FD_AFFINE_TOL * abs(J0) || return nothing
    end
    return (Float64.(g[1]), Float64.(g[2]), Float64.(g[3]))
end

# u ← (Mz ⊗ My ⊗ Mx) v, i fastest; buffers b1, b2 (same element type as the result)
@inline function fd_tensor!(out, Mx, My, Mz, v, b1, ngl::Int)
    n2 = ngl*ngl
    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl                    # x
        s = zero(eltype(b1)); base = (j-1)*ngl + (k-1)*n2
        for m = 1:ngl; s += Mx[i,m] * v[m + base]; end
        b1[i + base] = s
    end
    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl                    # y
        s = zero(eltype(out)); base = i + (k-1)*n2
        for m = 1:ngl; s += My[j,m] * b1[base + (m-1)*ngl]; end
        out[base + (j-1)*ngl] = s
    end
    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl                    # z (in place via b1)
        s = zero(eltype(b1)); base = i + (j-1)*ngl
        for m = 1:ngl; s += Mz[k,m] * out[base + (m-1)*n2]; end
        b1[base + (k-1)*n2] = s
    end
    copyto!(out, b1)
    return out
end

# y ← (L_ξ ⊕ L_η ⊕ L_ζ + diag κ) u   (real, tensor form; axis a uses |s_a| M_±)
function fd_apply_op!(y, u, wk::FDWork, κl, ngl::Int)
    n2 = ngl*ngl
    Mx = wk.AMx; My = wk.AMy; Mz = wk.AMz
    ax = wk.acx; ay = wk.acy; az = wk.acz
    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl
        p = i + (j-1)*ngl + (k-1)*n2
        sx = 0.0; sy = 0.0; sz = 0.0
        for m = 1:ngl
            sx += Mx[i,m] * u[m + (j-1)*ngl + (k-1)*n2]
            sy += My[j,m] * u[i + (m-1)*ngl + (k-1)*n2]
            sz += Mz[k,m] * u[i + (j-1)*ngl + (m-1)*n2]
        end
        y[p] = ax*sx + ay*sy + az*sz + κl[p]*u[p]
    end
    return y
end

# z ← P⁻¹ v by 1-D eigen-transforms (wk.PV*, set by fd_prepare!): P = L_ξ ⊕ L_η ⊕ L_ζ + κ̄ I (mean κ);
# den holds 1/(λ_ξ + λ_η + λ_ζ + κ̄)
function fd_precond!(z, v, den, wk::FDWork, ngl::Int)
    Vix = wk.PVix; Viy = wk.PViy; Viz = wk.PViz; Vx = wk.PVx; Vy = wk.PVy; Vz = wk.PVz
    c1 = wk.c1; c2 = wk.c2; c3 = wk.c3
    @inbounds for p in eachindex(v); c3[p] = v[p]; end
    fd_tensor!(c1, Vix, Viy, Viz, c3, c2, ngl)
    @inbounds for p in eachindex(c1); c3[p] = c1[p] * den[p]; end
    fd_tensor!(c1, Vx, Vy, Vz, c3, c2, ngl)
    @inbounds for p in eachindex(z); z[p] = real(c1[p]); end
    return z
end

"""
    fd_prepare!(mesh, metrics, dψ, ω, κ, iel, ngl, Ωx, Ωy, Ωz, m) -> s or nothing

Eligibility of one (element, direction) pair for the fast path, checked BEFORE any right-hand side is built:
affine element, inflow faces exactly the low/high ends and active, κ contrast within FD_KAPPA_SKIP. When eligible,
fills the work arrays (κl, W, den) and returns s = (Ω·∇ξ, Ω·∇η, Ω·∇ζ); otherwise nothing (counted as a skip).
"""
function fd_prepare!(mesh, metrics, dψ, ω, κ, iel::Int, ngl::Int, Ωx::Float64, Ωy::Float64, Ωz::Float64, m::UInt16)
    g0 = fd_affine(mesh, metrics, iel, ngl)
    g0 === nothing && (FD_STATS.n_skip += 1; FD_SKIP_WHY[FDW_AFFINE] += 1; return nothing)
    g = g0::NTuple{3,NTuple{3,Float64}}
    s = fd_svec(g, Ωx, Ωy, Ωz)
    why, mo = fd_axis_modes(s, m)
    why == 0 || (FD_STATS.n_skip += 1; FD_SKIP_WHY[why] += 1; return nothing)
    n = ngl^3; B = fd_base(dψ, ω, ngl)::FDBase; wk = fd_work(ngl)::FDWork
    wk.tan = mo[1] == 0 || mo[2] == 0 || mo[3] == 0
    wk.AMx = fd_axmat(B, mo[1]); wk.AMy = fd_axmat(B, mo[2]); wk.AMz = fd_axmat(B, mo[3])
    wk.acx = mo[1] == 0 ? s[1] : abs(s[1]); wk.acy = mo[2] == 0 ? s[2] : abs(s[2]); wk.acz = mo[3] == 0 ? s[3] : abs(s[3])
    κmin, κmax, κsum = fd_gather_kappa!(wk.κl, wk.W, κ, mesh.connijk, metrics.Je, ω, iel, ngl)::NTuple{3,Float64}   # function barrier (typed result)
    κ̄ = κsum / n
    smin = min(mo[1] == 0 ? Inf : abs(s[1]), mo[2] == 0 ? Inf : abs(s[2]), mo[3] == 0 ? Inf : abs(s[3]))
    (κmax - κmin) > FD_KAPPA_SKIP[] * (κ̄ + smin) && (FD_STATS.n_skip += 1; FD_SKIP_WHY[FDW_KAPPA] += 1; return nothing)
    # tangential axis: identity transform, zero spectrum (its tiny s_a D term is left to GMRES)
    wk.PVx  = fd_axsel(mo[1], B.Vp, B.Vm, B.Ic);   wk.PVy  = fd_axsel(mo[2], B.Vp, B.Vm, B.Ic);   wk.PVz  = fd_axsel(mo[3], B.Vp, B.Vm, B.Ic)
    wk.PVix = fd_axsel(mo[1], B.Vip, B.Vim, B.Ic); wk.PViy = fd_axsel(mo[2], B.Vip, B.Vim, B.Ic); wk.PViz = fd_axsel(mo[3], B.Vip, B.Vim, B.Ic)
    μx = fd_axsel(mo[1], B.μp, B.μm, B.μ0); μy = fd_axsel(mo[2], B.μp, B.μm, B.μ0); μz = fd_axsel(mo[3], B.μp, B.μm, B.μ0)
    ax, ay, az = abs(s[1]), abs(s[2]), abs(s[3]); den = wk.den
    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl
        den[i + (j-1)*ngl + (k-1)*ngl*ngl] = 1 / (ax*μx[i] + ay*μy[j] + az*μz[k] + κ̄)
    end
    return s
end

@inline fd_svec(g::NTuple{3,NTuple{3,Float64}}, Ωx::Float64, Ωy::Float64, Ωz::Float64) =
    (g[1][1]*Ωx + g[1][2]*Ωy + g[1][3]*Ωz, g[2][1]*Ωx + g[2][2]*Ωy + g[2][3]*Ωz, g[3][1]*Ωx + g[3][2]*Ωy + g[3][3]*Ωz)

# Reasons a pair goes to the dense solver (FD_SKIP_WHY counters, per process)
const FDW_AFFINE = 1      # element not affine
const FDW_S0     = 2      # active faces on an axis but s_a == 0 exactly
const FDW_ONESIDE = 3     # only one of an axis' two faces active
const FDW_INFLOW = 4      # the upwind face (by sign of s_a) is not marked inflow
const FDW_OUTFLOW = 5     # the downwind face is marked inflow
const FDW_ALLTAN = 6      # no active axis at all
const FDW_KAPPA  = 7      # κ contrast too large
const FD_SKIP_WHY = zeros(Int, 7)
const FDW_NAMES = ("not affine", "s_a = 0 on an active axis", "one face of an axis active", "upwind face not inflow",
                   "downwind face marked inflow", "no active axis", "κ contrast")

# Per axis: either both faces inactive (direction tangential to them: the dense path drops both face terms and
# keeps only the volume transport s_a D), or the inflow face is exactly the low end (s > 0) / high end (s < 0)
#   ξ: low = fid 5 (i=1), high = fid 6;  η: low = fid 2 (j=1), high = fid 1;  ζ: low = fid 4 (k=1), high = fid 3
# Returns (reason, modes): reason 0 = eligible, else an FDW_* code; mode per axis +1 (M+), -1 (M-), 0 (tangential).
function fd_axis_modes(s::NTuple{3,Float64}, m::UInt16)
    mo = (0, 0, 0)
    for (a, lo, hi) in ((1, 5, 6), (2, 2, 1), (3, 4, 3))
        alo = mask_active(m, lo); ahi = mask_active(m, hi)
        if !alo && !ahi
            mode = 0
        else
            alo && ahi || return FDW_ONESIDE, mo
            s[a] == 0 && return FDW_S0, mo
            fin, fout = s[a] > 0 ? (lo, hi) : (hi, lo)
            mask_inflow(m, fin) || return FDW_INFLOW, mo
            mask_inflow(m, fout) && return FDW_OUTFLOW, mo
            mode = s[a] > 0 ? 1 : -1
        end
        mo = a == 1 ? (mode, mo[2], mo[3]) : a == 2 ? (mo[1], mode, mo[3]) : (mo[1], mo[2], mode)
    end
    mo == (0, 0, 0) && return FDW_ALLTAN, mo
    return 0, mo
end

@inline fd_axmat(B::FDBase, mode::Int) = mode > 0 ? B.Mp : mode < 0 ? B.Mm : B.D
@inline fd_axsel(mode::Int, p, mm, z) = mode > 0 ? p : mode < 0 ? mm : z

function fd_gather_kappa!(κl, W, κ, connijk, Je, ω, iel::Int, ngl::Int)
    κmin = Inf; κmax = -Inf; κsum = 0.0
    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl
        p = i + (j-1)*ngl + (k-1)*ngl*ngl
        v = Float64(κ[connijk[iel,i,j,k]]); κl[p] = v; W[p] = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
        κmin = min(κmin, v); κmax = max(κmax, v); κsum += v
    end
    return κmin, κmax, κsum
end

"""
    fd_solve!(u, rhs, s, ngl) -> Bool

Fast solve of a pair prepared by fd_prepare! (same pair, nothing in between). u may alias rhs: rhs is read once
(into the work array f) before u is written, and u is left untouched when GMRES does not converge.
"""
function fd_solve!(u::AbstractVector{Float64}, rhs::AbstractVector{Float64}, s::NTuple{3,Float64}, ngl::Int)
    n = ngl^3; B = FD_BASE[]::FDBase; wk = FD_WORK[]::FDWork
    f = wk.f; W = wk.W
    @inbounds for p = 1:n; f[p] = rhs[p] / W[p]; end
    return fd_gmres!(u, f, s, B, wk, wk.κl, wk.den, n, ngl)
end

# both steps (benchmarks / tests); false (u untouched) when not eligible or not converged
function fd_try_solve!(u, rhs, mesh, metrics, dψ, ω, κ, iel::Int, ngl::Int, Ωx::Float64, Ωy::Float64, Ωz::Float64, m::UInt16)
    s = fd_prepare!(mesh, metrics, dψ, ω, κ, iel, ngl, Ωx, Ωy, Ωz, m)
    s === nothing && return false
    return fd_solve!(u, rhs, s, ngl)
end

# right-preconditioned GMRES from zero on the exact operator (typed inner function: the hot loops)
function fd_gmres!(u::AbstractVector{Float64}, f::Vector{Float64}, s::NTuple{3,Float64}, B::FDBase, wk::FDWork,
                   κl::Vector{Float64}, den::Vector{ComplexF64}, n::Int, ngl::Int)
    β = norm(f); β == 0 && (fill!(u, 0.0); FD_STATS.n_fd += 1; return true)
    V = wk.V; Z = wk.Z; H = wk.H; gv = wk.g; cs = wk.cs; sn = wk.sn; w = wk.w; maxit = wk.maxit
    fill!(H, 0.0); fill!(gv, 0.0); gv[1] = β
    @inbounds for p = 1:n; V[p,1] = f[p] / β; end
    its = 0; conv = false
    for jj = 1:maxit
        fd_precond!(view(Z, :, jj), view(V, :, jj), den, wk, ngl)
        fd_apply_op!(w, view(Z, :, jj), wk, κl, ngl)
        for i = 1:jj
            h = 0.0; @inbounds for p = 1:n; h += V[p,i] * w[p]; end
            H[i,jj] = h; @inbounds for p = 1:n; w[p] -= h * V[p,i]; end
        end
        hn = norm(w); H[jj+1,jj] = hn
        hn > 0 && (@inbounds for p = 1:n; V[p,jj+1] = w[p] / hn; end)
        for i = 1:jj-1
            t = cs[i]*H[i,jj] + sn[i]*H[i+1,jj]; H[i+1,jj] = -sn[i]*H[i,jj] + cs[i]*H[i+1,jj]; H[i,jj] = t
        end
        ρ = hypot(H[jj,jj], H[jj+1,jj]); cs[jj] = H[jj,jj]/ρ; sn[jj] = H[jj+1,jj]/ρ; H[jj,jj] = ρ; H[jj+1,jj] = 0.0
        gv[jj+1] = -sn[jj]*gv[jj]; gv[jj] = cs[jj]*gv[jj]; its = jj
        if abs(gv[jj+1]) <= FD_TOL[] * β; conv = true; break; end
        hn == 0 && break
    end
    conv || (FD_STATS.n_fallback += 1; wk.tan && (FD_STATS.n_fallback_tan += 1); return false)
    @inbounds for i = its:-1:1                         # back substitution in place: gv[1:its] ← H[1:its,1:its] \ gv
        t = gv[i]
        for l = i+1:its; t -= H[i,l] * gv[l]; end
        gv[i] = t / H[i,i]
    end
    fill!(u, 0.0)
    for i = 1:its; yi = gv[i]; @inbounds for p = 1:n; u[p] += yi * Z[p,i]; end; end
    FD_STATS.n_fd += 1; FD_STATS.its += its
    wk.tan && (FD_STATS.n_tan += 1; FD_STATS.its_tan += its)
    return true
end

function fd_stats_reset!()
    FD_STATS.n_fd = 0; FD_STATS.n_fallback = 0; FD_STATS.n_skip = 0; FD_STATS.its = 0
    FD_STATS.n_tan = 0; FD_STATS.its_tan = 0; FD_STATS.n_fallback_tan = 0
    fill!(FD_SKIP_WHY, 0)
    return nothing
end


# ---------------------------------------------------------------------------------------------------------------
# Per-element RT cost prediction (weights for repartitioning). Each (element, direction) pair is classified by the
# same tests fd_prepare! applies (affine element, face masks, κ contrast), without solving anything,
# and costed with per-class constants (µs per pair, fitted from SWEEP PROFILE runs). For fast pairs the GMRES
# iteration count is predicted from the contrast ratio r = (κmax - κmin)/(κ̄ + min|s_a|) ∈ [0, FD_KAPPA_SKIP].
# Classes: 1 fast, 2 fast tangential, 3 dense (κ contrast), 4 dense (masks / not affine)
# ---------------------------------------------------------------------------------------------------------------
# Constants fitted on BOMEX 32x32x30 6x6O4, 32 ranks (job 1383682, 30720 elements, 3 profiled sweeps; element R² 0.997,
# per-rank error <= 5.4%): dense 271.8 µs/pair; fast 43.2 + 14.84 µs per GMRES iteration; its = 3.96 + 7.18 r
const RTC_PAIR   = Ref(0.0)      # every solved pair (folded into the class costs below)
const RTC_FAST0  = Ref(43.2)     # fast pair: rhs flux terms, setup, upwind faces, source
const RTC_FASTIT = Ref(14.84)    # fast pair: per GMRES iteration
const RTC_ITS0   = Ref(3.96)     # predicted iterations: its0 + itsr * r
const RTC_ITSR   = Ref(7.18)
const RTC_DENSE  = Ref(271.8)    # dense pair: volume matrix + face terms + LU + upwind faces, source
const RTC_DENSE_HIT = Ref(56.0)  # dense pair solved from a cached factor (flux terms + LU solve); fitted on job 1384361
const RTC_NCLASS = 4

struct RTCostPrediction
    w::Vector{Float64}           # predicted µs per element per sweep
    cnt::Matrix{Int}             # nelem × RTC_NCLASS pairs per class
    its::Vector{Float64}         # predicted GMRES iterations summed over the element's fast pairs
end

"""
    rt_predict_costs(mesh, metrics, dψ, ω, κ, ngl, n_dirs, faceinfo) -> RTCostPrediction

No solves, no allocation per pair. `faceinfo` = the element cache's per-direction face masks.
"""
function rt_predict_costs(mesh, metrics, dψ, ω, κ, ngl::Int, n_dirs::Vector{Int}, faceinfo)
    nel = Int(mesh.nelem); n = ngl^3
    w = zeros(nel); cnt = zeros(Int, nel, RTC_NCLASS); its = zeros(nel)
    wk = fd_work(ngl)::FDWork
    fd_base(dψ, ω, ngl)
    for iel = 1:nel
        em = mesh.extra_mesh[iel]
        rt_predict_element!(w, cnt, its, iel, em.extra_coords, n_dirs[iel], faceinfo[iel]::Vector{UInt16},
                            fd_affine(mesh, metrics, iel, ngl)::FDAffine, wk, κ, mesh.connijk, metrics.Je, ω, ngl, n)
    end
    return RTCostPrediction(w, cnt, its)
end

function rt_predict_element!(w, cnt, its, iel::Int, ec, nd::Int, masks::Vector{UInt16}, g0::FDAffine, wk::FDWork,
                             κ, connijk, Je, ω, ngl::Int, n::Int)
    c = 0.0
    if g0 === nothing
        cnt[iel, 4] += nd; w[iel] = nd * (RTC_PAIR[] + RTC_DENSE[]); return nothing
    end
    g = g0::NTuple{3,NTuple{3,Float64}}
    κmin, κmax, κsum = fd_gather_kappa!(wk.κl, wk.W, κ, connijk, Je, ω, iel, ngl)::NTuple{3,Float64}
    κ̄ = κsum / n; Δκ = κmax - κmin
    for d = 1:nd
        θ = Float64(ec[1,d]); ϕ = Float64(ec[2,d])
        s = fd_svec(g, sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
        why, mo = fd_axis_modes(s, masks[d])
        c += RTC_PAIR[]
        if why != 0
            cnt[iel, 4] += 1; c += RTC_DENSE[]; continue
        end
        smin = min(mo[1] == 0 ? Inf : abs(s[1]), mo[2] == 0 ? Inf : abs(s[2]), mo[3] == 0 ? Inf : abs(s[3]))
        tan = mo[1] == 0 || mo[2] == 0 || mo[3] == 0
        if Δκ > FD_KAPPA_SKIP[] * (κ̄ + smin)
            cnt[iel, 3] += 1; c += RTC_DENSE[]; continue
        end
        cnt[iel, tan ? 2 : 1] += 1
        ip = RTC_ITS0[] + RTC_ITSR[] * Δκ / (κ̄ + smin)
        its[iel] += ip; c += RTC_FAST0[] + RTC_FASTIT[] * ip
    end
    w[iel] = c
    return nothing
end

