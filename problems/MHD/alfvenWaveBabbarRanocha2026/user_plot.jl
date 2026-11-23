#---------------------------------------------------------------------------------
# Convergence history of the smooth Alfvén wave (Babbar & Ranocha 2026,
# Fig. 13): the error at the final time against the number of degrees of
# freedom on log-log axes, drawn with the convergence figure of the smooth
# vortex (../smoothVortex/convergence_figure.jl), so that the two accuracy
# tests show their errors the same way: every order and stabilization on the
# same axes, dotted reference slopes, the measured rate in the legend.
#
# THE ERROR is the one the paper's runs report: Tenkai.jl's compute_error, read
# as "l2_error" by run/convergence_alfven_wave/convergence.jl of the
# reproducibility repository,
#
#     e(q) = sqrt( (1/|Ω|) ∫_Ω (q_h − q)² dΩ ),
#
# integrated with N + 4 Gauss-Legendre points per direction in every element
# (Tenkai: nd + 3 with nd = N + 1), the numerical solution interpolated there
# from its LGL nodes. Tenkai evaluates it for the FIRST conserved variable only,
# the density. The exact density is 1, so that number measures how far the
# scheme moves ρ off a constant (a by-product of the nonlinear coupling), not
# the wave itself. Every conserved variable is stored here, and so are the
# errors of the WAVE: the transverse field and velocity,
#
#     B⊥ = (B·e⊥, Bz),   u⊥ = (u·e⊥, w),   e⊥ = (−sin α, cos α),
#
# e(B⊥) = sqrt((1/|Ω|) ∫ |B⊥,h − B⊥|² dΩ) and the same for u⊥; the wave is
# 0.1 (sin, cos)(2π(x∥ + t)) in both. Each is also stored in L¹,
# (1/|Ω|) ∫ |q_h − q| dΩ (|B⊥,h − B⊥| for the wave), the norm of most MHD
# accuracy tests, and in L∞, the largest error over the Gauss points. The
# figures are a matrix of panels: rows B⊥, u⊥ (the wave's own error), ρ (the
# paper's measure) and the full state, sqrt(Σ_k e(q_k)²) in L² and L¹, the
# largest e(q_k) in L∞; columns L², L¹, L∞. The wave's L∞ is in every record
# written since 17 Nov 2026 and L¹ in every record written since 23 Nov 2026;
# an older record leaves those panels without its point.
#
# WHAT DSGS DID. Each record also stores the largest normalized residual r
# (before the cutoff) and the largest coefficient μ over the whole run
# (SGS.jl, dsgs_run_stats): μ = 0 with r below the cutoff means the cutoff,
# not the residual, switched the model off.
#
# At the final time every run stores its errors in
# `errors/nop<N>_nelx<M>_<visc>.dat` in this case directory, and the figure is
# drawn from every error stored there for the same final time: a sweep builds
# the whole figure and each run replaces only its own point.
# `rm -r problems/MHD/alfvenWaveBabbarRanocha2026/errors` starts afresh.
#
# #DOFs is the number of distinct unknowns per variable of this continuous
# Galerkin discretization on the doubly periodic box, (nelx N)². The paper's
# flux reconstruction has (nx (N + 1))² on nx² elements, so at equal #DOFs the
# curves of the two compare at equal cost, not at equal mesh.
#
# The hook is src/io/plotting/jeplots.jl (plot_triangulation, NSD_2D):
#   mesh     this rank's mesh (coordinates, connectivity, global extents)
#   q        flat npoin*nvar vector of the OUTPUT variables (`outvar`)
#   t        simulation time
#---------------------------------------------------------------------------------
const AW_ERR_DIR = joinpath(@__DIR__, "errors")
# the stem of every figure this case writes: alfven-convergence-it<n>.png, ...
const AW_FIG_NAME = "alfven-convergence"

# The convergence figure: the smooth vortex's, shared.
include(joinpath(@__DIR__, "..", "smoothVortex", "convergence_figure.jl"))
const AW_QVARS   = ("ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ")   # the state order of initialize.jl

# ATOMIC WRITES, as in smoothVortex: the runs of a sweep may redraw the same
# figure and read the same store at the same time. Write beside the target and
# rename (atomic on POSIX).
function _aw_atomic(path::AbstractString, write!::Function)
    mkpath(dirname(path))
    base, ext = splitext(path)          # keep the extension: Plots picks the format from it
    tmp = string(base, ".tmp", getpid(), ext)
    try
        write!(tmp)
        mv(tmp, path; force = true)
    catch err
        isfile(tmp) && (try; rm(tmp); catch; end)
        rethrow(err)
    end
    return nothing
end

#---------------------------------------------------------------------------------
# Quadrature. Gauss-Legendre by Golub-Welsch; the LGL nodes of the solution as
# ±1 and the zeros of P_N', the eigenvalues of the Jacobi matrix of P^(1,1)_{N-1}.
# Only the SET of nodes matters below, not the order Jexpresso numbers them in:
# the coordinates are interpolated with the same basis as the solution.
#---------------------------------------------------------------------------------
function _aw_gauss(n::Int)
    b = [k/sqrt(4.0*k*k - 1.0) for k = 1:n-1]
    E = LinearAlgebra.eigen(LinearAlgebra.SymTridiagonal(zeros(n), b))
    return E.values, 2.0 .* E.vectors[1, :].^2
end

function _aw_lgl(N::Int)
    b = [sqrt(k*(k + 2.0)/((2.0*k + 1.0)*(2.0*k + 3.0))) for k = 1:N-2]
    x = N >= 2 ? LinearAlgebra.eigvals(LinearAlgebra.SymTridiagonal(zeros(N - 1), b)) : Float64[]
    return vcat(-1.0, sort(x), 1.0)
end

# L[a,i] = ℓ_i(x_a) and D[a,i] = ℓ_i'(x_a), the Lagrange basis on the nodes ξ.
function _aw_lagrange(ξ, x)
    n = length(ξ); m = length(x)
    L = ones(m, n); D = zeros(m, n)
    for a = 1:m, i = 1:n
        for k = 1:n
            k == i && continue
            L[a,i] *= (x[a] - ξ[k])/(ξ[i] - ξ[k])
            d = 1.0/(ξ[i] - ξ[k])
            for j = 1:n
                (j == i || j == k) && continue
                d *= (x[a] - ξ[j])/(ξ[i] - ξ[j])
            end
            D[a,i] += d
        end
    end
    return L, D
end

#---------------------------------------------------------------------------------
# The error of every conserved variable at time t. Element by element, so each
# element is integrated once whatever the periodic node duplication or the rank
# count; the sums are reduced over the ranks. Collective: every rank must call it.
#---------------------------------------------------------------------------------
function _aw_errors(mesh, q, t, outvar)
    names = string.(outvar)
    ivn = [findfirst(==(n), names) for n in ("ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ")]
    any(isnothing, ivn) && return nothing
    iv = Int[i for i in ivn]
    npoin = mesh.npoin
    ngl   = mesh.ngl
    N     = ngl - 1
    γ     = γ_mhd
    Lx    = mesh.xmax - mesh.xmin      # global extents (Allreduced by the mesh reader)
    Ly    = mesh.ymax - mesh.ymin

    xg, wg = _aw_gauss(N + 4)
    L, D   = _aw_lagrange(_aw_lgl(N), xg)
    nq     = length(xg)

    nv  = length(AW_QVARS)
    s2  = zeros(nv)
    s1  = zeros(nv)
    sm  = zeros(nv)
    area = 0.0
    X = zeros(ngl, ngl); Y = zeros(ngl, ngl)
    U = zeros(ngl, ngl, nv)
    Xq = zeros(nq, nq); Yq = zeros(nq, nq)
    Xa = zeros(nq, nq); Xb = zeros(nq, nq); Ya = zeros(nq, nq); Yb = zeros(nq, nq)
    Uq = zeros(nq, nq, nv)
    # the wave: B·e⊥, u·e⊥ and w at the nodes (Bz is U[:,:,8]); sums of
    # |B⊥,h − B⊥|² and |u⊥,h − u⊥|², and their maxima
    sα, cα = sincos(AW_ALPHA)
    T  = zeros(ngl, ngl, 3)
    Tq = zeros(nq, nq, 3)
    sp2 = zeros(2)
    sp1 = zeros(2)
    spm = zeros(2)
    out(k, ip) = q[(iv[k] - 1)*npoin + ip]

    for iel = 1:mesh.nelem
        ip1 = mesh.connijk[iel, 1, 1]
        x1  = mesh.coords[1, ip1]; y1 = mesh.coords[2, ip1]
        for j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel, i, j]
            # unwrap across the periodic seam, in case a node of this element
            # carries the coordinates of its periodic twin
            dx = mesh.coords[1, ip] - x1; dy = mesh.coords[2, ip] - y1
            X[i,j] = x1 + dx - Lx*round(dx/Lx)
            Y[i,j] = y1 + dy - Ly*round(dy/Ly)
            ρ  = out(1, ip); u = out(2, ip); v = out(3, ip); w = out(4, ip); p = out(5, ip)
            Bx = out(6, ip); By = out(7, ip); Bz = out(8, ip); ψ = out(9, ip)
            U[i,j,1] = ρ
            U[i,j,2] = ρ*u
            U[i,j,3] = ρ*v
            U[i,j,4] = p/(γ - 1.0) + 0.5*ρ*(u*u + v*v + w*w) + 0.5*(Bx*Bx + By*By + Bz*Bz) + 0.5*ψ*ψ
            U[i,j,5] = ρ*w
            U[i,j,6] = Bx; U[i,j,7] = By; U[i,j,8] = Bz; U[i,j,9] = ψ
            T[i,j,1] = -sα*Bx + cα*By       # B·e⊥
            T[i,j,2] = -sα*u  + cα*v        # u·e⊥
            T[i,j,3] = w
        end
        # coordinates and their derivatives at the Gauss points (isoparametric)
        Xq .= L*X*L';  Yq .= L*Y*L'
        Xa .= D*X*L';  Xb .= L*X*D'
        Ya .= D*Y*L';  Yb .= L*Y*D'
        for k = 1:nv
            Uq[:,:,k] .= L*view(U, :, :, k)*L'
        end
        for k = 1:3
            Tq[:,:,k] .= L*view(T, :, :, k)*L'
        end
        for b = 1:nq, a = 1:nq
            wab = wg[a]*wg[b]*abs(Xa[a,b]*Yb[a,b] - Xb[a,b]*Ya[a,b])
            se  = aw_state(Xq[a,b], Yq[a,b], t, γ)
            for k = 1:nv
                e = Uq[a,b,k] - se[k]
                s2[k] += wab*e*e
                s1[k] += wab*abs(e)
                sm[k]  = max(sm[k], abs(e))
            end
            # the wave's transverse field and velocity against the exact ones
            eB2 = (Tq[a,b,1] - (-sα*se[6] + cα*se[7]))^2 + (Uq[a,b,8] - se[8])^2
            eU2 = (Tq[a,b,2] - (-sα*se[2] + cα*se[3])/se[1])^2 + (Tq[a,b,3] - se[5]/se[1])^2
            sp2[1] += wab*eB2;  sp1[1] += wab*sqrt(eB2);  spm[1] = max(spm[1], sqrt(eB2))
            sp2[2] += wab*eU2;  sp1[2] += wab*sqrt(eU2);  spm[2] = max(spm[2], sqrt(eU2))
            area += wab
        end
    end

    comm = get_mpi_comm()
    if MPI.Comm_size(comm) > 1
        r = MPI.Allreduce(vcat(s2, sp2, s1, sp1, area), MPI.SUM, comm)
        s2 .= r[1:nv]; sp2 .= r[nv+1:nv+2]; s1 .= r[nv+3:2nv+2]; sp1 .= r[2nv+3:2nv+4]
        area = r[end]
        sm  .= MPI.Allreduce(sm,  MPI.MAX, comm)
        spm .= MPI.Allreduce(spm, MPI.MAX, comm)
    end
    # the quadrature must integrate the box: a cheap check that the element
    # geometry (and the unwrapping above) is right
    abs(area - Lx*Ly) <= 1.0e-10*Lx*Ly ||
        @warn "alfvenWaveBabbarRanocha2026: the error quadrature integrates $(area), not the box area $(Lx*Ly)"
    return (l2 = sqrt.(s2 ./ area), l1 = s1 ./ area, linf = sm,
            perp_l2 = sqrt.(sp2 ./ area), perp_l1 = sp1 ./ area, perp_linf = spm,   # (B⊥, u⊥)
            area = area)
end

# The record's tag names what differs from the deck as committed: "dsgs" is
# DynSGS as the deck runs it (residual sensor, coefficient frozen per step, the
# smoothness cutoff AW_DSGS_CUTOFF), "dsgs_nocut" the same without the cutoff
# (AW_CUTOFF = 0 in user_inputs.jl) and "dsgs_cut<x>" with another one,
# "dsgs_perstage" the coefficient at every stage (JEXPRESSO_AW_FREEZE=0),
# "dsgs_legacy" the legacy sensor (never frozen). A time integrator other than
# the deck's for that order (_aw_default_solver in user_inputs.jl) adds its
# name, and a time step of its own (JEXPRESSO_AW_DT) adds "dt<Δt>", so that no
# variant replaces the record of another.
#
# AW_DSGS_CUTOFF is the cutoff of the deck as committed, AW_CUTOFF's value in
# user_inputs.jl; keep the two equal. It is a constant of its own so that a run
# with AW_CUTOFF edited to 0 is still stored as "dsgs_nocut", not as "dsgs".
# (Records stored before this convention, 16 Oct 2026, called the run without
# the cutoff "dsgs" and the one with it "dsgs_cut1e-03".)
const AW_DSGS_CUTOFF = 1.0e-3
function _aw_tag(inputs)
    s   = get(inputs, :ode_solver, nothing)
    nop = Int(get(inputs, :nop, 0))
    ode = (s === nothing || nameof(typeof(s)) == nameof(typeof(_aw_solver_named(_aw_default_solver(nop))))) ?
          "" : string("_", lowercase(string(nameof(typeof(s)))))
    tryparse(Float64, get(ENV, "JEXPRESSO_AW_DT", "")) === nothing ||
        (ode *= @sprintf("_dt%.3g", Float64(get(inputs, :Δt, 0.0))))
    # the diagnostics of initialize.jl: another GLM speed, another amplitude
    let ch = tryparse(Float64, get(ENV, "JEXPRESSO_AW_CH", ""))
        ch === nothing || (ode *= @sprintf("_ch%g", ch))
    end
    (isdefined(@__MODULE__, :AW_AMP) && AW_AMP != 0.1) && (ode *= @sprintf("_amp%g", AW_AMP))
    get(inputs, :lvisc, true) || return string("galerkin", ode)
    tag    = "dsgs"
    legacy = lowercase(string(get(inputs, :dsgs_sensor, "residual"))) == "legacy"
    legacy && (tag *= "_legacy")
    Float64(get(inputs, :dsgs_Cmin, 0.0)) > 0 && (tag *= "_cmin")
    cut = Float64(get(inputs, :dsgs_cutoff, 0.0))
    cut == AW_DSGS_CUTOFF || (tag *= cut > 0 ? @sprintf("_cut%.0e", cut) : "_nocut")
    !legacy && get(inputs, :dsgs_freeze_stage, false) != true && (tag *= "_perstage")
    tag *= _conv_bdf_suffix(inputs)      # "_bdf<k>" for :dsgs_bdf_order = k ≠ 2
    return string(tag, ode)
end

_aw_nelx_of(inputs) = (m = match(r"alfven_(\d+)x", string(get(inputs, :gmsh_filename, "")));
                       m === nothing ? Int(get(inputs, :nelx, 0)) : parse(Int, m.captures[1]))

# st = [max r, max μ] over the run (dsgs_run_stats; NaN, 0 for Galerkin)
function _aw_save(e, inputs, tnom, tnum, st)
    nop  = Int(get(inputs, :nop, 0))
    nelx = _aw_nelx_of(inputs)
    try
        # a final time other than the deck's 2 goes in the name, so that it does not
        # overwrite the t = 2 record of the same run (the loader reads t from the header)
        tsuf = abs(tnom - 2.0) > 1.0e-12 ? string("_t", tnom) : ""
        f = joinpath(AW_ERR_DIR, string("nop", nop, "_nelx", nelx, "_", _aw_tag(inputs), tsuf, ".dat"))
        _aw_atomic(f, tmp -> open(tmp, "w") do io
            println(io, "# smooth Alfvén wave: e(q) = sqrt((1/|Ω|) ∫ (q_h − q)² dΩ) against the exact solution,",
                        " N+4 Gauss points per direction; Linf at those points; L1 = (1/|Ω|) ∫ |q_h − q| dΩ.",
                        " Bperp = (B·e⊥, Bz), uperp = (u·e⊥, w): the wave, as vectors.",
                        " r_max, mu_max: the largest normalized residual (before the cutoff) and DSGS coefficient over the run")
            println(io, "# nop=", nop, " nelx=", nelx, " ndofs=", (nelx*nop)^2,
                        " t=", tnom, " tnum=", tnum, " visc=", _aw_tag(inputs),
                        " dt=", Float64(get(inputs, :Δt, 0.0)),
                        " Cmin=", Float64(get(inputs, :dsgs_Cmin, 0.0)),
                        " cutoff=", Float64(get(inputs, :dsgs_cutoff, 0.0)),
                        " solver=", nameof(typeof(get(inputs, :ode_solver, nothing))),
                        " area=", e.area,
                        " r_max=", st[1], " mu_max=", st[2])
            println(io, "# var L2 Linf L1")
            for (k, n) in enumerate(AW_QVARS)
                println(io, n, " ", e.l2[k], " ", e.linf[k], " ", e.l1[k])
            end
            println(io, "Bperp ", e.perp_l2[1], " ", e.perp_linf[1], " ", e.perp_l1[1])
            println(io, "uperp ", e.perp_l2[2], " ", e.perp_linf[2], " ", e.perp_l1[2])
        end)
    catch err
        @warn "alfvenWaveBabbarRanocha2026: could not store the error" exception=err
    end
    return nothing
end

function _aw_load(tnom)
    rows = NamedTuple[]
    isdir(AW_ERR_DIR) || return rows
    for fname in readdir(AW_ERR_DIR)
        (startswith(fname, "nop") && endswith(fname, ".dat")) || continue
        meta = Dict{String,String}(); l2 = Dict{String,Float64}(); linf = Dict{String,Float64}()
        l1 = Dict{String,Float64}()
        try
            for line in eachline(joinpath(AW_ERR_DIR, fname))
                tok = split(line)
                isempty(tok) && continue
                if tok[1] == "#"
                    for s in tok
                        occursin('=', s) || continue
                        k, v = split(s, '=', limit = 2); meta[k] = v
                    end
                elseif length(tok) >= 3
                    l2[tok[1]]   = parse(Float64, tok[2])
                    linf[tok[1]] = parse(Float64, tok[3])
                    length(tok) >= 4 && (l1[tok[1]] = parse(Float64, tok[4]))
                end
            end
        catch err
            @warn "alfvenWaveBabbarRanocha2026: skipping an unreadable error file" file=fname exception=err
            continue
        end
        tc = tryparse(Float64, get(meta, "t", ""))
        (tc === nothing || abs(tc - tnom) > 1.0e-8*max(1.0, abs(tnom))) && continue
        all(haskey(l2, n) for n in AW_QVARS) || continue
        nop  = tryparse(Int, get(meta, "nop", "")); nelx = tryparse(Int, get(meta, "nelx", ""))
        (nop === nothing || nelx === nothing) && continue
        # records stored before the wave's errors and the run maxima have NaN there
        push!(rows, (nop = nop, nelx = nelx, ndofs = (nelx*nop)^2, visc = get(meta, "visc", "dsgs"),
                     t = tc, dt = something(tryparse(Float64, get(meta, "dt", "")), NaN),
                     bperp = get(l2, "Bperp", NaN), uperp = get(l2, "uperp", NaN),
                     bperp_inf = get(linf, "Bperp", NaN), uperp_inf = get(linf, "uperp", NaN),
                     rho = l2["ρ"], state = sqrt(sum(abs2, l2[n] for n in AW_QVARS)),
                     state_inf = maximum(linf[n] for n in AW_QVARS),
                     rho_inf = linf["ρ"],
                     bperp_l1 = get(l1, "Bperp", NaN), uperp_l1 = get(l1, "uperp", NaN),
                     rho_l1 = get(l1, "ρ", NaN),
                     state_l1 = all(haskey(l1, n) for n in AW_QVARS) ?
                                sqrt(sum(abs2, l1[n] for n in AW_QVARS)) : NaN,
                     r_max  = something(tryparse(Float64, get(meta, "r_max",  "")), NaN),
                     mu_max = something(tryparse(Float64, get(meta, "mu_max", "")), NaN),
                     cutoff = something(tryparse(Float64, get(meta, "cutoff", "")), NaN)))
    end
    return sort!(rows, by = r -> (r.visc, r.nop, r.ndofs))
end

# Observed order between consecutive meshes: in 2D h ∝ #DOFs^(-1/2).
_aw_rate(d1, e1, d2, e2) = (e1 > 0 && e2 > 0) ? -2.0*log(e2/e1)/log(d2/d1) : NaN

function _aw_report(rows)
    println(" # alfvenWaveBabbarRanocha2026: L2 error at t = ", rows[1].t,
            " (B⊥, u⊥: the wave; rho: the paper's Fig. 13 measure, exact rho = 1;",
            " r, mu: largest normalized residual before the cutoff and largest DSGS coefficient over the run)")
    println(" #   visc         N  nelx    DOFs     L2(B⊥)  order     L1(B⊥)  order   Linf(B⊥)  order     L2(u⊥)  order    L2(rho)  order    L1(rho)  order       max r      max mu")
    _o(p) = isfinite(p) ? @sprintf("%.2f", p) : ""
    _g(x) = isfinite(x) ? @sprintf("%.3e", x) : "-"
    prev = nothing
    for r in rows
        same = prev !== nothing && prev.visc == r.visc && prev.nop == r.nop
        pb = same ? _aw_rate(prev.ndofs, prev.bperp, r.ndofs, r.bperp) : NaN
        pb1 = same ? _aw_rate(prev.ndofs, prev.bperp_l1, r.ndofs, r.bperp_l1) : NaN
        pr1 = same ? _aw_rate(prev.ndofs, prev.rho_l1, r.ndofs, r.rho_l1) : NaN
        pbi = same ? _aw_rate(prev.ndofs, prev.bperp_inf, r.ndofs, r.bperp_inf) : NaN
        pu = same ? _aw_rate(prev.ndofs, prev.uperp, r.ndofs, r.uperp) : NaN
        pr = same ? _aw_rate(prev.ndofs, prev.rho,   r.ndofs, r.rho)   : NaN
        println(@sprintf(" #   %-11s %2d  %4d  %6d  %10s  %5s %10s  %5s %10s  %5s %10s  %5s %10.3e  %5s %10s  %5s  %10s  %10s",
                         r.visc, r.nop, r.nelx, r.ndofs,
                         _g(r.bperp), _o(pb), _g(r.bperp_l1), _o(pb1), _g(r.bperp_inf), _o(pbi),
                         _g(r.uperp), _o(pu), r.rho, _o(pr), _g(r.rho_l1), _o(pr1),
                         _g(r.r_max), _g(r.mu_max)))
        prev = r
    end
    return nothing
end

# The figure (../smoothVortex/convergence_figure.jl). The line style of a tag:
# solid = DSGS as the deck runs it ("dsgs", or the first DSGS tag when there is
# none), dashed = Galerkin, and the DSGS variants (another cutoff, another
# integrator, :dsgs_freeze_stage, ...) the dash patterns left, in tag order;
# dotted is the reference slopes'.
const AW_VARIANT_LINES = (:dashdot, :dashdotdot)
function _aw_line(tag, tags)
    startswith(tag, "galerkin") && return :dash
    dsgs = sort([t for t in tags if !startswith(t, "galerkin")])
    base = "dsgs" in dsgs ? "dsgs" : first(dsgs)
    tag == base && return :solid
    variants = filter(!=(base), dsgs)
    return AW_VARIANT_LINES[mod1(something(findfirst(==(tag), variants), 1), length(AW_VARIANT_LINES))]
end
_aw_label(tag) = startswith(tag, "galerkin") ? replace(tag, "galerkin" => "Galerkin", "_" => " ") :
                 replace(tag, "dsgs" => "DSGS", "_nocut" => ", no cutoff",
                         r"_cut([0-9.e+-]+)" => s", cutoff \1", r"_bdf(\d+)" => s", BDF\1", "_" => " ")
# The error panels, a matrix: rows B⊥, u⊥, ρ and the state; columns L², L¹
# and L∞. The file key names the quantity and, but for L², the norm.
_aw_norm_tex(n) = n == "L2" ? "L^2" : n == "L1" ? "L^1" : "L^\\infty"
const AW_PANELS = [   # (field, ylabel, title, file)
    (Symbol(f, sfx),
     string("\\Vert", h, "\\Vert_{", _aw_norm_tex(n), "}"),
     string(_aw_norm_tex(n), "\\mathrm{-error\\ of", st ? "\\ the\\ state}" : "}\\ " * tt),
     n == "L2" ? key : string(key, "_", n == "L1" ? "L1" : "Linf"))
    for (f, h, tt, key, st) in (
        ("bperp", "\\mathbf{B}_{\\perp,h}-\\mathbf{B}_\\perp", "\\mathbf{B}_\\perp", "Bperp", false),
        ("uperp", "\\mathbf{u}_{\\perp,h}-\\mathbf{u}_\\perp", "\\mathbf{u}_\\perp", "uperp", false),
        ("rho",   "\\rho_h-\\rho",                             "\\rho",              "rho",   false),
        ("state", "\\mathbf{q}_h-\\mathbf{q}",                 "",                   "state", true))
    for (n, sfx) in (("L2", ""), ("L1", "_l1"), ("Linf", "_inf"))
]
# the figure of every variant side by side keeps the seven panels it had
const AW_ALL_PANELS = [p for p in AW_PANELS
                       if p[4] in ("Bperp", "Bperp_Linf", "uperp", "uperp_Linf", "rho", "state", "state_Linf")]
# the panels of tools/alfven_wave_plot.jl --panels=L2, L1, Linf: B⊥, u⊥ and ρ in that norm
const AW_PANEL_SETS = Dict(n => [p[4] for p in AW_PANELS
                                 if p[4] in (n == "L2" ? ("Bperp", "uperp", "rho") :
                                             (string("Bperp_", n), string("uperp_", n), string("rho_", n)))]
                           for n in ("L2", "L1", "Linf"))

#---------------------------------------------------------------------------------
# WHAT DSGS DID, against #DOFs, from the run maxima every record stores:
#
#   :r_max   the largest normalized residual r of the run, BEFORE the cutoff,
#            with the cutoff c drawn as a horizontal line. A point below the line
#            is a run in which μ = C_R Δ² max(0, r − c) was exactly zero
#            throughout: there DSGS is Galerkin because of the cutoff.
#   :mu_max  the largest coefficient μ of the run. μ = 0 has no place on a log
#            axis, so each curve's legend entry says on how many of its meshes μ
#            was exactly zero, and a curve that was zero on all of them has only
#            its entry. The dotted guide is O(h⁴), the scaling of the BDF2 floor
#            of the residual that sets μ without the cutoff (README, "What this
#            test found", item 4).
#
# One colour and marker per order and one line style per variant, as in the
# error panels. Galerkin records (r = NaN, μ = 0) are left out, and so are the
# records stored before the run maxima existed (r = NaN).
#---------------------------------------------------------------------------------
_aw_tex_num(c) = (e = log10(c); isapprox(e, round(e); atol = 1.0e-9) ?
                  string("10^{", round(Int, e), "}") : @sprintf("%.1e", c))

function _aw_diag_axes!(pl, xs, ys)
    x1, x2 = extrema(xs)
    xu  = sort(unique(xs))
    idx = length(xu) <= 4 ? eachindex(xu) : unique(round.(Int, range(1, length(xu); length = 4)))
    xt  = xu[idx]
    Plots.plot!(pl; xlims = (0.75*x1, 1.35*x2),
                xticks = (xt, [LaTeXStrings.latexstring(@sprintf("%d", round(Int, x))) for x in xt]))
    if isempty(ys)                     # μ panel with μ = 0 in every run drawn
        Plots.plot!(pl; ylims = (1.0e-8, 1.0e-2))
        Plots.annotate!(pl, sqrt(x1*x2), 1.0e-5,
                        Plots.text(LaTeXStrings.L"\mu = 0\ \mathrm{in\ every\ run}", CONV_FS_LEGEND, :gray25))
        return pl
    end
    ylo = 0.2*minimum(ys); yhi = 5.0*maximum(ys)
    decs = collect(ceil(Int, log10(ylo)):floor(Int, log10(yhi)))
    stp  = max(1, cld(length(decs), 7))
    Plots.plot!(pl; ylims = (ylo, yhi),
                yticks = ([10.0^e for e in decs[1:stp:end]],
                          [LaTeXStrings.latexstring(string("10^{", e, "}")) for e in decs[1:stp:end]]))
    return pl
end

function _aw_diag_panel(rows, fld::Symbol, ylab, title; tags, label, line, galerkin)
    dtags = [t for t in tags if !galerkin(t)]
    rs    = filter(r -> r.visc in dtags && isfinite(getfield(r, fld)), rows)
    isempty(rs) && return nothing
    lmu   = fld == :mu_max
    pl = Plots.plot(; xscale = :log10, yscale = :log10,
                    xlabel = LaTeXStrings.L"\#\mathrm{DOFs}",
                    ylabel = LaTeXStrings.latexstring(ylab),
                    framestyle = :box, grid = true, gridalpha = 0.25,
                    background_color_legend = Plots.RGBA(1, 1, 1, 0.72),
                    foreground_color_legend = :gray60,
                    legend = :outertopright, legendfontsize = CONV_FS_LEGEND,
                    titlefontsize = CONV_FS_TITLE, guidefontsize = CONV_FS_GUIDE,
                    tickfontsize = CONV_FS_TICK,
                    left_margin = 14Plots.mm, bottom_margin = 10Plots.mm,
                    top_margin = 4Plots.mm, right_margin = 6Plots.mm,
                    title = LaTeXStrings.latexstring(string(title, ",\\ t = ", rs[1].t)),
                    show = false)
    allx = Float64[]; ally = Float64[]
    base(t) = line(t, dtags) == :solid
    nops = sort(unique(r.nop for r in rs))
    curves = Dict{Tuple{Int,String},Tuple{Vector{Float64},Vector{Float64}}}()
    for nop in nops, tag in dtags
        g = sort(filter(r -> r.nop == nop && r.visc == tag, rs), by = r -> r.ndofs)
        isempty(g) || (curves[(nop, string(tag))] = ([Float64(r.ndofs) for r in g], [Float64(getfield(r, fld)) for r in g]))
    end
    attrs(nop, tag) = (c = _conv_style(nop);
                       (line = (c[1], CONV_LW, line(tag, dtags)),
                        marker = (c[2], base(tag) ? CONV_MS : CONV_MS + 2),
                        markerstrokecolor = c[1], markerstrokewidth = 1.0,
                        markercolor = base(tag) ? c[1] : :white, color = c[1]))
    # drawn: the variants first, as hollow markers, and the deck's DSGS last,
    # filled, so that where the two coincide it shows inside the ring
    for nop in nops, tag in vcat(filter(t -> !base(t), dtags), filter(base, dtags))
        haskey(curves, (nop, string(tag))) || continue
        xs, ys = curves[(nop, string(tag))]
        pos = ys .> 0
        any(pos) || continue
        Plots.plot!(pl, xs[pos], ys[pos]; label = "", attrs(nop, tag)...)
        append!(allx, xs[pos]); append!(ally, ys[pos])
    end
    # listed: the deck's DSGS first; an entry also for a curve that is zero on every mesh
    for nop in nops, tag in vcat(filter(base, dtags), filter(t -> !base(t), dtags))
        haskey(curves, (nop, string(tag))) || continue
        xs, ys = curves[(nop, string(tag))]
        pos  = ys .> 0
        at   = attrs(nop, tag)
        p    = count(pos) >= 2 ? _conv_rate(xs[pos], ys[pos]) : NaN
        nz   = count(!, pos)
        note = (!lmu || nz == 0) ? "" :
               nz == length(ys) ? ",\\ \\mu\\equiv 0\\ \\mathrm{on\\ every\\ mesh}" :
               string(",\\ \\mu\\equiv 0\\ \\mathrm{on\\ ", nz, "\\ of\\ ", length(ys), "\\ meshes}")
        lab  = LaTeXStrings.latexstring(string("\\mathbb{P}_", nop, "\\ \\mathrm{", _conv_mathrm(label(tag)), "}",
                                               isfinite(p) ? string("\\ (p=", round(p; digits = 2), ")") : "", note))
        Plots.plot!(pl, [NaN], [NaN]; label = lab, at...)
    end
    if !lmu
        # the cutoff of the runs drawn: r below the line means μ ≡ 0 for that run
        for c in sort(unique(r.cutoff for r in rs if isfinite(r.cutoff) && r.cutoff > 0))
            Plots.hline!(pl, [c]; line = (:gray25, 2.0, :dash),
                         label = LaTeXStrings.latexstring(string("\\mathrm{cutoff}\\ c = ", _aw_tex_num(c))))
            push!(ally, c)
        end
    elseif length(allx) >= 2
        # O(h⁴), drawn above every point: μ ≈ C_R Δ² (Δt²/3)|∂ₜ³q|/(ρ̄c̄) at a fixed Courant number
        A = 2.0*maximum(ally[i]*allx[i]^2 for i in eachindex(allx))
        x1, x2 = extrema(allx)
        Plots.plot!(pl, [x1, x2], [A/x1^2, A/x2^2]; line = (:gray40, 1.2, :dot),
                    label = LaTeXStrings.L"\mathcal{O}(h^4)")
    end
    _aw_diag_axes!(pl, isempty(allx) ? Float64[r.ndofs for r in rs] : allx, ally)
    return pl
end

#---------------------------------------------------------------------------------
# THE TWO COMPARISONS, each one pair of tags alone, with the reference slopes of
# the October 2026 figure: one dotted h^(N+1) guide PER ORDER, in the order's
# colour, over the order's widest curve, a factor 3 below its finest error and
# labelled at its end.
#   <name>-it<iout>.png               the deck's DSGS and Galerkin, and each
#                                     panel alone as <name>_<file>-it<iout>.png
#   <name>_dsgs_nocut-it<iout>.png    the deck's DSGS and DSGS without the cutoff
# The panels of AW_PANELS in a grid four wide, the last cell the legend. A
# comparison is drawn when each of its two tags has a curve of two meshes or
# more. The styles are those of the error panels: the deck's DSGS solid and
# filled, drawn last so that it shows inside a hollow marker it coincides with;
# Galerkin dashed with a larger hollow marker; the variant dash-dot.
#---------------------------------------------------------------------------------
const AW_COMPARISONS = [   # (tags, file suffix, title, each panel alone too)
    (["dsgs", "galerkin"],   "",           "DSGS (the deck) and Galerkin",                true),
    (["dsgs", "dsgs_nocut"], "dsgs_nocut", "DSGS (the deck) and DSGS without the cutoff", false),
]
const AW_INK = Plots.RGB(0x52/255, 0x51/255, 0x4e/255)   # the guides' labels
const AW_CMP_FS = (title = 15, guide = 15, tick = 13, legend = 14, note = 13)
const AW_CMP_CELL = (820, 640)

function _aw_curve(rows, fld, nop, tag)              # sorted by #DOFs, the positive finite errors
    g  = sort(filter(r -> r.nop == nop && r.visc == tag, rows), by = r -> r.ndofs)
    xs = [Float64(r.ndofs) for r in g]; ys = [Float64(getfield(r, fld)) for r in g]
    keep = isfinite.(ys) .& (ys .> 0)
    return xs[keep], ys[keep]
end

function _aw_cmp_attrs(nop, tag, tags)
    col, mk = _conv_style(nop)
    ls  = _aw_line(tag, tags)
    filled = ls == :solid
    return (line = (col, 2.0, ls),
            marker = (mk, startswith(tag, "galerkin") ? 9 : filled ? 5 : 7),
            markercolor = filled ? col : :white, markerstrokecolor = col,
            markerstrokewidth = 1.2, color = col)
end

function _aw_cmp_panel(rows, fld, ylab, title, tags)
    pl = Plots.plot(; xscale = :log10, yscale = :log10,
                    xlabel = "#DOFs = (nelx N)²", ylabel = LaTeXStrings.latexstring(ylab),
                    title = title, framestyle = :box, grid = true, gridalpha = 0.25,
                    legend = false, titlefontsize = AW_CMP_FS.title,
                    guidefontsize = AW_CMP_FS.guide, tickfontsize = AW_CMP_FS.tick, show = false)
    # the resolutions run, labelled as the paper labels them (32², 64², ...); with
    # too many to label, the powers of two
    dofs = sort(unique(Float64(r.ndofs) for r in rows))
    xlo, xhi = dofs[1]/1.25, dofs[end]*1.6           # room on the right for the guides' labels
    side = length(dofs) <= 8 ? sqrt.(dofs) :
           filter(s -> xlo <= s^2 <= xhi,
                  [2.0^k for k in floor(Int, log2(sqrt(dofs[1]))):ceil(Int, log2(sqrt(dofs[end])))])
    Plots.plot!(pl; xticks = (side.^2, [string(round(Int, s), "²") for s in side]),
                xlims = (xlo, xhi))
    drawn = sort(tags, by = t -> (_aw_line(t, tags) == :solid, t))   # the deck's DSGS last
    ally  = Float64[]
    for n in sort(unique(r.nop for r in rows))
        ends = Tuple{Vector{Float64},Float64}[]    # (#DOFs of a curve, its finest error)
        for tag in drawn
            xs, ys = _aw_curve(rows, fld, n, tag)
            isempty(xs) && continue
            append!(ally, ys)
            Plots.plot!(pl, xs, ys; label = "", _aw_cmp_attrs(n, tag, tags)...)
            length(xs) >= 2 && push!(ends, (xs, ys[end]))
        end
        isempty(ends) && continue
        # O(h^(N+1)) = #DOFs^(-(N+1)/2)
        xs = ends[argmax([length(e[1]) for e in ends])][1]
        sl = (n + 1)/2
        A  = minimum(e[2] for e in ends)/3.0*xs[end]^sl
        y2 = A*xs[end]^(-sl)
        Plots.plot!(pl, [xs[1], xs[end]], [A*xs[1]^(-sl), y2];
                    line = (_conv_style(n)[1], 1.0, :dot), label = "")
        Plots.annotate!(pl, xs[end]*1.08, y2,
                        Plots.text(LaTeXStrings.latexstring("h^{", n + 1, "}"), AW_INK, :left,
                                   AW_CMP_FS.note))
        push!(ally, y2)
    end
    if isempty(ally)              # nothing stored for this panel: L¹ before 23 Nov 2026
        Plots.plot!(pl; ylims = (1.0e-12, 1.0e-2))
        Plots.annotate!(pl, sqrt(xlo*xhi), 1.0e-7,
                        Plots.text(endswith(string(fld), "_l1") ?
                                   "no L¹ errors stored:\nruns from 23 Nov 2026 on store them" :
                                   "nothing stored", AW_INK, :center, AW_CMP_FS.note))
    end
    # decades (10^-7.5 is not a number anyone reports), set by the errors and the
    # guides' ends: a guide that starts above the errors is cut
    if !isempty(ally)
        lo = floor(Int, log10(minimum(ally))); hi = ceil(Int, log10(maximum(ally)))
        stp = max(1, cld(hi - lo + 1, 8))
        Plots.plot!(pl; ylims = (10.0^lo/1.5, 10.0^hi*1.5),
                    yticks = ([10.0^e for e in lo:stp:hi],
                              [LaTeXStrings.latexstring(string("10^{", e, "}")) for e in lo:stp:hi]))
    end
    return pl
end

# the legend entries of a comparison, as empty series (they draw nothing else)
function _aw_cmp_legend!(pl, rows, tags)
    listed = sort(tags, by = t -> (_aw_line(t, tags) != :solid, t))  # the deck's DSGS first
    for n in sort(unique(r.nop for r in rows)), tag in listed
        any(r -> r.nop == n && r.visc == tag, rows) || continue
        Plots.plot!(pl, [NaN], [NaN]; label = string("N = ", n, ", ", _aw_label(tag)),
                    _aw_cmp_attrs(n, tag, tags)...)
    end
    Plots.plot!(pl, [NaN], [NaN]; line = (AW_INK, 1.0, :dot),
                label = "h^(N+1), dotted, one per order")
    return pl
end
_aw_cmp_legend(rows, tags) =
    _aw_cmp_legend!(Plots.plot(; framestyle = :none, grid = false, xticks = false, yticks = false,
                               xlims = (0, 1), ylims = (0, 1), legend = :left,
                               legendfontsize = AW_CMP_FS.legend, foreground_color_legend = nothing,
                               background_color_legend = nothing, show = false), rows, tags)

# Returns the files written. `panels`, a selection of AW_PANELS, draws the
# comparisons with those panels alone, as <stem>_panels_<label>-it<iout>.png
# beside the full figures (tools/alfven_wave_plot.jl --panels); `label`
# defaults to the panels' file keys joined by "-".
function _aw_compare(rows, OUTPUT_DIR, iout; name = AW_FIG_NAME, panels = AW_PANELS,
                     label = join(last.(panels), "-"))
    full  = length(panels) == length(AW_PANELS)
    files = String[]
    for (tags, suffix, title, singles) in AW_COMPARISONS
        rs = filter(r -> r.visc in tags, rows)
        nops = unique(r.nop for r in rs)
        if !all(any(count(r -> r.visc == tag && r.nop == n, rs) >= 2 for n in nops) for tag in tags)
            println(" # no figure \"", title, "\": it needs ", join(tags, " and "),
                    " each on two meshes of one order")
            continue
        end
        stem = isempty(suffix) ? name : string(name, "_", suffix)
        full || (stem = string(stem, "_panels_", label))
        pls = Any[_aw_cmp_panel(rs, fld, ylab,
                                LaTeXStrings.latexstring(string("\\mathrm{(", 'a' + (i - 1), ")}\\ ", ttl)),
                                tags)
                  for (i, (fld, ylab, ttl, _)) in enumerate(panels)]
        if singles && full                             # each panel alone, with its legend
            for (i, (fld, ylab, ttl, fname)) in enumerate(AW_PANELS)
                p1 = _aw_cmp_legend!(_aw_cmp_panel(rs, fld, ylab, LaTeXStrings.latexstring(
                                         string(ttl, ",\\ t = ", rs[1].t)), tags), rs, tags)
                Plots.plot!(p1; legend = :outertopright, legendfontsize = AW_CMP_FS.legend)
                _conv_savefig(Plots.plot(p1; size = (1180, 700), left_margin = 8Plots.mm,
                                         bottom_margin = 6Plots.mm, show = false),
                              string(OUTPUT_DIR, "/", stem, "_", fname, "-it", iout, ".png"))
            end
        end
        # the panels in a grid (one row up to four; else three wide when they
        # come in threes, as the quantity × norm matrix does, four otherwise), and
        # the legend a narrower column on the right
        np = length(pls)
        nc = np <= 4 ? np : (np % 3 == 0 ? 3 : 4); nr = cld(np, nc)
        while length(pls) < nr*nc                      # an empty cell keeps the grid
            push!(pls, Plots.plot(; framestyle = :none, grid = false, show = false))
        end
        push!(pls, _aw_cmp_legend(rs, tags))
        wl  = 0.55                                     # the legend column, in panel widths
        lay = Plots.grid(1, 2; widths = [nc/(nc + wl), wl/(nc + wl)])
        lay[1, 1] = Plots.grid(nr, nc)
        sz  = (round(Int, (nc + wl)*AW_CMP_CELL[1]), nr*AW_CMP_CELL[2] + 60)
        plt = Plots.plot(pls...; layout = lay, size = sz,
                         plot_title = string(title, ", t = ", rs[1].t), plot_titlefontsize = 19,
                         left_margin = 12Plots.mm, bottom_margin = 10Plots.mm,
                         top_margin = 3Plots.mm, right_margin = 3Plots.mm, show = false)
        f = string(OUTPUT_DIR, "/", stem, "-it", iout, ".png")
        _conv_savefig(plt, f)
        push!(files, f)
    end
    return files
end

# <name>-it<iout>.png, the deck's DSGS and Galerkin in the seven error panels,
# and <name>_Bperp, _Bperp_Linf, _uperp, _uperp_Linf, _rho, _L2 and _Linf, one
# each; <name>_dsgs_nocut, the deck's DSGS and DSGS without the cutoff
# (_aw_compare); <name>_all, every variant stored in the seven panels side by
# side (_conv_figure). Then what DSGS did: <name>_rmax and
# _mumax (_aw_diag_panel), and <name>_paper, the three panels the paper's
# smooth-wave test needs side by side: the error of the wave (B⊥) for every
# variant, max r against the cutoff, and max μ. Returns the first file, or
# nothing when no order/tag has two meshes yet.
function _aw_plot(rows, OUTPUT_DIR, iout; name = AW_FIG_NAME)
    kw = (tags = sort(unique(r.visc for r in rows)), label = _aw_label, line = _aw_line,
          galerkin = t -> startswith(t, "galerkin"))
    # every variant stored, side by side (the shared convergence figure)
    fall = _conv_figure(rows, AW_ALL_PANELS, OUTPUT_DIR, iout; prefix = string(name, "_all"),
                        singles = false, kw...)
    fall === nothing && return nothing
    cmp = _aw_compare(rows, OUTPUT_DIR, iout; name = name)
    f = isempty(cmp) ? fall : first(cmp)
    pr = _aw_diag_panel(rows, :r_max,  "\\max_t\\, r",
                        "\\mathrm{largest\\ residual\\ before\\ the\\ cutoff}"; kw...)
    pm = _aw_diag_panel(rows, :mu_max, "\\max_t\\, \\mu",
                        "\\mathrm{largest\\ DSGS\\ coefficient}"; kw...)
    (pr === nothing || pm === nothing) && return f       # no DSGS record with the run maxima
    out(s) = string(OUTPUT_DIR, "/", name, "_", s, "-it", iout, ".png")
    _conv_savefig(Plots.plot(pr; size = CONV_PANEL_SIZE, show = false), out("rmax"))
    _conv_savefig(Plots.plot(pm; size = CONV_PANEL_SIZE, show = false), out("mumax"))
    pb = _conv_panel(rows, AW_PANELS[1][1], AW_PANELS[1][2], AW_PANELS[1][3]; kw...)
    _conv_savefig(Plots.plot(pb, pr, pm; layout = (1, 3),
                             size = (3*CONV_PANEL_SIZE[1], CONV_PANEL_SIZE[2]),
                             left_margin = 22Plots.mm, bottom_margin = 20Plots.mm,
                             top_margin = 13Plots.mm, right_margin = 6Plots.mm, show = false),
                  out("paper"))
    return f
end

function user_plot_2d(mesh, q, t, outvar, inputs, OUTPUT_DIR, iout; Minv = nothing)
    isfinite(t) && t > 0.0 || return nothing          # a final-time product
    # The time the run was ASKED for keys the store (it is the same for every
    # run of a sweep); the exact solution is placed at the time the numerical
    # one is AT, in case the two differ by a fraction of a step.
    tnom = Float64(get(inputs, :tend, t))
    Δt   = abs(Float64(get(inputs, :Δt, 0.0)))
    abs(t - tnom) <= max(1.0e-8*max(1.0, abs(t)), 0.51*Δt) || return nothing

    e = _aw_errors(mesh, q, t, outvar)                 # collective
    e === nothing && return nothing
    # what DSGS did up to here, the whole run (SGS.jl): [max r, max μ]; collective
    st = get(inputs, :lvisc, false) == true ? dsgs_run_stats(get_mpi_comm()) : [NaN, 0.0]
    if MPI.Comm_rank(get_mpi_comm()) == 0
        _aw_save(e, inputs, tnom, t, st)
        rows = _aw_load(tnom)
        isempty(rows) || (_aw_report(rows); _aw_plot(rows, OUTPUT_DIR, iout))
    end
    return nothing
end
