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
# 0.1 (sin, cos)(2π(x∥ + t)) in both. The figure has five panels: B⊥ and u⊥,
# the wave's own error; ρ, the paper's measure; the full state,
# sqrt(Σ_k e(q_k)²); and the L∞ error of the state.
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
                sm[k]  = max(sm[k], abs(e))
            end
            # the wave's transverse field and velocity against the exact ones
            eB2 = (Tq[a,b,1] - (-sα*se[6] + cα*se[7]))^2 + (Uq[a,b,8] - se[8])^2
            eU2 = (Tq[a,b,2] - (-sα*se[2] + cα*se[3])/se[1])^2 + (Tq[a,b,3] - se[5]/se[1])^2
            sp2[1] += wab*eB2;  spm[1] = max(spm[1], sqrt(eB2))
            sp2[2] += wab*eU2;  spm[2] = max(spm[2], sqrt(eU2))
            area += wab
        end
    end

    comm = get_mpi_comm()
    if MPI.Comm_size(comm) > 1
        r = MPI.Allreduce(vcat(s2, sp2, area), MPI.SUM, comm)
        s2 .= r[1:nv]; sp2 .= r[nv+1:nv+2]; area = r[end]
        sm  .= MPI.Allreduce(sm,  MPI.MAX, comm)
        spm .= MPI.Allreduce(spm, MPI.MAX, comm)
    end
    # the quadrature must integrate the box: a cheap check that the element
    # geometry (and the unwrapping above) is right
    abs(area - Lx*Ly) <= 1.0e-10*Lx*Ly ||
        @warn "alfvenWaveBabbarRanocha2026: the error quadrature integrates $(area), not the box area $(Lx*Ly)"
    return (l2 = sqrt.(s2 ./ area), linf = sm,
            perp_l2 = sqrt.(sp2 ./ area), perp_linf = spm,    # (B⊥, u⊥)
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
    get(inputs, :lvisc, true) || return string("galerkin", ode)
    tag    = "dsgs"
    legacy = lowercase(string(get(inputs, :dsgs_sensor, "residual"))) == "legacy"
    legacy && (tag *= "_legacy")
    Float64(get(inputs, :dsgs_Cmin, 0.0)) > 0 && (tag *= "_cmin")
    cut = Float64(get(inputs, :dsgs_cutoff, 0.0))
    cut == AW_DSGS_CUTOFF || (tag *= cut > 0 ? @sprintf("_cut%.0e", cut) : "_nocut")
    !legacy && get(inputs, :dsgs_freeze_stage, false) != true && (tag *= "_perstage")
    return string(tag, ode)
end

_aw_nelx_of(inputs) = (m = match(r"alfven_(\d+)x", string(get(inputs, :gmsh_filename, "")));
                       m === nothing ? Int(get(inputs, :nelx, 0)) : parse(Int, m.captures[1]))

# st = [max r, max μ] over the run (dsgs_run_stats; NaN, 0 for Galerkin)
function _aw_save(e, inputs, tnom, tnum, st)
    nop  = Int(get(inputs, :nop, 0))
    nelx = _aw_nelx_of(inputs)
    try
        f = joinpath(AW_ERR_DIR, string("nop", nop, "_nelx", nelx, "_", _aw_tag(inputs), ".dat"))
        _aw_atomic(f, tmp -> open(tmp, "w") do io
            println(io, "# smooth Alfvén wave: e(q) = sqrt((1/|Ω|) ∫ (q_h − q)² dΩ) against the exact solution,",
                        " N+4 Gauss points per direction; Linf at those points.",
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
            println(io, "# var L2 Linf")
            for (k, n) in enumerate(AW_QVARS)
                println(io, n, " ", e.l2[k], " ", e.linf[k])
            end
            println(io, "Bperp ", e.perp_l2[1], " ", e.perp_linf[1])
            println(io, "uperp ", e.perp_l2[2], " ", e.perp_linf[2])
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
                     rho = l2["ρ"], state = sqrt(sum(abs2, l2[n] for n in AW_QVARS)),
                     state_inf = maximum(linf[n] for n in AW_QVARS),
                     r_max  = something(tryparse(Float64, get(meta, "r_max",  "")), NaN),
                     mu_max = something(tryparse(Float64, get(meta, "mu_max", "")), NaN)))
    end
    return sort!(rows, by = r -> (r.visc, r.nop, r.ndofs))
end

# Observed order between consecutive meshes: in 2D h ∝ #DOFs^(-1/2).
_aw_rate(d1, e1, d2, e2) = (e1 > 0 && e2 > 0) ? -2.0*log(e2/e1)/log(d2/d1) : NaN

function _aw_report(rows)
    println(" # alfvenWaveBabbarRanocha2026: L2 error at t = ", rows[1].t,
            " (B⊥, u⊥: the wave; rho: the paper's Fig. 13 measure, exact rho = 1;",
            " r, mu: largest normalized residual before the cutoff and largest DSGS coefficient over the run)")
    println(" #   visc         N  nelx    DOFs     L2(B⊥)  order     L2(u⊥)  order    L2(rho)  order       max r      max mu")
    _o(p) = isfinite(p) ? @sprintf("%.2f", p) : ""
    _g(x) = isfinite(x) ? @sprintf("%.3e", x) : "-"
    prev = nothing
    for r in rows
        same = prev !== nothing && prev.visc == r.visc && prev.nop == r.nop
        pb = same ? _aw_rate(prev.ndofs, prev.bperp, r.ndofs, r.bperp) : NaN
        pu = same ? _aw_rate(prev.ndofs, prev.uperp, r.ndofs, r.uperp) : NaN
        pr = same ? _aw_rate(prev.ndofs, prev.rho,   r.ndofs, r.rho)   : NaN
        println(@sprintf(" #   %-11s %2d  %4d  %6d  %10s  %5s %10s  %5s %10.3e  %5s  %10s  %10s",
                         r.visc, r.nop, r.nelx, r.ndofs,
                         _g(r.bperp), _o(pb), _g(r.uperp), _o(pu), r.rho, _o(pr),
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
                         r"_cut([0-9.e+-]+)" => s", cutoff \1", "_" => " ")
const AW_PANELS = [   # (field, ylabel, title, file)
    (:bperp,     "\\Vert\\mathbf{B}_{\\perp,h}-\\mathbf{B}_\\perp\\Vert_{L^2}", "L^2\\mathrm{-error\\ of}\\ \\mathbf{B}_\\perp", "Bperp"),
    (:uperp,     "\\Vert\\mathbf{u}_{\\perp,h}-\\mathbf{u}_\\perp\\Vert_{L^2}", "L^2\\mathrm{-error\\ of}\\ \\mathbf{u}_\\perp", "uperp"),
    (:rho,       "\\Vert\\rho_h-\\rho\\Vert_{L^2}",                  "L^2\\mathrm{-error\\ of}\\ \\rho",               "rho"),
    (:state,     "\\Vert\\mathbf{q}_h-\\mathbf{q}\\Vert_{L^2}",      "L^2\\mathrm{-error\\ of\\ the\\ state}",       "L2"),
    (:state_inf, "\\Vert\\mathbf{q}_h-\\mathbf{q}\\Vert_{L^\\infty}", "L^\\infty\\mathrm{-error\\ of\\ the\\ state}", "Linf"),
]

# <name>-it<iout>.png, the five panels side by side, and <name>_Bperp, _uperp,
# _rho, _L2 and _Linf, one each. Returns the first, or nothing when no
# order/tag has two meshes yet.
_aw_plot(rows, OUTPUT_DIR, iout; name = "convergence") =
    _conv_figure(rows, AW_PANELS, OUTPUT_DIR, iout; prefix = name,
                 tags = sort(unique(r.visc for r in rows)), label = _aw_label, line = _aw_line,
                 galerkin = t -> startswith(t, "galerkin"))

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
