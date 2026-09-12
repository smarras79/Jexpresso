#---------------------------------------------------------------------------------
# Convergence history of the smooth MHD vortex, in the layout of
# Dao & Nazarov, J. Sci. Comput. 92:77 (2022), Fig. 1: the ABSOLUTE error of
# the velocity against the EXACT solution — ∫|u_h − u|dΩ, sqrt(∫|u_h − u|²dΩ),
# max|u_h − u| — against 1/sqrt(#DOFs) on log-log axes, one line per
# polynomial order, with slope guides at min(nop)+1 and max(nop)+1.
#
# At the final time every run stores its error in
# `errors/nop<N>_nelx<M>_<visc>.dat` in this case directory, and the figure is
# drawn from EVERY error stored there — so a sweep over orders and meshes
# builds the whole figure and each run replaces only its own point:
#
#     tools/smooth_vortex_mesh.sh     # the meshes, once
#     tools/smooth_vortex_scan.sh     # the sweep
#
# `rm -r problems/MHD/smoothVortex/errors` starts a fresh comparison. Only
# errors from the same final time are drawn together.
#
# The hook is src/io/plotting/jeplots.jl (plot_triangulation, NSD_2D):
#   mesh     this rank's mesh (coordinates, connectivity, extents)
#   q        flat npoin*nvar vector of the OUTPUT variables (`outvar`)
#   t        simulation time
#   Minv     the solver's assembled inverse lumped mass — the quadrature
#            weights the norms below are taken with
#---------------------------------------------------------------------------------
const SV_ERR_DIR = joinpath(@__DIR__, "errors")

# ATOMIC WRITES. A sweep runs several cases at the same time (SV_JOBS in
# tools/smooth_vortex_mpi_scan.sh), every one of them redraws the same figures
# from the shared store as it finishes, and two processes writing one PNG
# leave a torn file. Write beside the target and rename: on POSIX the rename
# is atomic, so a reader sees either the old file or the new one, never half
# of each. The same for the stored error, which the other cases are reading
# while it is being written.
function _sv_atomic(path::AbstractString, write!::Function)
    mkpath(dirname(path))
    # keep the extension: Plots picks the format from it
    base, ext = splitext(path)
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

_sv_savefig(plt, path) = _sv_atomic(path, f -> _savefig_silent(plt, f))

# One (colour, marker) per order, fixed so an order looks the same from one
# figure to the next.
const SV_STYLE = Dict(
    1 => (:darkorange, :dtriangle),
    2 => (:goldenrod, :utriangle),
    3 => (:purple,    :rect),
    4 => (:seagreen,  :diamond),
    5 => (:royalblue, :circle),
    6 => (:crimson,   :star5),
    7 => (:black,     :xcross),
)
_sv_style(nop) = get(SV_STYLE, nop, (:gray, :cross))

_sv_tag(inputs) = (get(inputs, :lvisc, true) ? "dsgs" : "galerkin")

#---------------------------------------------------------------------------------
# Error of the velocity (u, v) against the exact solution: the vortex of
# initialize.jl translated by v₀t and wrapped periodically. Integral norms are
# taken with the nodal quadrature weights of the mesh, which for LGL nodes is
# the mass-matrix lumping the solver itself uses.
#---------------------------------------------------------------------------------
function _sv_velocity_error(mesh, q, t, outvar, inputs, Minv)
    names = string.(outvar)
    iu = findfirst(==("u"), names)
    iv = findfirst(==("v"), names)
    (iu === nothing || iv === nothing) && return nothing
    npoin = mesh.npoin

    γ  = γ_mhd
    Lx = mesh.xmax - mesh.xmin
    Ly = mesh.ymax - mesh.ymin

    # Nodal quadrature weights: the SOLVER'S OWN assembled lumped mass,
    # w_i = M_ii = 1/Minv_i, handed down by the plotting hook. That is the
    # same quadrature the DSS and the mass-matrix inversion of the run use, so
    # the norms below are the solver's own integrals and not a rule invented
    # here. On this mesh they sum to the domain area to 14 digits.
    (Minv !== nothing && length(Minv) >= npoin) ||
        error("smoothVortex: the error norms need the solver's lumped mass (Minv) from the plotting hook")
    # One weight per DEGREE OF FREEDOM. On a doubly periodic mesh two local
    # nodes can be the same unknown — the right column of the box is the left
    # one — and they carry the same assembled mass, so summing over all npoin
    # counts the periodic edges twice (the weights integrate 105 instead of
    # the box's 100 at 4x4). The solver's own local-to-global map,
    # mesh.ip2gip, is what says which nodes are the same unknown; use it
    # rather than any geometric test of ours.
    #
    # ON MORE THAN ONE RANK the same unknown also lives on every rank that
    # touches it, and its assembled mass is the GLOBAL one, so a node counted
    # by two ranks is counted twice in the reduction below. mesh.gip2owner[ip]
    # is the rank that owns local node ip — the same map the DSS assembler
    # uses — so each unknown is weighed exactly once, by its owner.
    rank   = MPI.Comm_rank(get_mpi_comm())
    lowner = length(mesh.gip2owner) >= npoin
    seen   = Set{eltype(mesh.ip2gip)}()
    w      = zeros(npoin)
    for ip = 1:npoin
        lowner && mesh.gip2owner[ip] != rank && continue
        g = mesh.ip2gip[ip]
        g in seen && continue
        push!(seen, g)
        w[ip] = 1.0/Minv[ip]
    end

    s1 = 0.0; s2 = 0.0; si = 0.0
    r1 = 0.0; r2 = 0.0; ri = 0.0
    for ip = 1:npoin
        xr = sv_wrap(mesh.x[ip] - SV_XC - SV_U0*t, Lx)
        yr = sv_wrap(mesh.y[ip] - SV_YC - SV_V0*t, Ly)
        se = sv_state(xr, yr, γ)
        ue = se[2]/se[1]; ve = se[3]/se[1]
        uh = q[(iu - 1)*npoin + ip]
        vh = q[(iv - 1)*npoin + ip]
        e  = sqrt((uh - ue)^2 + (vh - ve)^2)
        m  = sqrt(ue*ue + ve*ve)
        s1 += w[ip]*e;   s2 += w[ip]*e*e;   si = max(si, e)
        r1 += w[ip]*m;   r2 += w[ip]*m*m;   ri = max(ri, m)
    end

    # The weights must integrate the domain: a cheap check that the lumped
    # mass handed down is this rank's and complete.
    if get(ENV, "JEXPRESSO_SV_DEBUG", "") == "1"
        @info "smoothVortex quadrature check" sum_w = sum(w) area = Lx*Ly npoin = npoin nunique = length(seen)
    end
    comm = get_mpi_comm()
    ndofs = length(seen)
    if MPI.Comm_size(comm) > 1
        sums  = MPI.Allreduce([s1, s2, r1, r2, Float64(ndofs)], MPI.SUM, comm)
        maxs  = MPI.Allreduce([si, ri], MPI.MAX, comm)
        s1, s2, r1, r2 = sums[1], sums[2], sums[3], sums[4]
        si, ri = maxs[1], maxs[2]
        ndofs  = round(Int, sums[5])
    end
    # ABSOLUTE norms, as in the paper's Fig. 1:
    #   L¹   = ∫|u_h − u|dΩ,   L² = sqrt(∫|u_h − u|²dΩ),   L∞ = max|u_h − u|.
    # The relative ones go in the stored header, since they cost nothing and
    # say how big the error is against the solution it is measured on.
    return (ndofs = ndofs,
            l1   = s1,
            l2   = sqrt(s2),
            linf = si,
            r1   = s1/max(r1, eps()),
            r2   = sqrt(s2)/max(sqrt(r2), eps()),
            rinf = si/max(ri, eps()))
end

function _sv_save_error(e, inputs, t)
    nop = Int(get(inputs, :nop, 0))
    # :nelx carries mod_inputs' placeholder for a gmsh case, so take the
    # element count from the mesh file name.
    m    = match(r"vortex_(\d+)x", string(get(inputs, :gmsh_filename, "")))
    nelx = m === nothing ? Int(get(inputs, :nelx, 0)) : parse(Int, m.captures[1])
    # The unique unknowns of the doubly periodic square are exactly (nelx·N)²,
    # so h = L/(nelx·N) ∝ 1/sqrt(#DOFs) exactly and the abscissa cannot move
    # with the number of ranks. The counted value is the fallback.
    ndofs = (nelx > 0 && nop > 0) ? (nelx*nop)^2 : e.ndofs
    try
        mkpath(SV_ERR_DIR)
        f = joinpath(SV_ERR_DIR, string("nop", nop, "_nelx", nelx, "_", _sv_tag(inputs), ".dat"))
        _sv_atomic(f, tmp -> open(tmp, "w") do io
            println(io, "# smooth MHD vortex: ABSOLUTE velocity error against the exact solution")
            println(io, "# nop=", nop, " nelx=", nelx, " ndofs=", ndofs,
                        " t=", t, " visc=", _sv_tag(inputs),
                        " norm=abs",            # absolute norms: see _sv_load_errors
                        " Cmin=", Float64(get(inputs, :dsgs_Cmin, 0.0)),
                        " rel=", Float64(get(inputs, :dsgs_rel, 1.0)),
                        " dt=", Float64(get(inputs, :Δt, 0.0)))
            println(io, "# relative, for reference: L1=", e.r1, " L2=", e.r2, " Linf=", e.rinf)
            println(io, "# L1 L2 Linf")
            println(io, e.l1, " ", e.l2, " ", e.linf)
        end)
    catch err
        @warn "smoothVortex: could not store the error" exception=err
    end
    return nothing
end

function _sv_load_errors(t)
    rows = NamedTuple[]
    isdir(SV_ERR_DIR) || return rows
    for fname in readdir(SV_ERR_DIR)
        (endswith(fname, ".dat") && startswith(fname, "nop")) || continue
        meta = Dict{String,String}(); vals = Float64[]
        try
            for line in eachline(joinpath(SV_ERR_DIR, fname))
                if startswith(line, "#")
                    for tok in split(line)
                        occursin('=', tok) || continue
                        k, v = split(tok, '=', limit = 2); meta[k] = v
                    end
                    continue
                end
                isempty(strip(line)) && continue
                append!(vals, parse.(Float64, split(line)))
            end
        catch err
            @warn "smoothVortex: skipping an unreadable error file" file=fname exception=err
            continue
        end
        length(vals) >= 3 || continue
        nop   = tryparse(Int, get(meta, "nop", ""));   nop   === nothing && continue
        ndofs = tryparse(Int, get(meta, "ndofs", "")); ndofs === nothing && continue
        tc    = tryparse(Float64, get(meta, "t", ""))
        (tc === nothing || abs(tc - t) > 1.0e-8*max(1.0, abs(t))) && continue
        # Errors written before the norms became absolute have no norm= key;
        # they are a different quantity and are not drawn with these.
        get(meta, "norm", "") == "abs" || continue
        push!(rows, (nop = nop, ndofs = ndofs, visc = get(meta, "visc", "dsgs"),
                     nelx = something(tryparse(Int, get(meta, "nelx", "")), 0),
                     t    = something(tc, NaN),
                     dt   = something(tryparse(Float64, get(meta, "dt", "")), NaN),
                     l1 = vals[1], l2 = vals[2], linf = vals[3]))
    end
    return sort!(rows, by = r -> (r.visc, r.nop, r.ndofs))
end

_sv_rate(xs, ys) = (length(xs) < 2 || ys[end-1] <= 0 || ys[end] <= 0) ? NaN :
                   log(ys[end-1]/ys[end])/log(xs[end-1]/xs[end])

#---------------------------------------------------------------------------------
# The figure, in the layout of Dao & Nazarov (2022), Fig. 1: the ABSOLUTE
# velocity error against 1/sqrt(#DOFs) ∝ h on log-log axes, one line per
# polynomial order, one panel per norm, one figure per stabilization (DynSGS /
# plain Galerkin — their two panels). Written both as the three-panel figure
# and as one file per norm, at publication sizes.
#
# The slope guides bracket the orders on the figure: the steeper one is
# max(nop) + 1, the nominal rate of the highest order present, and the
# shallower one min(nop) + 1.
#---------------------------------------------------------------------------------
const SV_FS_TITLE  = 22
const SV_FS_GUIDE  = 20
const SV_FS_TICK   = 17
const SV_FS_LEGEND = 15
const SV_LW        = 2.8
const SV_MS        = 9

function _sv_panel(sub, nops, fld, nm, tag)
    # The final time the errors were measured at, and the step they were taken
    # with, belong ON the figure: a sweep cut short for a pipeline check
    # (SV_TEND=0.05) produces a perfectly plausible-looking set of flat lines,
    # and nothing in the picture would otherwise say so.
    tt  = isempty(sub) ? NaN : sub[1].t
    dts = unique(r.dt for r in sub)
    stamp = string(",\\ t = ", isfinite(tt) ? tt : "?",
                   length(dts) == 1 && isfinite(dts[1]) ? string(",\\ \\Delta t = ", dts[1]) : "")
    allx = Float64[]; ally = Float64[]
    pl = Plots.plot(; xscale = :log10, yscale = :log10,
                    xlabel = LaTeXStrings.L"1/\sqrt{\#\mathrm{DOFs}}",
                    ylabel = LaTeXStrings.latexstring(string(
                        "\\|\\mathbf{u}_h-\\mathbf{u}_{exact}\\|_{", nm, "}")),
                    framestyle = :box, grid = true, gridalpha = 0.25,
                    legend = :bottomright, legendfontsize = SV_FS_LEGEND,
                    titlefontsize = SV_FS_TITLE, guidefontsize = SV_FS_GUIDE,
                    tickfontsize = SV_FS_TICK,
                    left_margin = 12Plots.mm, bottom_margin = 10Plots.mm,
                    top_margin = 4Plots.mm, right_margin = 6Plots.mm,
                    title = LaTeXStrings.latexstring(string(nm, "\\mathrm{-error},\\ \\mathrm{",
                             tag == "dsgs" ? "RV\\ (DynSGS)" : "Galerkin", "}", stamp)),
                    show = false)
    for nop in nops
        g  = sort(filter(r -> r.nop == nop, sub), by = r -> r.ndofs)
        xs = [1.0/sqrt(r.ndofs) for r in g]
        ys = [getfield(r, fld) for r in g]
        keep = isfinite.(ys) .& (ys .> 0)
        any(keep) || continue
        xs = xs[keep]; ys = ys[keep]
        append!(allx, xs); append!(ally, ys)
        col, mk = _sv_style(nop)
        p = _sv_rate(xs, ys)
        lab = isfinite(p) ?
              LaTeXStrings.latexstring(string("\\mathrm{nop}\\ ", nop, "\\ (p=", round(p; digits = 2), ")")) :
              LaTeXStrings.latexstring(string("\\mathrm{nop}\\ ", nop))
        Plots.plot!(pl, xs, ys; line = (col, SV_LW, :solid), marker = (mk, SV_MS),
                    markerstrokewidth = 1.0, color = col, label = lab)
    end
    # Slope guides: the nominal rates of the lowest and the highest order on
    # the figure, N+1. Each is drawn ALONGSIDE the curve it annotates — over a
    # 6x span of h a slope-7 line covers five decades, so a guide anchored at
    # the corner of the axes leaves the panel at once — and the axes are then
    # clipped to the data.
    if !isempty(allx)
        x1 = minimum(allx); x2 = maximum(allx)
        ymax = maximum(ally); ymin = minimum(ally)
        for (sl, nop, shift, col) in ((minimum(nops) + 1, minimum(nops), 3.0,   :gray40),
                                      (maximum(nops) + 1, maximum(nops), 1/4.0, :gray55))
            g = sort(filter(r -> r.nop == nop, sub), by = r -> r.ndofs)
            ys = [getfield(r, fld) for r in g]
            xs = [1.0/sqrt(r.ndofs) for r in g]
            keep = isfinite.(ys) .& (ys .> 0)
            any(keep) || continue
            xs = xs[keep]; ys = ys[keep]
            # anchor at the geometric middle of that order's own curve
            xm = exp(sum(log, xs)/length(xs)); ym = shift*exp(sum(log, ys)/length(ys))
            Plots.plot!(pl, [x1, x2], [ym*(xi/xm)^sl for xi in (x1, x2)];
                        line = (col, 2.0, :dash),
                        label = LaTeXStrings.latexstring(string("\\mathrm{slope}\\ ", sl)))
        end
        # Ticks at the resolutions actually run, thinned to four: 10^-1.8 says
        # nothing, and one label per (order, mesh) pair overlaps.
        xu = sort(unique(allx))
        idx = length(xu) <= 4 ? eachindex(xu) :
              unique(round.(Int, range(1, length(xu); length = 4)))
        xt  = xu[idx]
        Plots.plot!(pl; xlims = (0.85*x1, 1.18*x2), ylims = (0.25*ymin, 4.0*ymax),
                    xticks = (xt, [LaTeXStrings.latexstring(@sprintf("%.3g", x)) for x in xt]))
    end
    return pl
end

# The orders to draw on their own, beside the all-orders figure:
# JEXPRESSO_SV_PLOT_NOPS="1 3" gives the P1-vs-P3 comparison of the paper's
# Fig. 1 in its own file, convergence_<visc>_nop1-3-it<n>.png.
function _sv_plot_nops()
    v = strip(get(ENV, "JEXPRESSO_SV_PLOT_NOPS", ""))
    isempty(v) && return Int[]
    return sort(unique(filter(!isnothing, tryparse.(Int, split(v, r"[,\s]+")))))
end

function _sv_plot(rows, OUTPUT_DIR, iout; only::Vector{Int} = Int[], suffix::String = "")
    for tag in unique(r.visc for r in rows)
        sub  = filter(r -> r.visc == tag, rows)
        isempty(only) || (sub = filter(r -> r.nop in only, sub))
        nops = sort(unique(r.nop for r in sub))
        any(n -> count(r -> r.nop == n, sub) >= 2, nops) || continue

        panels = Plots.Plot[]
        for (fld, nm, fname) in ((:l1, "L^1", "L1"), (:l2, "L^2", "L2"), (:linf, "L^\\infty", "Linf"))
            pl = _sv_panel(sub, nops, fld, nm, tag)
            # one file per norm, for the paper
            plt1 = Plots.plot(pl; size = (900, 780), show = false)
            _sv_savefig(plt1, string(OUTPUT_DIR, "/convergence_", tag, suffix, "_", fname, "-it", iout, ".png"))
            push!(panels, pl)
        end
        plt = Plots.plot(panels...; layout = (1, 3), size = (2400, 800),
                         left_margin = 16Plots.mm, bottom_margin = 14Plots.mm, show = false)
        _sv_savefig(plt, string(OUTPUT_DIR, "/convergence_", tag, suffix, "-it", iout, ".png"))
    end
    return nothing
end

# Every figure a run writes: all the orders in the store, and — when
# JEXPRESSO_SV_PLOT_NOPS asks for it — the chosen subset on its own axes.
function _sv_plot_all(rows, OUTPUT_DIR, iout)
    _sv_plot(rows, OUTPUT_DIR, iout)
    sel = _sv_plot_nops()
    length(sel) >= 1 || return nothing
    _sv_plot(rows, OUTPUT_DIR, iout; only = sel,
             suffix = string("_nop", join(sel, "-")))
    return nothing
end

function _sv_report(rows)
    println(" # smoothVortex: ABSOLUTE velocity error against the exact solution")
    println(" #   visc      nop  nelx   DOFs           L1           L2         Linf")
    for r in rows
        println(@sprintf(" #   %-8s  %3d  %4d %6d   %10.3e   %10.3e   %10.3e",
                         r.visc, r.nop, r.nelx, r.ndofs, r.l1, r.l2, r.linf))
    end
    return nothing
end

function user_plot_2d(mesh, q, t, outvar, inputs, OUTPUT_DIR, iout; Minv = nothing)
    isfinite(t) && t > 0.0 || return nothing         # the error figure is a final-time product
    abs(t - Float64(get(inputs, :tend, t))) < 1.0e-8*max(1.0, abs(t)) || return nothing

    # Collective: every rank must enter it (the norms are Allreduced).
    e = _sv_velocity_error(mesh, q, t, outvar, inputs, Minv)
    e === nothing && return nothing

    # One writer. The norms are identical on every rank after the reduction,
    # so letting them all write the same file is a race, not redundancy.
    if MPI.Comm_rank(get_mpi_comm()) == 0
        _sv_save_error(e, inputs, t)
        rows = _sv_load_errors(t)
        isempty(rows) || (_sv_report(rows); _sv_plot_all(rows, OUTPUT_DIR, iout))
    end
    return nothing
end
