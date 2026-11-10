#---------------------------------------------------------------------------------
# The figure of Dao & Nazarov (2022), Fig. 2: density against the reference
# solution with four zoom boxes (the foot of the fast rarefaction, the
# compound wave, the contact, the slow shock), in place of the generic
# multi-panel figure of the 1D plotter. Written as density-it<n>.png at every
# output time; the reference and the insets appear at the final time only,
# when user_analytic.jl supplies the reference. Set :plot_user => false in
# user_inputs.jl to get the generic panels of the 1D plotter instead.
#
# Written at EVERY output time, beside the density figure. EVERY figure of
# this case is vector graphics (see _bw_fig_ext): .svg by default, .pdf with
# JEXPRESSO_BW_FIGFMT=pdf, .png only if raster is asked for.
#
#   <var>-it<n>.svg     ONE FILE PER SOLUTION QUANTITY -- rho, u, v, w, p, Bx,
#                       By, Bz -- each on the full canvas against the
#                       reference where user_analytic.jl has one, with the
#                       zoom boxes of BW_FIELD_INSETS on u, v, p and By where
#                       the DynSGS solution leaves the reference, and the four
#                       of the paper on rho
#   mu_dsgs-it<n>.svg   the structure of the DynSGS coefficient: ρ and By on
#                       top, ν per equation slot below, so that where the
#                       coefficient fires can be read against the waves
#
# ORDER COMPARISON AND CONVERGENCE HISTORY. At the final time the run stores
# its density profile and its DynSGS coefficient ν in
# `curves/nop<N>_dof<M>_<form>.dat` in this case directory, and the figures are drawn from EVERY curve stored there whose
# <form> matches this run's — the coefficient form (element or nodal) and the
# length-scale convention ("H" = Dao & Nazarov's Δ_K/k; none = the paper's
# Δ_K/(k+1) = h_e/(N+1), the default of this deck),
# which are what make one run a different method from another. Curves of any
# other form are reported on stdout and left off the figure:
#
#   density_dof<M>-it<n>.png      ONE FILE PER RESOLUTION: every order that ran
#                                 at ~M degrees of freedom, against the
#                                 reference, with the zoom boxes. A sweep
#                                 leaves density_dof150, density_dof300, ...
#   density-it<n>.png             the same for the finest resolution in the
#                                 store, under the plain name
#   rho_By_dof<M>-it<n>.png       THE COMPARISON FIGURE, in the layout of the
#   rho_By-it<n>.png              Brio-Wu figure of Dao, Nazarov & Tomas: ρ on
#                                 the left, B_y on the right, every order
#                                 against the reference, with zoom boxes; ONE
#                                 FIGURE PER RESOLUTION, and the finest under
#                                 the plain name (never two resolutions in one
#                                 B_y panel: too many curves to tell which one
#                                 overshoots)
#   nu_dof<M>-it<n>.png, nu-it<n>.png, nu_all-it<n>.png
#                                 the DynSGS coefficient ν alone, of the same
#                                 curves (nu_all: every resolution)
#   density_nu*-it<n>.png         ρ and ν side by side in place of ρ and B_y,
#                                 with JEXPRESSO_BW_RIGHT=nu (the earlier layout)
#   convergence-it<n>.png         Dao & Nazarov Fig. 1 layout: the L¹, L² and
#                                 L∞ error against 1/#DOFs on log-log axes, one
#                                 line per order, with slope guides
#   convergence_smooth-it<n>.png  the same restricted to BW_SMOOTH_WINDOW, the
#                                 one smooth non-constant part of this solution
#
# A sweep over orders and resolutions therefore builds the whole comparison,
# and each run replaces only its own (order, DOFs) point. The scan CLEARS the
# store first (BW_KEEP=1 to accumulate), so that a leftover sweep cannot show
# up on the figure of a new comparison before that comparison has run:
#
#     tools/brio_wu_order_scan.sh          # 4 orders x 4 resolutions
#
# or by hand,
#
#     for D in 150 300 600 1200; do
#       for N in 4 5 6 7; do
#         JEXPRESSO_BW_DOFS=$D JEXPRESSO_BW_NOP=$N \
#             julia --project=. src/Jexpresso.jl MHD brioWu1d
#       done
#     done
#
# `rm -r problems/MHD/brioWu1d/curves` starts a fresh comparison. Curves whose
# stored final time differs from the current one are ignored, so a change of
# `:tend` cannot silently mix solutions from different times.
#
# Signature required by src/io/plotting/jeplots.jl (plot_results, NSD_1D):
#   x        node coordinates (npoin)
#   q        npoin × nvar output variables (columns named by `outvar`)
#   qref     the reference of user_analytic_solution (NaN where absent) or nothing
#   μ_nodes  DynSGS coefficient at the nodes or nothing
#---------------------------------------------------------------------------------
const BW_INSETS = [   # (x-range, y-range, inset position as fractions of the axes from the bottom-left: x, y, w, h)
    ((0.30, 0.33), (0.94, 1.00), (0.07, 0.44, 0.28, 0.30)),   # foot of the fast rarefaction
    ((0.40, 0.46), (0.66, 0.74), (0.16, 0.06, 0.28, 0.30)),   # compound wave
    ((0.56, 0.59), (0.20, 0.35), (0.66, 0.50, 0.28, 0.26)),   # contact
    ((0.625, 0.675), (0.11, 0.26), (0.70, 0.16, 0.28, 0.26)), # slow shock
]

# A window inside the fast rarefaction: smooth, and not one of the constant
# states, so it is the one place on this solution where the order of the
# scheme can show. The errors are reported there as well as over the whole tube.
const BW_SMOOTH_WINDOW = (0.33, 0.41)

# OUTPUT FORMAT. A raster figure of a shock tube is unreadable exactly where
# it matters: the whole point of the zoom boxes is to look closely at two or
# three nodes across a shock, and at that magnification a PNG is pixels. The
# figures are therefore written as VECTOR graphics, default SVG (every viewer
# and every browser opens it, and one of these is ~100-300 kB because a 1D
# figure holds a few thousand path points and nothing else); "pdf" for direct
# \includegraphics into the paper, "png" to go back to raster.
#
#   JEXPRESSO_BW_FIGFMT=pdf   or   :plot_format => "pdf" in user_inputs.jl
function _bw_fig_ext(inputs)
    e = lowercase(strip(get(ENV, "JEXPRESSO_BW_FIGFMT",
                            string(get(inputs, :plot_format, "svg")))))
    if !(e in ("svg", "pdf", "png"))
        @warn "brioWu1d: unknown figure format $(e); using svg." _module=nothing _file=nothing
        return "svg"
    end
    return e
end

# Every figure goes through a temporary name in the same directory and is
# renamed into place, so that a figure is never half-written when another
# session reads or rewrites it: tools/brio_wu_order_scan.sh with BW_JOBS > 1
# runs several orders at once into the same output directory.
function _bw_savefig(plt, f)
    root, ext = splitext(f)
    tmp = string(root, ".tmp", getpid(), ext)
    _savefig_silent(plt, tmp)
    mv(tmp, f; force = true)
    return nothing
end

# ASCII file names: a figure called "ρ-it5.svg" is a nuisance to reference from
# a Makefile, a shell loop or a LaTeX \includegraphics.
const BW_SLUG = Dict("ρ" => "rho", "u" => "u", "v" => "v", "w" => "w", "p" => "p",
                     "Bx" => "Bx", "By" => "By", "Bz" => "Bz")
_bw_slug(name) = get(BW_SLUG, string(name), string(name))

# ZOOM BOXES, PER VARIABLE: (x-range, y-range, inset position as fractions of
# the axes from the bottom-left: x, y, w, h).
#
# Every box sits where the DynSGS solution leaves the reference, which on this
# problem is the neighbourhood of the slow shock at x = 0.64 and of the
# compound wave at x = 0.47 -- the two structures the coefficient has to
# resolve. The inset positions are chosen over the EMPTY part of each frame,
# which differs from variable to variable because the curves do; they are
# fractions, so they follow the axis limits computed from the data.
#
# ρ keeps the four boxes of the paper's Fig. 2 (BW_INSETS above).
const BW_FIELD_INSETS = Dict{String, Vector{Any}}(
    "u"  => [((0.600, 0.685), ( 0.520,  0.665), (0.11, 0.56, 0.29, 0.29)),   # top of the plateau into the slow shock
             ((0.605, 0.700), (-0.300, -0.175), (0.64, 0.56, 0.29, 0.29))],  # the post-shock state
    "v"  => [((0.595, 0.745), (-0.300,  0.040), (0.11, 0.38, 0.29, 0.29)),   # the jump back up across the slow shock
             ((0.555, 0.670), (-1.700, -1.470), (0.64, 0.09, 0.29, 0.29))],  # the deep plateau into the shock
    "p"  => [((0.520, 0.745), ( 0.040,  0.165), (0.64, 0.50, 0.29, 0.29)),   # the post-shock plateau
             ((0.430, 0.510), ( 0.420,  0.790), (0.11, 0.09, 0.29, 0.29))],  # the compound wave
    "By" => [((0.585, 0.700), (-1.030, -0.440), (0.64, 0.54, 0.29, 0.29)),   # the slow shock
             ((0.430, 0.520), (-0.640,  0.640), (0.11, 0.09, 0.29, 0.29))],  # the compound wave
)
_bw_field_insets(name) = get(BW_FIELD_INSETS, string(name), Any[])

# Where the per-(order, resolution) curves accumulate.
# JEXPRESSO_BW_STORE names a subdirectory of it, so that sweeps that are
# different methods (tools/brio_wu_order_scan.sh: one per C_min) keep their
# curves apart: curves/Cmin0.06, curves/Cmin0.03, ...
_bw_curve_dir() = joinpath(@__DIR__, "curves", strip(get(ENV, "JEXPRESSO_BW_STORE", "")))

# One (colour, line style, marker) per order, fixed so that a given order
# always looks the same from one figure to the next. EVERY numerical curve is
# broken (dash/dot); the reference is the only solid line, so the curves are
# told apart by line type and not by colour alone.
const BW_STYLE = Dict(
    1 => (:darkorange, :dash,     :star5),
    2 => (:teal,       :dashdot,  :cross),
    3 => (:magenta,    :dot,      :hexagon),
    4 => (:red,      :dash,       :circle),
    5 => (:blue,     :dashdot,    :rect),
    6 => (:seagreen, :dot,        :diamond),
    7 => (:purple,   :dashdotdot, :utriangle),
)
_bw_style(nop) = get(BW_STYLE, nop, (:darkorange, :dash, :xcross))

# Element- and nodal-form runs are different methods, so they are stored
# apart and only drawn together when the figure says which is which.
# With JEXPRESSO_BW_HSCALE=nop the deck folds Dao & Nazarov's Δ_K/k into the
# coefficients as a factor (k+1)/k on C_max, C_min. Divide it back out so the
# figure reports the method's C_min and not the per-order number that
# implements it. `_bw_hscale()` is the deck's own switch (user_inputs.jl):
# "ngl", the default, is the paper's h_e/(N+1) and needs no factor; nothing
# outside this case reads it.
function _bw_hfac(inputs)
    _bw_hscale() == "nop" || return 1.0
    N = Int(get(inputs, :nop, 0))
    return N > 0 ? (N + 1)/N : 1.0
end

# The key under which a run's curve is stored, and the only key drawn on one
# figure: the coefficient form AND the length-scale convention. A curve from a
# run with different settings is a different method and is never silently
# mixed into the comparison — the "H" suffix marks Δ_K/k.
function _bw_form(inputs)
    f = get(inputs, :ldsgs_nodal, false) ? "nodal" : "elem"
    _bw_hscale() == "nop" && (f = string(f, "H"))
    return f
end
_bw_curve_file(nop, ndofs, form) =
    joinpath(_bw_curve_dir(), string("nop", nop, "_dof", ndofs, "_", form, ".dat"))

#---------------------------------------------------------------------------------
# Store this run's density, coefficient and B_y profiles. One file per (order,
# DOFs), so a re-run at the same order and resolution replaces its own point
# and leaves the others.
#---------------------------------------------------------------------------------
function _bw_save_curve(xs, ρs, νs, Bys, inputs, t, iout)
    nop = Int(get(inputs, :nop, 0))
    nop > 0 || return nothing
    try
        mkpath(_bw_curve_dir())
        f   = _bw_curve_file(nop, length(xs), _bw_form(inputs))
        tmp = string(f, ".tmp", getpid())
        open(tmp, "w") do io
            println(io, "# Brio-Wu density profile, written by user_plot.jl for the order comparison.")
            println(io, "# nop=", nop,
                        " ndofs=", length(xs),
                        " nelx=", get(inputs, :nelx, 0),
                        " t=", t,
                        " dt=", Float64(get(inputs, :Δt, 0.0)),
                        " Cmin=", round(Float64(get(inputs, :dsgs_Cmin, 0.0))/_bw_hfac(inputs); digits = 6),
                        " CR=", Float64(get(inputs, :dsgs_CR, 1.0)),
                        " Cmax=", Float64(get(inputs, :dsgs_Cmax, 0.5)),
                        " sensor=", string(get(inputs, :dsgs_sensor, "residual")),
                        " iout=", iout,
                        " form=", _bw_form(inputs))
            # nu: the DynSGS coefficient at the node (the ρu slot, which in
            # this deck's conserved form is the one ν of every equation)
            println(io, "# x rho nu By")
            for i in eachindex(xs)
                println(io, xs[i], " ", ρs[i], " ", isempty(νs) ? NaN : νs[i],
                        " ", isempty(Bys) ? NaN : Bys[i])
            end
        end
        # Renamed into place: tools/brio_wu_order_scan.sh runs orders in
        # parallel, and a session drawing the store must never read half a file.
        mv(tmp, f; force = true)
    catch err
        @warn "brioWu1d: could not store the density curve for the order comparison" exception=err
    end
    return nothing
end

#---------------------------------------------------------------------------------
# Read back every stored curve whose final time matches this one.
#---------------------------------------------------------------------------------
function _bw_load_curves(t)
    curves = NamedTuple[]
    isdir(_bw_curve_dir()) || return curves
    for fname in readdir(_bw_curve_dir())
        (endswith(fname, ".dat") && startswith(fname, "nop")) || continue
        meta = Dict{String,String}()
        xs = Float64[]; ys = Float64[]; νs = Float64[]; bys = Float64[]
        try
            for line in eachline(joinpath(_bw_curve_dir(), fname))
                if startswith(line, "#")
                    for tok in split(line)
                        occursin('=', tok) || continue
                        k, v = split(tok, '=', limit = 2)
                        meta[k] = v
                    end
                    continue
                end
                isempty(strip(line)) && continue
                p = split(line)
                length(p) >= 2 || continue
                push!(xs, parse(Float64, p[1]))
                push!(ys, parse(Float64, p[2]))
                push!(νs,  length(p) >= 3 ? parse(Float64, p[3]) : NaN)
                push!(bys, length(p) >= 4 ? parse(Float64, p[4]) : NaN)
            end
        catch err
            @warn "brioWu1d: skipping an unreadable curve" file=fname exception=err
            continue
        end
        isempty(xs) && continue
        nop = tryparse(Int, get(meta, "nop", ""))
        nop === nothing && continue
        tc = tryparse(Float64, get(meta, "t", ""))
        # Only curves from the same final time belong on the same figure.
        (tc === nothing || abs(tc - t) > 1.0e-8*max(1.0, abs(t))) && continue
        push!(curves, (nop  = nop,
                       ndofs = length(xs),
                       Cmin = something(tryparse(Float64, get(meta, "Cmin", "")), 0.0),
                       form = get(meta, "form", "elem"),
                       mtime = mtime(joinpath(_bw_curve_dir(), fname)),
                       iout = something(tryparse(Int, get(meta, "iout", "")), 0),
                       x = xs, y = ys,
                       nu = all(isnan, νs) ? Float64[] : νs,    # empty: stored before ν was
                       By = all(isnan, bys) ? Float64[] : bys))  # empty: stored before B_y was
    end
    # One curve per (order, DOFs, form): the store has carried two file-naming
    # conventions, so keep the most recently written of any duplicates.
    sort!(curves, by = c -> (c.nop, c.ndofs, c.form, -c.mtime))
    out = eltype(curves)[]
    for c in curves
        isempty(out) && (push!(out, c); continue)
        l = out[end]
        (l.nop == c.nop && l.ndofs == c.ndofs && l.form == c.form) || push!(out, c)
    end
    return out
end

# The curves shown on the density figure: ONE per order, all at the same
# resolution — the comparison is only meaningful at equal degrees of freedom
# (Dao & Nazarov's Fig. 2 is captioned "under the same number of degrees of
# freedom"). Take the finest DOF count that every order in the store has;
# counts are bucketed to the nearest 50 because nelx = DOFs/nop is rounded
# (601, 601, 601, 603 points for orders 4, 5, 6, 7 at a target of 600). If the
# orders share no resolution, fall back to the finest of each.
_bw_bucket(ndofs) = 50*round(Int, ndofs/50)

function _bw_finest(curves)
    isempty(curves) && return curves
    nops    = sort(unique(c.nop for c in curves))
    buckets = [Set(_bw_bucket(c.ndofs) for c in curves if c.nop == n) for n in nops]
    common  = reduce(intersect, buckets)
    if !isempty(common)
        b = maximum(common)
        out = [first(sort(filter(c -> c.nop == n && _bw_bucket(c.ndofs) == b, curves),
                          by = c -> -c.ndofs)) for n in nops]
        return out
    end
    best = Dict{Int,Any}()
    for c in curves
        (!haskey(best, c.nop) || c.ndofs > best[c.nop].ndofs) && (best[c.nop] = c)
    end
    return sort!(collect(values(best)), by = c -> c.nop)
end

# Marker positions for curve `k` of `n`, staggered so that the orders do not
# stack their markers on the same nodes where the solutions coincide (they
# agree over most of the tube, which is the point of the comparison).
function _bw_marker_idx(npts, k, n, inputs)
    base = _sparse_marker_idx(npts, inputs)
    isempty(base) && return base
    stride = max(1, step(base))
    off    = (n <= 1) ? 0 : ((k - 1)*stride) ÷ n
    return (1 + off):stride:npts
end

#---------------------------------------------------------------------------------
# Errors against the reference solution.
#
# NOTE what the reference is (user_analytic.jl): a FIRST-ORDER finite-volume
# (HLL) solution on 10 000 cells, sampled at 2000 points. Its own error is
# O(Δx) — about 1e-4 in the smooth fan, but a few cells of smearing at every
# discontinuity, where it is far larger than the difference between two
# spectral-element orders. A norm over the whole tube is therefore dominated
# by the jumps and says little about the order; the same norm over
# BW_SMOOTH_WINDOW is the one that can.
#---------------------------------------------------------------------------------
_bw_trapz(x, f) = sum(0.5*(f[i] + f[i+1])*(x[i+1] - x[i]) for i = 1:length(x)-1; init = 0.0)

function _bw_ref_at(xq)
    tab = _bw_read_reference()
    xr  = view(tab, :, 1)
    yr  = view(tab, :, 2)          # ρ
    out = similar(xq)
    for (i, xi) in enumerate(xq)
        xc = clamp(xi, xr[1], xr[end])
        j  = clamp(searchsortedlast(xr, xc), 1, length(xr) - 1)
        θ  = (xc - xr[j])/(xr[j+1] - xr[j])
        out[i] = (1 - θ)*yr[j] + θ*yr[j+1]
    end
    return out
end

# Relative L¹, L² (integral norms) and L∞ of ρ_h − ρ_ref over the window `w`.
function _bw_norms(xs, e, ρr, w)
    sel = (xs .>= w[1]) .& (xs .<= w[2])
    count(sel) >= 2 || return (NaN, NaN, NaN)
    x  = xs[sel]
    l1 = _bw_trapz(x, abs.(e[sel]))    / max(_bw_trapz(x, abs.(ρr[sel])), eps())
    l2 = sqrt(_bw_trapz(x, e[sel].^2)) / max(sqrt(_bw_trapz(x, ρr[sel].^2)), eps())
    li = maximum(abs, e[sel])          / max(maximum(abs, ρr[sel]), eps())
    return (l1, l2, li)
end

function _bw_error_table(curves)
    rows = NamedTuple[]
    for c in curves
        ρr = _bw_ref_at(c.x)
        e  = c.y .- ρr
        a1, a2, ai = _bw_norms(c.x, e, ρr, (0.0, 1.0))
        s1, s2, si = _bw_norms(c.x, e, ρr, BW_SMOOTH_WINDOW)
        push!(rows, (nop = c.nop, ndofs = c.ndofs, Cmin = c.Cmin,
                     l1_all = a1, l2_all = a2, linf_all = ai,
                     l1_sm  = s1, l2_sm  = s2, linf_sm  = si))
    end
    return sort!(rows, by = r -> (r.nop, r.ndofs))
end

function _bw_report_errors(rows)
    println(" # brioWu1d: density error against the reference (relative norms)")
    println(" #   nop  DOFs     L1 (0,1)     L2 (0,1)   Linf (0,1)     L1 (fan)     L2 (fan)   Linf (fan)")
    for r in rows
        println(@sprintf(" #   %3d %5d   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e",
                         r.nop, r.ndofs, r.l1_all, r.l2_all, r.linf_all,
                         r.l1_sm, r.l2_sm, r.linf_sm))
    end
    try
        mkpath(_bw_curve_dir())
        open(joinpath(_bw_curve_dir(), "errors.dat"), "w") do io
            println(io, "# density error against reference_hll.dat, relative norms")
            println(io, "# smooth window = ", BW_SMOOTH_WINDOW)
            println(io, "# nop ndofs Cmin L1_all L2_all Linf_all L1_sm L2_sm Linf_sm")
            for r in rows
                println(io, r.nop, " ", r.ndofs, " ", r.Cmin, " ",
                        r.l1_all, " ", r.l2_all, " ", r.linf_all, " ",
                        r.l1_sm, " ", r.l2_sm, " ", r.linf_sm)
            end
        end
    catch err
        @warn "brioWu1d: could not write the error table" exception=err
    end
    return nothing
end

# Observed rate from the two finest resolutions of one order: error ~ h^p with
# h ∝ 1/DOFs.
function _bw_rate(xs, ys)
    length(xs) >= 2 || return NaN
    (ys[end-1] > 0 && ys[end] > 0) || return NaN
    return log(ys[end-1]/ys[end])/log(xs[end-1]/xs[end])
end

#---------------------------------------------------------------------------------
# Convergence history in the layout of Dao & Nazarov (2022), Fig. 1: the error
# against 1/#DOFs on log-log axes, one line per polynomial order, with dashed
# slope guides. (Their abscissa is 1/sqrt(#DOFs) because the vortex problem is
# 2D; in 1D the mesh size is h ∝ 1/#DOFs, so that is the abscissa here, and it
# plays the same role.) The legend carries the rate measured between the two
# finest resolutions of each order.
#---------------------------------------------------------------------------------
function _bw_plot_convergence(rows, inputs, OUTPUT_DIR, iout; smooth::Bool)
    fields = smooth ? (:l1_sm, :l2_sm, :linf_sm) : (:l1_all, :l2_all, :linf_all)
    names  = ("L^1", "L^2", "L^\\infty")
    nops   = sort(unique(r.nop for r in rows))
    # Nothing to show unless at least one order has two resolutions.
    any(nop -> count(r -> r.nop == nop, rows) >= 2, nops) || return nothing

    panels = Plots.Plot[]
    for (fld, nm) in zip(fields, names)
        allx = Float64[]; ally = Float64[]
        pl = Plots.plot(; xscale = :log10, yscale = :log10,
                        xlabel = LaTeXStrings.L"1/\#\mathrm{DOFs}",
                        ylabel = LaTeXStrings.latexstring(string(
                            "\\Vert\\rho_h-\\rho_{ref}\\Vert_{", nm, "}\\ /\\ \\Vert\\rho_{ref}\\Vert_{", nm, "}")),
                        framestyle = :box, grid = true,
                        legend = :bottomright, legendfontsize = 8,
                        titlefontsize = 13, guidefontsize = 11, tickfontsize = 10,
                        title = LaTeXStrings.latexstring(string(nm, "\\mathrm{-error},\\ ",
                                 smooth ? "\\mathrm{fast\\ rarefaction}" : "\\mathrm{whole\\ tube}")),
                        show = false)
        for nop in nops
            sub = sort(filter(r -> r.nop == nop, rows), by = r -> r.ndofs)
            xs  = [1.0/r.ndofs for r in sub]
            ys  = [getfield(r, fld) for r in sub]
            keep = isfinite.(ys) .& (ys .> 0)
            any(keep) || continue
            xs = xs[keep]; ys = ys[keep]
            append!(allx, xs); append!(ally, ys)
            col, ls, mk = _bw_style(nop)
            p   = _bw_rate(xs, ys)
            lab = isfinite(p) ?
                  LaTeXStrings.latexstring(string("\\mathrm{nop}\\ ", nop, "\\ (p=", round(p; digits = 2), ")")) :
                  LaTeXStrings.latexstring(string("\\mathrm{nop}\\ ", nop))
            Plots.plot!(pl, xs, ys; line = (col, 1.8, :solid), marker = (mk, 5),
                        markerstrokewidth = 0.8, color = col, label = lab)
        end
        # Slope guides bracketing the data, as in the paper.
        if !isempty(allx)
            x2   = maximum(allx)
            ymax = maximum(ally); ymin = minimum(ally)
            for (sl, col, anchor) in ((1, :gray40, 1.8*ymax), (2, :gray70, 0.55*ymin))
                xg = [minimum(allx), x2]
                yg = [anchor*(xi/x2)^sl for xi in xg]
                Plots.plot!(pl, xg, yg; line = (col, 1.4, :dash),
                            label = LaTeXStrings.latexstring(string("\\mathrm{slope}\\ ", sl)))
            end
        end
        push!(panels, pl)
    end
    plt = Plots.plot(panels...; layout = (1, 3), size = (1500, 450),
                     left_margin = 9Plots.mm, bottom_margin = 8Plots.mm, show = false)
    _bw_savefig(plt, string(OUTPUT_DIR, "/convergence", smooth ? "_smooth" : "", "-it", iout, ".", _bw_fig_ext(inputs)))
    return nothing
end

# Legend entry: "nop 4, 601 DOFs" — with C_min appended only when the curves
# on the figure disagree about it (when they agree it goes in the title, once).
function _bw_label(c, lshow_cmin)
    s = string("\\mathrm{nop}\\ ", c.nop, ",\\ ", c.ndofs, "\\ \\mathrm{DOFs}")
    lshow_cmin && (s = string(s, ",\\ C_{min} = ", c.Cmin))
    return LaTeXStrings.latexstring(s)
end
# the reference's entry, in the same (LaTeX) font as the curves'
const BW_REF_LABEL = LaTeXStrings.L"\mathrm{Reference\ solution}"

#---------------------------------------------------------------------------------
# The density figure itself: the curves handed to it (one per order, ALL at the
# same resolution — the comparison is only meaningful at equal degrees of
# freedom), the reference, and the four zoom boxes. `tag` goes into the file
# name, so one call per resolution leaves one file per resolution behind.
#---------------------------------------------------------------------------------
function _bw_plot_density(curves, xs, href, inputs, OUTPUT_DIR, iout, tag)
    isempty(curves) && return nothing
    # C_min in the title when every curve used the same one, in the legend
    # entries when they differ.
    lcommon    = length(unique(c -> round(c.Cmin; digits = 12), curves)) == 1
    lshow_cmin = !lcommon
    fname = string("\\ (", startswith(curves[1].form, "nodal") ? "\\mathrm{nodal}" : "\\mathrm{element}",
                   "\\ \\nu,\\ \\Delta_K/", endswith(curves[1].form, "H") ? "k" : "(k{+}1)", ")")
    ttl = lcommon ? LaTeXStrings.latexstring(string("\\mathrm{Density},\\ C_{min} = ", curves[1].Cmin, fname)) :
                    LaTeXStrings.latexstring(string("\\mathrm{Density}", fname))

    xl = (0.0, 1.0); yl = (0.1, 1.0)
    plt = Plots.plot(; title = ttl,
                     xlabel = LaTeXStrings.L"x", ylabel = LaTeXStrings.L"\rho",
                     xlims = xl, ylims = yl,
                     xticks = 0.0:0.2:1.0, yticks = 0.1:0.1:1.0,
                     framestyle = :box, grid = false,
                     legend = :topright, legendfontsize = 11,
                     titlefontsize = 20, guidefontsize = 18, tickfontsize = 13,
                     size = (900, 720), left_margin = 6Plots.mm, bottom_margin = 5Plots.mm,
                     show = false)

    for (k, c) in enumerate(curves)
        col, ls, mk = _bw_style(c.nop)
        Plots.plot!(plt, c.x, c.y; line = (col, 1.6, ls), label = _bw_label(c, lshow_cmin))
        midx = _bw_marker_idx(length(c.x), k, length(curves), inputs)
        isempty(midx) || Plots.scatter!(plt, c.x[midx], c.y[midx];
                                        marker = (mk, 3.5), markerstrokewidth = 0,
                                        color = col, label = "")
    end

    if href !== nothing
        # The only solid line on the figure.
        Plots.plot!(plt, xs, href; line = (:black, 1.8, :solid), label = BW_REF_LABEL)

        # Zoom boxes: the paper's four (_bw_add_insets! draws the frame, the
        # connector and the inset; this fills it with every stored curve).
        _bw_add_insets!(plt, BW_INSETS, xl, yl, (sp, xr, yr) -> begin
            sel = (xs .>= xr[1]) .& (xs .<= xr[2])
            Plots.plot!(sp, xs[sel], href[sel]; line = (:black, 1.8, :solid), label = "")
            for c in curves
                col, ls, _ = _bw_style(c.nop)
                q = (c.x .>= xr[1]) .& (c.x .<= xr[2])
                any(q) && Plots.plot!(sp, c.x[q], c.y[q]; line = (col, 1.6, ls), label = "")
            end
        end)
    end
    _bw_savefig(plt, string(OUTPUT_DIR, "/density", tag, "-it", iout, ".", _bw_fig_ext(inputs)))
    return nothing
end

#---------------------------------------------------------------------------------
# DENSITY AND COEFFICIENT, side by side: the layout of Dao & Nazarov's order
# comparison (2024, Fig. 2, "first-order solutions to the Brio-Wu problem and
# the corresponding amounts of artificial viscosity"). Left, the density of
# every curve handed to it against the reference, with the zoom boxes of the
# density figure; right, the DynSGS coefficient ν of the same runs at the same
# time. Written as density_nu<tag>-it<n>, beside the density figure.
#
# ν is the coefficient of the ρu slot, which in this deck's conserved form is
# THE coefficient (one kinematic ν on every equation). The element form makes
# it constant per element, so through the nodes it is a staircase. Its axis is
# linear, as in the paper, unless it spans more than two decades (C_min = 0,
# where it falls to the residual's floor between the waves): then log10.
#
# `style(c)` gives a curve its (colour, line style, marker), by default the
# order's own; `markers = false` drops the markers and `legend_nu_only` keeps
# the curves' legend in the ν panel, where a crowded figure has room for it.
#---------------------------------------------------------------------------------
function _bw_plot_density_nu(curves, xs, href, inputs, OUTPUT_DIR, iout, tag;
                             style = c -> _bw_style(c.nop), markers = true, legend_nu_only = false)
    withnu = filter(c -> !isempty(c.nu), curves)
    if isempty(withnu)
        println(" #   brioWu1d: no coefficient stored with these curves (written before ν was ",
                "stored: rerun them); no density_nu", tag, " figure.")
        return nothing
    end
    lcommon    = length(unique(c -> round(c.Cmin; digits = 12), withnu)) == 1
    lshow_cmin = !lcommon
    form  = withnu[1].form
    fdesc = string(startswith(form, "nodal") ? "\\mathrm{nodal}" : "\\mathrm{element}",
                   "\\ \\nu,\\ \\Delta_K/", endswith(form, "H") ? "k" : "(k{+}1)")
    ttl_ρ = LaTeXStrings.latexstring(string("\\mathrm{Density}",
                                            lcommon ? string(",\\ C_{min} = ", withnu[1].Cmin) : "",
                                            "\\ (", fdesc, ")"))
    νs   = reduce(vcat, [filter(isfinite, c.nu) for c in withnu])
    νmax = maximum(νs; init = 0.0)
    νpos = minimum(filter(>(0), νs); init = νmax)
    llog = νmax > 0 && νmax/max(νpos, floatmin(Float64)) > 100.0

    xl = (0.0, 1.0); yl = (0.1, 1.0)
    font = (legendfontsize = 11, titlefontsize = 16, guidefontsize = 18, tickfontsize = 13)
    plt = Plots.plot(; layout = (1, 2), size = (1800, 720), show = false,
                     left_margin = 9Plots.mm, bottom_margin = 10Plots.mm,
                     top_margin = 2Plots.mm, right_margin = 4Plots.mm)
    Plots.plot!(plt[1]; title = ttl_ρ, xlabel = LaTeXStrings.L"x", ylabel = LaTeXStrings.L"\rho_h",
                xlims = xl, ylims = yl, xticks = 0.0:0.2:1.0, yticks = 0.1:0.1:1.0,
                framestyle = :box, grid = false, legend = :topright, font...)
    Plots.plot!(plt[2]; title = LaTeXStrings.L"\mathrm{DynSGS\ coefficient}\ \nu",
                xlabel = LaTeXStrings.L"x", ylabel = LaTeXStrings.L"\nu",
                xlims = xl, xticks = 0.0:0.2:1.0, framestyle = :box, grid = false,
                # the left state carries the smallest ν: the top-left corner is free
                legend = legend_nu_only ? :topleft : :topright, yscale = llog ? :log10 : :identity,
                ylims = llog ? :auto : (0.0, 1.08*νmax), font...)

    for (k, c) in enumerate(withnu)
        col, ls, mk = style(c)
        lab  = _bw_label(c, lshow_cmin)
        midx = markers ? _bw_marker_idx(length(c.x), k, length(withnu), inputs) : (1:0)
        ν    = llog ? [v > 0 ? v : NaN for v in c.nu] : c.nu
        for (sp, y) in ((1, c.y), (2, ν))
            Plots.plot!(plt, c.x, y; subplot = sp, line = (col, 1.6, ls),
                        label = (legend_nu_only && sp == 1) ? "" : lab)
            isempty(midx) || Plots.scatter!(plt, c.x[midx], y[midx]; subplot = sp,
                                            marker = (mk, 3.5), markerstrokewidth = 0,
                                            color = col, label = "")
        end
    end
    if href !== nothing
        Plots.plot!(plt, xs, href; subplot = 1, line = (:black, 1.8, :solid), label = BW_REF_LABEL)
        # the insets of the left panel become subplots 3, 4, ...: 1 and 2 are the panels
        _bw_add_insets!(plt, BW_INSETS, xl, yl, (sp, xr, yr) -> begin
            sel = (xs .>= xr[1]) .& (xs .<= xr[2])
            Plots.plot!(sp, xs[sel], href[sel]; line = (:black, 1.8, :solid), label = "")
            for c in withnu
                col, ls, _ = style(c)
                q = (c.x .>= xr[1]) .& (c.x .<= xr[2])
                any(q) && Plots.plot!(sp, c.x[q], c.y[q]; line = (col, 1.6, ls), label = "")
            end
        end; parent = 1, first = 3)
    end
    _bw_savefig(plt, string(OUTPUT_DIR, "/density_nu", tag, "-it", iout, ".", _bw_fig_ext(inputs)))
    return nothing
end

#---------------------------------------------------------------------------------
# ρ AND B_y SIDE BY SIDE, the layout of the Brio-Wu figure of Dao, Nazarov &
# Tomas (arXiv:2310.18467, Figs. 1-2), THE DEFAULT comparison figure: every
# curve handed to it against the reference, ρ with the four zoom boxes of the
# density figure, B_y with the three of theirs and a fourth at the foot of the
# slow shock, which theirs leave out (BW_PAIR_BY_INSETS). B_y falls from 1 on
# the left to -1 on the right, so its panel is empty to the left, where two of
# its boxes go, and to the right, where the other two go, one above the other;
# the legend (every curve and the reference) goes to the top right of the ρ
# panel, empty above its boxes, as in their figure. Written as
# rho_By<tag>-it<n>. `style` and `markers` as for _bw_plot_density_nu.
#---------------------------------------------------------------------------------
const BW_PAIR_BY_INSETS = [   # (x-range, y-range, inset position: x, y, w, h of the axes, connector[, x ticks])
    # Dao, Nazarov & Tomas's three, read off their Fig. 2: the same windows, the
    # first two insets where theirs are (moved in a little, for their tick
    # labels to clear the frame, and the first one smaller, to clear the fast
    # rarefaction), and their connectors, from the bottom-left corner of the
    # box to the bottom-left corner of the inset
    ((0.416, 0.476), ( 0.355,  0.615), (0.085, 0.58,  0.27, 0.29), (:bottomleft, :bottomleft)), # top of the compound wave
    ((0.560, 0.650), (-0.630, -0.510), (0.115, 0.065, 0.29, 0.30), (:bottomleft, :bottomleft)), # the plateau up to the slow shock
    # the foot of the slow shock, where an undershoot below the plateau sits:
    # the inset at the top right, its left edge on the right edge of the box,
    # so that the connector runs straight up through empty space, past the
    # rarefaction inset, and no tick at the window's ends, where a label would
    # sit on the connector
    ((0.632, 0.698), (-0.950, -0.850), (0.698, 0.66,  0.24, 0.28), (:topright,   :bottomleft), [0.64, 0.66, 0.68]), # the foot of the slow shock
    # the right fast rarefaction, below the slow-shock inset rather than at the
    # height of theirs, a little smaller, and with their ticks
    ((0.763, 0.913), (-1.010, -0.880), (0.72,  0.30,  0.26, 0.28), (:bottomleft, :bottomleft), [0.80, 0.85, 0.90]), # the right fast rarefaction
]

# "element ν, Δ_K/(k+1)": the method of a store form, for the titles
_bw_form_desc(form) = string(startswith(form, "nodal") ? "\\mathrm{nodal}" : "\\mathrm{element}",
                             "\\ \\nu,\\ \\Delta_K/", endswith(form, "H") ? "k" : "(k{+}1)")

# the reference and every curve's `field` inside one zoom box
function _bw_fill_inset!(sp, xr, xs, yref, curves, field, style)
    sel = (xs .>= xr[1]) .& (xs .<= xr[2])
    any(sel) && Plots.plot!(sp, xs[sel], yref[sel]; line = (:black, 1.8, :solid), label = "")
    for c in curves
        col, ls, _ = style(c)
        q = (c.x .>= xr[1]) .& (c.x .<= xr[2])
        any(q) && Plots.plot!(sp, c.x[q], field(c)[q]; line = (col, 1.6, ls), label = "")
    end
    return nothing
end

function _bw_plot_rho_By(curves, xs, href, hBy, inputs, OUTPUT_DIR, iout, tag;
                         style = c -> _bw_style(c.nop), markers = true, legendfontsize = 11)
    withBy = filter(c -> !isempty(c.By), curves)
    if isempty(withBy)
        println(" #   brioWu1d: no B_y stored with these curves (written before B_y was ",
                "stored: rerun them); no rho_By", tag, " figure.")
        return nothing
    end
    lcommon    = length(unique(c -> round(c.Cmin; digits = 12), withBy)) == 1
    lshow_cmin = !lcommon
    ttl_ρ = LaTeXStrings.latexstring(string("\\rho",
                                            lcommon ? string(",\\ C_{min} = ", withBy[1].Cmin) : "",
                                            "\\ (", _bw_form_desc(withBy[1].form), ")"))
    # B_y's frame: the data and the reference, padded by 2% as in Dao, Nazarov &
    # Tomas's figure, which the positions of the insets are taken from
    lo = minimum(minimum(c.By) for c in withBy); hi = maximum(maximum(c.By) for c in withBy)
    hBy === nothing || (lo = min(lo, minimum(hBy)); hi = max(hi, maximum(hBy)))
    pad = 0.02*(hi - lo)
    xl = (0.0, 1.0); yl = (0.1, 1.0); ylB = (lo - pad, hi + pad)
    font = (titlefontsize = 16, guidefontsize = 18, tickfontsize = 13)
    plt = Plots.plot(; layout = (1, 2), size = (1800, 720), show = false,
                     left_margin = 9Plots.mm, bottom_margin = 10Plots.mm,
                     top_margin = 2Plots.mm, right_margin = 4Plots.mm)
    Plots.plot!(plt[1]; title = ttl_ρ, xlabel = LaTeXStrings.L"x", ylabel = LaTeXStrings.L"\rho_h",
                xlims = xl, ylims = yl, xticks = 0.0:0.2:1.0, yticks = 0.1:0.1:1.0,
                framestyle = :box, grid = false, legend = :topright,
                legendfontsize = legendfontsize, font...)
    Plots.plot!(plt[2]; title = LaTeXStrings.L"B_y", xlabel = LaTeXStrings.L"x",
                ylabel = LaTeXStrings.L"B_{y,h}", xlims = xl, ylims = ylB, xticks = 0.0:0.2:1.0,
                framestyle = :box, grid = false, legend = false, font...)
    for (k, c) in enumerate(withBy)
        col, ls, mk = style(c)
        midx = markers ? _bw_marker_idx(length(c.x), k, length(withBy), inputs) : (1:0)
        for (sp, y) in ((1, c.y), (2, c.By))
            Plots.plot!(plt, c.x, y; subplot = sp, line = (col, 1.6, ls),
                        label = sp == 1 ? _bw_label(c, lshow_cmin) : "")
            isempty(midx) || Plots.scatter!(plt, c.x[midx], y[midx]; subplot = sp,
                                            marker = (mk, 3.5), markerstrokewidth = 0,
                                            color = col, label = "")
        end
    end
    nsub = 3     # the panels are subplots 1 and 2; the zoom boxes follow
    if href !== nothing
        Plots.plot!(plt, xs, href; subplot = 1, line = (:black, 1.8, :solid), label = BW_REF_LABEL)
        _bw_add_insets!(plt, BW_INSETS, xl, yl,
                        (sp, xr, yr) -> _bw_fill_inset!(sp, xr, xs, href, withBy, c -> c.y, style);
                        parent = 1, first = nsub)
        nsub += length(BW_INSETS)
    end
    if hBy !== nothing
        Plots.plot!(plt, xs, hBy; subplot = 2, line = (:black, 1.8, :solid), label = "")
        _bw_add_insets!(plt, BW_PAIR_BY_INSETS, xl, ylB,
                        (sp, xr, yr) -> _bw_fill_inset!(sp, xr, xs, hBy, withBy, c -> c.By, style);
                        parent = 2, first = nsub)
    end
    _bw_savefig(plt, string(OUTPUT_DIR, "/rho_By", tag, "-it", iout, ".", _bw_fig_ext(inputs)))
    return nothing
end

#---------------------------------------------------------------------------------
# THE DYNSGS COEFFICIENT ALONE, for the same curves: ν of the ρu slot (see
# _bw_plot_density_nu) against x, on a linear axis unless it spans more than
# two decades (C_min = 0: then log10). The legend sits top left, over the left
# state, where ν is smallest. Written as nu<tag>-it<n>.
#---------------------------------------------------------------------------------
function _bw_plot_nu(curves, inputs, OUTPUT_DIR, iout, tag;
                     style = c -> _bw_style(c.nop), markers = true, legendfontsize = 11)
    withnu = filter(c -> !isempty(c.nu), curves)
    if isempty(withnu)
        println(" #   brioWu1d: no coefficient stored with these curves (written before ν was ",
                "stored: rerun them); no nu", tag, " figure.")
        return nothing
    end
    lcommon    = length(unique(c -> round(c.Cmin; digits = 12), withnu)) == 1
    lshow_cmin = !lcommon
    ttl  = LaTeXStrings.latexstring(string("\\mathrm{DynSGS\\ coefficient}\\ \\nu",
                                           lcommon ? string(",\\ C_{min} = ", withnu[1].Cmin) : "",
                                           "\\ (", _bw_form_desc(withnu[1].form), ")"))
    νs   = reduce(vcat, [filter(isfinite, c.nu) for c in withnu])
    νmax = maximum(νs; init = 0.0)
    νpos = minimum(filter(>(0), νs); init = νmax)
    llog = νmax > 0 && νmax/max(νpos, floatmin(Float64)) > 100.0
    plt = Plots.plot(; title = ttl, xlabel = LaTeXStrings.L"x", ylabel = LaTeXStrings.L"\nu",
                     xlims = (0.0, 1.0), xticks = 0.0:0.2:1.0, framestyle = :box, grid = false,
                     legend = :topleft, legendfontsize = legendfontsize,
                     yscale = llog ? :log10 : :identity, ylims = llog ? :auto : (0.0, 1.08*νmax),
                     titlefontsize = 16, guidefontsize = 18, tickfontsize = 13,
                     size = (1000, 650), left_margin = 8Plots.mm, bottom_margin = 6Plots.mm,
                     right_margin = 5Plots.mm, show = false)   # room for the last x tick
    for (k, c) in enumerate(withnu)
        col, ls, mk = style(c)
        ν    = llog ? [v > 0 ? v : NaN for v in c.nu] : c.nu
        midx = markers ? _bw_marker_idx(length(c.x), k, length(withnu), inputs) : (1:0)
        Plots.plot!(plt, c.x, ν; line = (col, 1.6, ls), label = _bw_label(c, lshow_cmin))
        isempty(midx) || Plots.scatter!(plt, c.x[midx], ν[midx]; marker = (mk, 3.5),
                                        markerstrokewidth = 0, color = col, label = "")
    end
    _bw_savefig(plt, string(OUTPUT_DIR, "/nu", tag, "-it", iout, ".", _bw_fig_ext(inputs)))
    return nothing
end

#---------------------------------------------------------------------------------
# EVERY RESOLUTION ON ONE FIGURE, as Dao & Nazarov's Fig. 2 puts its two: the
# colour is the order, the line type the resolution (dashed the coarsest, then
# dotted, dash-dot, dash-dot-dot), no markers. Written as nu_all-it<n> (and
# density_nu_all-it<n> with JEXPRESSO_BW_RIGHT=nu) when the store holds two
# resolutions or more; with four, that is sixteen curves. Not ρ and B_y: with
# every resolution in one B_y panel there are too many curves to tell which
# one overshoots, so those stay one figure per resolution (rho_By_dof<M>).
#---------------------------------------------------------------------------------
const BW_RES_LINES = [:dash, :dot, :dashdot, :dashdotdot]

function _bw_plot_all(stored, xs, href, inputs, OUTPUT_DIR, iout, right)
    buckets = sort(unique(_bw_bucket(c.ndofs) for c in stored))
    length(buckets) >= 2 || return nothing
    curves = reduce(vcat, [_bw_finest(filter(c -> _bw_bucket(c.ndofs) == b, stored)) for b in buckets])
    rline(c) = BW_RES_LINES[mod1(findfirst(==(_bw_bucket(c.ndofs)), buckets), length(BW_RES_LINES))]
    style = c -> (_bw_style(c.nop)[1], rline(c), :none)
    right == "nu" && _bw_plot_density_nu(curves, xs, href, inputs, OUTPUT_DIR, iout, "_all";
                                         style = style, markers = false, legend_nu_only = true)
    _bw_plot_nu(curves, inputs, OUTPUT_DIR, iout, "_all"; style = style, markers = false, legendfontsize = 9)
    return nothing
end

#---------------------------------------------------------------------------------
# Every comparison figure of the store at the final time: per resolution the
# density figure, the ρ-and-B_y figure (with JEXPRESSO_BW_RIGHT=nu the
# density-and-coefficient one instead) and the coefficient alone, the finest of
# each under the plain name, the coefficient of every resolution together, the
# error table and the two convergence histories. Called by the run
# (user_plot_1d) and by tools/brio_wu_plot.jl, which redraws them from the
# store alone at the end of a sweep. `hBy` is the B_y reference on the
# abscissae `xs` of `href`.
#---------------------------------------------------------------------------------
function _bw_draw_store(stored, xs, href, inputs, OUTPUT_DIR, iout;
                        hBy = nothing, lerrors = href !== nothing)
    isempty(stored) && return nothing
    right = lowercase(strip(get(ENV, "JEXPRESSO_BW_RIGHT", "by")))   # the right panel: "by" or "nu"
    # Curves stored before B_y was have none to draw: when no curve has it, the
    # ρ-and-ν figures stand in for the ρ-and-B_y ones.
    noBy = count(c -> isempty(c.By), stored)
    if right != "nu" && noBy > 0
        if noBy == length(stored)
            println(" #   brioWu1d: no B_y stored with these curves (written before it was): ",
                    "density_nu* (ρ and ν) instead of rho_By*; rerun the sweep for ρ and B_y.")
            right = "nu"
        else
            println(" #   brioWu1d: ", noBy, " curve(s) stored without B_y (written before it was) ",
                    "are left out of rho_By*; rerun them.")
        end
    end
    figs = Tuple{Any,String}[]
    # ONE FIGURE PER RESOLUTION. Orders may only be compared at equal
    # degrees of freedom, so each DOF bucket gets its own file and a sweep
    # leaves density_dof150, density_dof300, ... behind, each holding every
    # order that ran at that resolution.
    for b in sort(unique(_bw_bucket(c.ndofs) for c in stored))
        push!(figs, (_bw_finest(filter(c -> _bw_bucket(c.ndofs) == b, stored)), string("_dof", b)))
    end
    # The finest resolution keeps the plain name as the headline figure.
    push!(figs, (_bw_finest(stored), ""))
    for (sub, tag) in figs
        _bw_plot_density(sub, xs, href, inputs, OUTPUT_DIR, iout, tag)
        try
            right == "nu" ? _bw_plot_density_nu(sub, xs, href, inputs, OUTPUT_DIR, iout, tag) :
                            _bw_plot_rho_By(sub, xs, href, hBy, inputs, OUTPUT_DIR, iout, tag)
            _bw_plot_nu(sub, inputs, OUTPUT_DIR, iout, tag)
        catch err
            @warn "brioWu1d: a comparison figure failed; the density figure is unaffected." exception=err
        end
    end
    # Every order at every resolution on one figure. ρ and B_y no longer get
    # one: an rho_By_all left by an earlier redraw would pass for current.
    isdir(OUTPUT_DIR) && foreach(f -> occursin(r"^rho_By_all-it\d+\.(svg|pdf|png)$", f) &&
                                      rm(joinpath(OUTPUT_DIR, f); force = true), readdir(OUTPUT_DIR))
    try
        _bw_plot_all(stored, xs, href, inputs, OUTPUT_DIR, iout, right)
    catch err
        @warn "brioWu1d: the all-resolutions figures failed; the other figures are unaffected." exception=err
    end
    # Error table and the convergence history, from every stored curve.
    lerrors || return nothing
    rows = _bw_error_table(stored)
    _bw_report_errors(rows)
    _bw_plot_convergence(rows, inputs, OUTPUT_DIR, iout; smooth = false)
    _bw_plot_convergence(rows, inputs, OUTPUT_DIR, iout; smooth = true)
    return nothing
end

# ---------------------------------------------------------------------------
# ZOOM BOXES. Shared by the density figure and the per-variable figures: the
# gray frame on the main axes, the dashed connector from the corner of the box
# nearest the inset to the corner of the inset nearest the box (or between the
# two corners a box names in a fourth entry, (box corner, inset corner), each
# :bottomleft, :bottomright, :topleft or :topright), and the inset subplot
# itself, on subplot `parent` (the insets become subplots `first`, `first`+1,
# ...), with the x ticks of a fifth entry when there is one (a narrow inset
# would run its automatic ones together). `draw!(sp, xr, yr)` fills the inset
# -- the two callers have different curves to put in it (a whole store of
# orders, or this run and its reference), which is the only part that differs.
# ---------------------------------------------------------------------------
_bw_corner(xr, yr, c) = (c in (:bottomleft, :topleft) ? xr[1] : xr[2],
                         c in (:bottomleft, :bottomright) ? yr[1] : yr[2])

function _bw_add_insets!(plt, boxes, xl, yl, draw!; parent = 1, first = 2)
    for (k, box) in enumerate(boxes)
        xr, yr, pos = box[1], box[2], box[3]
        bx = [xr[1], xr[2], xr[2], xr[1], xr[1]]
        by = [yr[1], yr[1], yr[2], yr[2], yr[1]]
        Plots.plot!(plt, bx, by; subplot = parent, line = (:gray, 1.0), label = "")
        ix  = (xl[1] + pos[1]*(xl[2] - xl[1]), xl[1] + (pos[1] + pos[3])*(xl[2] - xl[1]))
        iy  = (yl[1] + pos[2]*(yl[2] - yl[1]), yl[1] + (pos[2] + pos[4])*(yl[2] - yl[1]))
        if length(box) >= 4
            cornerx, cornery = _bw_corner(xr, yr, box[4][1])
            insx, insy       = _bw_corner(ix, iy, box[4][2])
        else
            bcx = 0.5*(xr[1] + xr[2]); bcy = 0.5*(yr[1] + yr[2])
            icx = 0.5*(ix[1] + ix[2]); icy = 0.5*(iy[1] + iy[2])
            cornerx = (icx < bcx) ? xr[1] : xr[2]
            cornery = (icy < bcy) ? yr[1] : yr[2]
            insx    = (icx < bcx) ? ix[2] : ix[1]
            insy    = (icy < bcy) ? iy[2] : iy[1]
        end
        Plots.plot!(plt, [cornerx, insx], [cornery, insy]; subplot = parent,
                    line = (:black, 1.0, :dash), label = "")
        Plots.plot!(plt; inset = (parent, Plots.bbox(pos[1], pos[2], pos[3], pos[4], :bottom, :left)),
                    subplot = first + k - 1)
        sp = plt[first + k - 1]
        Plots.plot!(sp; xlims = xr, ylims = yr, framestyle = :box, grid = true,
                    xticks = length(box) >= 5 ? box[5] : :auto,
                    tickfontsize = 8, legend = false, background_color_inside = :white)
        draw!(sp, xr, yr)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# ONE FIGURE PER SOLUTION QUANTITY: <name>-it<n>.<svg|pdf|png>
#
# rho-it5.svg, u-it5.svg, v-it5.svg, w-it5.svg, p-it5.svg, Bx-it5.svg,
# By-it5.svg, Bz-it5.svg -- each the full canvas rather than one cell of a
# panel matrix, so a figure can go into a paper or a slide on its own and the
# zoom boxes have room to be read.
#
# The reference of user_analytic.jl covers rho, u, v, w, p, By, Bz; Bx has
# none because it is a constant of the 1D system, and a quantity that is
# constant by construction is centred on its own value with the spread in the
# title, so the panel reads as the check it is rather than as an arbitrary
# auto-scaled window.
#
# Boxes come from BW_FIELD_INSETS (rho: the paper's four), and are drawn only
# where there is a reference to compare against -- at the earlier output times
# there is none, and a zoom on a single curve says nothing.
# ---------------------------------------------------------------------------
function _bw_plot_field(xs, y, yref, name, inputs, OUTPUT_DIR, iout, t, ext)
    lo, hi = extrema(y)
    if yref !== nothing
        rlo, rhi = extrema(yref)
        lo = min(lo, rlo); hi = max(hi, rhi)
    end
    mid  = 0.5*(lo + hi)
    ttl  = string(name)
    lconst = hi - lo <= 1.0e-8*max(1.0, abs(mid))
    if lconst
        d  = max(hi - lo, 1.0e-12*max(1.0, abs(mid)))
        yl = (mid - 5d, mid + 5d)
        ttl = string(name, @sprintf("  (const, spread %.1e)", hi - lo))
    else
        pad = 0.06*(hi - lo)
        yl  = (lo - pad, hi + pad)
    end
    # ρ keeps the fixed frame of the paper's density figure, because its four
    # inset positions were tuned against exactly these limits.
    string(name) == "ρ" && (yl = (0.1, 1.0))
    xl = (minimum(xs), maximum(xs))

    plt = Plots.plot(; title = LaTeXStrings.latexstring(string("\\mathrm{Brio-Wu},\\ t = ",
                                                            @sprintf("%.3f", t === nothing ? NaN : t))),
                     xlabel = LaTeXStrings.L"x", ylabel = ttl,
                     xlims = xl, ylims = yl, xticks = 0.0:0.2:1.0,
                     framestyle = :box, grid = false,
                     # OUTSIDE the frame: with two entries and up to four zoom
                     # insets, any in-frame corner is one the insets want, and
                     # a legend under an inset is a legend nobody can read.
                     legend = :outertop, legend_columns = yref === nothing ? 1 : 2,
                     legendfontsize = 12,
                     titlefontsize = 20, guidefontsize = 18, tickfontsize = 13,
                     size = (900, 720), left_margin = 6Plots.mm, bottom_margin = 5Plots.mm,
                     show = false)

    midx = _bw_marker_idx(length(xs), 1, 1, inputs)
    Plots.plot!(plt, xs, y; line = (:blue, 1.8, :dash), label = "Jexpresso")
    isempty(midx) || Plots.scatter!(plt, xs[midx], y[midx];
                                    marker = (:circle, 3.5, :blue),
                                    markerstrokewidth = 0, label = "")
    yref === nothing || Plots.plot!(plt, xs, yref; line = (:black, 1.8, :solid),
                                    label = "Reference solution")

    if yref !== nothing && !lconst
        _bw_add_insets!(plt, _bw_field_insets(name), xl, yl, (sp, xr, yr) -> begin
            sel = (xs .>= xr[1]) .& (xs .<= xr[2])
            if any(sel)
                Plots.plot!(sp, xs[sel], yref[sel]; line = (:black, 1.8, :solid), label = "")
                Plots.plot!(sp, xs[sel], y[sel];    line = (:blue,  1.8, :dash),  label = "")
                Plots.scatter!(sp, xs[sel], y[sel]; marker = (:circle, 3.0, :blue),
                               markerstrokewidth = 0, label = "")
            end
        end)
    end
    _bw_savefig(plt, string(OUTPUT_DIR, "/", _bw_slug(name), "-it", iout, ".", ext))
    return nothing
end

function _bw_plot_fields(xs, qs, qrefs, outvar, inputs, OUTPUT_DIR, iout, t)
    names = string.(outvar)
    ext   = _bw_fig_ext(inputs)
    for ivar = 1:length(names)
        yref = (qrefs !== nothing && any(isfinite, @view(qrefs[:, ivar]))) ?
               collect(@view qrefs[:, ivar]) : nothing
        _bw_plot_field(xs, collect(@view qs[:, ivar]), yref,
                       names[ivar], inputs, OUTPUT_DIR, iout, t, ext)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# THE STRUCTURE OF THE DynSGS COEFFICIENT: mu_dsgs-it<n>.png
#
# Two panels on one x axis. The top one carries rho and By, which is where
# this solution's five waves are; the bottom one carries nu, so the figure
# answers the question the coefficient exists to answer -- WHERE it fires and
# by how much -- rather than showing a curve with nothing to read it against.
#
# nu is drawn per equation slot. This deck runs the conserved form
# (:dsgs_conserved => true) with :mu = ones(8), so every slot carries the same
# kinematic nu and the eight curves coincide; slots whose curve is identical
# to one already drawn are therefore folded into a single entry rather than
# overplotted eight deep. A deck that scales the slots apart (:mu != ones, or
# Nazarov's kappa on the energy slot) shows them as separate curves with no
# change here.
#
# The y axis is log10 when nu spans more than two decades, which it does at
# the cap-to-floor contrast of a shock tube; nu = 0 is then dropped rather
# than clipped, and the caption says how many nodes that was. The element
# form (:ldsgs_nodal => false, the default) gives the staircase of a
# per-element coefficient broadcast to its nodes; the nodal form gives a
# continuous curve. Both are drawn as they are, markers at the nodes, so the
# figure never smooths over which one produced it.
# ---------------------------------------------------------------------------
const BW_MU_EQNAMES = ["ρ", "ρu", "ρv", "E", "ρw", "Bx", "By", "Bz", "ψ"]

function _bw_plot_mu(xs, qs, μs, outvar, inputs, OUTPUT_DIR, iout, t)
    names = string.(outvar)
    nμ    = size(μs, 2)

    # Fold slots whose coefficient is the same curve into one legend entry.
    groups = Tuple{Vector{Int}, Vector{Float64}}[]
    for ieq = 1:nμ
        col = collect(@view μs[:, ieq])
        all(iszero, col) && continue
        k = findfirst(g -> g[2] == col, groups)
        k === nothing ? push!(groups, ([ieq], col)) : push!(groups[k][1], ieq)
    end

    μmax = maximum(maximum(g[2]) for g in groups; init = 0.0)
    if isempty(groups) || μmax <= 0.0
        println(" #   brioWu1d: the DynSGS coefficient is identically zero at this output time; no mu_dsgs figure.")
        return nothing
    end
    μpos  = minimum(minimum(filter(>(0), g[2]); init = μmax) for g in groups)
    llog  = μmax/max(μpos, eps()) > 100.0
    nzero = count(iszero, groups[1][2])

    # Top: where the waves are.
    top = Plots.plot(; xlabel = "", ylabel = "", legend = :best,
                     titlefontsize = 13, guidefontsize = 10, tickfontsize = 9,
                     legendfontsize = 8, show = false,
                     title = @sprintf("Brio-Wu, t = %.4f", t === nothing ? NaN : t))
    for (nm, st) in (("ρ", (:black, 1.6, :solid)), ("By", (:steelblue, 1.6, :dash)))
        j = findfirst(==(nm), names)
        j === nothing || Plots.plot!(top, xs, @view(qs[:, j]); line = st, label = nm)
    end

    # Bottom: the coefficient itself.
    form = _bw_form(inputs)
    ylab = llog ? "log10 ν" : "ν"
    bot  = Plots.plot(; xlabel = "x", ylabel = ylab, legend = :best,
                      guidefontsize = 10, tickfontsize = 9, legendfontsize = 8,
                      titlefontsize = 12, show = false,
                      title = string("DynSGS coefficient (", form, " form), max ν = ",
                                     @sprintf("%.3e", μmax)))
    cols = [:red, :darkorange, :seagreen, :purple, :teal, :magenta, :brown, :navy, :gray]
    midx = _bw_marker_idx(length(xs), 1, 1, inputs)
    for (k, (eqs, col)) in enumerate(groups)
        lab = join(getindex.(Ref(BW_MU_EQNAMES), eqs), ", ")
        c   = cols[mod1(k, length(cols))]
        y   = llog ? [v > 0 ? log10(v) : NaN for v in col] : col
        Plots.plot!(bot, xs, y; line = (c, 1.6), label = lab)
        isempty(midx) || Plots.scatter!(bot, xs[midx], y[midx];
                                        marker = (:circle, 2.5, c),
                                        markerstrokewidth = 0, label = "")
    end
    llog && nzero > 0 &&
        println(" #   brioWu1d: mu_dsgs on a log axis; ν = 0 at ", nzero, " of ",
                length(xs), " nodes (not drawn).")

    plt = Plots.plot(top, bot; layout = (2, 1), size = (900, 700), link = :x,
                     left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, show = false)
    _bw_savefig(plt, string(OUTPUT_DIR, "/mu_dsgs-it", iout, ".", _bw_fig_ext(inputs)))
    return nothing
end

function user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)
    iρ = findfirst(==("ρ"), string.(outvar))
    iρ === nothing && error("user_plot_1d (brioWu1d): no ρ among the output variables")
    idx = sortperm(x)
    xs  = x[idx]
    ρs  = q[idx, iρ]
    href = (qref !== nothing && any(isfinite, @view(qref[:, iρ]))) ? qref[idx, iρ] : nothing

    # Every solution quantity, and the structure of the coefficient that
    # stabilized them, at EVERY output time. The density figure below is the
    # paper's comparison and stays the headline; these two are the state of
    # the run. They are independent of the curve store, so they are drawn
    # before it and a failure in either cannot cost the density figure.
    qs = q[idx, :]
    try
        _bw_plot_fields(xs, qs, qref === nothing ? nothing : qref[idx, :],
                        outvar, inputs, OUTPUT_DIR, iout, t)
    catch err
        @warn "brioWu1d: the fields figure failed; the density figure is unaffected." exception=err
    end
    if μ_nodes !== nothing && size(μ_nodes, 1) == length(x)
        try
            _bw_plot_mu(xs, qs, μ_nodes[idx, :], outvar, inputs, OUTPUT_DIR, iout, t)
        catch err
            @warn "brioWu1d: the DynSGS coefficient figure failed." exception=err
        end
    end

    nop   = Int(get(inputs, :nop, 0))
    Cmin  = Float64(get(inputs, :dsgs_Cmin, 0.0))
    # the coefficient stored with the curve: the ρu slot (see _bw_plot_density_nu)
    νs = (μ_nodes !== nothing && size(μ_nodes, 1) == length(x)) ?
         Float64.(μ_nodes[idx, min(2, size(μ_nodes, 2))]) : Float64[]
    # and B_y, the other quantity of Dao & Nazarov's figure, with its reference
    iBy = findfirst(==("By"), string.(outvar))
    Bys = iBy === nothing ? Float64[] : Float64.(q[idx, iBy])
    hBy = (iBy !== nothing && qref !== nothing && any(isfinite, @view(qref[:, iBy]))) ?
          qref[idx, iBy] : nothing

    # The reference is supplied at the final time only; that is the figure the
    # orders are compared on, so that is where the curve is stored and where
    # every stored curve is drawn.
    lfinal = href !== nothing
    stored = NamedTuple[]
    form   = _bw_form(inputs)
    if lfinal
        _bw_save_curve(xs, ρs, νs, Bys, inputs, t, iout)
        # Inside a sweep (tools/brio_wu_order_scan.sh sets JEXPRESSO_BW_DRAWSTORE=0)
        # the run only stores its curve: the sweep draws every comparison figure
        # once, from the complete store, when its last run is done.
        get(ENV, "JEXPRESSO_BW_DRAWSTORE", "1") == "0" && return nothing
        # Only this run's form: the element and the nodal coefficient, and the
        # two length-scale conventions, are different methods and do not belong
        # on one comparison. Say so when the store holds curves from another —
        # a leftover sweep silently appearing on the figure is the one way this
        # comparison can lie.
        every  = _bw_load_curves(t)
        stored = filter(c -> c.form == form, every)
        other  = sort(unique(c.form for c in every if c.form != form))
        isempty(other) || println(" #   brioWu1d: ", length(every) - length(stored),
                                  " stored curve(s) of another method (", join(other, ", "),
                                  ") are NOT drawn; this figure is the \"", form,
                                  "\" comparison. rm -r ", _bw_curve_dir(), " clears the store.")
    end
    if isempty(stored)
        # Not the final time (or an empty store): this run's own curve only.
        _bw_plot_density([(nop = nop, ndofs = length(xs), Cmin = Cmin, form = form, x = xs, y = ρs)],
                         xs, href, inputs, OUTPUT_DIR, iout, "")
    else
        # every comparison figure, the error table and the convergence history
        _bw_draw_store(stored, xs, href, inputs, OUTPUT_DIR, iout; hBy = hBy)
        println(" #   (", form, " form of the DynSGS coefficient)")
    end
    return nothing
end
