#---------------------------------------------------------------------------------
# The figure of Dao & Nazarov (2022), Fig. 2: density against the reference
# solution with four zoom boxes (the foot of the fast rarefaction, the
# compound wave, the contact, the slow shock), in place of the generic
# multi-panel figure of the 1D plotter. Written as density-it<n>.png at every
# output time; the reference and the insets appear at the final time only,
# when user_analytic.jl supplies the reference. Set :plot_user => false in
# user_inputs.jl to get the generic fields-it<n>.png with every output
# variable and the DynSGS coefficient panel instead.
#
# ORDER COMPARISON AND CONVERGENCE HISTORY. At the final time the run stores
# its density profile in `curves/nop<N>_dof<M>_<form>.dat` in this case
# directory, and the figures are drawn from EVERY curve stored there whose
# <form> matches this run's — the coefficient form (element or nodal) and the
# length-scale convention ("H" = the paper's Δ_K/k, the default of this deck),
# which are what make one run a different method from another. Curves of any
# other form are reported on stdout and left off the figure:
#
#   density_dof<M>-it<n>.png      ONE FILE PER RESOLUTION: every order that ran
#                                 at ~M degrees of freedom, against the
#                                 reference, with the zoom boxes. A sweep
#                                 leaves density_dof150, density_dof300, ...
#   density-it<n>.png             the same for the finest resolution in the
#                                 store, under the plain name
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

# Where the per-(order, resolution) curves accumulate.
const BW_CURVE_DIR = joinpath(@__DIR__, "curves")

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
# The deck folds the Δ_K/k convention into the coefficients as a factor
# (k+1)/k on C_max, C_min. Divide it back out so the figure reports the
# method's C_min and not the per-order number that implements it.
# `_bw_hscale()` is the deck's own switch (user_inputs.jl); it defaults to the
# paper's "nop" here and to nothing at all anywhere else in the code.
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
    joinpath(BW_CURVE_DIR, string("nop", nop, "_dof", ndofs, "_", form, ".dat"))

#---------------------------------------------------------------------------------
# Store this run's density profile. One file per (order, DOFs), so a re-run at
# the same order and resolution replaces its own point and leaves the others.
#---------------------------------------------------------------------------------
function _bw_save_curve(xs, ρs, inputs, t)
    nop = Int(get(inputs, :nop, 0))
    nop > 0 || return nothing
    try
        mkpath(BW_CURVE_DIR)
        open(_bw_curve_file(nop, length(xs), _bw_form(inputs)), "w") do io
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
                        " form=", _bw_form(inputs))
            println(io, "# x rho")
            for i in eachindex(xs)
                println(io, xs[i], " ", ρs[i])
            end
        end
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
    isdir(BW_CURVE_DIR) || return curves
    for fname in readdir(BW_CURVE_DIR)
        (endswith(fname, ".dat") && startswith(fname, "nop")) || continue
        meta = Dict{String,String}()
        xs = Float64[]; ys = Float64[]
        try
            for line in eachline(joinpath(BW_CURVE_DIR, fname))
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
                       mtime = mtime(joinpath(BW_CURVE_DIR, fname)),
                       x = xs, y = ys))
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
        mkpath(BW_CURVE_DIR)
        open(joinpath(BW_CURVE_DIR, "errors.dat"), "w") do io
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
function _bw_plot_convergence(rows, OUTPUT_DIR, iout; smooth::Bool)
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
                            "\\|\\rho_h-\\rho_{ref}\\|_{", nm, "}\\ /\\ \\|\\rho_{ref}\\|_{", nm, "}")),
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
    _savefig_silent(plt, string(OUTPUT_DIR, "/convergence", smooth ? "_smooth" : "", "-it", iout, ".png"))
    return nothing
end

# Legend entry: "nop 4, 601 DOFs" — with C_min appended only when the curves
# on the figure disagree about it (when they agree it goes in the title, once).
function _bw_label(c, lshow_cmin)
    s = string("\\mathrm{nop}\\ ", c.nop, ",\\ ", c.ndofs, "\\ \\mathrm{DOFs}")
    lshow_cmin && (s = string(s, ",\\ C_{min} = ", c.Cmin))
    return LaTeXStrings.latexstring(s)
end

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
        Plots.plot!(plt, xs, href; line = (:black, 1.8, :solid), label = "Reference solution")

        # Zoom boxes: a gray frame on the main axes, a dashed connector to the
        # inset, and the inset itself (every curve, tick labels only).
        for (k, (xr, yr, pos)) in enumerate(BW_INSETS)
            bx = [xr[1], xr[2], xr[2], xr[1], xr[1]]
            by = [yr[1], yr[1], yr[2], yr[2], yr[1]]
            Plots.plot!(plt, bx, by; line = (:gray, 1.0), label = "")
            # The inset rectangle in data coordinates (pos are fractions of
            # the axes span from the bottom-left); the dashed connector joins
            # the box corner nearest the inset to the inset corner nearest
            # the box, as in the paper.
            ix = (xl[1] + pos[1]*(xl[2] - xl[1]), xl[1] + (pos[1] + pos[3])*(xl[2] - xl[1]))
            iy = (yl[1] + pos[2]*(yl[2] - yl[1]), yl[1] + (pos[2] + pos[4])*(yl[2] - yl[1]))
            bcx = 0.5*(xr[1] + xr[2]); bcy = 0.5*(yr[1] + yr[2])
            icx = 0.5*(ix[1] + ix[2]); icy = 0.5*(iy[1] + iy[2])
            cornerx = (icx < bcx) ? xr[1] : xr[2]
            cornery = (icy < bcy) ? yr[1] : yr[2]
            insx    = (icx < bcx) ? ix[2] : ix[1]
            insy    = (icy < bcy) ? iy[2] : iy[1]
            Plots.plot!(plt, [cornerx, insx], [cornery, insy]; line = (:black, 1.0, :dash), label = "")
            Plots.plot!(plt; inset = (1, Plots.bbox(pos[1], pos[2], pos[3], pos[4], :bottom, :left)),
                        subplot = k + 1)
            sel = (xs .>= xr[1]) .& (xs .<= xr[2])
            Plots.plot!(plt[k + 1], xs[sel], href[sel]; line = (:black, 1.8, :solid), label = "",
                        xlims = xr, ylims = yr, framestyle = :box, grid = true,
                        tickfontsize = 9, background_color_inside = :white)
            for c in curves
                col, ls, _ = _bw_style(c.nop)
                s = (c.x .>= xr[1]) .& (c.x .<= xr[2])
                any(s) && Plots.plot!(plt[k + 1], c.x[s], c.y[s]; line = (col, 1.6, ls), label = "")
            end
        end
    end
    _savefig_silent(plt, string(OUTPUT_DIR, "/density", tag, "-it", iout, ".png"))
    return nothing
end

function user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)
    iρ = findfirst(==("ρ"), string.(outvar))
    iρ === nothing && error("user_plot_1d (brioWu1d): no ρ among the output variables")
    idx = sortperm(x)
    xs  = x[idx]
    ρs  = q[idx, iρ]
    href = (qref !== nothing && any(isfinite, @view(qref[:, iρ]))) ? qref[idx, iρ] : nothing

    nop   = Int(get(inputs, :nop, 0))
    Cmin  = Float64(get(inputs, :dsgs_Cmin, 0.0))

    # The reference is supplied at the final time only; that is the figure the
    # orders are compared on, so that is where the curve is stored and where
    # every stored curve is drawn.
    lfinal = href !== nothing
    stored = NamedTuple[]
    form   = _bw_form(inputs)
    if lfinal
        _bw_save_curve(xs, ρs, inputs, t)
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
                                  "\" comparison. rm problems/MHD/brioWu1d/curves to clear the store.")
    end
    if isempty(stored)
        # Not the final time (or an empty store): this run's own curve only.
        _bw_plot_density([(nop = nop, ndofs = length(xs), Cmin = Cmin, form = form, x = xs, y = ρs)],
                         xs, href, inputs, OUTPUT_DIR, iout, "")
    else
        # ONE FIGURE PER RESOLUTION. Orders may only be compared at equal
        # degrees of freedom, so each DOF bucket gets its own file and a sweep
        # leaves density_dof150, density_dof300, ... behind, each holding every
        # order that ran at that resolution.
        for b in sort(unique(_bw_bucket(c.ndofs) for c in stored))
            sub = _bw_finest(filter(c -> _bw_bucket(c.ndofs) == b, stored))
            _bw_plot_density(sub, xs, href, inputs, OUTPUT_DIR, iout, string("_dof", b))
        end
        # The finest resolution keeps the plain name as the headline figure.
        _bw_plot_density(_bw_finest(stored), xs, href, inputs, OUTPUT_DIR, iout, "")
    end

    # Error table and the convergence history, from every stored curve.
    if lfinal && !isempty(stored)
        rows = _bw_error_table(stored)
        _bw_report_errors(rows)
        _bw_plot_convergence(rows, OUTPUT_DIR, iout; smooth = false)
        _bw_plot_convergence(rows, OUTPUT_DIR, iout; smooth = true)
        println(" #   (", form, " form of the DynSGS coefficient)")
    end
    return nothing
end
