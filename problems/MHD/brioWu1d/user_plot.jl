#---------------------------------------------------------------------------------
# The figure of Dao & Nazarov (2022), Fig. 2: density against the reference
# solution with three zoom boxes (the foot of the fast rarefaction, the
# compound wave, the contact), in place of the generic multi-panel figure
# of the 1D plotter. Written as density-it<n>.png at every output time; the
# reference and the insets appear at the final time only, when
# user_analytic.jl supplies the reference. Set :plot_user => false in
# user_inputs.jl to get the generic fields-it<n>.png with every output
# variable and the DynSGS coefficient panel instead.
#
# MULTI-ORDER OVERLAY. At the final time the run also writes its own density
# profile to `curves/nop<N>.dat` in this case directory, and the figure is
# then drawn from EVERY curve stored there, not only from the run that is
# writing it. Running the case at several polynomial orders therefore builds
# one figure with all of them superimposed, each re-run replacing its own
# order's curve:
#
#     for N in 4 5 6 7; do
#         JEXPRESSO_BW_NOP=$N julia --project=. src/Jexpresso.jl MHD brioWu1d
#     done
#
# `JEXPRESSO_BW_NOP` also picks the element count that keeps the number of
# DOFs at ~600 (user_inputs.jl), so the orders are compared at equal cost.
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
    ((0.56, 0.59), (0.20, 0.35), (0.66, 0.28, 0.28, 0.30)),   # contact
]

# Where the per-order curves accumulate, one file per polynomial order.
const BW_CURVE_DIR = joinpath(@__DIR__, "curves")

# One (colour, line style, marker) per order, fixed so that a given order
# always looks the same from one figure to the next. EVERY numerical curve is
# broken (dash/dot); the reference is the only solid line, so the curves are
# told apart by line type and not by colour alone.
const BW_STYLE = Dict(
    4 => (:red,      :dash,       :circle),
    5 => (:blue,     :dashdot,    :rect),
    6 => (:seagreen, :dot,        :diamond),
    7 => (:purple,   :dashdotdot, :utriangle),
)
_bw_style(nop) = get(BW_STYLE, nop, (:darkorange, :dash, :xcross))

_bw_curve_file(nop) = joinpath(BW_CURVE_DIR, string("nop", nop, ".dat"))

#---------------------------------------------------------------------------------
# Store this run's density profile for the overlay. One file per order, so a
# re-run at the same order replaces its own curve and leaves the others.
#---------------------------------------------------------------------------------
function _bw_save_curve(xs, ρs, inputs, t)
    nop = Int(get(inputs, :nop, 0))
    nop > 0 || return nothing
    try
        mkpath(BW_CURVE_DIR)
        open(_bw_curve_file(nop), "w") do io
            println(io, "# Brio-Wu density profile, written by user_plot.jl for the multi-order overlay.")
            println(io, "# nop=", nop,
                        " ndofs=", length(xs),
                        " nelx=", get(inputs, :nelx, 0),
                        " t=", t,
                        " Cmin=", Float64(get(inputs, :dsgs_Cmin, 0.0)),
                        " CR=", Float64(get(inputs, :dsgs_CR, 1.0)),
                        " Cmax=", Float64(get(inputs, :dsgs_Cmax, 0.5)),
                        " sensor=", string(get(inputs, :dsgs_sensor, "residual")),
                        " nodal=", string(get(inputs, :ldsgs_nodal, false)))
            println(io, "# x rho")
            for i in eachindex(xs)
                println(io, xs[i], " ", ρs[i])
            end
        end
    catch err
        @warn "brioWu1d: could not store the density curve for the multi-order overlay" exception=err
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
        endswith(fname, ".dat") || continue
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
            @warn "brioWu1d: skipping an unreadable overlay curve" file=fname exception=err
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
                       x = xs, y = ys))
    end
    sort!(curves, by = c -> c.nop)
    return curves
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

# Legend entry: "nop 4, 601 DOFs" — with C_min appended only when the curves
# on the figure disagree about it (when they agree it goes in the title, once).
function _bw_label(c, lshow_cmin)
    s = string("\\mathrm{nop}\\ ", c.nop, ",\\ ", c.ndofs, "\\ \\mathrm{DOFs}")
    lshow_cmin && (s = string(s, ",\\ C_{min} = ", c.Cmin))
    return LaTeXStrings.latexstring(s)
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
    if lfinal
        _bw_save_curve(xs, ρs, inputs, t)
        curves = _bw_load_curves(t)
    else
        curves = NamedTuple[]
    end
    if isempty(curves)   # intermediate time, or the store could not be read
        curves = [(nop = nop, ndofs = length(xs), Cmin = Cmin, x = xs, y = ρs)]
    end

    # C_min in the title when every curve used the same one, in the legend
    # entries when they differ.
    cmins       = unique(c -> round(c.Cmin; digits = 12), curves)
    lcommon     = length(cmins) == 1
    lshow_cmin  = !lcommon
    ttl = lcommon ? LaTeXStrings.latexstring(string("\\mathrm{Density},\\ C_{min} = ", curves[1].Cmin)) :
                    LaTeXStrings.latexstring("\\mathrm{Density}")

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
    _savefig_silent(plt, string(OUTPUT_DIR, "/density-it", iout, ".png"))
    return nothing
end
