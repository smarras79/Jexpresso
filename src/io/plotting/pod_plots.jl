#---------------------------------------------------------------------------------
# pod_plots.jl — the standard POD figures, drawn on an equirectangular map.
#
# Three pictures, which together are what a POD is reported as:
#
#   pod_<field>_modes.png         the leading spatial modes φ_i(λ,φ), one panel
#                                 each, on a diverging colour scale centred on
#                                 zero and SYMMETRIC — a mode has no preferred
#                                 sign (pod_core.jl only fixes it by convention),
#                                 so a scale that is not symmetric about zero
#                                 invents structure that is not there.
#   pod_<field>_spectrum.png      the energy spectrum λ_i/Σλ on a logarithmic
#                                 axis, and the cumulative energy beside it. This
#                                 is the plot that says whether a reduced-order
#                                 model is possible at all: a spectrum that falls
#                                 off a cliff after a handful of modes means a
#                                 handful of ODEs can carry the flow, and one
#                                 that decays slowly means it cannot.
#   pod_<field>_coefficients.png  the temporal coefficients a_i(t), with the
#                                 (a_1,a_2) phase portrait beside them. The
#                                 portrait is there because the leading modes of
#                                 a TRAVELLING structure — which is what a
#                                 barotropically unstable jet produces — come in
#                                 near-degenerate PAIRS with λ_1 ≈ λ_2, whose
#                                 coefficients are in quadrature; the pair traces
#                                 a circle, and that circle is the propagation.
#                                 Two modes of a standing structure trace a line.
#
# The maps are equirectangular (see equirectangular.jl). The rasterizer renders
# the element tiling rather than interpolating scattered points, so what is drawn
# is the discrete mode itself: the colour scale is the mode's own range, and no
# small-scale structure has been smoothed away on the way to the picture.
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

using Plots

export pod_plot, plot_pod_modes, plot_pod_spectrum, plot_pod_coefficients
export plot_sphere_field


"""
    pod_plot(P, mesh, set, OUTPUT_DIR; verbose = true)

All three figures of one decomposition, plus the temporal mean when it was
subtracted. Serial only: the raster needs the whole sphere on one rank.
"""
function pod_plot(P::St_pod, mesh, set::St_pod_settings, OUTPUT_DIR::String;
                  verbose::Bool = true)

    isdir(OUTPUT_DIR) || mkpath(OUTPUT_DIR)

    plot_pod_modes(P, mesh, OUTPUT_DIR;
                   nmodes = set.nmodes_plot, nlon = set.nlon, nlat = set.nlat,
                   cmap = set.cmap)
    plot_pod_spectrum(P, OUTPUT_DIR)
    plot_pod_coefficients(P, OUTPUT_DIR;
                          nmodes = set.nmodes_plot, tscale = set.tscale, tlabel = set.tlabel)

    if P.lmean
        for c = 1:P.ncomp
            cn = P.ncomp == 1 ? "" : string("_", P.comps[c])
            plot_sphere_field(view(P.q̄, :, c), mesh,
                              joinpath(OUTPUT_DIR, string("pod_", P.name, "_mean", cn, ".png"));
                              title = string("temporal mean — ", P.comps[c]),
                              nlon = set.nlon, nlat = set.nlat,
                              cmap = :viridis, lsymmetric = false)
        end
    end

    verbose && @printf(" #     %s{_modes,_spectrum,_coefficients}.png\n",
                       joinpath(abspath(OUTPUT_DIR), string("pod_", P.name)))
    return nothing
end


"""
    plot_pod_modes(P, mesh, OUTPUT_DIR; kwargs...)

The leading `nmodes` spatial modes as equirectangular maps: one multi-panel
figure per component (`pod_<field>_modes.png`), and one PNG per mode beside it
for when a single mode has to go into a paper.
"""
function plot_pod_modes(P::St_pod, mesh, OUTPUT_DIR::String;
                        nmodes::Int = 6, nlon::Int = 720, nlat::Int = 360,
                        cmap::Symbol = :balance, ncols::Int = 2,
                        lindividual::Bool = true)

    r = min(nmodes <= 0 ? length(P.λ) : nmodes, length(P.λ))
    r >= 1 || return nothing
    grad = _pod_cgrad(cmap)

    for c = 1:P.ncomp
        cn    = P.ncomp == 1 ? "" : string("_", P.comps[c])
        plts  = Any[]

        for i = 1:r
            λg, φg, F = equirectangular_raster(view(P.Φ, :, c, i), mesh; nlon = nlon, nlat = nlat)
            m    = _robust_extreme(F)
            ttl  = @sprintf("mode %d — E = %.2f %%", i, 100*P.energy[i])
            plt  = _equirect_panel(λg, φg, F, ttl, grad, (-m, m))
            push!(plts, plt)

            if lindividual
                single = _equirect_panel(λg, φg, F,
                                         string(P.name, cn, "  ", ttl), grad, (-m, m);
                                         wide = true)
                _savefig_silent(single,
                                joinpath(OUTPUT_DIR, @sprintf("pod_%s%s_mode_%03d.png", P.name, cn, i)))
            end
        end

        nrows = cld(r, ncols)
        fig   = Plots.plot(plts...;
                           layout = (nrows, ncols),
                           size   = (700*ncols, 380*nrows),
                           plot_title = string("POD modes — ", P.name,
                                               P.ncomp == 1 ? "" : string(" (", P.comps[c], ")")),
                           plot_titlefontsize = 16,
                           # the panels carry axis labels of their own, and the
                           # default margins of a grid layout clip them
                           left_margin = 8Plots.mm, bottom_margin = 6Plots.mm,
                           top_margin = 3Plots.mm)
        _savefig_silent(fig, joinpath(OUTPUT_DIR, string("pod_", P.name, cn, "_modes.png")))
    end
    return nothing
end


"""
    plot_pod_spectrum(P, OUTPUT_DIR)

Energy per mode (log axis) and cumulative energy, side by side — the two curves
a truncation is chosen from. The 90 % and 99 % lines are drawn because those are
the thresholds a ROM order is normally quoted against.
"""
function plot_pod_spectrum(P::St_pod, OUTPUT_DIR::String)

    r = length(P.λ)
    r >= 1 || return nothing
    idx = collect(1:r)

    # A mode at exactly zero energy cannot be drawn on a log axis; it is also not
    # a mode. Clip to the smallest positive value rather than dropping points, so
    # the index axis stays the mode number.
    E    = 100 .* P.energy
    Epos = filter(>(0), E)
    floorE = isempty(Epos) ? 1.0e-16 : minimum(Epos)
    Eplot  = max.(E, floorE)

    pa = Plots.plot(idx, Eplot;
                    seriestype = :line, marker = (:circle, 5), line = (:solid, 2),
                    color = :black, yscale = :log10, legend = false,
                    xlabel = "mode index  i", ylabel = "E_i = λ_i / Σλ  [%]",
                    title  = "energy spectrum", framestyle = :box,
                    xlims  = (0, r+1), titlefontsize = 13,
                    guidefontsize = 11, tickfontsize = 10)

    pb = Plots.plot(idx, 100 .* P.cumenergy;
                    seriestype = :line, marker = (:circle, 5), line = (:solid, 2),
                    color = :darkblue, legend = false,
                    xlabel = "modes retained  r", ylabel = "Σ_{i≤r} E_i  [%]",
                    title  = "cumulative energy", framestyle = :box,
                    xlims  = (0, r+1), ylims = (0, 102), titlefontsize = 13,
                    guidefontsize = 11, tickfontsize = 10)
    # The labels go BELOW their line: the cumulative curve reaches the top of the
    # axis on the right, which is exactly where a label above the 99 % line would
    # sit.
    for (lev, col) in ((90.0, :orange), (99.0, :red))
        Plots.hline!(pb, [lev]; color = col, linestyle = :dash, label = "")
        n = pod_rank_for_energy(P, lev/100)
        Plots.annotate!(pb, 0.60*r, lev - 5.0,
                        Plots.text(@sprintf("%.0f %% : r = %d", lev, n), 9, col, :left))
    end

    fig = Plots.plot(pa, pb;
                     layout = (1, 2), size = (1100, 440),
                     plot_title = string("POD spectrum — ", P.name,
                                         @sprintf("  (%d snapshots)", P.nsnap)),
                     plot_titlefontsize = 15,
                     bottom_margin = 6Plots.mm, left_margin = 6Plots.mm)
    _savefig_silent(fig, joinpath(OUTPUT_DIR, string("pod_", P.name, "_spectrum.png")))
    return nothing
end


"""
    plot_pod_coefficients(P, OUTPUT_DIR; nmodes, tscale, tlabel)

`a_i(t)` for the leading modes, and the `(a_1,a_2)` phase portrait. See the
header on why the portrait is part of the standard set.
"""
function plot_pod_coefficients(P::St_pod, OUTPUT_DIR::String;
                               nmodes::Int = 6, tscale::Float64 = 1.0,
                               tlabel::String = "t")

    r = min(nmodes <= 0 ? length(P.λ) : nmodes, length(P.λ))
    r >= 1 || return nothing
    tt = P.t .* tscale

    pa = Plots.plot(; xlabel = tlabel, ylabel = "a_i(t)", framestyle = :box,
                    title = "temporal coefficients", titlefontsize = 13,
                    guidefontsize = 11, tickfontsize = 10, legendfontsize = 9,
                    legend = :outertopright)
    for i = 1:r
        Plots.plot!(pa, tt, P.a[:, i];
                    line = (:solid, 2), label = @sprintf("a%d  (%.1f %%)", i, 100*P.energy[i]))
    end

    fig = if r >= 2
        pb = Plots.plot(P.a[:,1], P.a[:,2];
                        line = (:solid, 2), marker = (:circle, 3), color = :black,
                        legend = false, aspect_ratio = :equal, framestyle = :box,
                        xlabel = "a₁", ylabel = "a₂",
                        title = "phase portrait (a₁, a₂)", titlefontsize = 13,
                        guidefontsize = 11, tickfontsize = 10)
        Plots.scatter!(pb, [P.a[1,1]], [P.a[1,2]]; marker = (:star5, 9), color = :red, label = "")
        Plots.plot(pa, pb; layout = (1, 2), size = (1200, 450))
    else
        Plots.plot(pa; size = (800, 450))
    end

    fig = Plots.plot(fig;
                     plot_title = string("POD coefficients — ", P.name),
                     plot_titlefontsize = 15,
                     bottom_margin = 6Plots.mm, left_margin = 6Plots.mm)
    _savefig_silent(fig, joinpath(OUTPUT_DIR, string("pod_", P.name, "_coefficients.png")))
    return nothing
end


"""
    plot_sphere_field(f, mesh, fout_name; kwargs...)

One nodal field on an equirectangular map. Used here for the POD mean, and
usable on its own for any field on the shell — the vorticity at a given output
time, a reconstruction error, the difference between two runs.

`lsymmetric = true` centres the colour scale on zero, which is what a signed
anomaly wants and what a positive-definite field (a depth, a speed) does not.
"""
function plot_sphere_field(f::AbstractVector, mesh, fout_name::String;
                           title::String = "", nlon::Int = 720, nlat::Int = 360,
                           cmap::Symbol = :viridis, lsymmetric::Bool = true,
                           clims = nothing)

    λg, φg, F = equirectangular_raster(f, mesh; nlon = nlon, nlat = nlat)
    grad = _pod_cgrad(cmap)

    cl = if clims !== nothing
        (Float64(clims[1]), Float64(clims[2]))
    elseif lsymmetric
        m = _robust_extreme(F); (-m, m)
    else
        finite = filter(isfinite, F)
        lo, hi = isempty(finite) ? (0.0, 1.0) : (minimum(finite), maximum(finite))
        hi > lo ? (lo, hi) : (lo - 0.5, hi + 0.5)
    end

    plt = _equirect_panel(λg, φg, F, title, grad, cl; wide = true)
    _savefig_silent(plt, fout_name)
    return nothing
end


#---------------------------------------------------------------------------------
# helpers
#---------------------------------------------------------------------------------
#
# One map panel. Ticks every 60° in longitude and 30° in latitude, equal aspect
# so the canvas is the 2:1 rectangle a plate-carrée map is, and the data clamped
# to the colour range so that a value outside it takes the end colour instead of
# being drawn as a hole (the same convention as the flat-case plotter).
#
function _equirect_panel(λg, φg, F, ttl::String, grad, clims; wide::Bool = false)
    Fc = clamp.(F, clims[1], clims[2])
    return Plots.contourf(λg, φg, Fc';
                          color = grad, clims = clims, levels = 31, linewidth = 0,
                          colorbar = true, legend = false,
                          aspect_ratio = :equal,
                          xlims = (-180, 180), ylims = (-90, 90),
                          xticks = -180:60:180, yticks = -90:30:90,
                          framestyle = :box,
                          xlabel = "longitude [deg]", ylabel = "latitude [deg]",
                          title = ttl, titlefontsize = wide ? 14 : 12,
                          guidefontsize = 10, tickfontsize = 9,
                          size = wide ? (960, 480) : (700, 380),
                          # A standalone panel pays for its own margins; one that
                          # goes into a grid gets them from the outer plot. Without
                          # the right margin the colour-bar tick labels of a field
                          # with small values run off the canvas.
                          left_margin   = wide ? 5Plots.mm : 0Plots.mm,
                          right_margin  = wide ? 10Plots.mm : 0Plots.mm,
                          bottom_margin = wide ? 5Plots.mm : 0Plots.mm,
                          show = false)
end

#
# The colour range of a mode. maximum(abs, ·) would hand the whole scale to a
# single node whenever one exists — and the POD of a high-order field regularly
# has one, since the modes are as spiky as the solution is. The 99.8th percentile
# keeps the structure visible; values beyond it are clamped, not hidden.
#
function _robust_extreme(F::AbstractArray; q::Float64 = 0.998)
    v = Float64[abs(x) for x in F if isfinite(x)]
    isempty(v) && return 1.0
    sort!(v)
    m = v[clamp(ceil(Int, q*length(v)), 1, length(v))]
    m > 0 || (m = maximum(v))
    m > 0 || (m = 1.0)
    return m
end

#
# Colour scheme by name, with a fallback: :balance is cmocean's diverging map and
# ships with ColorSchemes, but a Plots backend that does not know a name throws,
# and a missing colour map is not a reason to lose the figure.
#
function _pod_cgrad(sym::Symbol)
    try
        return Plots.cgrad(sym)
    catch
        @warn string("POD: colour map :", sym, " is not available; falling back to :RdBu.")
        return Plots.cgrad(:RdBu, rev = true)
    end
end
