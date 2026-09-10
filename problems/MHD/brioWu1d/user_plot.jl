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

function user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)
    iρ = findfirst(==("ρ"), string.(outvar))
    iρ === nothing && error("user_plot_1d (brioWu1d): no ρ among the output variables")
    idx = sortperm(x)
    xs  = x[idx]
    ρs  = q[idx, iρ]
    href = (qref !== nothing && any(isfinite, @view(qref[:, iρ]))) ? qref[idx, iρ] : nothing
    npoin = length(x)
    nop   = get(inputs, :nop, 0)
    lbl   = string("\\mathbb{P}_{", nop, "},\\ ", npoin, "\\ \\mathrm{DOFs}")

    xl = (0.0, 1.0); yl = (0.1, 1.0)
    plt = Plots.plot(xs, ρs;
                     line = (:red, 1.5), label = LaTeXStrings.latexstring(lbl),
                     title = "Density", xlabel = LaTeXStrings.L"x", ylabel = LaTeXStrings.L"\rho",
                     xlims = xl, ylims = yl,
                     xticks = 0.0:0.2:1.0, yticks = 0.1:0.1:1.0,
                     framestyle = :box, grid = false,
                     legend = :topright, legendfontsize = 12,
                     titlefontsize = 20, guidefontsize = 18, tickfontsize = 13,
                     size = (900, 720), left_margin = 6Plots.mm, bottom_margin = 5Plots.mm,
                     show = false)
    if href !== nothing
        Plots.plot!(plt, xs, href; line = (:black, 1.5), label = "Reference solution")

        # Zoom boxes: a gray frame on the main axes, a dashed connector to the
        # inset, and the inset itself (both curves, tick labels only).
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
            Plots.plot!(plt[k + 1], xs[sel], href[sel]; line = (:black, 1.5), label = "",
                        xlims = xr, ylims = yr, framestyle = :box, grid = true,
                        tickfontsize = 9, background_color_inside = :white)
            Plots.plot!(plt[k + 1], xs[sel], ρs[sel]; line = (:red, 1.5), label = "")
        end
    end
    _savefig_silent(plt, string(OUTPUT_DIR, "/density-it", iout, ".png"))
    return nothing
end
