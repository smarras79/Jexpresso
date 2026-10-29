#---------------------------------------------------------------------------------
# Figure for the soliton test, in place of the generic 1D panels.
#
# Called by plot_results (src/io/plotting/jeplots.jl) at every output time as
#
#     user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)
#
# with q the numerical (u, v, w) and qref the exact KdV soliton and its first
# and second derivatives from user_analytic.jl. Writes kdv-it<n>.png: u, v, w
# against the exact curves, and the pointwise error u - η.
#
# It also prints the max-norm error of each field, relative to the field's own
# peak, so that convergence in τ, Δt or the mesh can be read off stdout without
# opening a figure. At t = 0 the error is the interpolation error only.
#---------------------------------------------------------------------------------
function user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)

    idx = sortperm(x)
    xs  = x[idx]
    tt  = @sprintf "t = %.4f" t

    labels = ("u ≈ η", "v ≈ ∂ₓη", "w ≈ ∂ₓ²η")
    plts   = []
    errs   = zeros(3)

    for ivar = 1:3
        plt = Plots.plot(xs, q[idx, ivar];
                         line = (:blue, 2), label = "Jexpresso",
                         title = string(labels[ivar], "   ", tt), xlabel = "x",
                         titlefontsize = 16, guidefontsize = 14,
                         legendfontsize = 11, tickfontsize = 11,
                         legend = :topright, show = false)
        if qref !== nothing
            Plots.plot!(plt, xs, qref[idx, ivar];
                        line = (:black, 2, :dash), label = "exact KdV")
            errs[ivar] = maximum(abs.(q[:, ivar] .- qref[:, ivar])) /
                         maximum(abs.(qref[:, ivar]))
        end
        push!(plts, plt)
    end

    if qref !== nothing
        plt = Plots.plot(xs, q[idx, 1] .- qref[idx, 1];
                         line = (:red, 2), legend = false,
                         title = string("u - η   ", tt), xlabel = "x",
                         titlefontsize = 16, guidefontsize = 14, tickfontsize = 11,
                         show = false)
        push!(plts, plt)

        @printf " # KdV soliton, t = %8.4f: max|err|/max|exact|  u: %.3e  v: %.3e  w: %.3e\n" t errs[1] errs[2] errs[3]
    end

    fig = Plots.plot(plts...; layout = (2, 2), size = (1400, 900), margin = 5Plots.mm)
    _savefig_silent(fig, string(OUTPUT_DIR, "/kdv-it", iout, ".png"))

    return nothing
end
