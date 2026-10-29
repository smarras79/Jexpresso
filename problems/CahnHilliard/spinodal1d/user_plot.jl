#---------------------------------------------------------------------------------
# Figure for the spinodal test, in place of the generic 1D panels.
#
# Called by plot_results (src/io/plotting/jeplots.jl) at every output time as
#
#     user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)
#
# with q the numerical (c, μ) and qref the linear solution of user_analytic.jl
# (all NaN once the linear phase is over). Writes ch-it<n>.png with c and μ,
# against the linear solution while it applies.
#
# Prints, at every output time, the amplitude max|c| and the growth rate it
# implies, ln(max|c|/δ)/t, next to the Cahn-Hilliard rate σ(k), and, while
# the linear solution applies, the max-norm error of c and μ relative to their
# peaks.
#---------------------------------------------------------------------------------
function user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)

    p   = ch_parameters()
    σ   = ch_growth_rate(p.k)
    idx = sortperm(x)
    xs  = x[idx]
    tt  = @sprintf "t = %.4f" t

    llinear = qref !== nothing && all(isfinite, qref)

    labels = ("c", "μ")
    plts   = []
    errs   = zeros(2)

    for ivar = 1:2
        plt = Plots.plot(xs, q[idx, ivar];
                         line = (:blue, 2), label = "Jexpresso",
                         title = string(labels[ivar], "   ", tt), xlabel = "x",
                         titlefontsize = 16, guidefontsize = 14,
                         legendfontsize = 11, tickfontsize = 11,
                         legend = :topright, show = false)
        if llinear
            Plots.plot!(plt, xs, qref[idx, ivar];
                        line = (:black, 2, :dash), label = "linear CH")
            errs[ivar] = maximum(abs.(q[:, ivar] .- qref[:, ivar])) /
                         maximum(abs.(qref[:, ivar]))
        end
        push!(plts, plt)
    end

    A = maximum(abs.(q[:, 1]))
    if t > 0.0
        @printf " # Cahn-Hilliard, t = %7.4f: max|c| = %.4e, ln(max|c|/δ)/t = %.5f (σ = %.5f)" t A log(A/p.δ)/t σ
    else
        @printf " # Cahn-Hilliard, t = %7.4f: max|c| = %.4e" t A
    end
    if llinear
        @printf ", max|err|/max|exact|  c: %.3e  μ: %.3e" errs[1] errs[2]
    end
    @printf "\n"

    fig = Plots.plot(plts...; layout = (2, 1), size = (1000, 900), margin = 5Plots.mm)
    _savefig_silent(fig, string(OUTPUT_DIR, "/ch-it", iout, ".png"))

    return nothing
end
