#---------------------------------------------------------------------------------
# Figures for the Kassam & Trefethen run, in place of the generic 1D panels.
#
# Called by plot_results (src/io/plotting/jeplots.jl) at every output time as
#
#     user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)
#
# with q the numerical (u, w). Writes
#
#   ks-it<n>.png       u and w ≈ ∂ₓ²u at this output time
#   ks-spacetime.png   u(x,t) from t = 0 to now, the picture of Kassam &
#                      Trefethen's Fig. 4.1, redrawn at every output time
#
# The frames of the space-time diagram are kept in ks_history (reset by
# initialize.jl), so set :diagnostics_at_times densely enough for it.
#---------------------------------------------------------------------------------
if !@isdefined(ks_history)
    const ks_history = (t = Float64[], x = Ref(Float64[]), u = Vector{Vector{Float64}}())
end

function ks_reset_history!()
    empty!(ks_history.t)
    empty!(ks_history.u)
    ks_history.x[] = Float64[]
    return nothing
end

function user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)

    idx = sortperm(x)
    xs  = x[idx]
    tt  = @sprintf "t = %.2f" t

    #-----------------------------------------------------------------------------
    # This output time.
    #-----------------------------------------------------------------------------
    labels = ("u", "w ≈ ∂ₓ²u")
    plts   = []
    for ivar = 1:2
        plt = Plots.plot(xs, q[idx, ivar];
                         line = (:blue, 2), legend = false,
                         title = string(labels[ivar], "   ", tt), xlabel = "x",
                         titlefontsize = 16, guidefontsize = 14, tickfontsize = 11,
                         show = false)
        push!(plts, plt)
    end
    fig = Plots.plot(plts...; layout = (2, 1), size = (1200, 800), margin = 5Plots.mm)
    _savefig_silent(fig, string(OUTPUT_DIR, "/ks-it", iout, ".png"))

    #-----------------------------------------------------------------------------
    # Space-time diagram so far.
    #-----------------------------------------------------------------------------
    if isempty(ks_history.t) || t > ks_history.t[end]
        ks_history.x[] = xs
        push!(ks_history.t, t)
        push!(ks_history.u, q[idx, 1])
    end
    if length(ks_history.t) >= 2
        U = reduce(hcat, ks_history.u)          # npoin × ntimes
        hm = Plots.heatmap(ks_history.x[], ks_history.t, permutedims(U);
                           c = :viridis, xlabel = "x", ylabel = "t",
                           title = "Kuramoto-Sivashinsky u(x,t)",
                           titlefontsize = 16, guidefontsize = 14, tickfontsize = 11,
                           size = (900, 900), show = false)
        _savefig_silent(hm, string(OUTPUT_DIR, "/ks-spacetime.png"))
    end

    @printf " # Kuramoto-Sivashinsky, t = %8.3f: max|u| = %.4f, max|w| = %.4f\n" t maximum(abs.(q[:, 1])) maximum(abs.(q[:, 2]))

    return nothing
end
