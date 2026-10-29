#---------------------------------------------------------------------------------
# Figures for the Kassam & Trefethen run, in place of the generic 1D panels.
#
# Called by plot_results (src/io/plotting/jeplots.jl) at every output time as
#
#     user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)
#
# with q the numerical (u, w). Writes
#
#   ks-it<n>.png       u and w ≈ ∂ₓ²u at this output time, every
#                      ks_frame_every()-th output
#   ks-spacetime.png   u(x,t) over the whole run, the picture of Kassam &
#                      Trefethen's Fig. 4.1, drawn once, at tend. Redrawing
#                      the growing diagram during the run cost 4.7 M of the
#                      5.6 M allocations of this hook (measured, 301 outputs)
#
# The space-time diagram is drawn on the NATIVE grid, not resampled: every
# LGL node is one column of cells, bounded by the midpoints to its neighbours
# (so the columns have the non-uniform LGL widths), and every output time is
# one row, bounded by the midpoints to the neighbouring output times. On the
# periodic domain the node at x = 0 is repeated at x = L, so the gap between
# the last node and L belongs to it.
#
# Colours: inferno, saturated. :ks_clims (default (-1.9, 1.0)) clips the
# colour range well inside max|u| ≈ 2.6, so crests saturate to yellow and
# troughs to black, and u = 0 maps to orange.
#
# The frames of the space-time diagram are kept in ks_history (reset by
# initialize.jl), so set :diagnostics_at_times densely enough for it.
#---------------------------------------------------------------------------------
#
# A method and not a `const` for the frame interval, so that re-including this
# file in a session where ks_history already exists (from an earlier run of
# this case) still defines it; the `if` below skips everything inside it then.
#---------------------------------------------------------------------------------
@inline ks_frame_every() = 10

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

    lfinal = t >= inputs[:tend] - 1.0e-10
    ldraw  = mod(iout - 1, ks_frame_every()) == 0 || lfinal

    #-----------------------------------------------------------------------------
    # This output time.
    #-----------------------------------------------------------------------------
    if ldraw
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
    end

    #-----------------------------------------------------------------------------
    # Space-time diagram so far.
    #-----------------------------------------------------------------------------
    if isempty(ks_history.t) || t > ks_history.t[end]
        ks_history.x[] = xs
        push!(ks_history.t, t)
        push!(ks_history.u, q[idx, 1])
    end
    if lfinal && length(ks_history.t) >= 2
        xn = ks_history.x[]
        tn = ks_history.t
        U  = permutedims(reduce(hcat, ks_history.u))     # ntimes × npoin
        if inputs[:lperiodic_1d]
            xn = [xn; inputs[:xmax]]
            U  = hcat(U, U[:, 1])
        end
        # cell edges: the nodes' midpoints, closed by the first and last node
        xe = [xn[1]; 0.5 .* (xn[1:end-1] .+ xn[2:end]); xn[end]]
        te = [tn[1]; 0.5 .* (tn[1:end-1] .+ tn[2:end]); tn[end]]

        hm = Plots.heatmap(xe, te, U;
                           c = :inferno, clims = get(inputs, :ks_clims, (-1.9, 1.0)),
                           xlabel = "x", ylabel = "t",
                           title = "Kuramoto-Sivashinsky u(x,t)",
                           titlefontsize = 16, guidefontsize = 14, tickfontsize = 11,
                           size = (900, 900), dpi = 150, show = false)
        _savefig_silent(hm, string(OUTPUT_DIR, "/ks-spacetime.png"))
    end

    @printf " # Kuramoto-Sivashinsky, t = %8.3f: max|u| = %.4f, max|w| = %.4f\n" t maximum(abs, @view(q[:, 1])) maximum(abs, @view(q[:, 2]))

    return nothing
end
