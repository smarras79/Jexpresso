#---------------------------------------------------------------------------------
# Figures for the Kassam & Trefethen run, in place of the generic 1D panels.
#
# Called by plot_results (src/io/plotting/jeplots.jl) at every output time as
#
#     user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)
#
# with q the numerical u. Writes
#
#   ks-it<n>.png       u at this output time, every :ks_frame_every-th output
#                      (default 10, and the last one; 0 = none)
#   ks-spacetime.png   u(x,t) over the whole run, the picture of Kassam &
#                      Trefethen's Fig. 4.1, drawn once, at tend
#
# Cost: the time loop with output off takes 0.02 s; a frame takes ~30 ms
# (GR rasterizing and PNG-encoding 1200×450 pixels), and the space-time
# diagram ~0.12 s. The default 31 frames are therefore most of a run with
# output on (1.2 s): :ks_frame_every => 0 keeps the diagram alone (0.14 s).
# The other outputs only store u for the diagram.
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
# The rows of the space-time diagram are kept in ks_history (reset by
# initialize.jl), so set :diagnostics_at_times densely enough for it.
#
# :ks_plot => false (user_inputs.jl) switches this hook off: it returns at once,
# with no figures, no history and no diagnostic line. (:plot_user => false
# would instead hand the output to the generic 1D panels, which also plot;
# :outformat => "none" switches off all output.)
#---------------------------------------------------------------------------------
#
# The history is created only once per session: re-including this file (every
# run of this case re-includes it) must not redefine the constants.
# ks_order is the order of the nodes by x, computed at the first output
# of a run (the nodes do not move).
#---------------------------------------------------------------------------------
if !@isdefined(ks_history)
    const ks_history = (t = Float64[], x = Ref(Float64[]), u = Vector{Vector{Float64}}())
end
if !@isdefined(ks_order)
    const ks_order = Ref(Int[])
end

function ks_reset_history!()
    empty!(ks_history.t)
    empty!(ks_history.u)
    ks_history.x[] = Float64[]
    ks_order[] = Int[]
    return nothing
end

function user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)

    get(inputs, :ks_plot, true)::Bool || return nothing

    if length(ks_order[]) != length(x)
        ks_order[]     = sortperm(x)
        ks_history.x[] = x[ks_order[]]
    end
    idx = ks_order[]
    xs  = ks_history.x[]
    us  = q[idx, 1]

    nframe = get(inputs, :ks_frame_every, 10)::Int
    lfinal = t >= inputs[:tend] - 1.0e-10
    ldraw  = nframe > 0 && (mod(iout - 1, nframe) == 0 || lfinal)

    #-----------------------------------------------------------------------------
    # This output time.
    #-----------------------------------------------------------------------------
    if ldraw
        fig = Plots.plot(xs, us;
                         line = (:blue, 2), legend = false,
                         title = @sprintf("u   t = %.2f", t), xlabel = "x",
                         titlefontsize = 16, guidefontsize = 14, tickfontsize = 11,
                         size = (1200, 450), margin = 5Plots.mm, show = false)
        _savefig_silent(fig, string(OUTPUT_DIR, "/ks-it", iout, ".png"))
    end

    #-----------------------------------------------------------------------------
    # Space-time diagram.
    #-----------------------------------------------------------------------------
    if isempty(ks_history.t) || t > ks_history.t[end]
        push!(ks_history.t, t)
        push!(ks_history.u, us)
    end
    if lfinal && length(ks_history.t) >= 2
        xn = xs
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

    @printf " # Kuramoto-Sivashinsky, t = %8.3f: max|u| = %.4f\n" t maximum(abs, us)

    return nothing
end
