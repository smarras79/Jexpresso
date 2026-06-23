using Dierckx
using Plots
using LaTeXStrings
using ColorSchemes
#using Interpolations

#
# Curves (1D) or Contours (2D) with Plots.jl
#

#
# Write a figure to file without touching the screen. With the GR
# backend, the png export of a figure whose :overwrite_figure attribute
# is true (the Plots default) first draws to the ACTIVE workstation "to
# set the canvas viewport" -- that is what flashes one gksqt window per
# variable. With :overwrite_figure = false the export goes through a
# dedicated file workstation and never appears on screen. NOTE: this
# path also closes the GKS session, which takes any open gksqt window
# with it -- use it only when no live window is wanted (see
# render_plot_matrix).
#
function _savefig_silent(plt, fout_name)
    plt[:overwrite_figure] = false
    Plots.savefig(plt, string(fout_name))
    return nothing
end

function plot_initial(SD::NSD_1D, x, q, ivar, OUTPUT_DIR::String)

    npoin = length(q)
    plt = Plots.scatter(x[1:npoin], q[1:npoin];
                        markersize = 5,
                        color = :blue,
                        xlabel = "x",
                        ylabel = "q(x)",
                        title = "u",
                        titlefontsize = 24,
                        guidefontsize = 18,
                        legendfontsize = 14,
                        tickfontsize = 14,
                        legend = false,
                        size = (800, 600))

    fout_name = string(OUTPUT_DIR, "/INIT-", ivar, ".png")

    _savefig_silent(plt, fout_name)
    plt
end

function plot_results(SD::NSD_1D, mesh::St_mesh, q, title::String, OUTPUT_DIR::String, outvar, inputs; iout=1, nvar=1, PT=nothing, μ_nodes=nothing)

    epsi = 1.1
    npoin = mesh.npoin

    lmatrix = get(inputs, :plot_matrix, true)

    qout = reshape(q, npoin, nvar)   # 2D view (3x4) - NO allocation
    x_coords = mesh.coords[1:npoin, 1]
    sort_idx = sortperm(x_coords)
    plts = []
    for ivar=1:nvar

        idx = (ivar - 1)*npoin

        plt = Plots.plot(x_coords[sort_idx], qout[sort_idx, ivar];
                        line = (:blue, 2),
                        marker = (:circle, 5, :blue),
                        title = string(outvar[ivar], "  ", title),
                        xlabel = "x",
                        titlefontsize = 22,
                        guidefontsize = 18,
                        legendfontsize = 14,
                        tickfontsize = 14,
                        legend = false,
                        show = false,
                        size = (600, 400))

        vlines = inputs[:plot_vlines]
        hlines = inputs[:plot_hlines]
        axis = inputs[:plot_axis]
        if !(vlines == "empty")
            for i=1:size(vlines,1)
                Plots.vline!(plt, [vlines[i]]; color = :red, linestyle = :solid, label = "")
            end
        end
        if !(hlines == "empty")
            for i=1:size(hlines,1)
                Plots.hline!(plt, [hlines[i]]; color = :red, linestyle = :solid, label = "")
            end
        end
        if !(axis == "empty")
            idx = (ivar-1)*2
            Plots.ylims!(plt, axis[1+idx], axis[2+idx])
        end
        if !lmatrix
            fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")
            _savefig_silent(plt, fout_name)
        end
        push!(plts, plt)
    end

    # DSGS runs: show the per-element viscosity staircase as one more
    # panel of the same output time.
    if μ_nodes !== nothing
        ieq = min(2, size(μ_nodes, 2))
        plt_μ = Plots.plot(x_coords[sort_idx], μ_nodes[sort_idx, ieq];
                           line = (:red, 2),
                           marker = (:circle, 3, :red),
                           title = string("μ_dsgs  ", title),
                           xlabel = "x",
                           titlefontsize = 22,
                           guidefontsize = 18,
                           tickfontsize = 14,
                           legend = false,
                           show = false,
                           size = (600, 400))
        if !lmatrix
            _savefig_silent(plt_μ, string(OUTPUT_DIR, "/mu_dsgs-it", iout, ".png"))
        end
        push!(plts, plt_μ)
    end

    render_plot_matrix(lmatrix, plts, OUTPUT_DIR, iout; wfig=600, hfig=400)
end

#
# Combine the per-variable figures of one output time into a single
# plot-matrix figure and render it ONCE with a plain savefig. With the
# GR backend a savefig of a figure whose :overwrite_figure attribute is
# true (the Plots default) paints the active workstation in place
# (clearws/draw/updatews -- the GKS session stays open, so the gksqt
# window is replaced on the fly, never closed and reopened) and then
# prints the very same canvas to fields-it<iout>.png. In a headless run
# (GKSwstype=100/nul) the screen workstation is inert and only the file
# is produced. This is the only GR-friendly way to have BOTH a live
# window and file output without flicker; per-variable files instead
# require the silent export path, which closes the GKS session and with
# it the window -- that is what :plot_matrix => false selects.
#
function render_plot_matrix(lmatrix, plts, OUTPUT_DIR, iout; wfig=600, hfig=400)
    lmatrix || return nothing
    nplt = length(plts)
    nplt == 0 && return nothing
    comm    = get_mpi_comm()
    mpisize = MPI.Comm_size(comm)
    piece   = mpisize > 1 ? string("-rank", MPI.Comm_rank(comm)) : ""
    try
        ncols = ceil(Int, sqrt(nplt))
        nrows = ceil(Int, nplt/ncols)
        figm  = Plots.plot(plts...;
                           layout = (nrows, ncols),
                           show = false,
                           size = (ncols*wfig, nrows*hfig))
        Plots.savefig(figm, string(OUTPUT_DIR, "/fields", piece, "-it", iout, ".png"))
    catch
    end
    return nothing
end


function plot_results!(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, outvar, inputs; iout=1, nvar=1, fig=nothing, color ="blue", p=[], marker = :circle, PT=nothing)
    
    epsi = 1.1
    npoin = mesh.npoin

    for ivar=1:1
        idx = (ivar - 1)*npoin

        if fig === nothing
            fig = Plots.plot(xlabel = "x",
                           title = string(outvar[ivar]),
                           titlefontsize = 18,
                           guidefontsize = 14,
                           legend = false,
                           show = false)
        end

        if !(p==[])
            # Add to existing plot without decorations
            Plots.scatter!(fig, mesh.x[1:mesh.npoin_original],
                         q[idx+1:(ivar-1)*npoin+mesh.npoin_original];
                         marker = marker,
                         markersize = 5,
                         color = color,
                         label = "")
        else
            Plots.scatter!(fig, mesh.x[1:mesh.npoin_original],
                         q[idx+1:(ivar-1)*npoin+mesh.npoin_original];
                         marker = marker,
                         markersize = 5,
                         color = color,
                         label = "")
        end

        Plots.ylims!(fig, -0.03, 0.03)
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".eps")
        _savefig_silent(fig, fout_name)
        fig
    end
end

#
# Plot the per-element DSGS viscosity as a piecewise-constant staircase
# against x.  μ_dsgs[1:nelem, 1:neqs] is filled by compute_dsgs_viscosity!
# (one column per equation).  By default we draw column ieq=2 (the
# momentum equation) since for 1D E-form Marras gives a single μ shared
# by every equation. Every node of an element gets that element's
# value so the staircase is rendered cleanly.
#
function plot_dsgs_1d(mesh::St_mesh, μ_dsgs::AbstractMatrix, t, OUTPUT_DIR::String, inputs;
                      iout = 1, varname = "μ_dsgs", ieq = min(2, size(μ_dsgs, 2)))

    nelem = mesh.nelem
    ngl   = mesh.ngl

    xs = Vector{Float64}(undef, nelem*ngl)
    ys = Vector{Float64}(undef, nelem*ngl)
    @inbounds for ie = 1:nelem
        for i = 1:ngl
            ip = mesh.connijk[ie, i, 1, 1]
            xs[(ie-1)*ngl + i] = mesh.coords[ip, 1]
            ys[(ie-1)*ngl + i] = μ_dsgs[ie, ieq]
        end
    end
    sort_idx = sortperm(xs)

    plt = Plots.plot(xs[sort_idx], ys[sort_idx];
                     line = (:red, 2),
                     marker = (:circle, 3, :red),
                     title = string(varname, " (DSGS)  t = ", round(t, digits=4)),
                     xlabel = "x",
                     ylabel = varname,
                     titlefontsize = 18,
                     guidefontsize = 14,
                     legendfontsize = 12,
                     tickfontsize = 12,
                     legend = false,
                     show = false,
                     size = (600, 400))

    fout_name = string(OUTPUT_DIR, "/mu_dsgs-it", iout, ".png")
    _savefig_silent(plt, fout_name)
    plt
end

function plot_1d_grid(mesh::St_mesh)

    plt = Plots.plot() #Clear plot
    for i=1:mesh.npoin
        display(Plots.scatter(mesh.x[1:mesh.npoin], zeros(mesh.npoin),
                             markersize = 4,
                             color = :blue,
                             legend = false))
    end
end


function plot_initial(SD::NSD_2D, x::Array, q::Array, ivar, OUTPUT_DIR::String)
    nothing
end

#
# Nearest-neighbour rasterization of scattered nodal data onto a regular
# grid, used by plot_triangulation to draw filled contours. The nodes are
# binned into coarse cells once; every pixel then only searches its own
# and the surrounding cells. Unlike a global spline fit (cf. plot_surf3d)
# this cannot overshoot at solution kinks such as a shallow water wet/dry
# front.
#
function _grid_nearest(x, y, v, xg, yg)

    npoin = length(x)
    xmin, xmax = first(xg), last(xg)
    ymin, ymax = first(yg), last(yg)

    nbx = max(1, floor(Int, sqrt(npoin/2)))
    nby = nbx
    fx  = nbx/(xmax - xmin + eps(xmax - xmin))
    fy  = nby/(ymax - ymin + eps(ymax - ymin))

    bins = [Int[] for _ in 1:nbx, _ in 1:nby]
    for ip in 1:npoin
        i = clamp(1 + floor(Int, (x[ip] - xmin)*fx), 1, nbx)
        j = clamp(1 + floor(Int, (y[ip] - ymin)*fy), 1, nby)
        push!(bins[i,j], ip)
    end

    z = Matrix{Float64}(undef, length(xg), length(yg))
    for (jj, yp) in enumerate(yg), (ii, xp) in enumerate(xg)
        i0 = clamp(1 + floor(Int, (xp - xmin)*fx), 1, nbx)
        j0 = clamp(1 + floor(Int, (yp - ymin)*fy), 1, nby)
        best  = 0
        bestd = Inf
        ring  = 1
        while best == 0 && ring <= max(nbx, nby)
            for j in max(1,j0-ring):min(nby,j0+ring), i in max(1,i0-ring):min(nbx,i0+ring)
                for ip in bins[i,j]
                    d = (x[ip] - xp)^2 + (y[ip] - yp)^2
                    if d < bestd
                        bestd = d
                        best  = ip
                    end
                end
            end
            ring += 1
        end
        z[ii,jj] = best == 0 ? NaN : v[best]
    end
    return z
end

function plot_triangulation(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, inputs; iout=1, nvar=1, varnames=nothing)

    """
        Plot arbitrarily gridded unstructured 2D nodal data as filled
        contours. By default (:plot_matrix => true) all variables of one
        output time are rendered as a single plot-matrix figure that
        updates the interactive window (gksqt) in place and is written
        to fields-it<iout>.png. With :plot_matrix => false one silent
        PNG per variable is written instead (<var>-it<iout>.png) and no
        window is opened (see render_plot_matrix for why these two modes
        are mutually exclusive). On multi-rank runs each rank writes its
        own piece with a -rankN suffix.

        Inputs honoured: :plot_matrix, :plot_colormap (default :balance —
        a desaturated diverging map that brings the waves out),
        :plot_vlines/:plot_hlines.
    """

    comm    = get_mpi_comm()
    rank    = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    piece   = mpisize > 1 ? string("-rank", rank) : ""

    npoin = mesh.npoin
    xn = view(mesh.x, 1:npoin)
    yn = view(mesh.y, 1:npoin)
    xmin, xmax = extrema(xn)
    ymin, ymax = extrema(yn)
    Lx = xmax - xmin
    Ly = ymax - ymin

    # raster resolution proportional to the domain aspect ratio
    nmax = 400
    if Lx >= Ly
        nxi = nmax
        nyi = max(64, round(Int, nmax*Ly/Lx))
    else
        nyi = nmax
        nxi = max(64, round(Int, nmax*Lx/Ly))
    end
    xg = LinRange(xmin, xmax, nxi)
    yg = LinRange(ymin, ymax, nyi)

    # figure size that matches the domain aspect (axes flush with the
    # data, no dead white space); extra width for the colorbar
    hfig = 500
    wfig = clamp(round(Int, hfig*Lx/Ly) + 150, 350, 1300)

    cmap = Plots.cgrad(Symbol(get(inputs, :plot_colormap, :balance)))

    lmatrix = get(inputs, :plot_matrix, true)

    plts = []
    for ivar=1:nvar
        idx  = (ivar - 1)*npoin
        var  = (varnames === nothing || length(varnames) < ivar) ?
                   string("ivar", ivar) : string(varnames[ivar])
        fout_name = string(OUTPUT_DIR, "/", var, piece, "-it", iout, ".png")

        qvar = @view q[idx+1:idx+npoin]

        # Color range: explicit, and padded when the field is uniform (a
        # degenerate range breaks the colorbar).
        minq = minimum(qvar)
        maxq = maximum(qvar)
        clims = maxq > minq ? (minq, maxq) : (minq - 0.5, maxq + 0.5)

        zg = _grid_nearest(xn, yn, qvar, xg, yg)

        # Filled contours (no contour lines) of the rasterized field
        plt = Plots.contourf(xg, yg, zg';
                            color = cmap,
                            clims = clims,
                            levels = 30,
                            linewidth = 0,
                            colorbar = true,
                            legend = false,
                            aspect_ratio = :equal,
                            xlims = (xmin, xmax),
                            ylims = (ymin, ymax),
                            framestyle = :box,
                            xlabel = "x",
                            ylabel = "y",
                            title = string(var, "  ", title),
                            show = false,
                            size = (wfig, hfig))

        vlines = inputs[:plot_vlines]
        hlines = inputs[:plot_hlines]
        if !(vlines == "empty")
            for i=1:size(vlines,1)
                Plots.vline!(plt, [vlines[i]]; color = :red, linestyle = :dash, label = "")
            end
        end
        if !(hlines == "empty")
            for i=1:size(hlines,1)
                Plots.hline!(plt, [hlines[i]]; color = :red, linestyle = :dash, label = "")
            end
        end

        if !lmatrix
            _savefig_silent(plt, fout_name)
        end
        push!(plts, plt)
    end

    render_plot_matrix(lmatrix, plts, OUTPUT_DIR, iout; wfig=wfig, hfig=hfig)
end
function plot_triangulation(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, inputs; nvar=1) nothing end
function plot_triangulation(SD::NSD_3D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, inputs; nvar=1) nothing end

function plot_surf3d(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1, smoothing_factor=1e-3, varnames=nothing)

    xmin = minimum(mesh.x); xmax = maximum(mesh.x);
    ymin = minimum(mesh.y); ymax = maximum(mesh.y);

    comm    = get_mpi_comm()
    rank    = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    piece   = mpisize > 1 ? string("-rank", rank) : ""

    nxi = 500
    nyi = 500
    npoin = mesh.npoin
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        var  = (varnames === nothing || length(varnames) < ivar) ?
                   string("ivar", ivar) : string(varnames[ivar])
        fout_name = string(OUTPUT_DIR, "/", var, piece, "-it", iout, ".png")

        #Spline2d
        spl = Spline2D(mesh.x[1:npoin], mesh.y[1:npoin], q[idx+1:idx+npoin]; kx=4, ky=4, s=smoothing_factor)
        xg = LinRange(xmin, xmax, nxi); yg = LinRange(ymin, ymax, nyi);
        zspl = evalgrid(spl, xg, yg);
        #End spline2d

        #figure:
        plt = Plots.surface(xg, yg, zspl;
                           color = :viridis,
                           camera = (0, 90),  # Top-down view similar to Makie's azimuth/elevation
                           colorbar = true,
                           legend = false,
                           title = string(var, "  ", title),
                           show = false,
                           size = (1200, 400))

        _savefig_silent(plt, fout_name)
        plt
    end

end
