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

function plot_results(SD::NSD_1D, mesh::St_mesh, q, title::String, OUTPUT_DIR::String, outvar, inputs; iout=1, nvar=1, PT=nothing, μ_nodes=nothing, t=nothing)

    epsi = 1.1
    npoin = mesh.npoin

    lmatrix = get(inputs, :plot_matrix, true)

    qout = reshape(q, npoin, nvar)   # 2D view (3x4) - NO allocation
    x_coords = mesh.coords[1, 1:npoin]
    sort_idx = sortperm(x_coords)
    plts = []

    # Optional per-case reference/analytic solution.
    #
    # This 1D plotter is shared by every 1D case, and most have no closed
    # form to compare against, so nothing is hardcoded here. Instead a case
    # may define
    #
    #     user_analytic_solution(x, t, outvar, inputs) -> npoin × nvar matrix
    #
    # in its own problem directory (see problems/CompEuler/sod1d). Case files
    # are included into this module, so the method simply exists for the case
    # being run and is absent for every other one. Entries the case cannot
    # supply should be NaN — Plots skips them, so a case can provide an exact
    # solution for some variables and not others. Failures are reported once
    # and never abort the run: a broken reference solution must not take the
    # simulation output with it.
    qref = nothing
    if t !== nothing && get(inputs, :_has_analytic, false) &&
        isdefined(@__MODULE__, :user_analytic_solution)
        try
            qref = user_analytic_solution(x_coords, t, outvar, inputs)
            if qref !== nothing && size(qref) != (npoin, nvar)
                @warn "user_analytic_solution returned $(size(qref)), expected $((npoin, nvar)); ignoring."
                qref = nothing
            end
        catch err
            @warn "user_analytic_solution failed; plotting numerical solution only." exception=err
            qref = nothing
        end
    end
    # Optional per-case figure. A case may ship a user_plot.jl defining
    #
    #     user_plot_1d(x, q, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)
    #
    # (x: node coordinates, q: npoin × nvar output variables, qref: the
    # reference of user_analytic_solution or nothing, μ_nodes: the DynSGS
    # coefficient at the nodes or nothing) that renders and saves its own
    # figure — e.g. the single density plot of a paper — in place of the
    # generic panels below (see problems/MHD/brioWu1d). :plot_user => false
    # in the inputs falls back to the generic panels.
    if get(inputs, :_has_user_plot, false) && get(inputs, :plot_user, true) &&
        isdefined(@__MODULE__, :user_plot_1d)
        try
            user_plot_1d(x_coords, qout, qref, μ_nodes, t, outvar, inputs, OUTPUT_DIR, iout)
            return nothing
        catch err
            @warn "user_plot_1d failed; falling back to the generic panels." exception=err
        end
    end

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
                        label = "Jexpresso",
                        show = false,
                        size = (600, 400))

        if qref !== nothing && any(isfinite, @view(qref[:, ivar]))
            Plots.plot!(plt, x_coords[sort_idx], qref[sort_idx, ivar];
                        line = (:black, 2, :dash),
                        marker = :none,
                        label = "exact",
                        legend = :best)
        end

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
function render_plot_matrix(lmatrix, plts, OUTPUT_DIR, iout; wfig=600, hfig=400, piece=nothing)
    lmatrix || return nothing
    nplt = length(plts)
    nplt == 0 && return nothing
    if piece === nothing
        comm    = get_mpi_comm()
        mpisize = MPI.Comm_size(comm)
        piece   = mpisize > 1 ? string("-rank", MPI.Comm_rank(comm)) : ""
    end
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
            xs[(ie-1)*ngl + i] = mesh.coords[1, ip]
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

#
# Concatenate what MPI.gather hands back on the root (a vector of the
# per-rank arrays, or already one flat array depending on the MPI.jl
# version) into one flat vector.
#
_flat_gathered(g::AbstractVector{<:AbstractArray}) = reduce(vcat, g)
_flat_gathered(g) = g

#
# In-plane vector potential A on a raster grid from the rasterized field
# components (Bx, By):  Bx = ∂A/∂y, By = -∂A/∂x, i.e.
#
#   A(x, y) = ∫₀ʸ Bx(x, y') dy'  -  ∫₀ˣ By(x', y₀) dx'
#
# by cumulative trapezoids. Its isocontours are the magnetic field lines
# (exact for a divergence-free field, and a faithful rendering otherwise).
#
function _vector_potential(Bxr, Byr, xg, yg)
    nx, ny = size(Bxr)
    A = zeros(nx, ny)
    # bottom row: -∫ By dx
    for i = 2:nx
        A[i,1] = A[i-1,1] - 0.5*(Byr[i,1] + Byr[i-1,1])*(xg[i] - xg[i-1])
    end
    # columns: +∫ Bx dy
    for i = 1:nx, j = 2:ny
        A[i,j] = A[i,j-1] + 0.5*(Bxr[i,j] + Bxr[i,j-1])*(yg[j] - yg[j-1])
    end
    return A
end

#
# Isolines of the raster field z(xg, yg) at the given levels by marching
# squares, returned as one polyline pair with NaN separators so that a
# single Plots.plot! draws them all. Used instead of Plots.contour! for
# overlays: under GR a contour series added on top of a filled-contour
# panel takes the PANEL's color limits as the range of its level set, so
# the requested levels of a field with a different range are never drawn.
#
function _isolines(xg, yg, z, levels)
    nx, ny = size(z)
    xs = Float64[]; ys = Float64[]
    interp(p1, p2, v1, v2, lev) = p1 + (lev - v1)/(v2 - v1)*(p2 - p1)
    pts = NTuple{2,Float64}[]
    for lev in levels
        for j = 1:ny-1, i = 1:nx-1
            v1 = z[i,j]; v2 = z[i+1,j]; v3 = z[i+1,j+1]; v4 = z[i,j+1]
            (isfinite(v1) && isfinite(v2) && isfinite(v3) && isfinite(v4)) || continue
            empty!(pts)
            # edge crossings: bottom (1-2), right (2-3), top (4-3), left (1-4)
            if (v1 < lev) != (v2 < lev); push!(pts, (interp(xg[i], xg[i+1], v1, v2, lev), yg[j]));   end
            if (v2 < lev) != (v3 < lev); push!(pts, (xg[i+1], interp(yg[j], yg[j+1], v2, v3, lev))); end
            if (v4 < lev) != (v3 < lev); push!(pts, (interp(xg[i], xg[i+1], v4, v3, lev), yg[j+1])); end
            if (v1 < lev) != (v4 < lev); push!(pts, (xg[i], interp(yg[j], yg[j+1], v1, v4, lev)));   end
            np = length(pts)
            if np == 2
                push!(xs, pts[1][1], pts[2][1], NaN); push!(ys, pts[1][2], pts[2][2], NaN)
            elseif np == 4
                # saddle cell: pair the crossings by the cell-center value
                vc = 0.25*(v1 + v2 + v3 + v4)
                if (vc < lev) == (v1 < lev)
                    push!(xs, pts[1][1], pts[2][1], NaN, pts[3][1], pts[4][1], NaN)
                    push!(ys, pts[1][2], pts[2][2], NaN, pts[3][2], pts[4][2], NaN)
                else
                    push!(xs, pts[1][1], pts[4][1], NaN, pts[2][1], pts[3][1], NaN)
                    push!(ys, pts[1][2], pts[4][2], NaN, pts[2][2], pts[3][2], NaN)
                end
            end
        end
    end
    return xs, ys
end

function plot_triangulation(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, inputs; iout=1, nvar=1, varnames=nothing, μ_nodes=nothing, μ_names=nothing)

    """
        Plot arbitrarily gridded unstructured 2D nodal data as filled
        contours. By default (:plot_matrix => true) all variables of one
        output time are rendered as a single plot-matrix figure that
        updates the interactive window (gksqt) in place and is written
        to fields-it<iout>.png. With :plot_matrix => false one silent
        PNG per variable is written instead (<var>-it<iout>.png) and no
        window is opened (see render_plot_matrix for why these two modes
        are mutually exclusive). Under MPI the nodal data of all ranks is
        gathered on rank 0, which renders the whole domain; the other
        ranks return at once.

        Inputs honoured (all optional, see mod_inputs.jl for the defaults):
          :plot_matrix, :plot_colormap (default :balance — a desaturated
          diverging map that brings the waves out), :plot_vlines/:plot_hlines,
          :plot_xlabel/:plot_ylabel, :plot_raster_nmax,
          :plot_vars             names of the variables to render (default: all)
          :plot_log10            names rendered as log10(var)
          :plot_clims            Dict(name => (lo, hi)) fixed color range, data clamped to it
          :plot_fieldlines       (Bx_name, By_name): overlay the isocontours of the
                                 vector potential of that in-plane field (black lines,
                                 :plot_fieldlines_levels of them)
          :plot_vectors          (u_name, v_name): overlay a velocity-vector field
                                 (white arrows, :plot_vectors_n = (nx, ny) arrows,
                                 reference arrow of speed :plot_vectors_ref)
          :plot_overlay_on       names of the panels that get the overlays (default: all)
          :plot_user             1D: use the case's user_plot_1d (user_plot.jl) figure when
                                 the case ships one (default: true); false gives the generic panels
          :plot_dsgs             render the μ_dsgs panels of a DynSGS run (default: true)
          :plot_dsgs_vars        names of the damped variables whose μ_dsgs panel is
                                 written (default: all slots with a non-zero coefficient)
          :plot_dsgs_log10       render log₁₀ μ_dsgs floored at :plot_dsgs_floor
                                 (default: false; floor 1e-6), file log10_μ_dsgs_<var>-it<n>.png
          :plot_profile_x        x at which a vertical profile figure profile-it<iout>.png
                                 of :plot_profile_vars is written (nodes on that line,
                                 or the nearest raster column); :plot_profile_log10,
                                 :plot_profile_ylims (Dict(name => (lo, hi))) and
                                 :plot_profile_vlines (heights marked) style it.
    """

    comm    = get_mpi_comm()
    rank    = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)

    npoin = mesh.npoin
    names = [(varnames === nothing || length(varnames) < ivar) ?
                 string("ivar", ivar) : string(varnames[ivar]) for ivar = 1:nvar]
    nμ    = μ_nodes === nothing ? 0 : size(μ_nodes, 2)
    μnames = [(μ_names === nothing || length(μ_names) < ieq) ?
                  string("μ_dsgs_", ieq) : string("μ_dsgs_", μ_names[ieq]) for ieq = 1:nμ]

    #
    # Nodal data, gathered on rank 0 under MPI. Points shared by two
    # partitions appear twice, which the nearest-neighbour raster below does
    # not mind.
    #
    xn = collect(view(mesh.x, 1:npoin))
    yn = collect(view(mesh.y, 1:npoin))
    qv = [collect(view(q, (ivar - 1)*npoin + 1:ivar*npoin)) for ivar = 1:nvar]
    μv = [collect(view(μ_nodes, 1:npoin, ieq)) for ieq = 1:nμ]
    if mpisize > 1
        xg_ = MPI.gather(xn, comm)
        yg_ = MPI.gather(yn, comm)
        qg_ = [MPI.gather(qv[ivar], comm) for ivar = 1:nvar]
        μg_ = [MPI.gather(μv[ieq], comm) for ieq = 1:nμ]
        rank == 0 || return nothing
        xn = _flat_gathered(xg_)
        yn = _flat_gathered(yg_)
        qv = [_flat_gathered(qg_[ivar]) for ivar = 1:nvar]
        μv = [_flat_gathered(μg_[ieq]) for ieq = 1:nμ]
    end
    npts = length(xn)

    xmin, xmax = extrema(xn)
    ymin, ymax = extrema(yn)
    Lx = xmax - xmin
    Ly = ymax - ymin

    # raster resolution proportional to the domain aspect ratio
    nmax = get(inputs, :plot_raster_nmax, 400)
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

    lmatrix    = get(inputs, :plot_matrix, true)
    xlab       = string(get(inputs, :plot_xlabel, "x"))
    ylab       = string(get(inputs, :plot_ylabel, "y"))
    plot_vars  = get(inputs, :plot_vars, nothing)
    logvars    = get(inputs, :plot_log10, String[])
    clims_d    = get(inputs, :plot_clims, Dict{String,Any}())
    fieldlines = get(inputs, :plot_fieldlines, nothing)
    vectors    = get(inputs, :plot_vectors, nothing)
    overlay_on = get(inputs, :plot_overlay_on, nothing)
    vlines     = get(inputs, :plot_vlines, "empty")
    hlines     = get(inputs, :plot_hlines, "empty")

    findvar(name) = findfirst(==(string(name)), names)

    # Rasterize on demand, once per variable
    raster = Dict{Int, Matrix{Float64}}()
    getraster(ivar) = get!(raster, ivar) do
        _grid_nearest(xn, yn, qv[ivar], xg, yg)
    end

    #
    # Overlays: magnetic field lines (isocontours of the vector potential)
    # and a velocity-vector field, both built on the raster.
    #
    Araster = nothing
    if fieldlines !== nothing
        ib = findvar(fieldlines[1]); jb = findvar(fieldlines[2])
        if ib === nothing || jb === nothing
            @warn " plot_triangulation: :plot_fieldlines => $(fieldlines) names a variable that is not in the output set $(names); no field lines drawn."
        else
            Araster = _vector_potential(getraster(ib), getraster(jb), xg, yg)
        end
    end
    quiv = nothing
    if vectors !== nothing
        iu = findvar(vectors[1]); iv = findvar(vectors[2])
        if iu === nothing || iv === nothing
            @warn " plot_triangulation: :plot_vectors => $(vectors) names a variable that is not in the output set $(names); no vectors drawn."
        else
            Ur = getraster(iu); Vr = getraster(iv)
            nqx, nqy = get(inputs, :plot_vectors_n, (30, 13))
            ix  = unique(round.(Int, range(1, nxi, length=nqx + 2)[2:end-1]))
            jy  = unique(round.(Int, range(1, nyi, length=nqy + 2)[2:end-1]))
            ref = get(inputs, :plot_vectors_ref, nothing)
            if ref === nothing
                ref = maximum(sqrt.(Ur.^2 .+ Vr.^2))
                ref = ref > 0 ? ref : 1.0
            end
            Lref  = 0.09*Lx                 # drawn length of the reference speed
            scale = Lref/ref
            xs = Float64[]; ys = Float64[]; us = Float64[]; vs = Float64[]
            vmin = 0.01*ref          # vectors below 1% of the reference speed are not drawn
            for j in jy, i in ix
                sqrt(Ur[i,j]^2 + Vr[i,j]^2) >= vmin || continue
                push!(xs, xg[i]); push!(ys, yg[j])
                push!(us, scale*Ur[i,j]); push!(vs, scale*Vr[i,j])
            end
            quiv = (xs, ys, us, vs, Lref, ref)
        end
    end

    flines = nothing
    if Araster !== nothing
        amin, amax = extrema(Araster)
        nlev = get(inputs, :plot_fieldlines_levels, 40)
        if amax > amin
            lev    = collect(range(amin, amax, length=nlev + 2)[2:end-1])
            flines = _isolines(xg, yg, Araster, lev)
        end
    end

    # Arrow heads: the small open head, explicitly — GR's default closed
    # head is drawn at full size even for a zero-length vector.
    arrowhead = Plots.arrow(:simple, :head, 0.1, 0.1)

    function _overlay!(plt)
        if flines !== nothing
            Plots.plot!(plt, flines[1], flines[2]; color = :black, linewidth = 0.8, label = "")
        end
        if quiv !== nothing
            xs, ys, us, vs, Lref, ref = quiv
            Plots.quiver!(plt, xs, ys; quiver = (us, vs), color = :white, linewidth = 0.6, arrow = arrowhead)
            # reference arrow, bottom-left corner
            x0 = xmin + 0.02*Lx
            y0 = ymin + 0.05*Ly
            Plots.quiver!(plt, [x0], [y0]; quiver = ([Lref], [0.0]), color = :white, linewidth = 1.2, arrow = arrowhead)
            Plots.annotate!(plt, x0 + Lref + 0.01*Lx, y0, Plots.text(string("= ", ref), 8, :white, :left))
        end
        return plt
    end

    function _add_lines!(plt)
        if !(vlines == "empty")
            for i = 1:size(vlines, 1)
                Plots.vline!(plt, [vlines[i]]; color = :red, linestyle = :dash, label = "")
            end
        end
        if !(hlines == "empty")
            for i = 1:size(hlines, 1)
                Plots.hline!(plt, [hlines[i]]; color = :red, linestyle = :dash, label = "")
            end
        end
        return plt
    end

    plts = []
    for ivar = 1:nvar
        var = names[ivar]
        (plot_vars === nothing || var in plot_vars) || continue
        fout_name = string(OUTPUT_DIR, "/", var, "-it", iout, ".png")

        zg    = copy(getraster(ivar))
        label = var
        if var in logvars
            zg    = log10.(max.(zg, 1e-300))
            label = string("log10(", var, ")")
        end

        # Color range: fixed by the user (data clamped to it so that the
        # out-of-range values take the end colors, as a colorbar with
        # "extend" would show them), otherwise the data range, padded when
        # the field is uniform (a degenerate range breaks the colorbar).
        if haskey(clims_d, var)
            clims = (Float64(clims_d[var][1]), Float64(clims_d[var][2]))
            zg    = clamp.(zg, clims[1], clims[2])
        else
            finite = filter(isfinite, zg)
            minq = isempty(finite) ? 0.0 : minimum(finite)
            maxq = isempty(finite) ? 0.0 : maximum(finite)
            clims = maxq > minq ? (minq, maxq) : (minq - 0.5, maxq + 0.5)
        end

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
                            xlabel = xlab,
                            ylabel = ylab,
                            title = string(label, "  ", title),
                            show = false,
                            size = (wfig, hfig),
                            bottom_margin = 5Plots.mm, left_margin = 3Plots.mm)
        _add_lines!(plt)
        if overlay_on === nothing || var in overlay_on
            _overlay!(plt)
        end

        if !lmatrix
            _savefig_silent(plt, fout_name)
        end
        push!(plts, plt)
    end

    # DSGS runs: one extra filled-contour panel per equation showing the
    # eddy viscosity actually applied, mirroring what the 1D plotter does
    # with its μ_dsgs staircase. Slots that are identically zero (mass
    # everywhere; ρw and Bz in a 2D MHD run) carry no information as a
    # contour plot and are skipped rather than drawn as a blank map.
    #
    # NOTE the slots are not in a common unit — momentum/energy carry the
    # dynamic ρ̄μ, magnetic and ψ slots the kinematic μ — so each panel is
    # scaled to its own range and should be read against itself over time,
    # not against a neighbouring panel.
    #
    # :plot_dsgs_vars selects the slots by the name of the variable they damp
    # (a DynSGS-MHD run in its conserved form gives every slot the same
    # kinematic coefficient, so one panel says it all); :plot_dsgs_log10
    # renders log₁₀ μ floored at :plot_dsgs_floor, since a residual-based
    # coefficient spans several decades between the quiet flow and the cap.
    μsel   = get(inputs, :plot_dsgs_vars, nothing)
    μlog   = get(inputs, :plot_dsgs_log10, false)
    μfloor = get(inputs, :plot_dsgs_floor, 1.0e-6)
    if nμ > 0 && get(inputs, :plot_dsgs, true)
        for ieq = 1:nμ
            if μsel !== nothing && !(μ_names !== nothing && length(μ_names) >= ieq && string(μ_names[ieq]) in μsel)
                continue
            end
            μvar = μv[ieq]
            μmax = maximum(μvar)
            μmax > 0 || continue

            name = μnames[ieq]
            if μlog
                μvar = log10.(max.(μvar, μfloor))
                μmax = maximum(μvar)
                name = string("log10_", name)
            end
            zgμ  = _grid_nearest(xn, yn, μvar, xg, yg)
            pltμ = Plots.contourf(xg, yg, zgμ';
                                  color = cmap,
                                  clims = (minimum(μvar), μmax),
                                  levels = 30,
                                  linewidth = 0,
                                  colorbar = true,
                                  legend = false,
                                  aspect_ratio = :equal,
                                  xlims = (xmin, xmax),
                                  ylims = (ymin, ymax),
                                  framestyle = :box,
                                  xlabel = xlab,
                                  ylabel = ylab,
                                  title = string(name, "  ", title),
                                  show = false,
                                  size = (wfig, hfig),
                                  bottom_margin = 5Plots.mm, left_margin = 3Plots.mm)
            if !lmatrix
                _savefig_silent(pltμ, string(OUTPUT_DIR, "/", name, "-it", iout, ".png"))
            end
            push!(plts, pltμ)
        end
    end

    #
    # Vertical profiles at x = :plot_profile_x (e.g. the centerline of a
    # rising loop): the nodes lying on that line if the mesh has any,
    # otherwise the nearest raster column.
    #
    xp = get(inputs, :plot_profile_x, nothing)
    if xp !== nothing
        pvars = get(inputs, :plot_profile_vars, nothing)
        pvars = pvars === nothing ? names : pvars
        plog  = get(inputs, :plot_profile_log10, String[])
        pyl   = get(inputs, :plot_profile_ylims, Dict{String,Any}())
        pvl   = get(inputs, :plot_profile_vlines, Float64[])

        sel   = findall(ip -> abs(xn[ip] - xp) <= 1e-6*Lx, 1:npts)
        use_nodes = length(sel) >= 3
        if use_nodes
            order = sortperm(yn[sel])
            zs    = yn[sel][order]
        else
            ii = argmin(abs.(xg .- xp))
            zs = collect(yg)
        end

        pl = []
        for var in pvars
            ivp = findvar(var)
            ivp === nothing && continue
            vals = use_nodes ? qv[ivp][sel][order] : getraster(ivp)[ii, :]
            lab  = string(var)
            if var in plog
                vals = log10.(max.(vals, 1e-300))
                lab  = string("log10(", var, ")")
            end
            p = Plots.scatter(zs, vals;
                              markersize = 2.0, markerstrokewidth = 0, color = :steelblue,
                              xlabel = ylab, ylabel = lab,
                              xlims = (ymin, ymax),
                              title = string(lab, " at ", xlab, " = ", xp, "  ", title),
                              titlefontsize = 9,
                              legend = false, framestyle = :box, show = false)
            if haskey(pyl, var)
                Plots.ylims!(p, (Float64(pyl[var][1]), Float64(pyl[var][2])))
            end
            for v in pvl
                Plots.vline!(p, [v]; color = :black, linestyle = :dashdot, label = "")
            end
            push!(pl, p)
        end
        if !isempty(pl)
            np  = length(pl)
            fig = Plots.plot(pl...; layout = (1, np), size = (420*np, 360),
                             left_margin = 6Plots.mm, bottom_margin = 6Plots.mm, show = false)
            _savefig_silent(fig, string(OUTPUT_DIR, "/profile-it", iout, ".png"))
        end
    end

    render_plot_matrix(lmatrix, plts, OUTPUT_DIR, iout; wfig=wfig, hfig=hfig, piece="")

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
