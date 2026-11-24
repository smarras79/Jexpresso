#!/usr/bin/env julia
#---------------------------------------------------------------------------------
# Figures for problems/MHD/orszagTangBormanis2024 (the Orszag-Tang vortex).
#
# A normal run of that case writes VTK only (`:outformat => "vtk"`), so the
# figures that appear in README.md are NOT produced by the solver — this script
# produces them, off-line, from the `.pvtu` files a finished run leaves in
# `output/MHD/orszagTangBormanis2024/<run>/`.
#
# It renders three 2x2 panels, each at four output times (0.7, 0.8, 0.9, 1.0 by
# default), on one shared color scale per figure so the panels are comparable:
#
#   assets/MHD_OT_rho.png        density, on the [0.1, 0.4] scale of the "True"
#                                row of Fig. 3 of Bormanis et al. (2024)
#   assets/MHD_OT_By.png         the y magnetic field B_y, on a symmetric
#                                diverging scale, for comparison against the
#                                same paper's magnetic-field figures
#   assets/MHD_OT_mu_dsgs_By.png the DynSGS eddy viscosity applied to the B_y
#                                equation, with B_y isolines drawn on top: the
#                                viscosity should sit where the isolines bunch,
#                                i.e. on the large gradients of B_y
#
# The viscosity field is only present in the VTK output when the run used
# `:visc_model => DSGS_MHD()` (the shipped setting); `write_vtk` writes it as
# point data named `mu_dsgs_<var>`, one field per equation.
#
# Run (from the repository root, after the case has been run):
#
#     julia --project=. tools/plot_orszag_tang.jl
#     julia --project=. tools/plot_orszag_tang.jl --dir output/MHD/orszagTangBormanis2024/run-1
#     julia --project=. tools/plot_orszag_tang.jl --times 0.4,0.6,0.8,1.0 --out /tmp/figs
#
# Options:
#   --dir DIR      run directory holding iter_*.pvtu. Default: the most recently
#                  modified subdirectory of output/MHD/orszagTangBormanis2024/.
#   --out DIR      where to write the PNGs. Default: assets/
#   --times a,b,.. output times to draw. Default: 0.7,0.8,0.9,1.0
#   --t0 T         time of the FIRST snapshot, for runs whose VTK files predate
#                  the TimeValue field data. Default: 0.0
#   --dtout DT     interval between snapshots, same fallback. Default: 0.05
#                  (`:diagnostics_at_times => 0.0:0.05:1.0` in user_inputs.jl)
#   --rho-clims a,b  density color limits. Default: 0.1,0.4
#---------------------------------------------------------------------------------

using Printf
using Plots

include(joinpath(@__DIR__, "vtu_reader.jl"))

const CASE_OUTPUT = joinpath("output", "MHD", "orszagTangBormanis2024")

#---------------------------------------------------------------------------------
# Command line
#---------------------------------------------------------------------------------

function parse_args(argv)
    opts = Dict{String,String}(
        "dir"       => "",
        "out"       => "assets",
        "times"     => "0.7,0.8,0.9,1.0",
        "t0"        => "0.0",
        "dtout"     => "0.05",
        "rho-clims" => "0.1,0.4",
    )
    i = 1
    while i <= length(argv)
        a = argv[i]
        startswith(a, "--") || error("unrecognized argument \"$a\" (see the header of this file)")
        key = a[3:end]
        haskey(opts, key) || error("unknown option \"--$key\" (see the header of this file)")
        i + 1 <= length(argv) || error("option --$key needs a value")
        opts[key] = argv[i+1]
        i += 2
    end
    return opts
end

_floats(s) = [parse(Float64, strip(t)) for t in split(s, ",")]

#---------------------------------------------------------------------------------
# Locating the snapshots
#
# write_vtk names its files iter_<n>, n = 1, 2, ... in output order, and the
# simulation time is recorded as TimeValue field data. Files written before
# that was added carry no time at all, hence the --t0/--dtout fallback: with
# `:lwrite_initial => true` and `:diagnostics_at_times => 0.0:0.05:1.0`,
# snapshot n is at t = t0 + (n-1)*dtout.
#---------------------------------------------------------------------------------

function find_run_dir(explicit::AbstractString)
    if !isempty(explicit)
        isdir(explicit) || error("no such directory: $explicit")
        return explicit
    end
    isdir(CASE_OUTPUT) ||
        error("$CASE_OUTPUT does not exist — run the case first:\n" *
              "    julia --project=. src/Jexpresso.jl MHD orszagTangBormanis2024")

    cands = String[]
    for e in readdir(CASE_OUTPUT; join=true)
        isdir(e) || continue
        isempty(filter(f -> endswith(f, ".pvtu") || endswith(f, ".vtu"), readdir(e))) && continue
        push!(cands, e)
    end
    # The case itself may be the run directory when :loverwrite_output is on.
    isempty(cands) && !isempty(filter(f -> endswith(f, ".pvtu"), readdir(CASE_OUTPUT))) &&
        return CASE_OUTPUT
    isempty(cands) && error("found no VTK output under $CASE_OUTPUT")

    return cands[argmax(mtime.(cands))]
end

"All snapshots in `dir`, as (index, path, time), sorted by index."
function list_snapshots(dir::AbstractString, t0::Float64, dtout::Float64)
    files = filter(f -> occursin(r"^iter_\d+\.p?vtu$", f), readdir(dir))
    isempty(files) && error("no iter_*.pvtu files in $dir")

    # Prefer the .pvtu master when both it and a bare .vtu are present.
    byidx = Dict{Int,String}()
    for f in files
        n = parse(Int, match(r"^iter_(\d+)\.", f).captures[1])
        if !haskey(byidx, n) || endswith(f, ".pvtu")
            byidx[n] = f
        end
    end

    out = NamedTuple{(:index, :path, :time),Tuple{Int,String,Float64}}[]
    for n in sort(collect(keys(byidx)))
        path = joinpath(dir, byidx[n])
        t = try
            g = read_vtk_grid(path)
            g.time
        catch err
            @warn "could not read $path" err
            nothing
        end
        push!(out, (index = n, path = path,
                    time = t === nothing ? t0 + (n - 1)*dtout : t))
    end
    return out
end

function pick(snaps, t::Float64)
    i = argmin(abs.(getfield.(snaps, :time) .- t))
    s = snaps[i]
    abs(s.time - t) > 1.0e-6 &&
        @info @sprintf("requested t = %.4f, using the nearest snapshot %s at t = %.4f",
                       t, basename(s.path), s.time)
    return s
end

#---------------------------------------------------------------------------------
# Rasterization
#
# Nearest-neighbour, same choice as src/io/plotting/jeplots.jl: the LGL nodes
# are non-uniformly spaced and the solution has shocks, so a global spline fit
# would overshoot at exactly the features these figures are about.
#---------------------------------------------------------------------------------

function raster(x, y, v, xg, yg)
    npoin      = length(x)
    xmin, xmax = first(xg), last(xg)
    ymin, ymax = first(yg), last(yg)

    nbx = max(1, floor(Int, sqrt(npoin/2)))
    nby = nbx
    fx  = nbx/(xmax - xmin + eps(xmax - xmin))
    fy  = nby/(ymax - ymin + eps(ymax - ymin))

    bins = [Int[] for _ = 1:nbx, _ = 1:nby]
    for ip = 1:npoin
        i = clamp(1 + floor(Int, (x[ip] - xmin)*fx), 1, nbx)
        j = clamp(1 + floor(Int, (y[ip] - ymin)*fy), 1, nby)
        push!(bins[i,j], ip)
    end

    z = Matrix{Float64}(undef, length(xg), length(yg))
    for (jj, yp) in enumerate(yg), (ii, xp) in enumerate(xg)
        i0 = clamp(1 + floor(Int, (xp - xmin)*fx), 1, nbx)
        j0 = clamp(1 + floor(Int, (yp - ymin)*fy), 1, nby)
        best, bestd, ring = 0, Inf, 1
        while best == 0 && ring <= max(nbx, nby)
            for j in max(1,j0-ring):min(nby,j0+ring), i in max(1,i0-ring):min(nbx,i0+ring)
                for ip in bins[i,j]
                    d = (x[ip] - xp)^2 + (y[ip] - yp)^2
                    if d < bestd
                        bestd, best = d, ip
                    end
                end
            end
            ring += 1
        end
        z[ii,jj] = best == 0 ? NaN : v[best]
    end
    return z
end

#---------------------------------------------------------------------------------
# Panels
#---------------------------------------------------------------------------------

const NRASTER = 320

function grids(g)
    xmin, xmax = extrema(g.x)
    ymin, ymax = extrema(g.y)
    return LinRange(xmin, xmax, NRASTER), LinRange(ymin, ymax, NRASTER)
end

"Fetch a point-data field, with a clear error listing what the file does have."
function field(g, name::AbstractString)
    haskey(g.data, name) && return g.data[name]
    error("field \"$name\" is not in this output. Available: " *
          join(sort(collect(keys(g.data))), ", "))
end

function panel(g, name, t, clims, cmap; title = name, isolines = nothing,
               isolevels = 8, logscale = false)
    xg, yg = grids(g)
    z = raster(g.x, g.y, field(g, name), xg, yg)

    zplot, cl = if logscale
        floor_ = max(clims[1], 1.0e-14)
        (log10.(max.(z, floor_)), (log10(floor_), log10(max(clims[2], 10*floor_))))
    else
        (z, clims)
    end

    plt = Plots.heatmap(xg, yg, zplot';
                        color        = cmap,
                        clims        = cl,
                        aspect_ratio = :equal,
                        xlims        = (first(xg), last(xg)),
                        ylims        = (first(yg), last(yg)),
                        framestyle   = :box,
                        title        = @sprintf("%s   t = %.2f", title, t),
                        titlefontsize = 10,
                        tickfontsize  = 7,
                        colorbar      = true,
                        colorbar_tickfontsize = 7,
                        legend        = false,
                        show          = false)

    if isolines !== nothing
        ziso = raster(g.x, g.y, field(g, isolines), xg, yg)
        a    = maximum(abs, filter(isfinite, ziso))
        if a > 0
            # The overlaid field lives on a completely different scale from the
            # background (B_y is O(0.5), log10 μ is O(-4)), and `clims` is a
            # SUBPLOT attribute: GR clips the contour series against it too and
            # drops every line. Map the isoline field affinely onto the color
            # range instead. An affine map of z does not move its level sets, so
            # as long as the levels are mapped the same way the drawn curves are
            # exactly the isolines of the original field — and they are black, so
            # the rescaling is invisible.
            lo, hi = cl
            sc(v)  = lo + (v + a)/(2a)*(hi - lo)
            lv     = collect(LinRange(-a, a, isolevels + 2))[2:end-1]  # skip the
                                                                       # degenerate
                                                                       # end levels
            # `colorbar_entry = false`, not `colorbar = false`: the latter is a
            # subplot attribute and would drop the heatmap's colorbar with it.
            Plots.contour!(plt, xg, yg, sc.(ziso)';
                           levels         = sc.(lv),
                           linecolor      = :black,
                           linewidth      = 0.7,
                           alpha          = 0.8,
                           colorbar_entry = false,
                           legend         = false)
        end
    end
    return plt
end

#---------------------------------------------------------------------------------
# Does the viscosity follow the large gradients of B_y?
#
# The figure invites that reading, and the eye is easy to fool with a blocky
# log-scaled map, so the same question is answered numerically on every run.
#
# μ is piecewise constant per element, so both quantities are reduced to the
# element grid before they are compared: μ_e directly, and |∇B_y| as the maximum
# over the element of a central difference on the raster of unique nodes. A
# Spearman (rank) correlation is used rather than Pearson because both fields
# are heavy-tailed and only the ordering is meaningful.
#---------------------------------------------------------------------------------

_rank(v) = (p = sortperm(v); r = similar(v, Float64); r[p] = collect(1:length(v)); r)

function _spearman(a, b)
    ra, rb = _rank(collect(a)), _rank(collect(b))
    ma, mb = sum(ra)/length(ra), sum(rb)/length(rb)
    da, db = ra .- ma, rb .- mb
    return sum(da .* db)/sqrt(sum(da .^ 2)*sum(db .^ 2))
end

"Reduce nodal data to a per-element grid of `nel × nel` by binning on (x, y)."
function _element_max(x, y, v, nel)
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    fx = nel/(xmax - xmin + eps(xmax - xmin))
    fy = nel/(ymax - ymin + eps(ymax - ymin))
    o  = zeros(nel, nel)
    for ip in eachindex(x)
        i = clamp(1 + floor(Int, (x[ip] - xmin)*fx), 1, nel)
        j = clamp(1 + floor(Int, (y[ip] - ymin)*fy), 1, nel)
        o[i,j] = max(o[i,j], v[ip])
    end
    return o
end

# `nel` is the element count per direction of the mesh the run used — 32 for the
# `OT_32x32_periodic.msh` that ships with the case. It only has to be right for
# the binning to land one element per bin; a wrong value still produces a
# correlation, just of a coarser or finer field than μ actually is.
function mu_gradient_stats(gs, chosen, muname, isoname; nel = 32, nras = 128)
    println()
    println("  Does mu follow |grad ", isoname, "|?  (per-element, ", nel, "x", nel, " elements)")
    @printf("  %6s %11s %11s %9s %14s %12s\n",
            "t", "min mu", "max mu", "decades", "frac>max/10", "spearman")
    for (g, s) in zip(gs, chosen)
        xg = LinRange(extrema(g.x)..., nras)
        yg = LinRange(extrema(g.y)..., nras)
        z  = raster(g.x, g.y, field(g, isoname), xg, yg)
        gm = [hypot(z[i == nras ? 1 : i+1, j] - z[i == 1 ? nras : i-1, j],
                    z[i, j == nras ? 1 : j+1] - z[i, j == 1 ? nras : j-1])
              for i = 1:nras, j = 1:nras]

        gel = zeros(nel, nel)
        for i = 1:nras, j = 1:nras
            ii = clamp(1 + floor(Int, (i-1)/(nras-1)*nel), 1, nel)
            jj = clamp(1 + floor(Int, (j-1)/(nras-1)*nel), 1, nel)
            gel[ii,jj] = max(gel[ii,jj], gm[i,j])
        end

        μel = vec(_element_max(g.x, g.y, field(g, muname), nel))
        lo, hi = minimum(μel), maximum(μel)
        @printf("  %6.2f %11.2e %11.2e %9.2f %14.2f %12.2f\n",
                s.time, lo, hi, lo > 0 ? log10(hi/lo) : NaN,
                count(>=(hi/10), μel)/length(μel), _spearman(μel, vec(gel)))
    end
    println()
end

function save_matrix(plts, path)
    mkpath(dirname(path))
    fig = Plots.plot(plts...; layout = (2, 2), size = (1100, 1000))
    fig[:overwrite_figure] = false
    Plots.savefig(fig, path)
    println("  wrote ", path)
end

#---------------------------------------------------------------------------------

function main(argv)
    opts   = parse_args(argv)
    outdir = opts["out"]
    times  = _floats(opts["times"])
    length(times) == 4 || error("--times needs exactly four values (the figures are 2x2 panels)")

    dir   = find_run_dir(opts["dir"])
    snaps = list_snapshots(dir, parse(Float64, opts["t0"]), parse(Float64, opts["dtout"]))
    println("Reading ", dir, " (", length(snaps), " snapshots, t = ",
            @sprintf("%.3f", first(snaps).time), " … ",
            @sprintf("%.3f", last(snaps).time), ")")

    chosen = [pick(snaps, t) for t in times]
    gs     = [read_vtk_grid(s.path) for s in chosen]

    Plots.gr()
    Plots.default(fontfamily = "Computer Modern")

    # --- density, on the paper's fixed [0.1, 0.4] scale ------------------------
    rho_clims = Tuple(_floats(opts["rho-clims"]))
    save_matrix([panel(g, "ρ", s.time, rho_clims, :viridis; title = "ρ")
                 for (g, s) in zip(gs, chosen)],
                joinpath(outdir, "MHD_OT_rho.png"))

    # --- By, on a symmetric diverging scale shared by the four panels ---------
    bmax = maximum(maximum(abs, field(g, "By")) for g in gs)
    by_clims = (-bmax, bmax)
    save_matrix([panel(g, "By", s.time, by_clims, :balance; title = "B_y")
                 for (g, s) in zip(gs, chosen)],
                joinpath(outdir, "MHD_OT_By.png"))

    # --- DynSGS eddy viscosity of the By equation, with By isolines on top ----
    #
    # Drawn in log10 because μ can span decades between the smooth regions and
    # the shocks. The color range is the ACTUAL range over the four panels, not
    # a fixed number of decades below the maximum: how many decades the field
    # spans is part of what the figure is reporting, so padding it would hide
    # the answer. Read the colorbar, not the colors.
    muname = "mu_dsgs_By"
    if any(g -> !haskey(g.data, muname), gs)
        @warn """no "$muname" field in this output — the DynSGS viscosity figure is skipped.
                 It is written only when the run used :visc_model => DSGS_MHD()."""
    else
        μmax = maximum(maximum(field(g, muname)) for g in gs)
        μmin = minimum(minimum(filter(>(0), field(g, muname)); init = μmax) for g in gs)
        μmax > 0 || error("$muname is identically zero in every snapshot")
        mu_clims = (μmin, μmax)
        save_matrix([panel(g, muname, s.time, mu_clims, :inferno;
                           title = "log10 mu_DynSGS (By eq.) + By isolines",
                           isolines = "By", logscale = true)
                     for (g, s) in zip(gs, chosen)],
                    joinpath(outdir, "MHD_OT_mu_dsgs_By.png"))

        mu_gradient_stats(gs, chosen, muname, "By")
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
