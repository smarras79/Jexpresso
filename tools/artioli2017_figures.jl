#!/usr/bin/env julia
#---------------------------------------------------------------------------------
# tools/artioli2017_figures.jl — the plots and tables of
#
#   F. Artioli, L. Beirão da Veiga, C. Lovadina, E. Sacco, "Arbitrary order 2D
#   virtual elements for polygonal meshes: part II, inelastic problem",
#   Comput. Mech. 60 (2017), arXiv:1701.06676 — Sec. 4.2 and 4.3,
#
# redrawn from Jexpresso's VTK output (the iter_<k>.pvtu files a run writes,
# read with tools/vtu_reader.jl) and laid over the paper's own curves
# (tools/artioli2017_data: read from the vector graphics of the arXiv PDF).
#
#   julia --project=. tools/artioli2017_figures.jl                 run the cases, then plot
#   julia --project=. tools/artioli2017_figures.jl --no-run        plot the runs of an earlier call
#   julia --project=. tools/artioli2017_figures.jl --only=plate    (or --only=cylinder)
#   julia --project=. tools/artioli2017_figures.jl --no-run \
#         --plate-dir=output/Elasticity/perforated_plate2d/output  plot a run of your own
#
# Options
#   --no-run            do not run anything: read <outdir>/runs/<label>, or the
#                       directory given with --<label>-dir (a case with neither
#                       is skipped)
#   --only=WHAT         all (default) | cylinder (Fig. 5.2) | plate (Fig. 5.4, Table 1)
#   --outdir=DIR        where the figures go [output/artioli2017]
#   --ve1-dir=DIR, --ve2-dir=DIR, --plate-dir=DIR
#                       read that run's output directory (the one holding
#                       simulation.pvd) instead of running the case
#   --scale=S           displacement magnification of the deformed plate [5]
#   --no-anim           skip the animated gif of the plate
#
# The runs (each case at the paper's settings, through the ENV overrides its
# user_inputs.jl documents; the output is copied to <outdir>/runs/<label>):
#   ve1, ve2  problems/Elasticity/viscoelastic_cylinder2d, Prony sets 1 and 2,
#             nop = 4, Δt = 1 and the load entering in the FIRST increment —
#             FEAP's convention, which is what Fig. 5.2 shows. VTK at t = 0..20.
#             About a minute each.
#   plate     problems/Elasticity/perforated_plate2d, nop = 4, 400 increments of
#             the top displacement, VTK every 4 increments. Several minutes.
#
# What comes out, in <outdir>:
#   fig5_2a_ve1.png, fig5_2b_ve2.png   radial displacement of A (r = a) and B (r = b)
#                                      against the integration step
#   fig5_4_reaction.png                reaction sum on the top edge against the load step
#   plate_plastic_strain.png           ēᵖ on the deformed quarter plate, end of loading
#   plate_deforming.gif                the plate deforming, beside its reaction curve
#   table1.md                          u_A, v_B at the end of the loading (Table 1) and
#                                      the Newton iterations per increment (Table 2)
#   fig5_2_sem.csv, fig5_4_sem.csv     the Jexpresso data plotted
# The plate run also leaves the solver's own deformed-*.png and deformed.gif in
# <outdir>/runs/plate.
#
# Not reproducible here: Figs. 5.5–5.6 (the shape-memory-alloy device, Sec. 4.4)
# — that constitutive model is not implemented.
#---------------------------------------------------------------------------------

get!(ENV, "GKSwstype", "100")                   # GR without a display
using Printf
using DelimitedFiles
using Plots

include(joinpath(@__DIR__, "vtu_reader.jl"))
using .VTUReader

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DATA = joinpath(@__DIR__, "artioli2017_data")

# Colours: the two categorical slots of a validated palette (blue, orange: CVD
# ΔE 24.7, both ≥ 3:1 on white) for the two points A and B; the paper's curves
# in neutral grey ink; one sequential blue ramp, light → dark, for ēᵖ.
const C_A    = "#2a78d6"
const C_B    = "#eb6834"
const C_REF  = "#52514e"
const C_BAND = "#c9c8c3"
const SEQ    = cgrad(["#e8f1fc", "#b7d3f6", "#86b6ef", "#5598e7", "#2a78d6",
                      "#1c5cab", "#104281", "#0d366b"])
const STYLE  = (; dpi = 150, size = (720, 470), legendfontsize = 8, guidefontsize = 10,
                  tickfontsize = 8, titlefontsize = 10, gridalpha = 0.15,
                  foreground_color_legend = nothing, background_color_legend = RGBA(1, 1, 1, 0.85),
                  left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)

#---------------------------------------------------------------------------------
# Cases and their settings
#---------------------------------------------------------------------------------
const VE_ENV = ["JEXPRESSO_VE_DT" => "1.0", "JEXPRESSO_VE_LOAD" => "first_increment",
                "JEXPRESSO_VE_INTEGRATOR" => "exponential", "JEXPRESSO_VE_NOP" => "4",
                "JEXPRESSO_VE_OUTFORMAT" => "vtk", "JEXPRESSO_VE_PLOT" => "false",
                "JEXPRESSO_VE_HISTORY" => "nodal"]
const CASES = Dict(
    "ve1"   => (case = "viscoelastic_cylinder2d", env = ["JEXPRESSO_VE_SET" => "ve1"; VE_ENV]),
    "ve2"   => (case = "viscoelastic_cylinder2d", env = ["JEXPRESSO_VE_SET" => "ve2"; VE_ENV]),
    "plate" => (case = "perforated_plate2d",
                env = ["JEXPRESSO_PP_STEPS" => "400", "JEXPRESSO_PP_NOP" => "4",
                       "JEXPRESSO_PP_OUTFORMAT" => "vtk", "JEXPRESSO_PP_PLOT" => "true"]),
)

function parse_options(args)
    o = Dict{String,Any}("run" => true, "only" => "all", "anim" => true, "scale" => 5.0,
                         "outdir" => joinpath(ROOT, "output", "artioli2017"),
                         "ve1-dir" => nothing, "ve2-dir" => nothing, "plate-dir" => nothing)
    for a in args
        if a in ("-h", "--help")
            for l in eachline(@__FILE__)
                startswith(l, "#") || break
                println(l)
            end
            exit(0)
        elseif a == "--no-run"
            o["run"] = false
        elseif a == "--no-anim"
            o["anim"] = false
        elseif startswith(a, "--") && occursin('=', a)
            k, v = split(a[3:end], '='; limit = 2)
            haskey(o, k) || error("unknown option --$k (see --help)")
            o[k] = k == "scale" ? parse(Float64, v) : k == "only" ? String(v) : abspath(expanduser(v))
        else
            error("unknown argument \"$a\" (see --help)")
        end
    end
    o["only"] in ("all", "cylinder", "plate") || error("--only=all|cylinder|plate")
    return o
end

labels_for(only) = only == "cylinder" ? ["ve1", "ve2"] : only == "plate" ? ["plate"] : ["ve1", "ve2", "plate"]

#
# Run one case with its settings and copy what THIS run wrote (files newer
# than its start: the case's output directory may hold older runs) to dest.
#
function run_and_collect(label, dest)
    c   = CASES[label]
    src = joinpath(ROOT, "output", "Elasticity", c.case, "output")
    println("\n=== running ", label, " (problems/Elasticity/", c.case, ") ===")
    t0 = time()
    cd(ROOT) do
        withenv(c.env...) do
            Base.invokelatest(getfield(Main, :Jexpresso).run_case, "Elasticity", c.case)
        end
    end
    rm(dest; force = true, recursive = true)
    mkpath(dest)
    for (dir, _, files) in walkdir(src), f in files
        p = joinpath(dir, f)
        mtime(p) >= t0 - 1 || continue
        q = joinpath(dest, relpath(p, src))
        mkpath(dirname(q))
        cp(p, q; force = true)
    end
    @printf("    %s finished in %.0f s; output copied to %s\n", label, time() - t0, dest)
    return dest
end

#---------------------------------------------------------------------------------
# Reading a run: simulation.pvd lists (time, file) for every VTK output.
#---------------------------------------------------------------------------------
function read_series(dir)
    pvd = joinpath(dir, "simulation.pvd")
    isfile(pvd) || error("no simulation.pvd in $dir: run the case with :outformat => \"vtk\" ",
                         "(or point --<label>-dir at the directory that holds it).")
    snaps = Tuple{Float64,VTKGrid}[]
    for m in eachmatch(r"<DataSet\b[^>]*>", read(pvd, String))
        t = match(r"timestep\s*=\s*\"([^\"]*)\"", m.match)
        f = match(r"file\s*=\s*\"([^\"]*)\"", m.match)
        (t === nothing || f === nothing) && continue
        push!(snaps, (parse(Float64, t.captures[1]), read_vtk_grid(joinpath(dir, f.captures[1]))))
    end
    isempty(snaps) && error("$pvd lists no output")
    return sort!(snaps; by = first)
end

nearest(g, x0, y0) = argmin((g.x .- x0) .^ 2 .+ (g.y .- y0) .^ 2)

function field(g, name)
    haskey(g.data, name) && return g.data[name]
    error("the VTK files have no field \"$name\" (they have: ", join(sort!(collect(keys(g.data))), ", "),
          "). Re-run the case with this version of its initialize.jl.")
end

function read_reference(path)
    lines = filter(l -> !startswith(l, "#") && !isempty(strip(l)), readlines(path))
    hdr   = split(lines[1], ",")
    vals  = reduce(vcat, (permutedims(parse.(Float64, split(l, ","))) for l in lines[2:end]))
    return Dict(String(h) => vals[:, i] for (i, h) in enumerate(hdr))
end

#---------------------------------------------------------------------------------
# Fig. 5.2: u_A and u_B against the integration step, one Prony set per panel
#---------------------------------------------------------------------------------
function figure_5_2(series, set, panel, ref, outdir)
    g0 = series[1][2]
    # A and B on the x axis: the inner and the outer radius of the quarter annulus
    r  = hypot.(g0.x, g0.y)
    a, b = minimum(r), maximum(r)
    iA = nearest(g0, a, 0.0); iB = nearest(g0, b, 0.0)
    t  = first.(series)
    uA = [field(g, "ux")[iA] for (_, g) in series]
    uB = [field(g, "ux")[iB] for (_, g) in series]

    μ  = set == "ve1" ? "(0.01, 0.99)" : "(0.3, 0.7)"
    ytop = max(set == "ve1" ? 0.7 : 0.12, 1.04*max(maximum(uA), maximum(ref["$(set)_uA"])))   # the paper's axes
    plt = plot(; xlabel = "integration step  [−]   (Δt = λ = 1)", ylabel = "displacement  [mm]",
               title = "Fig. 5.2($panel)  viscoelastic cylinder, (μ₀, μ₁) = $μ",
               legend = :topleft, xlims = (0, 20.5), ylims = (0, ytop), xticks = 0:2:20, STYLE...)
    plot!(plt, t, uA; color = C_A, lw = 2, label = "u_A  (r = a)   Jexpresso SEM, nop = 4")
    plot!(plt, t, uB; color = C_B, lw = 2, label = "u_B  (r = b)   Jexpresso SEM, nop = 4")
    scatter!(plt, ref["step"], ref["$(set)_uA"]; marker = :circle, ms = 5, mc = :white, msc = C_A,
             msw = 1.5, label = "u_A   paper, Ref (FEAP)")
    scatter!(plt, ref["step"], ref["$(set)_uB"]; marker = :utriangle, ms = 5, mc = :white, msc = C_B,
             msw = 1.5, label = "u_B   paper, Ref (FEAP)")
    annotate!(plt, t[end] - 0.3, uA[end], text("A", 9, C_REF, :right, :bottom))
    annotate!(plt, t[end] - 0.3, uB[end], text("B", 9, C_REF, :right, :bottom))
    f = joinpath(outdir, "fig5_2$(panel)_$(set).png")
    savefig(plt, f)

    # the gap at the paper's steps
    dA = 0.0; dB = 0.0
    for (k, s) in enumerate(ref["step"])
        j = findfirst(tt -> abs(tt - s) < 1e-9, t)
        j === nothing && continue
        dA = max(dA, abs(uA[j] - ref["$(set)_uA"][k])); dB = max(dB, abs(uB[j] - ref["$(set)_uB"][k]))
    end
    @printf("  Fig. 5.2(%s) %s: max |u − u_Ref| over the steps: A %.2e, B %.2e mm   -> %s\n",
            panel, set, dA, dB, f)
    return (; t, uA, uB, dA, dB)
end

#---------------------------------------------------------------------------------
# The perforated plate: geometry read off the grid
#---------------------------------------------------------------------------------
function plate_geometry(g)
    L = maximum(g.x); H = maximum(g.y); R = minimum(hypot.(g.x, g.y))
    top = findall(y -> abs(y - H) <= 1e-8*H, g.y)
    return (; L, H, R, top)
end

# load step from the prescribed top displacement δ(t): 400 increments to 2 mm
const PLATE_DELTA = 2.0
const PLATE_STEPS = 400

# The paper's "reaction sum" is TEN times the force Jexpresso computes, which is
# per unit thickness (plane strain, unit depth): over the elastic part of the
# curve the ratio is 10.00 for every one of the paper's five methods
# (elastic_ratio below measures and prints it), while the displacements of
# Table 1 agree — the same plate, its force reported for a 10 mm thickness the
# paper does not state. Jexpresso's curve is multiplied by this factor to land
# on the paper's axes; fig5_4_sem.csv has it both ways.
const PAPER_FORCE_SCALE = 10.0
const METHODS = ["VEM_Quad", "VEM_Tri", "VEM_Voronoi", "FEM_Q9", "FEM_T6"]
const METHOD_NAMES = Dict("VEM_Quad" => "VEM Quad", "VEM_Tri" => "VEM Tri", "VEM_Voronoi" => "VEM Voronoi",
                          "FEM_Q9" => "FEM Q9", "FEM_T6" => "FEM T6")
# for the tables: the paper's quadratic methods — VEM of order k = 2, Q9 and T6 FEM
paper_method(name) = startswith(name, "VEM") ? "$name, k = 2" : name

# linear interpolation of (xs, ys) at x (xs increasing)
function interp(xs, ys, x)
    j = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    w = (x - xs[j])/(xs[j + 1] - xs[j])
    return (1 - w)*ys[j] + w*ys[j + 1]
end

# paper / Jexpresso over the elastic steps (the curve is straight up to ~45)
function elastic_ratio(curve, ref; steps = 4:4:40)
    r = [ref[m][findfirst(==(s), ref["step"])]/interp(curve.step, curve.Rsum, s) for m in METHODS, s in steps]
    return (; mean = sum(r)/length(r), lo = minimum(r), hi = maximum(r))
end

function plate_curve(series)
    geo  = plate_geometry(series[1][2])
    δ    = [sum(field(g, "uy")[geo.top])/length(geo.top) for (_, g) in series]
    Rsum = [sum(field(g, "ry")[geo.top]) for (_, g) in series]
    return (; geo, step = δ ./ PLATE_DELTA .* PLATE_STEPS, δ, Rsum)
end

function figure_5_4(series, ref, outdir)
    c   = plate_curve(series)
    er  = elastic_ratio(c, ref)
    Rp  = PAPER_FORCE_SCALE .* c.Rsum                  # on the paper's scale
    lo  = [minimum(ref[m][i] for m in METHODS) for i in eachindex(ref["step"])]
    hi  = [maximum(ref[m][i] for m in METHODS) for i in eachindex(ref["step"])]
    lbl = @sprintf("Jexpresso SEM, nop = 4  (force × %g)", PAPER_FORCE_SCALE)

    # (a) the paper's figure, on its axes (reaction in units of 10⁴ kg)
    pa = plot(; xlabel = "load step  [-]   (400 steps to δ = 2 mm)", ylabel = "reaction sum  [10⁴ kg]",
              title = "Fig. 5.4  perforated plate: reaction sum on the top edge",
              legend = :bottomright, xlims = (0, 405), ylims = (0, 2.1), yticks = 0:0.2:2,
              STYLE...)
    for (k, m) in enumerate(METHODS)
        plot!(pa, ref["step"], ref[m] ./ 1e4; color = C_BAND, lw = 1.5,
              label = k == 1 ? "paper: VEM Quad, Tri, Voronoi; FEM Q9, T6 (k = 2)" : "")
    end
    plot!(pa, c.step, Rp ./ 1e4; color = C_A, lw = 2, label = lbl)

    # (b) the plateau, magnified: the five paper curves told apart by their
    # end labels (nudged apart where two ends nearly coincide), Jexpresso in colour
    zlo = 150.0
    sel = ref["step"] .>= zlo; selc = c.step .>= zlo - 1e-9
    ylo = min(minimum(lo[sel]), minimum(Rp[selc])); yhi = max(maximum(hi[sel]), maximum(Rp[selc]))
    pad = 0.08*(yhi - ylo)
    pb = plot(; xlabel = "load step  [-]", ylabel = "reaction sum  [kg]",
              title = "the plateau, magnified", legend = false, xlims = (zlo, 475),
              xticks = zlo:50:400, ylims = (ylo - pad, yhi + pad), STYLE...)
    ends = [(ref[m][end], METHOD_NAMES[m], C_REF) for m in METHODS]
    push!(ends, (Rp[end], "Jexpresso", C_A))
    sort!(ends; by = first)
    ylab = [e[1] for e in ends]
    gap  = 0.045*(yhi - ylo + 2pad)
    for _ = 1:50, k = 2:length(ylab)                 # push labels apart, keep the order
        d = ylab[k] - ylab[k-1]
        d < gap && (ylab[k-1] -= (gap - d)/2; ylab[k] += (gap - d)/2)
    end
    for m in METHODS
        plot!(pb, ref["step"][sel], ref[m][sel]; color = C_BAND, lw = 1.5)
    end
    plot!(pb, c.step[selc], Rp[selc]; color = C_A, lw = 2)
    for (k, (yend, name, col)) in enumerate(ends)
        plot!(pb, [400, 404], [yend, ylab[k]]; color = col, lw = 0.6)
        annotate!(pb, 406, ylab[k], text(name, 7, col, :left))
    end

    plt = plot(pa, pb; layout = (1, 2), size = (1240, 480), dpi = 150,
               left_margin = 7Plots.mm, bottom_margin = 7Plots.mm)
    f = joinpath(outdir, "fig5_4_reaction.png")
    savefig(plt, f)

    @printf("  Fig. 5.4: the paper's reaction sums are %.3f x Jexpresso's force per unit thickness over the elastic steps (%.3f-%.3f over its five methods)\n",
            er.mean, er.lo, er.hi)
    @printf("            end of the loading: Jexpresso %.1f kg (x %g = %.0f); the paper's quadratic methods %.0f-%.0f kg   -> %s\n",
            c.Rsum[end], PAPER_FORCE_SCALE, Rp[end], lo[end], hi[end], f)
    return merge(c, (; er))
end

#---------------------------------------------------------------------------------
# A field on the deformed plate. The VTK cells are the sub-quads between
# neighbouring LGL nodes; each is split into two triangles and the nodal values
# interpolated linearly onto a pixel grid (pixels outside the body: NaN).
#---------------------------------------------------------------------------------
function rasterize(x, y, f, cells; npix = 520, box = nothing)
    xmin, xmax, ymin, ymax = box === nothing ? (extrema(x)..., extrema(y)...) : box
    nx = npix; ny = max(2, round(Int, npix*(ymax - ymin)/(xmax - xmin)))
    gx = range(xmin, xmax; length = nx); gy = range(ymin, ymax; length = ny)
    img = fill(NaN, ny, nx)
    for c in cells, (p, q, r) in ((c[1], c[2], c[3]), (c[1], c[3], c[4]))
        x1, y1, x2, y2, x3, y3 = x[p], y[p], x[q], y[q], x[r], y[r]
        det = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
        abs(det) <= eps()*max(abs(x1), abs(y1), 1.0)^2 && continue
        for i in searchsortedfirst(gx, min(x1, x2, x3)):searchsortedlast(gx, max(x1, x2, x3)),
            j in searchsortedfirst(gy, min(y1, y2, y3)):searchsortedlast(gy, max(y1, y2, y3))
            l1 = ((y2 - y3)*(gx[i] - x3) + (x3 - x2)*(gy[j] - y3))/det
            l2 = ((y3 - y1)*(gx[i] - x3) + (x1 - x3)*(gy[j] - y3))/det
            l3 = 1 - l1 - l2
            (l1 >= -1e-10 && l2 >= -1e-10 && l3 >= -1e-10) || continue
            img[j, i] = l1*f[p] + l2*f[q] + l3*f[r]
        end
    end
    return gx, gy, img
end

# boundary of the sub-quad mesh: the edges that belong to one cell only
function boundary_segments(cells)
    count = Dict{Tuple{Int,Int},Int}()
    for c in cells, k in eachindex(c)
        e = minmax(c[k], c[mod1(k + 1, length(c))])
        count[e] = get(count, e, 0) + 1
    end
    return [e for (e, n) in count if n == 1]
end

function segments_xy(segs, x, y)
    X = Float64[]; Y = Float64[]
    for (p, q) in segs
        append!(X, (x[p], x[q], NaN)); append!(Y, (y[p], y[q], NaN))
    end
    return X, Y
end

function plate_panel(g, geo, scale, cmax, segs; title = "", box = nothing)
    xd = g.x .+ scale .* field(g, "ux"); yd = g.y .+ scale .* field(g, "uy")
    gx, gy, img = rasterize(xd, yd, field(g, "ep_bar"), g.cells; box = box)
    plt = heatmap(gx, gy, img; c = SEQ, clims = (0, cmax), aspect_ratio = :equal,
                  colorbar = :right, colorbar_title = "ēᵖ",
                  xlabel = "x  [mm]", ylabel = "y  [mm]", title = title,
                  xlims = (box === nothing ? extrema(gx) : box[1:2]),
                  ylims = (box === nothing ? extrema(gy) : box[3:4]), STYLE..., size = (560, 640),
                  legend = false, right_margin = 4Plots.mm)
    # the undeformed outline (dashed) and the deformed one
    θ  = range(0, π/2; length = 60)
    ox = [geo.R; geo.L; geo.L; 0.0; 0.0; geo.R .* sin.(θ)]
    oy = [0.0; 0.0; geo.H; geo.H; geo.R; geo.R .* cos.(θ)]
    plot!(plt, ox, oy; color = C_REF, lw = 1, ls = :dash, label = "")
    X, Y = segments_xy(segs, xd, yd)
    plot!(plt, X, Y; color = "#0b0b0b", lw = 1, label = "")
    return plt
end

function figure_plate_field(series, outdir, scale)
    g   = series[end][2]
    geo = plate_geometry(g)
    isempty(g.cells) && error("the VTK files carry no cells")
    segs = boundary_segments(g.cells)
    cmax = maximum(field(g, "ep_bar"))
    plt  = plate_panel(g, geo, scale, cmax, segs;
                       title = @sprintf("accumulated plastic strain ēᵖ\nt = %.3g: δ = %.2f mm, displacement × %g",
                                        series[end][1], sum(field(g, "uy")[geo.top])/length(geo.top), scale))
    f = joinpath(outdir, "plate_plastic_strain.png")
    savefig(plt, f)
    println("  plastic zone at the end of the loading   -> ", f)
    return nothing
end

function animate_plate(series, curve, ref, outdir, scale)
    g1   = series[end][2]
    geo  = plate_geometry(g1)
    segs = boundary_segments(g1.cells)
    cmax = max(maximum(field(g1, "ep_bar")), eps())
    # one frame box for all: the extent of the most deformed state
    xd = g1.x .+ scale .* field(g1, "ux"); yd = g1.y .+ scale .* field(g1, "uy")
    box = (min(0.0, minimum(xd)), max(geo.L, maximum(xd)), min(0.0, minimum(yd)), max(geo.H, maximum(yd)) + 2)
    Rp   = PAPER_FORCE_SCALE .* curve.Rsum ./ 1e4
    anim = Animation()
    for (k, (t, g)) in enumerate(series)
        p1 = plate_panel(g, geo, scale, cmax, segs; box = box,
                         title = @sprintf("plastic strain ēᵖ\nδ = %.3f mm  (step %d), displacement × %g", curve.δ[k],
                                          round(Int, curve.step[k]), scale))
        p2 = plot(; xlabel = "load step", ylabel = "reaction sum on the top edge  [10⁴ kg]", xlims = (0, 405),
                  ylims = (0, 2.1), yticks = 0:0.2:2, legend = :bottomright, STYLE..., size = (520, 640))
        for (i, m) in enumerate(METHODS)
            plot!(p2, ref["step"], ref[m] ./ 1e4; color = C_BAND, lw = 1.5,
                  label = i == 1 ? "paper (5 quadratic methods)" : "")
        end
        plot!(p2, curve.step[1:k], Rp[1:k]; color = C_A, lw = 2,
              label = @sprintf("Jexpresso SEM, nop = 4 (× %g)", PAPER_FORCE_SCALE))
        scatter!(p2, [curve.step[k]], [Rp[k]]; color = C_A, ms = 5, msw = 0, label = "")
        frame(anim, plot(p1, p2; layout = (1, 2), size = (1040, 640), dpi = 100))
    end
    f = joinpath(outdir, "plate_deforming.gif")
    gif(anim, f; fps = 12, show_msg = false)
    println("  animation                                 -> ", f)
    return nothing
end

#---------------------------------------------------------------------------------
# Table 1 (displacements at the end of the loading) and Table 2 (Newton)
#---------------------------------------------------------------------------------
function table_1(series, curve, ref, rundir, outdir)
    t, g = series[end]
    geo  = plate_geometry(g)
    iA   = nearest(g, geo.R, 0.0); iB = nearest(g, 0.0, geo.R)
    uA   = field(g, "ux")[iA]; vB = field(g, "uy")[iB]
    paper = [("Ref (overkill mixed u-p-εᵛ FEM)", 2.741, 1.859),
             ("VEM Quad, k = 2", 2.720, 1.851), ("VEM Tri, k = 2", 2.719, 1.852),
             ("VEM Voronoi, k = 2", 2.714, 1.852), ("FEM Q9", 2.733, 1.853), ("FEM T6", 2.719, 1.853)]
    newton_paper = [("VEM Quad", 6.43), ("VEM Tri", 6.11), ("VEM Voronoi", 6.21), ("FEM Q9", 6.70), ("FEM T6", 6.11)]

    io = IOBuffer()
    println(io, "# Perforated plastic plate — Artioli et al. (2017), Sec. 4.3\n")
    @printf(io, "End of the loading (t = %.4g, δ = %.4f mm). A = (R, 0), B = (0, R), R = %g mm.\n\n",
            t, sum(field(g, "uy")[geo.top])/length(geo.top), geo.R)
    println(io, "## Table 1: displacements [mm]\n")
    println(io, "| method | u_A | v_B | u_A − Ref | v_B − Ref |")
    println(io, "|---|---:|---:|---:|---:|")
    @printf(io, "| **Jexpresso SEM, nop = 4, 132 elements** | **%.4f** | **%.4f** | %+.2f%% | %+.2f%% |\n",
            uA, vB, 100*(uA/2.741 - 1), 100*(vB/1.859 - 1))
    for (m, a, b) in paper
        @printf(io, "| paper: %s | %.3f | %.3f | %+.2f%% | %+.2f%% |\n", m, a, b, 100*(a/2.741 - 1), 100*(b/1.859 - 1))
    end

    println(io, "\n## Fig. 5.4: reaction sum on the top edge at the end of the loading\n")
    println(io, "| method | reaction sum [kg] |")
    println(io, "|---|---:|")
    @printf(io, "| **Jexpresso SEM, nop = 4** (per mm of thickness × %g) | **%.0f** |\n",
            PAPER_FORCE_SCALE, PAPER_FORCE_SCALE*curve.Rsum[end])
    for m in METHODS
        @printf(io, "| paper: %s | %.0f |\n", paper_method(METHOD_NAMES[m]), ref[m][end])
    end
    @printf(io, "\nJexpresso computes the force per unit thickness, %.2f kg/mm. Over the elastic steps the ",
            curve.Rsum[end])
    @printf(io, "paper's sums are %.3f times it (%.3f–%.3f across its five methods) while the displacements ",
            curve.er.mean, curve.er.lo, curve.er.hi)
    println(io, "above agree, so the paper's force is that of a plate 10 mm thick (a thickness it does not state); ",
            "the figure multiplies Jexpresso's by ", PAPER_FORCE_SCALE, ".")

    csv = joinpath(rundir, "solid_history.csv")
    if isfile(csv)
        lines = readlines(csv)
        hdr   = split(lines[1], ",")
        ji    = findfirst(==("newton_iters"), hdr); jc = findfirst(==("newton_cuts"), hdr)
        if ji !== nothing && length(lines) > 2
            its  = [parse(Float64, split(l, ",")[ji]) for l in lines[3:end]]
            cuts = jc === nothing ? zeros(length(its)) : [parse(Float64, split(l, ",")[jc]) for l in lines[3:end]]
            println(io, "\n## Table 2: Newton iterations per load increment\n")
            println(io, "| method | average |")
            println(io, "|---|---:|")
            @printf(io, "| **Jexpresso SEM, nop = 4** | **%.2f** |\n", sum(its)/length(its))
            for (m, v) in newton_paper
                @printf(io, "| paper: %s | %.2f |\n", paper_method(m), v)
            end
            @printf(io, "\nJexpresso: %d increments, at most %d iterations in one, %d increment(s) cut in halves ",
                    length(its), maximum(its), count(>(0), cuts))
            println(io, "(the iterations of the halves are counted against the increment).",
                    " Convergence: ‖r‖ ≤ 1e-10 ‖f‖ on the free DOFs (solid_history.csv).")
        end
    end
    s = String(take!(io))
    f = joinpath(outdir, "table1.md")
    write(f, s)
    @printf("  Table 1: u_A = %.4f mm (Ref 2.741, %+.2f%%), v_B = %.4f mm (Ref 1.859, %+.2f%%)   -> %s\n",
            uA, 100*(uA/2.741 - 1), vB, 100*(vB/1.859 - 1), f)
    return (; uA, vB)
end

#---------------------------------------------------------------------------------
function main(args)
    o = parse_options(args)
    outdir = mkpath(o["outdir"])
    labels = labels_for(o["only"])

    # where each run's output is: given, run now, or from an earlier call
    # (with --no-run, a case that has neither is skipped)
    dirs = Dict{String,String}()
    torun = String[]
    for l in copy(labels)
        if o["$l-dir"] !== nothing
            dirs[l] = o["$l-dir"]
        elseif o["run"]
            push!(torun, l)
            dirs[l] = joinpath(outdir, "runs", l)
        elseif isfile(joinpath(outdir, "runs", l, "simulation.pvd"))
            dirs[l] = joinpath(outdir, "runs", l)
        else
            println("  (no run of $l in $(joinpath(outdir, "runs", l)) and no --$l-dir: skipped)")
            filter!(!=(l), labels)
        end
    end
    isempty(labels) && error("--no-run, and no output to read: drop --no-run, or give --ve1-dir=, --ve2-dir=, --plate-dir=.")
    if !isempty(torun)
        @eval Main using Jexpresso
        for l in torun
            run_and_collect(l, dirs[l])
        end
    end

    println("\n=== figures in ", outdir, " ===")
    if any(in(labels), ("ve1", "ve2"))
        ref = read_reference(joinpath(DATA, "fig5_2_ref.csv"))
        rows = Vector{Vector{Float64}}()
        res  = Dict{String,Any}()
        for (set, panel) in (("ve1", "a"), ("ve2", "b"))
            set in labels || continue
            res[set] = figure_5_2(read_series(dirs[set]), set, panel, ref, outdir)
        end
        open(joinpath(outdir, "fig5_2_sem.csv"), "w") do io
            println(io, "# Jexpresso SEM (nop = 4) radial displacement at A (r = a) and B (r = b) [mm]; Δt = 1, FEAP load convention")
            println(io, "step,", join(["$(s)_uA,$(s)_uB" for s in sort!(collect(keys(res)))], ","))
            n  = minimum(length(r.t) for r in values(res))
            t1 = first(values(res)).t
            for k = 1:n
                vals = [t1[k]]
                for s in sort!(collect(keys(res)))
                    push!(vals, res[s].uA[k], res[s].uB[k])
                end
                println(io, join((@sprintf("%.10g", v) for v in vals), ","))
            end
        end
    end
    if "plate" in labels
        series = read_series(dirs["plate"])
        ref    = read_reference(joinpath(DATA, "fig5_4.csv"))
        curve  = figure_5_4(series, ref, outdir)
        open(joinpath(outdir, "fig5_4_sem.csv"), "w") do io
            println(io, "# Jexpresso SEM (nop = 4): reaction sum on the top edge against the load step. ",
                    "reaction_per_mm: per unit thickness [kg/mm]; reaction_paper_scale: x ", PAPER_FORCE_SCALE,
                    ", the scale of the paper's Fig. 5.4 [kg]")
            println(io, "step,delta_mm,reaction_per_mm,reaction_paper_scale")
            for k in eachindex(curve.step)
                @printf(io, "%.6g,%.10g,%.10g,%.10g\n", curve.step[k], curve.δ[k], curve.Rsum[k],
                        PAPER_FORCE_SCALE*curve.Rsum[k])
            end
        end
        table_1(series, curve, ref, dirs["plate"], outdir)
        figure_plate_field(series, outdir, o["scale"])
        o["anim"] && animate_plate(series, curve, ref, outdir, o["scale"])
    end
    return nothing
end

main(ARGS)
