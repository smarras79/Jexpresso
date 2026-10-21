#---------------------------------------------------------------------------------
# THE CONVERGENCE FIGURE of the accuracy tests, in the layout of the
# convergence figure of Dao & Nazarov (2022), Fig. 1: an error at the final
# time against the NUMBER OF DEGREES OF FREEDOM on log-log axes — their
# abscissa, and the honest one for comparing orders, since it compares them at
# equal cost rather than at equal h — with every polynomial order and every
# stabilization ON THE SAME AXES (one colour and marker per order, one line
# style per stabilization), dotted reference slopes, and the measured rate p of
# each curve in its legend. Written as the panels side by side and as one file
# per panel, at publication sizes.
#
# Drawn by this case (user_plot.jl) and by alfvenWaveBabbarRanocha2026, which
# includes this file, so that the two accuracy tests show their errors the same
# way. Plain Julia on Plots and LaTeXStrings, so that the tools that redraw a
# figure from a store (tools/smooth_vortex_plot.jl, tools/alfven_wave_plot.jl)
# load it without Jexpresso; `_savefig_silent` comes from
# src/io/plotting/jeplots.jl in a run and from those tools outside one.
#
# The records are NamedTuples with nop, ndofs, visc (the record's tag), t
# (the final time), dt (the time step, NaN when unknown) and the error fields
# the panels name. In two dimensions #DOFs ∝ h^(-2), so an O(h^(N+1)) method
# falls as #DOFs^(-(N+1)/2); the guides carry that slope and are labelled
# O(h^(N+1)).
#---------------------------------------------------------------------------------
const CONV_FS_TITLE  = 22
const CONV_FS_GUIDE  = 20
const CONV_FS_TICK   = 17
const CONV_FS_LEGEND = 15
const CONV_LW        = 1.4
const CONV_MS        = 7

# One (colour, marker) per order, fixed so an order looks the same from one
# figure to the next, and from one test to the other.
const CONV_STYLE = Dict(
    1 => (:darkorange, :dtriangle),
    2 => (:goldenrod, :utriangle),
    3 => (:purple,    :rect),
    4 => (:seagreen,  :diamond),
    5 => (:royalblue, :circle),
    6 => (:crimson,   :star5),
    7 => (:black,     :xcross),
)
_conv_style(nop) = get(CONV_STYLE, nop, (:gray, :cross))

# The measured order of accuracy from the last two points. The abscissa is
# the DOF COUNT, and in two dimensions h ∝ #DOFs^(-1/2), so an error ∝ h^p
# falls as #DOFs^(-p/2): the order is minus twice the log-log slope.
_conv_rate(xs, ys) = (length(xs) < 2 || ys[end-1] <= 0 || ys[end] <= 0) ? NaN :
                     -2.0*log(ys[end-1]/ys[end])/log(xs[end-1]/xs[end])

# A name set inside \mathrm, where a space must be "\ ".
_conv_mathrm(s) = replace(string(s), " " => "\\ ")

# ATOMIC WRITES: the runs of a sweep finish at different times and every one
# redraws the same figures from the shared store, so a figure is written
# beside its target and renamed (atomic on POSIX): a reader sees either the
# old file or the new one, never half of each.
function _conv_savefig(plt, path::AbstractString)
    mkpath(dirname(path))
    base, ext = splitext(path)                   # keep the extension: Plots picks the format from it
    tmp = string(base, ".tmp", getpid(), ext)
    try
        _savefig_silent(plt, tmp)
        mv(tmp, path; force = true)
    catch err
        isfile(tmp) && (try; rm(tmp); catch; end)
        rethrow(err)
    end
    return nothing
end

#---------------------------------------------------------------------------------
# One panel: every (order, tag) curve of `rows` with the error field `fld`.
#
#   ylab, title  LaTeX for the axis and for the title; the final time, and the
#                time step when every record shares one, are appended to the
#                title — a sweep cut short for a pipeline check draws a
#                perfectly plausible set of flat lines, and nothing in the
#                picture would otherwise say so
#   tags         the tags to draw, in legend order among the variants
#   label(tag)   the name of a tag in the legend, plain text
#   line(tag, tags)  its line style: solid for the stabilization the case runs
#                by default, dashed for plain Galerkin, the other dash patterns
#                for the variants. DOTTED is reserved for the reference slopes,
#                so that no guide can be mistaken for a curve.
#   galerkin(tag)  whether a tag is the unstabilized run
#
# Points that COINCIDE nest instead of hiding one another: the smoothness
# cutoff makes ν exactly zero on a resolved smooth solution, so there DynSGS IS
# Galerkin, point for point. Galerkin is drawn first, as the largest hollow
# marker, the variants over it as smaller hollow ones, and the default
# stabilization last, as a filled marker in the middle. (A hollow marker has a
# white inside: GR draws neither a :transparent nor an alpha-0 fill as one.)
# The legend is added separately, so that it lists the default stabilization
# first, then the variants, then Galerkin, whatever order they are drawn in.
#---------------------------------------------------------------------------------
function _conv_panel(rows, fld::Symbol, ylab, title; tags = sort(unique(r.visc for r in rows)),
                     label = string, line = (tag, tags) -> :solid,
                     galerkin = tag -> startswith(string(tag), "galerkin"))
    nops = sort(unique(r.nop for r in rows))
    tt   = isempty(rows) ? NaN : rows[1].t
    dts  = unique(r.dt for r in rows)
    stamp = string(",\\ t = ", isfinite(tt) ? tt : "?",
                   length(dts) == 1 && isfinite(dts[1]) ? string(",\\ \\Delta t = ", dts[1]) : "")
    allx = Float64[]; ally = Float64[]
    pl = Plots.plot(; xscale = :log10, yscale = :log10,
                    xlabel = LaTeXStrings.L"\#\mathrm{DOFs}",
                    ylabel = LaTeXStrings.latexstring(ylab),
                    framestyle = :box, grid = true, gridalpha = 0.25,
                    background_color_legend = Plots.RGBA(1, 1, 1, 0.72),
                    foreground_color_legend = :gray60,
                    legend = :outertopright, legendfontsize = CONV_FS_LEGEND,
                    titlefontsize = CONV_FS_TITLE, guidefontsize = CONV_FS_GUIDE,
                    tickfontsize = CONV_FS_TICK,
                    left_margin = 14Plots.mm, bottom_margin = 10Plots.mm,
                    top_margin = 4Plots.mm, right_margin = 6Plots.mm,
                    title = LaTeXStrings.latexstring(string(title, stamp)),
                    show = false)

    base(t)  = !galerkin(t) && line(t, tags) == :solid
    drawn    = vcat(filter(galerkin, tags), filter(t -> !galerkin(t) && !base(t), tags), filter(base, tags))
    listed   = vcat(filter(base, tags), filter(t -> !galerkin(t) && !base(t), tags), filter(galerkin, tags))
    attrs(nop, tag) = (c = _conv_style(nop);
                       (line = (c[1], CONV_LW, line(tag, tags)),
                        marker = (c[2], galerkin(tag) ? CONV_MS + 4 : base(tag) ? CONV_MS : CONV_MS + 2),
                        markerstrokecolor = c[1], markerstrokewidth = 1.0,
                        markercolor = base(tag) ? c[1] : :white, color = c[1]))
    curves = Dict{Tuple{Int,String},Tuple{Vector{Float64},Vector{Float64}}}()
    for nop in nops, tag in tags
        g = sort(filter(r -> r.nop == nop && r.visc == tag, rows), by = r -> r.ndofs)
        isempty(g) && continue
        xs = [Float64(r.ndofs) for r in g]
        ys = [Float64(getfield(r, fld)) for r in g]
        keep = isfinite.(ys) .& (ys .> 0)
        any(keep) || continue
        curves[(nop, string(tag))] = (xs[keep], ys[keep])
        append!(allx, xs[keep]); append!(ally, ys[keep])
    end
    for nop in nops, tag in drawn
        haskey(curves, (nop, string(tag))) || continue
        xs, ys = curves[(nop, string(tag))]
        Plots.plot!(pl, xs, ys; label = "", attrs(nop, tag)...)
    end
    # the legend: per curve an empty (NaN) series, which draws nothing but its entry
    for nop in nops, tag in listed
        haskey(curves, (nop, string(tag))) || continue
        xs, ys = curves[(nop, string(tag))]
        p   = _conv_rate(xs, ys)
        lab = LaTeXStrings.latexstring(string("\\mathbb{P}_", nop, "\\ \\mathrm{", _conv_mathrm(label(tag)), "}",
                                              isfinite(p) ? string("\\ (p=", round(p; digits = 2), ")") : ""))
        Plots.plot!(pl, [NaN], [NaN]; label = lab, attrs(nop, tag)...)
    end

    # Reference slopes, for the lowest and the highest order on the figure (one
    # when there is one order): an error ∝ h^(N+1) falls as #DOFs^(-(N+1)/2).
    # They BRACKET the data rather than run through it: the shallow one is drawn
    # ABOVE every curve and the steep one BELOW every curve, so neither crosses
    # the lines it is there to measure. For y = A x^s the offset A is read off
    # the data itself — the largest y_i x_i^(-s) puts the line through the
    # topmost point, the smallest through the bottom one — and a factor of two
    # moves it clear. The axes stay clipped to the data, so a steep guide simply
    # leaves the panel at its ends.
    if !isempty(allx)
        x1 = minimum(allx); x2 = maximum(allx)
        ymax = maximum(ally); ymin = minimum(ally)
        guides = minimum(nops) == maximum(nops) ? ((minimum(nops), true, :gray40),) :
                 ((minimum(nops), true, :gray40), (maximum(nops), false, :gray55))
        for (nop, above, col) in guides
            sl = -(nop + 1)/2
            r  = [ally[i]*allx[i]^(-sl) for i in eachindex(allx)]
            A  = above ? 2.0*maximum(r) : 0.5*minimum(r)
            Plots.plot!(pl, [x1, x2], [A*xi^sl for xi in (x1, x2)];
                        line = (col, 1.2, :dot),
                        label = LaTeXStrings.latexstring(string("\\mathcal{O}(h^{", nop + 1, "})")))
        end

        # Decades on the error axis, and the resolutions actually run on the
        # DOF axis: 10^{-4.5} is not a number anyone reports. Only the decades
        # inside the limits, which span more than one, so there is always one:
        # GR would draw the others' labels outside the frame.
        ylo = 0.2*ymin; yhi = 5.0*ymax
        lo = ceil(Int, log10(ylo)); hi = floor(Int, log10(yhi))
        decs = collect(lo:hi)
        stp  = max(1, cld(length(decs), 7))
        yt   = [10.0^e for e in decs[1:stp:end]]
        xu   = sort(unique(allx))
        idx  = length(xu) <= 4 ? eachindex(xu) :
               unique(round.(Int, range(1, length(xu); length = 4)))
        xt   = xu[idx]
        Plots.plot!(pl; xlims = (0.75*x1, 1.35*x2), ylims = (ylo, yhi),
                    xticks = (xt, [LaTeXStrings.latexstring(@sprintf("%d", round(Int, x))) for x in xt]),
                    yticks = (yt, [LaTeXStrings.latexstring(string("10^{", e, "}"))
                                   for e in decs[1:stp:end]]))
    end
    return pl
end

#---------------------------------------------------------------------------------
# The figure: each panel on its own at publication size in
# <prefix>_<file>-it<iout>.png, and all of them side by side, each at that same
# size and with its own legend (the rates in it are that panel's), in
# <prefix>-it<iout>.png. `panels` is a list of (field, ylabel, title, file); the
# keywords are those of _conv_panel. Returns the side-by-side file, or nothing
# when no (order, tag) has two points to join yet.
#---------------------------------------------------------------------------------
const CONV_PANEL_SIZE = (1180, 780)

function _conv_figure(rows, panels, OUTPUT_DIR, iout; prefix = "convergence", kw...)
    any(count(r -> r.nop == n && r.visc == v, rows) >= 2
        for n in unique(r.nop for r in rows), v in unique(r.visc for r in rows)) || return nothing
    pls = Plots.Plot[]
    for (fld, ylab, title, fname) in panels
        pl = _conv_panel(rows, fld, ylab, title; kw...)
        _conv_savefig(Plots.plot(pl; size = CONV_PANEL_SIZE, show = false),
                      string(OUTPUT_DIR, "/", prefix, "_", fname, "-it", iout, ".png"))
        push!(pls, pl)
    end
    # Wider margins than a panel's own: Plots' GR backend sizes the room for the
    # title, the axis labels and the x tick labels by the figure's HEIGHT, but GR
    # scales text by its longer side, so on a canvas this wide they would be
    # clipped.
    plt = Plots.plot(pls...; layout = (1, length(pls)),
                     size = (CONV_PANEL_SIZE[1]*length(pls), CONV_PANEL_SIZE[2]),
                     left_margin = 22Plots.mm, bottom_margin = 20Plots.mm,
                     top_margin = 13Plots.mm, right_margin = 6Plots.mm, show = false)
    f = string(OUTPUT_DIR, "/", prefix, "-it", iout, ".png")
    _conv_savefig(plt, f)
    return f
end
