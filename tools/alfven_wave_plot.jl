#!/usr/bin/env julia
#
# Redraw the smooth-Alfvén-wave convergence figure (Babbar & Ranocha 2026,
# Fig. 13) and its table from the errors already stored, without running a
# simulation: at the end of tools/alfven_wave_scan.sh, after a sweep split
# over several sessions or machines, or to plot a subset of the orders.
#
#   julia --project=. tools/alfven_wave_plot.jl                  # everything stored
#   julia --project=. tools/alfven_wave_plot.jl 3,5,7            # N = 3, 5, 7 alone
#   julia --project=. tools/alfven_wave_plot.jl --visc=dsgs,galerkin --out=figs
#   julia --project=. tools/alfven_wave_plot.jl --t=1.0          # another final time
#   julia --project=. tools/alfven_wave_plot.jl --panels         # L² of B⊥, u⊥, ρ, one legend
#   julia --project=. tools/alfven_wave_plot.jl --panels=L1      # the same in L¹
#
# It reads problems/MHD/alfvenWaveBabbarRanocha2026/errors/*.dat (what every
# run writes) and uses the case's own plotting code, so the figure is the one
# a run draws. Written in output/MHD/alfvenWaveBabbarRanocha2026/output unless
# --out says otherwise (alfven-convergence_nop3-5-7... for a subset of the orders), all
# -it0.png:
#
#   alfven-convergence         the deck's DSGS and Galerkin alone, nine error
#                       panels (B⊥, u⊥, ρ; in L², L¹, L∞) in a grid, one dotted
#                       h^(N+1) guide per order
#   alfven-convergence_Bperp, _Bperp_L1, _Bperp_Linf, _uperp, ..., _rho, _rho_L1,
#   _rho_Linf           each of those panels alone, with its legend
#   alfven-convergence_dsgs_nocut   the deck's DSGS and DSGS without the cutoff alone,
#                       the same way
#   alfven-convergence_all     every variant stored (or selected with --visc), five
#                       panels side by side (B⊥ and u⊥ in L² and L∞, and ρ),
#                       the rate of each curve in the legend
#   alfven-convergence_rmax    the largest residual r of every DSGS run, before the
#                       cutoff, against the cutoff c: a point below c is a run
#                       in which μ was exactly zero throughout
#   alfven-convergence_mumax   the largest DSGS coefficient μ of every run; the legend
#                       says on how many meshes it was exactly zero
#   alfven-convergence_paper   the error of B⊥ | max r | max μ, side by side
#
# A comparison is drawn only when both of its variants are stored; a
# Galerkin-only selection has no rmax, mumax or paper figure.
#
#   --visc=TAGS   only these records, by the tag of their file name: dsgs (the
#                 deck as committed), dsgs_nocut (cutoff 0), galerkin, dsgs_perstage, ...
#                 (comma list)
#   --panels[=P]  instead of the figures above, the two comparisons with only
#                 the panels P, in a row (or a grid) with one legend for all:
#                 L2, L1 or Linf for B⊥, u⊥ and ρ in that norm (without "=P":
#                 L2), or panel keys (Bperp, Bperp_L1, Bperp_Linf, uperp, ...,
#                 rho, ..., and state, state_L1, state_Linf: the state q, which
#                 the other figures leave out), comma separated, mixed at will:
#                 --panels=L2,L1 gives B⊥, u⊥, ρ in L² over L¹. Written as
#                 alfven-convergence_panels_<P>-it0.png and
#                 alfven-convergence_dsgs_nocut_panels_<P>-it0.png
using Plots, LaTeXStrings, Printf, LinearAlgebra

const CASE = joinpath(@__DIR__, "..", "problems", "MHD", "alfvenWaveBabbarRanocha2026")

# The case's hook calls this; here it is the whole of the output side.
function _savefig_silent(plt, f)
    mkpath(dirname(f))
    Plots.savefig(plt, f)
    # the case writes through a temporary name and renames; report the destination
    println(" # wrote ", replace(String(f), r"\.tmp\d+(?=\.)" => ""))
end

include(joinpath(CASE, "user_plot.jl"))

function main(args)
    nops  = Int[]
    viscs = String[]
    t     = nothing
    panels = nothing
    out   = joinpath(@__DIR__, "..", "output", "MHD", "alfvenWaveBabbarRanocha2026", "output")
    for a in args
        if startswith(a, "--t=")
            t = parse(Float64, split(a, '='; limit = 2)[2])
        elseif startswith(a, "--out=")
            out = split(a, '='; limit = 2)[2]
        elseif a == "--panels"
            panels = ["L2"]
        elseif startswith(a, "--panels=")
            panels = String.(filter(!isempty, split(split(a, '='; limit = 2)[2], r"[,\s]+")))
        elseif startswith(a, "--visc=")
            append!(viscs, filter(!isempty, split(split(a, '='; limit = 2)[2], r"[,\s]+")))
        elseif startswith(a, "--")
            error("unknown option $a")
        else
            append!(nops, parse.(Int, filter(!isempty, split(a, r"[,\s]+"))))
        end
    end

    isdir(AW_ERR_DIR) || error("no error store at $(AW_ERR_DIR): run the case first")
    # the final times present in the store; without --t, the largest
    times = Float64[]
    for f in readdir(AW_ERR_DIR)
        (startswith(f, "nop") && endswith(f, ".dat")) || continue
        for line in eachline(joinpath(AW_ERR_DIR, f))
            startswith(line, "#") || continue
            for tok in split(line)
                startswith(tok, "t=") || continue
                v = tryparse(Float64, split(tok, '='; limit = 2)[2])
                v === nothing || push!(times, v)
            end
        end
    end
    isempty(times) && error("no stored errors under $(AW_ERR_DIR)")
    t === nothing && (t = maximum(times))

    rows = _aw_load(t)
    isempty(nops)  || filter!(r -> r.nop in nops, rows)
    isempty(viscs) || filter!(r -> r.visc in viscs, rows)
    isempty(rows) && error("no stored errors at t = $t for these orders/tags (times present: $(sort(unique(times))))")
    _aw_report(rows)
    name = isempty(nops) ? AW_FIG_NAME : string(AW_FIG_NAME, "_nop", join(sort(unique(nops)), "-"))
    if panels !== nothing
        # L2, L1, Linf: B⊥, u⊥ and ρ in that norm; else panel keys
        keys = [p[4] for p in AW_PANEL_MATRIX]
        sets = Dict(lowercase(k) => v for (k, v) in AW_PANEL_SETS)
        sel  = vcat(map(panels) do k
            haskey(sets, lowercase(k)) && return [AW_PANEL_MATRIX[findfirst(==(x), keys)] for x in sets[lowercase(k)]]
            i = findfirst(x -> lowercase(x) == lowercase(k), keys)
            i === nothing && error("unknown panel \"$k\": use L2, L1 or Linf (B⊥, u⊥ and ρ in ",
                                   "that norm) or panels among $(join(keys, ", "))")
            [AW_PANEL_MATRIX[i]]
        end...)
        isempty(_aw_compare(rows, out, 0; name = name, panels = sel, label = join(panels, "-"))) &&
            println(" # no figure: no comparison has both of its variants on two meshes")
        return nothing
    end
    _aw_plot(rows, out, 0; name = name) === nothing &&
        println(" # no figure: every order/tag has fewer than two meshes")
    return nothing
end

main(ARGS)
