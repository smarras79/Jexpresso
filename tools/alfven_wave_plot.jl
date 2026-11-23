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
#
# It reads problems/MHD/alfvenWaveBabbarRanocha2026/errors/*.dat (what every
# run writes) and uses the case's own plotting code, so the figure is the one
# a run draws. Written in output/MHD/alfvenWaveBabbarRanocha2026/output unless
# --out says otherwise (convergence_nop3-5-7... for a subset of the orders), all
# -it0.png:
#
#   convergence         the deck's DSGS and Galerkin alone, the seven error
#                       panels in a grid, one dotted h^(N+1) guide per order
#   convergence_Bperp, _Bperp_Linf, _uperp, _uperp_Linf, _rho, _L2, _Linf
#                       each of those panels alone, with its legend
#   convergence_dsgs_nocut   the deck's DSGS and DSGS without the cutoff alone,
#                       the same way
#   convergence_all     every variant stored (or selected with --visc), the
#                       seven panels side by side, the rate of each curve in
#                       the legend
#   convergence_rmax    the largest residual r of every DSGS run, before the
#                       cutoff, against the cutoff c: a point below c is a run
#                       in which μ was exactly zero throughout
#   convergence_mumax   the largest DSGS coefficient μ of every run; the legend
#                       says on how many meshes it was exactly zero
#   convergence_paper   the error of B⊥ | max r | max μ, side by side
#
# A comparison is drawn only when both of its variants are stored; a
# Galerkin-only selection has no rmax, mumax or paper figure.
#
#   --visc=TAGS   only these records, by the tag of their file name: dsgs (the
#                 deck as committed), dsgs_nocut (cutoff 0), galerkin, dsgs_perstage, ...
#                 (comma list)
#   --panels[=P]  instead of the figures above, the two comparisons with only
#                 the panels P (comma list of Bperp, Bperp_Linf, uperp,
#                 uperp_Linf, rho, L2, Linf; without "=P": Bperp,uperp,rho, the
#                 L² errors of B⊥, u⊥ and ρ) in one row and one legend for all,
#                 as convergence_panels_Bperp-uperp-rho-it0.png and
#                 convergence_dsgs_nocut_panels_Bperp-uperp-rho-it0.png
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
            panels = ["Bperp", "uperp", "rho"]
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
    name = isempty(nops) ? "convergence" : string("convergence_nop", join(sort(unique(nops)), "-"))
    if panels !== nothing
        keys = [p[4] for p in AW_PANELS]
        sel  = map(panels) do k
            i = findfirst(x -> lowercase(x) == lowercase(k), keys)
            i === nothing && error("unknown panel \"$k\": the panels are $(join(keys, ", "))")
            AW_PANELS[i]
        end
        isempty(_aw_compare(rows, out, 0; name = name, panels = sel)) &&
            println(" # no figure: no comparison has both of its variants on two meshes")
        return nothing
    end
    _aw_plot(rows, out, 0; name = name) === nothing &&
        println(" # no figure: every order/tag has fewer than two meshes")
    return nothing
end

main(ARGS)
