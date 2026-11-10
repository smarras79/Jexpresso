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
#
# It reads problems/MHD/alfvenWaveBabbarRanocha2026/errors/*.dat (what every
# run writes) and uses the case's own plotting code, so the figure is the one
# a run draws. It is written as convergence-it0.png (convergence_nop3-5-7-it0.png
# for a subset) in output/MHD/alfvenWaveBabbarRanocha2026/output unless --out
# says otherwise.
#
#   --visc=TAGS   only these records, by the tag of their file name: dsgs (the
#                 deck as committed), galerkin, dsgs_perstage, ... (comma list)
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
    out   = joinpath(@__DIR__, "..", "output", "MHD", "alfvenWaveBabbarRanocha2026", "output")
    for a in args
        if startswith(a, "--t=")
            t = parse(Float64, split(a, '='; limit = 2)[2])
        elseif startswith(a, "--out=")
            out = split(a, '='; limit = 2)[2]
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
    _aw_plot(rows, out, 0; name = name) === nothing &&
        println(" # no figure: every order/tag has fewer than two meshes")
    return nothing
end

main(ARGS)
