#!/usr/bin/env julia
#
# Redraw the smooth-vortex convergence figures from the errors already stored,
# without running a simulation — to plot a SUBSET of the orders on their own
# axes (the P1-vs-P3 comparison of Dao & Nazarov (2022), Fig. 1, say), or to
# restyle a figure after a sweep that took hours.
#
#   julia --project=. tools/smooth_vortex_plot.jl              # every order
#   julia --project=. tools/smooth_vortex_plot.jl 1,3          # P1 and P3 alone
#   julia --project=. tools/smooth_vortex_plot.jl 1,3 --t=1.0 --out=figs
#
# It reads problems/MHD/smoothVortex/errors/*.dat (what every run writes) and
# uses the case's own plotting code, so the figures are identical to the ones
# a run produces. Files are written as convergence_<visc>[_nop1-3]-it0.png.
using Plots, LaTeXStrings, Printf

const CASE = joinpath(@__DIR__, "..", "problems", "MHD", "smoothVortex")

# The case's hook calls this; here it is the whole of the output side.
function _savefig_silent(plt, f)
    mkpath(dirname(f))
    Plots.savefig(plt, f)
    # the case writes through a temporary name and renames (atomic, so that
    # concurrent cases of a sweep cannot tear a figure); report the destination
    println(" # wrote ", replace(String(f), r"\.tmp\d+(?=\.)" => ""))
end

include(joinpath(CASE, "user_plot.jl"))

function main(args)
    nops = Int[]
    t    = nothing
    out  = joinpath(@__DIR__, "..", "output", "MHD", "smoothVortex", "output")
    for a in args
        if startswith(a, "--t=")
            t = parse(Float64, split(a, '=')[2])
        elseif startswith(a, "--out=")
            out = split(a, '='; limit = 2)[2]
        elseif startswith(a, "--")
            error("unknown option $a")
        else
            append!(nops, parse.(Int, split(a, r"[,\s]+")))
        end
    end

    isdir(SV_ERR_DIR) || error("no error store at $(SV_ERR_DIR) — run the case first")

    # Final times present in the store; without --t, take the largest.
    times = Float64[]
    for f in readdir(SV_ERR_DIR)
        endswith(f, ".dat") || continue
        for line in eachline(joinpath(SV_ERR_DIR, f))
            startswith(line, "#") || continue
            for tok in split(line)
                startswith(tok, "t=") || continue
                v = tryparse(Float64, split(tok, '=')[2])
                v === nothing || push!(times, v)
            end
        end
    end
    isempty(times) && error("no stored errors under $(SV_ERR_DIR)")
    t === nothing && (t = maximum(times))

    rows = _sv_load_errors(t)
    isempty(rows) && error("no stored errors at t = $t (present: $(sort(unique(times))))")
    _sv_report(rows)
    if isempty(nops)
        _sv_plot(rows, out, 0)
    else
        _sv_plot(rows, out, 0; only = sort(unique(nops)),
                 suffix = string("_nop", join(sort(unique(nops)), "-")))
    end
    return nothing
end

main(ARGS)
