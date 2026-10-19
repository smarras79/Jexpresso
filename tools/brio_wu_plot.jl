#!/usr/bin/env julia
#
# Redraw the Brio-Wu comparison figures from the curves the runs stored,
# without running the case: at the end of tools/brio_wu_order_scan.sh (whose
# orders may run in parallel, so that no single run sees the whole sweep), or
# to redraw a subset.
#
#   julia --project=. tools/brio_wu_plot.jl                  # everything stored
#   julia --project=. tools/brio_wu_plot.jl 4,7              # N = 4 and 7 alone
#   julia --project=. tools/brio_wu_plot.jl --dofs=600       # one resolution
#   julia --project=. tools/brio_wu_plot.jl --fmt=pdf --out=figs
#
# It reads problems/MHD/brioWu1d/curves/*.dat and draws with the case's own
# plotting code (user_plot.jl), so the figures are the ones a run draws, in
# output/MHD/brioWu1d/output unless --out says otherwise:
#
#   density_nu_dof<M>-it<n>   the density of every order at ~M DOFs with the
#                             zoom boxes, and their DynSGS coefficient ν: the
#                             layout of Dao & Nazarov (2024), Fig. 2
#   density_nu-it<n>          the same at the finest resolution
#   density_dof<M>-it<n>, density-it<n>         the density alone
#   convergence-it<n>, convergence_smooth-it<n> error against 1/#DOFs
#
#   --form=F   the method to draw: elem (the deck: element ν, Δ_K/(k+1)),
#              elemH (Δ_K/k), nodal or nodalH; default the form of the newest
#              stored curve
#   --t=T      the final time; default the latest in the store
#   --fmt=F    svg (default), pdf or png
using Plots, LaTeXStrings, Printf

const CASE = joinpath(@__DIR__, "..", "problems", "MHD", "brioWu1d")

# What src/io/plotting/jeplots.jl gives the case's hook during a run.
function _savefig_silent(plt, f)
    mkpath(dirname(f))
    Plots.savefig(plt, f)
    # the case writes through a temporary name and renames; report the destination
    println(" # wrote ", replace(String(f), r"\.tmp\d+(?=\.)" => ""))
end
function _sparse_marker_idx(n::Int, inputs)
    nmark = get(inputs, :plot_markers, 20)
    (nmark <= 0 || n <= nmark) && return nmark <= 0 ? (1:0) : (1:n)
    stride = max(1, n ÷ nmark)
    return 1:stride:n
end

include(joinpath(CASE, "user_analytic.jl"))   # the reference, reference_hll.dat
include(joinpath(CASE, "user_plot.jl"))

function main(args)
    nops = Int[]; dofs = Int[]
    form = nothing; t = nothing
    out  = joinpath(@__DIR__, "..", "output", "MHD", "brioWu1d", "output")
    inputs = Dict{Symbol,Any}(:plot_markers => 20)
    for a in args
        k, v = occursin('=', a) ? split(a, '='; limit = 2) : (a, "")
        if k == "--t"
            t = parse(Float64, v)
        elseif k == "--out"
            out = v
        elseif k == "--form"
            form = v
        elseif k == "--fmt"
            inputs[:plot_format] = v
        elseif k == "--dofs"
            append!(dofs, parse.(Int, filter(!isempty, split(v, r"[,\s]+"))))
        elseif startswith(a, "--")
            error("unknown option $a")
        else
            append!(nops, parse.(Int, filter(!isempty, split(a, r"[,\s]+"))))
        end
    end

    isdir(BW_CURVE_DIR) || error("no curve store at $(BW_CURVE_DIR): run the case first")
    # the final times present in the store; without --t, the latest
    times = Float64[]
    for f in readdir(BW_CURVE_DIR)
        (startswith(f, "nop") && endswith(f, ".dat")) || continue
        for line in eachline(joinpath(BW_CURVE_DIR, f))
            startswith(line, "#") || break
            m = match(r"\bt=(\S+)", line)
            m === nothing && continue
            v = tryparse(Float64, m.captures[1])
            v === nothing || push!(times, v)
        end
    end
    isempty(times) && error("no stored curves under $(BW_CURVE_DIR)")
    t === nothing && (t = maximum(times))

    every = _bw_load_curves(t)
    isempty(every) && error("no stored curves at t = $t (times present: $(sort(unique(times))))")
    form === nothing && (form = every[argmax([c.mtime for c in every])].form)
    stored = filter(c -> c.form == form, every)
    other  = sort(unique(c.form for c in every if c.form != form))
    isempty(other) || println(" # ", length(every) - length(stored), " stored curve(s) of another method (",
                              join(other, ", "), ") are not drawn; --form= selects one")
    isempty(nops) || filter!(c -> c.nop in nops, stored)
    isempty(dofs) || filter!(c -> _bw_bucket(c.ndofs) in _bw_bucket.(dofs), stored)
    isempty(stored) && error("no stored curves of form $form for these orders/resolutions")

    tab  = _bw_read_reference()
    xs   = collect(tab[:, 1])
    href = abs(t - BW_REF_T) < 1.0e-6 ? collect(tab[:, 2]) : nothing
    href === nothing && println(" # t = $t: the reference is for t = $(BW_REF_T); figures without it, no errors")
    iout = maximum(c.iout for c in stored)
    iout > 0 || (iout = 5)    # curves stored before the output index was: the deck's last output
    println(" # Brio-Wu, t = ", t, ", form ", form, ": ",
            join([string("N = ", c.nop, " (", c.ndofs, " DOFs)") for c in sort(stored, by = c -> (c.ndofs, c.nop))], ", "))
    _bw_draw_store(stored, xs, href, inputs, out, iout)
    return nothing
end

main(ARGS)
