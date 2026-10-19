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
#   julia --project=. tools/brio_wu_plot.jl --store=Cmin0.03   # the sweep at C_min = 0.03
#
# It reads problems/MHD/brioWu1d/curves/*.dat (curves/<store>/*.dat with
# --store; tools/brio_wu_order_scan.sh keeps one store per C_min, Cmin<value>)
# and draws with the case's own plotting code (user_plot.jl), so the figures
# are the ones a run draws, in output/MHD/brioWu1d/output (output/<store> with
# --store) unless --out says otherwise:
#
#   rho_By_dof<M>-it<n>       ρ (left) and B_y (right) of every order at ~M DOFs
#                             against the reference, with zoom boxes, as in Dao
#                             & Nazarov's Brio-Wu figure
#   rho_By-it<n>              the same at the finest resolution
#   rho_By_all-it<n>          EVERY order at EVERY resolution; colour = order,
#                             line type = resolution
#   nu_dof<M>-it<n>, nu-it<n>, nu_all-it<n>   the DynSGS coefficient ν alone,
#                             of the same curves
#   density_dof<M>-it<n>, density-it<n>         the density alone
#   convergence-it<n>, convergence_smooth-it<n> error against 1/#DOFs
#
#   --form=F   the method to draw: elem (the deck: element ν, Δ_K/(k+1)),
#              elemH (Δ_K/k), nodal or nodalH; default the form of the newest
#              stored curve
#   --t=T      the final time; default the latest in the store
#   --fmt=F    svg (default), pdf or png
#   --right=nu ρ and ν side by side instead of ρ and B_y (density_nu_*-it<n>,
#              the earlier layout; JEXPRESSO_BW_RIGHT=nu does the same)
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
    form = nothing; t = nothing; out = nothing
    inputs = Dict{Symbol,Any}(:plot_markers => 20)
    for a in args
        k, v = occursin('=', a) ? split(a, '='; limit = 2) : (a, "")
        if k == "--t"
            t = parse(Float64, v)
        elseif k == "--out"
            out = v
        elseif k == "--store"
            ENV["JEXPRESSO_BW_STORE"] = v      # read by _bw_curve_dir()
        elseif k == "--right"
            ENV["JEXPRESSO_BW_RIGHT"] = v      # read by _bw_draw_store
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

    store = strip(get(ENV, "JEXPRESSO_BW_STORE", ""))
    out === nothing && (out = joinpath(@__DIR__, "..", "output", "MHD", "brioWu1d", "output", store))
    isdir(_bw_curve_dir()) || error("no curve store at $(_bw_curve_dir()): run the case first")
    # the final times present in the store; without --t, the latest
    times = Float64[]
    for f in readdir(_bw_curve_dir())
        (startswith(f, "nop") && endswith(f, ".dat")) || continue
        for line in eachline(joinpath(_bw_curve_dir(), f))
            startswith(line, "#") || break
            m = match(r"\bt=(\S+)", line)
            m === nothing && continue
            v = tryparse(Float64, m.captures[1])
            v === nothing || push!(times, v)
        end
    end
    isempty(times) && error("no stored curves under $(_bw_curve_dir())")
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
    hBy  = abs(t - BW_REF_T) < 1.0e-6 ? collect(tab[:, 7]) : nothing   # x ρ u v w p By Bz
    href === nothing && println(" # t = $t: the reference is for t = $(BW_REF_T); figures without it, no errors")
    iout = maximum(c.iout for c in stored)
    iout > 0 || (iout = 5)    # curves stored before the output index was: the deck's last output
    println(" # Brio-Wu, t = ", t, ", form ", form, ": ",
            join([string("N = ", c.nop, " (", c.ndofs, " DOFs)") for c in sort(stored, by = c -> (c.ndofs, c.nop))], ", "))
    _bw_draw_store(stored, xs, href, inputs, out, iout; hBy = hBy)
    return nothing
end

main(ARGS)
