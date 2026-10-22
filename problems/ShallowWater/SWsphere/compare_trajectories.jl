#!/usr/bin/env julia
#-------------------------------------------------------------------------------
# compare_trajectories.jl — how different are the SWsphere trajectories?
#
# Reads the VTK output of every <dir>/trajectory_X/ written by
# run_trajectories.jl and prints, for each trajectory and each output time,
# max|f_X - f_ref| and rms(f_X - f_ref) against a reference trajectory, for
# the height h and the zonal wind u.
#
#   julia --project=. problems/ShallowWater/SWsphere/compare_trajectories.jl [dir] [ref]
#
#   dir   parent directory of the trajectory_X/   (default: ./trajectories)
#   ref   reference trajectory number              (default: the one with
#                                                   ĥ = 120 m, φ₂ = 45°, from
#                                                   trajectories.csv; else the first)
#
# dir is the OUTDIR the sweep was submitted with; the submit script printed it
# ("output: …"), and every trajectory log has it on its "Output written to" line.
#
# WHAT TO EXPECT AT t = 0 (file sphere_grid_ho). The trajectories differ only
# in the height bump h' = ĥ cos φ exp[-(λ/α)²] exp[-((φ₂-φ)/β)²], so:
#   * h:  of the order of the difference in ĥ·cos φ₂ (e.g. 40 vs 120 m at
#         45° ≈ 57 m) between amplitudes, and up to ~ĥ between latitudes whose
#         bumps do not overlap. Exactly 0 means the perturbation did NOT reach
#         the initial condition.
#   * u:  identical (0) at t = 0 for every trajectory — the Galewsky
#         perturbation is in h only; the wind only picks it up as the flow
#         adjusts. Differences in u then GROW with time if the jet is unstable
#         at this resolution and damping, and stay small if it is not.
#
# Needs no Jexpresso: it reads the files with tools/vtu_reader.jl.
#-------------------------------------------------------------------------------
using Printf

include(joinpath(@__DIR__, "..", "..", "..", "tools", "vtu_reader.jl"))
using .VTUReader

dir    = abspath(length(ARGS) >= 1 ? ARGS[1] : "./trajectories")
ref_arg = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : nothing
FIELDS = ("h", "u")

# Trajectory parameters, from the table run_trajectories.jl wrote.
params  = Dict{Int,String}()
galewsky = nothing          # the published case, ĥ = 120 m and φ₂ = 45°
csv = joinpath(dir, "trajectories.csv")
if isfile(csv)
    for line in Iterators.drop(eachline(csv), 1)
        c = split(line, ',')
        id, h, φ = parse(Int, c[1]), parse(Float64, c[5]), parse(Float64, c[6])
        params[id] = @sprintf("ĥ=%.3f m, φ₂=%.4f°", h, φ)
        isapprox(h, 120.0) && isapprox(φ, 45.0) && (global galewsky = id)
    end
end

isdir(dir) || error("""$dir does not exist. Pass the OUTDIR the sweep wrote to as the first argument;
       find it with:   grep -h "Output written to" logs/SWtraj.*.out | sort -u""")
trajdirs = filter(d -> occursin(r"^trajectory_\d+$", d) && isdir(joinpath(dir, d)), readdir(dir))
ids = sort!([parse(Int, d[12:end]) for d in trajdirs])
isempty(ids) && error("no trajectory_X directories in $dir")
ref_id = something(ref_arg, galewsky, first(ids))
ref_id in ids || error("reference trajectory_$ref_id not found in $dir (have: $(join(ids, ", ")))")

# Output files of one trajectory, in time order: the initial condition first.
function output_files(tdir)
    fs  = readdir(tdir)
    ext = any(endswith(".pvtu"), fs) ? ".pvtu" : ".vtu"
    ic  = filter(f -> f == "sphere_grid_ho" * ext, fs)
    out = sort!(filter(f -> occursin(Regex("^sphere_\\d{4}\\" * ext * "\$"), f), fs))
    return [joinpath(tdir, f) for f in vcat(ic, out)]
end

ref_files = output_files(joinpath(dir, "trajectory_$ref_id"))
isempty(ref_files) && error("no sphere_*.vtu/.pvtu output in trajectory_$ref_id")
ref = [read_vtk_grid(f) for f in ref_files]

println(" # reference: trajectory_", ref_id, " ", get(params, ref_id, ""), " — ", length(ref_files), " output times")
for id in ids
    id == ref_id && continue
    files = output_files(joinpath(dir, "trajectory_$id"))
    println()
    println(" # trajectory_", id, " ", get(params, id, ""))
    @printf(" #   %-22s %10s", "file", "t [h]")
    for f in FIELDS
        @printf("   %12s %12s", "max|Δ$f|", "rms Δ$f")
    end
    println()
    for (k, fpath) in enumerate(files)
        k <= length(ref) || break
        g = read_vtk_grid(fpath)
        length(g.x) == length(ref[k].x) ||
            (println(" #   ", basename(fpath), ": different node count, skipped"); continue)
        t = g.time === nothing ? NaN : g.time / 3600
        @printf(" #   %-22s %10.2f", basename(fpath), t)
        for f in FIELDS
            if haskey(g.data, f) && haskey(ref[k].data, f)
                d = g.data[f] .- ref[k].data[f]
                @printf("   %12.4e %12.4e", maximum(abs, d), sqrt(sum(abs2, d) / length(d)))
            else
                @printf("   %12s %12s", "n/a", "n/a")
            end
        end
        println()
    end
end
