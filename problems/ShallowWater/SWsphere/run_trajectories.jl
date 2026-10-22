#!/usr/bin/env julia
#-------------------------------------------------------------------------------
# run_trajectories.jl — ensemble of SWsphere (Galewsky jet) runs for ROM training.
#
# Each TRAJECTORY is one SWsphere run, identical to the shipped deck
# (user_inputs.jl) except for the height perturbation and the final time,
# which is 3 days (:tend = 259200 s) by default instead of the deck's 6
# (change it with --tend <seconds>). Output is written every 6 hours
# (--output-interval), i.e. :ndiagnostics_outputs = tend / 6 h = 12 dumps.
# The perturbation is
#
#     h'(λ,φ) = ĥ cos φ exp[-(λ/α)²] exp[-((φ₂-φ)/β)²]
#
# The amplitude ĥ and latitude φ₂ are each varied over a symmetric family
# around the published case (ĥ = 120 m, φ₂ = 45°), and EVERY combination of
# the two is run. By default each family has 5 members within ±5%:
#
#   ĥ  = 120 m × (1 + [-5, -2.5, 0, +2.5, +5] %) = 114, 117, 120, 123, 126 m
#   φ₂ =  45°  × (1 + [-5, -2.5, 0, +2.5, +5] %) = 42.75, 43.875, 45, 46.125, 47.25°
#
# i.e. 5 × 5 = 25 trajectories, numbered with the latitude varying fastest:
#
#   trajectory_1  … trajectory_5    ĥ = 114 m, φ₂ = 42.75° … 47.25°
#   trajectory_6  … trajectory_10   ĥ = 117 m, φ₂ = 42.75° … 47.25°
#   …
#   trajectory_21 … trajectory_25   ĥ = 126 m, φ₂ = 42.75° … 47.25°
#
# (trajectory_13 is the reference Galewsky run.)
# Trajectory X writes everything (VTK, POD, the deck copy,
# input_overrides.txt) to
#
#     <outdir>/trajectory_X/
#
# and <outdir>/trajectories.csv lists X → (ĥ, φ₂) for the whole ensemble.
# All runs share ONE Julia session, so compilation is paid only once.
#
# USAGE (from the repository root):
#
#   julia --project=. problems/ShallowWater/SWsphere/run_trajectories.jl
#
#   # parallel: every rank runs every trajectory together
#   mpiexec -n 4 julia --project=. problems/ShallowWater/SWsphere/run_trajectories.jl
#
#   # wider/denser grid (7 × 5 = 35 runs), shorter runs:
#   julia --project=. problems/ShallowWater/SWsphere/run_trajectories.jl \
#         --amp-variation 0.10 --namp 7 --lat-variation 0.05 --nlat 5 --tend 86400
#
#   # exactly the values you want (3 × 2 = 6 runs):
#   julia --project=. problems/ShallowWater/SWsphere/run_trajectories.jl \
#         --amplitudes 100,120,140 --latitudes 40,50
#
#   # change ANY other deck input, in every trajectory (repeat --set as needed;
#   # the value is a Julia expression, the key a user_inputs.jl key):
#   julia --project=. problems/ShallowWater/SWsphere/run_trajectories.jl \
#         --set μ=1.0e5 --set filter_alpha=0.02 --set lvisc=false
#
# Latitudes are in DEGREES. The jet occupies φ₀ = π/7 ≈ 25.7° to
# φ₁ = π/2 - π/7 ≈ 64.3°.
#
# A trajectory whose directory already holds a COMPLETED marker is skipped,
# so an interrupted sweep can simply be relaunched (use --force to redo).
# Use --dry-run to print the trajectory table without running anything.
#-------------------------------------------------------------------------------
using ArgParse
using Printf

function parse_sweep_args()
    s = ArgParseSettings(description = "Run SWsphere for a sweep of Galewsky perturbations.")
    @add_arg_table! s begin
        "--hhat0"
            help     = "reference perturbation amplitude ĥ [m]"
            arg_type = Float64
            default  = 120.0
        "--lat0"
            help     = "reference perturbation latitude φ₂ [deg]"
            arg_type = Float64
            default  = 45.0
        "--amp-variation"
            help     = "relative half-width of the amplitude family (0.05 means plus or minus 5 percent)"
            arg_type = Float64
            default  = 0.05
        "--lat-variation"
            help     = "relative half-width of the latitude family (0.05 means plus or minus 5 percent)"
            arg_type = Float64
            default  = 0.05
        "--namp"
            help     = "number of amplitude values (evenly spaced over ±variation)"
            arg_type = Int
            default  = 5
        "--nlat"
            help     = "number of latitude values (evenly spaced over ±variation)"
            arg_type = Int
            default  = 5
        "--amplitudes"
            help    = "explicit amplitudes ĥ [m], e.g. \"100,120,140\"; replaces --hhat0/--amp-variation/--namp"
            default = ""
        "--latitudes"
            help    = "explicit latitudes φ₂ [deg], e.g. \"40,45,50\"; replaces --lat0/--lat-variation/--nlat"
            default = ""
        "--set"
            help   = "override any deck input in every trajectory, KEY=VALUE with VALUE a Julia expression, e.g. --set μ=1.0e5 (repeatable)"
            action = :append_arg
        "--outdir", "-o"
            help    = "parent directory of the trajectory_X directories"
            default = "./trajectories"
        "--tend"
            help    = "final time :tend of every trajectory [s] (default: 3 days)"
            arg_type = Float64
            default  = 6*24*3600.0
        "--output-interval"
            help     = "time between outputs [s]; sets :ndiagnostics_outputs = tend/interval (default: 6 h)"
            arg_type = Float64
            default  = 6*3600.0
        "--first"
            help    = "number of the first trajectory (directory trajectory_<first>)"
            arg_type = Int
            default  = 1
        "--only"
            help    = "run only these trajectory numbers, e.g. \"3,7\" (numbering unchanged)"
            default = ""
        "--force"
            help   = "re-run trajectories that are already COMPLETED"
            action = :store_true
        "--stop-on-error"
            help   = "abort the sweep on the first failed trajectory (default: log it and go on)"
            action = :store_true
        "--dry-run"
            help   = "print the trajectory table and exit"
            action = :store_true
    end
    return parse_args(s)
end

const EQS  = "ShallowWater"
const CASE = "SWsphere"

args   = parse_sweep_args()
outdir = abspath(args["outdir"])
only   = isempty(strip(args["only"])) ? nothing : Set(parse.(Int, split(args["only"], ',')))

family(x0, δ, n) = n == 1 ? [x0] : collect(x0 .* (1 .+ range(-δ, δ; length = n)))

parse_list(str) = parse.(Float64, split(str, ','))

hhat0, lat0 = args["hhat0"], args["lat0"]
amps = isempty(strip(args["amplitudes"])) ? family(hhat0, args["amp-variation"], args["namp"]) :
                                            parse_list(args["amplitudes"])
lats = isempty(strip(args["latitudes"]))  ? family(lat0,  args["lat-variation"], args["nlat"]) :
                                            parse_list(args["latitudes"])

# --set KEY=VALUE: parsed now (so a typo fails before anything is submitted),
# evaluated after Jexpresso is loaded, so values like CarpenterKennedy2N54() work.
const SWEEP_KEYS = (:galewsky_hhat, :galewsky_phi2)
user_sets = map(args["set"]) do kv
    i = findfirst('=', kv)
    i === nothing && error("--set expects KEY=VALUE, got \"$kv\"")
    key = Symbol(lstrip(strip(kv[1:prevind(kv, i)]), ':'))
    key in SWEEP_KEYS && error("--set $key: the sweep sets it; use --amplitudes/--latitudes instead")
    key => Meta.parse(strip(kv[nextind(kv, i):end]))
end
# Every (amplitude, latitude) pair; the latitude varies fastest.
params = [(hhat = h, lat_deg = φ, iamp = ia, ilat = il)
          for (ia, h) in enumerate(amps) for (il, φ) in enumerate(lats)]
trajectories = [(; id = args["first"] + i - 1, p...) for (i, p) in enumerate(params)]

let n = args["tend"] / args["output-interval"]
    isapprox(n, round(n)) ||
        @warn @sprintf("tend = %g s is not a multiple of --output-interval = %g s; outputs will be every %g s",
                       args["tend"], args["output-interval"], args["tend"] / max(1, round(Int, n)))
end

φ0_deg, φ1_deg = rad2deg(π/7), rad2deg(π/2 - π/7)
for t in trajectories
    (φ0_deg < t.lat_deg < φ1_deg) ||
        @warn @sprintf("trajectory_%d: φ₂ = %.2f° lies outside the jet (%.1f°–%.1f°)",
                       t.id, t.lat_deg, φ0_deg, φ1_deg)
end

function print_table(io, trajectories)
    println(io, "trajectory,directory,amp_index,lat_index,galewsky_hhat_m,galewsky_phi2_deg,galewsky_phi2_rad")
    for t in trajectories
        @printf(io, "%d,trajectory_%d,%d,%d,%.10g,%.10g,%.16g\n",
                t.id, t.id, t.iamp, t.ilat, t.hhat, t.lat_deg, deg2rad(t.lat_deg))
    end
end

if args["dry-run"]
    print_table(stdout, trajectories)
    for (k, ex) in user_sets
        println(stderr, " # --set :", k, " => ", ex, "  (every trajectory)")
    end
    exit(0)
end

using Jexpresso
using MPI
MPI.Initialized() || MPI.Init()
const RANK = MPI.Comm_rank(MPI.COMM_WORLD)

user_overrides = Dict{Symbol,Any}(k => Core.eval(Jexpresso, ex) for (k, ex) in user_sets)

if RANK == 0
    mkpath(outdir)
    # Write-then-rename: under a Slurm job array every task writes this same
    # file at the same moment, and a rename is atomic where a write is not.
    tmp = joinpath(outdir, "trajectories.csv.$(getpid()).$(gethostname()).tmp")
    open(io -> print_table(io, trajectories), tmp, "w")
    mv(tmp, joinpath(outdir, "trajectories.csv"); force = true)
    println(" # SWsphere sweep: ", length(trajectories), " trajectories → ", outdir)
end
MPI.Barrier(MPI.COMM_WORLD)

failed = Int[]
for t in trajectories
    only !== nothing && !(t.id in only) && continue

    tdir   = joinpath(outdir, "trajectory_$(t.id)")
    marker = joinpath(tdir, "COMPLETED")
    if isfile(marker) && !args["force"]
        RANK == 0 && println(" # trajectory_$(t.id): already COMPLETED, skipping")
        continue
    end

    overrides = Dict{Symbol,Any}(
        :lgalewsky_perturbation => true,
        :galewsky_hhat          => t.hhat,
        :galewsky_phi2          => deg2rad(t.lat_deg),
        :loverwrite_output      => true,
    )
    overrides[:tend] = args["tend"]
    # The sphere solver writes every tend/:ndiagnostics_outputs seconds.
    overrides[:ndiagnostics_outputs] = max(1, round(Int, args["tend"] / args["output-interval"]))
    # --set wins over everything except the swept ĥ, φ₂ (rejected above).
    merge!(overrides, user_overrides)

    if RANK == 0
        @printf(" # ===== trajectory_%d of %d: ĥ = %g m, φ₂ = %g° =====\n",
                t.id, trajectories[end].id, t.hhat, t.lat_deg)
        flush(stdout)
    end

    t0 = time()
    try
        Jexpresso.run_case(EQS, CASE; overrides = overrides, output_dir = tdir)
        MPI.Barrier(MPI.COMM_WORLD)
        if RANK == 0
            open(marker, "w") do io
                @printf(io, "galewsky_hhat = %.10g\ngalewsky_phi2_deg = %.10g\nwall_seconds = %.1f\n",
                        t.hhat, t.lat_deg, time() - t0)
            end
        end
    catch err
        push!(failed, t.id)
        if RANK == 0
            @error "trajectory_$(t.id) failed" exception = (err, catch_backtrace())
            mkpath(tdir)
            open(io -> showerror(io, err), joinpath(tdir, "FAILED"), "w")
        end
        args["stop-on-error"] && rethrow()
    end
end

if RANK == 0
    if isempty(failed)
        println(" # SWsphere sweep finished: all trajectories completed.")
    else
        println(" # SWsphere sweep finished; FAILED trajectories: ", join(failed, ", "))
    end
end
