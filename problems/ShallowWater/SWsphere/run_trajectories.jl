#!/usr/bin/env julia
#-------------------------------------------------------------------------------
# run_trajectories.jl — ensemble of SWsphere (Galewsky jet) runs for ROM training.
#
# Each TRAJECTORY is one SWsphere run, identical to the shipped deck
# (user_inputs.jl) except for the height perturbation
#
#     h'(λ,φ) = ĥ cos φ exp[-(λ/α)²] exp[-((φ₂-φ)/β)²]
#
# Every trajectory runs to 20 days (--tend). Output is written every hour for
# the first 2 days (--early-interval, --early-until), when the wave grows and
# breaks, and every 6 hours from then to the end (--output-interval):
#
#   t = 1, 2, …, 48 h, then 54, 60, …, 480 h   → 48 + 72 = 120 outputs
#
# passed to the solver as :diagnostics_at_times.
#
# The amplitude ĥ and latitude φ₂ are each varied over an evenly spaced range
# around the published case (ĥ = 120 m, φ₂ = 45°), and EVERY combination of
# the two is run. By default:
#
#   ĥ  = 120 m ± 80 m  (--amp-halfwidth) = 40, 80, 120, 160, 200 m
#   φ₂ =  45°  ± 15°   (--lat-halfwidth) = 30, 37.5, 45, 52.5, 60°
#
# By default (--design subset) only 10 of the 5 × 5 = 25 combinations are
# run: two "Latin" diagonals of the grid, so every amplitude and every
# latitude is used exactly twice and the two parameters are not tied to each
# other. They are numbered 1–10:
#
#   trajectory   1   2    3    4    5   6   7    8    9   10
#   ĥ [m]       40  40   80   80  120 120 160  160  200  200
#   φ₂ [°]      30  45 37.5 52.5   45  60  30 52.5 37.5   60
#
# (trajectory_5 is the reference Galewsky run.) --design grid runs all 25,
# numbered 1–25 with the latitude varying fastest; that numbering differs
# from the subset's, so give it its own --outdir. Narrower ranges (e.g. ±5%)
# give trajectories that are hard to tell apart after the wave breaks: the
# jet's own fastest-growing mode takes over whatever the details of the bump.
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
#         --amp-halfwidth 90 --namp 7 --lat-halfwidth 10 --nlat 5 --tend 86400
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
        "--amp-halfwidth"
            help     = "amplitudes span hhat0 ± this [m]"
            arg_type = Float64
            default  = 80.0
        "--lat-halfwidth"
            help     = "latitudes span lat0 ± this [deg]"
            arg_type = Float64
            default  = 15.0
        "--namp"
            help     = "number of amplitude values (evenly spaced over hhat0 ± amp-halfwidth)"
            arg_type = Int
            default  = 5
        "--nlat"
            help     = "number of latitude values (evenly spaced over lat0 ± lat-halfwidth)"
            arg_type = Int
            default  = 5
        "--design"
            help    = "grid = every (ĥ, φ₂) combination; subset = 2 Latin diagonals of that grid, i.e. 2n runs for n × n values (default)"
            default = "subset"
        "--amplitudes"
            help    = "explicit amplitudes ĥ [m], e.g. \"100,120,140\"; replaces --hhat0/--amp-halfwidth/--namp"
            default = ""
        "--latitudes"
            help    = "explicit latitudes φ₂ [deg], e.g. \"40,45,50\"; replaces --lat0/--lat-halfwidth/--nlat"
            default = ""
        "--set"
            help   = "override any deck input in every trajectory, KEY=VALUE with VALUE a Julia expression, e.g. --set μ=1.0e5 (repeatable)"
            action = :append_arg
        "--outdir", "-o"
            help    = "parent directory of the trajectory_X directories"
            default = "/scratch/smarras/trajectories"
        "--tend"
            help    = "final time :tend of every trajectory [s] (default: 20 days)"
            arg_type = Float64
            default  = 100*24*3600.0
        "--early-interval"
            help     = "time between outputs up to --early-until [s] (default: 1 h)"
            arg_type = Float64
            default  = 3600.0
        "--early-until"
            help     = "end of the densely written early phase [s] (default: 4 days; 0 = none)"
            arg_type = Float64
            default  = 4*24*3600.0
        "--output-interval"
            help     = "time between outputs after --early-until [s] (default: 6 h)"
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

family(x0, w, n) = n == 1 ? [x0] : collect(range(x0 - w, x0 + w; length = n))

parse_list(str) = parse.(Float64, split(str, ','))

hhat0, lat0 = args["hhat0"], args["lat0"]
amps = isempty(strip(args["amplitudes"])) ? family(hhat0, args["amp-halfwidth"], args["namp"]) :
                                            parse_list(args["amplitudes"])
lats = isempty(strip(args["latitudes"]))  ? family(lat0,  args["lat-halfwidth"], args["nlat"]) :
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

# --design subset: the cells (i, i) and (i, i + ⌊n/2⌋ mod n) of the n × n grid,
# so each amplitude and each latitude appears exactly twice.
if args["design"] == "subset"
    n = length(amps)
    length(lats) == n ||
        error("--design subset needs as many amplitudes as latitudes (have $(length(amps)) and $(length(lats))); use --design grid")
    keep = Set(vcat([(i, i) for i in 1:n], [(i, mod(i - 1 + n ÷ 2, n) + 1) for i in 1:n]))
    filter!(p -> (p.iamp, p.ilat) in keep, params)
elseif args["design"] != "grid"
    error("--design must be \"grid\" or \"subset\", got \"$(args["design"])\"")
end
# Contiguous numbering of whatever was kept, in grid order.
trajectories = [(; id = args["first"] + i - 1, p...) for (i, p) in enumerate(params)]

# Output times: every --early-interval up to --early-until, then every
# --output-interval to tend (always ending exactly at tend).
function output_times(tend, dt_early, t_early, dt_late)
    (dt_early > 0 && dt_late > 0) || error("output intervals must be positive")
    t_early = min(t_early, tend)
    ts = collect(dt_early:dt_early:t_early)
    t0 = isempty(ts) ? 0.0 : ts[end]
    append!(ts, (t0 + dt_late):dt_late:tend)
    (isempty(ts) || !isapprox(ts[end], tend)) && push!(ts, tend)
    return ts
end
const OUTPUT_TIMES = output_times(args["tend"], args["early-interval"],
                                  args["early-until"], args["output-interval"])

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
    @printf(stderr, " # tend = %g h; %d outputs at t = %s h\n", args["tend"]/3600, length(OUTPUT_TIMES),
            join((@sprintf("%g", x/3600) for x in OUTPUT_TIMES), ", "))
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
    # An explicit output schedule: the sphere solver writes at these times
    # when :ndiagnostics_outputs is 0.
    overrides[:diagnostics_at_times] = OUTPUT_TIMES
    overrides[:ndiagnostics_outputs] = 0
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
