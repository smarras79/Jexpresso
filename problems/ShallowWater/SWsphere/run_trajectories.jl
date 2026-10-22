#!/usr/bin/env julia
#-------------------------------------------------------------------------------
# run_trajectories.jl — ensemble of SWsphere (Galewsky jet) runs for ROM training.
#
# Each TRAJECTORY is one full SWsphere run, identical to the shipped deck
# (user_inputs.jl) except for the height perturbation
#
#     h'(λ,φ) = ĥ cos φ exp[-(λ/α)²] exp[-((φ₂-φ)/β)²]
#
# The default ensemble is 10 trajectories, varying one parameter at a time
# around the published case (ĥ = 120 m, φ₂ = 45°):
#
#   trajectory_1  … trajectory_5   amplitude ĥ  = 120 m × (1 + [-5, -2.5, 0, +2.5, +5] %),  φ₂ = 45°
#   trajectory_6  … trajectory_10  latitude  φ₂ =  45°  × (1 + [-5, -2.5, 0, +2.5, +5] %),  ĥ = 120 m
#
# (the 0 % member of each group is the reference Galewsky run; it appears in
# both groups so each group is a complete, symmetric 5-point family.)
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
#   # wider/denser families, a different reference, shorter runs:
#   julia --project=. problems/ShallowWater/SWsphere/run_trajectories.jl \
#         --amp-variation 0.10 --namp 7 --lat-variation 0.05 --nlat 5 --tend 86400
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
            help     = "number of amplitude trajectories (evenly spaced over ±variation)"
            arg_type = Int
            default  = 5
        "--nlat"
            help     = "number of latitude trajectories (evenly spaced over ±variation)"
            arg_type = Int
            default  = 5
        "--outdir", "-o"
            help    = "parent directory of the trajectory_X directories"
            default = "./trajectories"
        "--tend"
            help    = "override :tend [s] (default: the deck's value)"
            arg_type = Float64
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

hhat0, lat0 = args["hhat0"], args["lat0"]
params = vcat([(hhat = h, lat_deg = lat0, family = "amplitude")
               for h in family(hhat0, args["amp-variation"], args["namp"])],
              [(hhat = hhat0, lat_deg = φ, family = "latitude")
               for φ in family(lat0, args["lat-variation"], args["nlat"])])
trajectories = [(; id = args["first"] + i - 1, p...) for (i, p) in enumerate(params)]

φ0_deg, φ1_deg = rad2deg(π/7), rad2deg(π/2 - π/7)
for t in trajectories
    (φ0_deg < t.lat_deg < φ1_deg) ||
        @warn @sprintf("trajectory_%d: φ₂ = %.2f° lies outside the jet (%.1f°–%.1f°)",
                       t.id, t.lat_deg, φ0_deg, φ1_deg)
end

function print_table(io, trajectories)
    println(io, "trajectory,directory,family,galewsky_hhat_m,galewsky_phi2_deg,galewsky_phi2_rad")
    for t in trajectories
        @printf(io, "%d,trajectory_%d,%s,%.10g,%.10g,%.16g\n",
                t.id, t.id, t.family, t.hhat, t.lat_deg, deg2rad(t.lat_deg))
    end
end

if args["dry-run"]
    print_table(stdout, trajectories)
    exit(0)
end

using Jexpresso
using MPI
MPI.Initialized() || MPI.Init()
const RANK = MPI.Comm_rank(MPI.COMM_WORLD)

if RANK == 0
    mkpath(outdir)
    open(io -> print_table(io, trajectories), joinpath(outdir, "trajectories.csv"), "w")
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
    args["tend"] === nothing || (overrides[:tend] = args["tend"])

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
