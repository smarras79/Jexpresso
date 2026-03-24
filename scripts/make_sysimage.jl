# scripts/make_sysimage.jl
#
# Build a PackageCompiler system image that contains precompiled native code
# for Jexpresso's full execution path.  After building once, every subsequent
# run starts in seconds instead of minutes.
#
# USAGE
#   julia --project=. scripts/make_sysimage.jl <equations> <case>
#
# EXAMPLE
#   julia --project=. scripts/make_sysimage.jl CompEuler thetaAlya
#
# OUTPUT
#   jexpresso.so  (in the project root — add --sysimage jexpresso.so to julia)
#
# TWO BUILD MODES (chosen automatically based on available RAM)
#
#   Full  (≥ 4 GiB available)
#     Runs one warmup timestep; every JIT-compiled method is baked in.
#     Fastest possible cold start — package loading AND first timestep JIT
#     are both eliminated.  Requires ~4–6 GiB peak RSS from LLVM.
#     Override: set JEXPRESSO_SYSIMAGE_FULL=1 to request full mode.
#              (Still aborts if RAM is below the hard minimum of 3.5 GiB.)
#
#   Light (< 4 GiB available, automatically chosen)
#     No warmup execution; only each package's own @compile_workload
#     blocks are compiled.  Package-loading latency is eliminated; the
#     first-timestep JIT still runs at runtime.  Requires ~2 GiB.
#
#   Skip (< 2 GiB available)
#     Not enough RAM to build any useful sysimage.  Exits with code 0
#     without creating jexpresso.so so the caller can run Julia directly.
#     Close other applications and retry to get a sysimage.
#
# WHY --heap-size-hint DOES NOT HELP
#   --heap-size-hint only controls Julia's GC heap aggressiveness.  LLVM
#   (which does the actual native code generation) allocates memory directly
#   from the OS, completely outside Julia's GC.  The only way to reduce peak
#   LLVM RSS is to reduce the number of methods it has to compile.
#
# REBUILD WHEN
#   Rebuild after significant source changes (new functions, changed type
#   signatures).  Minor bug fixes usually don't require a rebuild.

using PackageCompiler

# ── Argument handling ──────────────────────────────────────────────────────────
if length(ARGS) < 2
    println("""
    Usage: julia --project=. scripts/make_sysimage.jl <equations> <case>

    Example:
        julia --project=. scripts/make_sysimage.jl CompEuler thetaAlya
    """)
    exit(1)
end

equations     = ARGS[1]
case          = ARGS[2]
root          = dirname(@__DIR__)           # project root
sysimage_path = joinpath(root, "jexpresso.so")

# ── Available-memory helper ────────────────────────────────────────────────────
function available_memory_gib()
    if Sys.islinux()
        for line in eachline("/proc/meminfo")
            if startswith(line, "MemAvailable:")
                return parse(Int, split(line)[2]) / 2^20   # KiB → GiB
            end
        end
    elseif Sys.isapple()
        try
            out = readchomp(`vm_stat`)
            pm  = match(r"page size of (\d+) bytes", out)
            pagesize = pm === nothing ? 4096 : parse(Int, pm[1])
            pages = Dict{String,Int}()
            for line in split(out, '\n')
                m = match(r"^Pages\s+(.+):\s+(\d+)\.", line)
                m === nothing && continue
                pages[strip(String(m[1]))] = parse(Int, m[2])
            end
            avail = (get(pages, "free", 0) + get(pages, "inactive", 0) +
                     get(pages, "speculative", 0)) * pagesize / 2^30
            avail > 0 && return avail
        catch
        end
    end
    return Sys.free_memory() / 2^30
end

# ── Memory thresholds ──────────────────────────────────────────────────────────
# LLVM code generation runs OUTSIDE Julia's GC heap.  --heap-size-hint does
# not limit LLVM memory.  The only way to reduce peak RSS is to reduce the
# number of methods LLVM compiles.  These thresholds reflect measured peaks:
#
#   Full  build: ~4–6 GiB LLVM  (tens of thousands of solver specialisations)
#   Light build: ~2–3 GiB LLVM  (package @compile_workload blocks only)
#
const FULL_MIN_GiB  = 3.5   # hard minimum for full build — abort below this
const FULL_WANT_GiB = 4.0   # auto-select threshold: use full only if ≥ this
const LIGHT_MIN_GiB = 2.0   # hard minimum for light build — skip below this

avail_gib   = available_memory_gib()
total_gib   = Sys.total_memory() / 2^30
force_full  = haskey(ENV, "JEXPRESSO_SYSIMAGE_FULL")

println("Memory: $(round(avail_gib, digits=1)) GiB available / $(round(total_gib, digits=1)) GiB total")

# ── Decide: full / light / skip ───────────────────────────────────────────────
if force_full && avail_gib < FULL_MIN_GiB
    # Hard abort: user asked for full but there is no way it can succeed.
    # Do this BEFORE any expensive work so we don't waste 13+ minutes.
    println(stderr, """

    ERROR: Full sysimage build requires at least $(FULL_MIN_GiB) GiB of free RAM.
           Available right now: $(round(avail_gib, digits=1)) GiB / $(round(total_gib, digits=1)) GiB total.

    To fix this, choose one of:
      1. Close other applications to free ~$(FULL_MIN_GiB) GiB, then retry:
             JEXPRESSO_SYSIMAGE_FULL=1 REBUILD_SYSIMAGE=1 ./run_coupled.sh 2 2

      2. Drop JEXPRESSO_SYSIMAGE_FULL=1 to let the script pick the lighter
         build mode (needs ≥ $(LIGHT_MIN_GiB) GiB free):
             REBUILD_SYSIMAGE=1 ./run_coupled.sh 2 2

      3. Skip the sysimage entirely and run Julia directly (slower startup,
         but works on any machine).  Just run without REBUILD_SYSIMAGE=1
         and without a pre-existing jexpresso.so.
    """)
    exit(1)

elseif avail_gib < LIGHT_MIN_GiB
    # Not enough RAM for even the light build.  Exit gracefully (code 0) so
    # run_coupled.sh can continue without a sysimage rather than aborting.
    println(stderr, """

    WARNING: Not enough free RAM ($(round(avail_gib, digits=1)) GiB) to build any sysimage.
             Minimum needed for light build: $(LIGHT_MIN_GiB) GiB.
             Skipping sysimage — Julia will run without --sysimage (slower startup).

    To build a sysimage later, close other applications to free ≥ $(LIGHT_MIN_GiB) GiB and run:
        REBUILD_SYSIMAGE=1 ./run_coupled.sh 2 2
    """)
    exit(0)   # caller checks for jexpresso.so to decide whether to use it
end

use_full_mode = force_full || avail_gib >= FULL_WANT_GiB

if !use_full_mode
    println("""
    Low memory ($(round(avail_gib, digits=1)) GiB available) — using LIGHT build mode.
      Package-loading latency eliminated; first-timestep JIT still runs at runtime.
      Set JEXPRESSO_SYSIMAGE_FULL=1 to request a full build (needs ≥ $(FULL_MIN_GiB) GiB free).
    """)
end

# ── GC heap hint ───────────────────────────────────────────────────────────────
# Applies only to Julia's GC inside the subprocess, not LLVM.
heap_gib = clamp(floor(Int, avail_gib * 0.70), 2, 32)

# ── CPU target ─────────────────────────────────────────────────────────────────
# On aarch64 (Apple Silicon, AWS Graviton) default_app_cpu_target() returns
# "generic", which cannot lower the `vscale` LLVM intrinsic emitted by
# HostCPUFeatures.jl → signal 6 abort.
cpu_target = Sys.ARCH === :aarch64 ? "native" : PackageCompiler.default_app_cpu_target()

@info "Building Jexpresso sysimage" equations case output=sysimage_path mode=(use_full_mode ? "full" : "light") avail_gib=round(avail_gib, digits=1) cpu_target

# ── Warmup script (full mode only) ────────────────────────────────────────────
warmup_script = nothing
if use_full_mode
    @info "Step 1/2: Running warmup execution (1 timestep) — expect a few minutes..."
    warmup_script = joinpath(root, "scripts", "_warmup_tmp.jl")
    write(warmup_script, """
# Auto-generated by make_sysimage.jl — do not edit; deleted after use.
import Pkg
Pkg.activate("$(root)")

# JEXPRESSO_WARMUP=1 caps the run to 1 timestep and suppresses output.
# Do NOT set JEXPRESSO_SYSIMAGE_BUILD here — that guard prevents
# jexpresso_main() from running, which defeats the entire purpose of a
# warmup: we need the full execution path (driver → sem_setup →
# mod_mesh_mesh_driver → build_Interpolation_basis! → time_loop! → …)
# to actually execute so PackageCompiler can bake those specialisations
# into the sysimage.  Without this, the sysimage only covers package
# loading and the first-timestep JIT still hits at runtime.
ENV["JEXPRESSO_WARMUP"] = "1"  # 1 timestep, no output files
push!(empty!(ARGS), "$(equations)", "$(case)")
include(joinpath("$(root)", "src", "Jexpresso.jl"))
""")
else
    @info "Step 1/2: Skipped warmup (light mode) — compiling package workloads only..."
end

# ── Build the sysimage ────────────────────────────────────────────────────────
# Signal to Jexpresso.jl that it should skip __precompile__(false).
# This env var is inherited by every subprocess PackageCompiler spawns,
# including the --output-o compilation subprocess that loads `using Jexpresso`.
ENV["JEXPRESSO_BUILDING_SYSIMAGE"] = "1"

try
    kwargs = Dict{Symbol,Any}(
        :sysimage_path       => sysimage_path,
        :cpu_target          => cpu_target,
        :sysimage_build_args => `--heap-size-hint=$(heap_gib)G`,
    )
    use_full_mode && (kwargs[:precompile_execution_file] = warmup_script)

    create_sysimage(["Jexpresso"]; kwargs...)

    mode_note = use_full_mode ?
        "Full sysimage — package loading AND first-timestep JIT eliminated." :
        "Light sysimage — package loading eliminated; first-timestep JIT runs at runtime.\n  │  Rebuild with JEXPRESSO_SYSIMAGE_FULL=1 on a machine with ≥ $(FULL_MIN_GiB) GiB free for best performance."

    @info "Step 2/2: Sysimage written to $sysimage_path"
    @info """\

    ┌─────────────────────────────────────────────────────────────────────┐
    │  Sysimage ready!  Add this flag to every julia invocation:          │
    │                                                                     │
    │    julia --sysimage $(relpath(sysimage_path)) ...                   │
    │                                                                     │
    │  $(mode_note)
    │                                                                     │
    │  Rebuild after significant source changes:                          │
    │    julia --project=. scripts/make_sysimage.jl $(equations) $(case) │
    └─────────────────────────────────────────────────────────────────────┘
    """
catch e
    @error "Sysimage build failed" exception=(e, catch_backtrace())
    rethrow()
finally
    warmup_script !== nothing && rm(warmup_script, force=true)
end
