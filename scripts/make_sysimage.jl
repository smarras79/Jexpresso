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
# THREE BUILD MODES (chosen automatically based on available RAM)
#
#   Full  (≥ 4 GiB available, or JEXPRESSO_SYSIMAGE_FULL=1 with ≥ 3.5 GiB)
#     Warmup runs one timestep; every JIT-compiled method is baked in at
#     default (-O2) optimisation.  Fastest possible cold start — package
#     loading AND first-timestep JIT are both eliminated.
#     Peak LLVM RSS: ~4–6 GiB.
#
#   Lean  (≥ 1.5 GiB available, automatically chosen on constrained machines)
#     Same warmup execution as Full mode, but LLVM compiles at -O1 instead
#     of the default -O2.  This roughly halves LLVM's memory footprint while
#     still eliminating all first-timestep JIT.  Runtime performance is
#     slightly less than Full but far better than no sysimage.
#     Ideal for laptops with 8 GiB RAM where the OS leaves 1.5–3 GiB free.
#     Force lean mode: JEXPRESSO_SYSIMAGE_LEAN=1
#     Peak LLVM RSS: ~2–3 GiB.
#
#   Skip  (< 1.5 GiB available)
#     Not enough RAM to build any useful sysimage.  Exits with code 0
#     without creating jexpresso.so so the caller can run Julia directly.
#     Close other applications and retry to get a sysimage.
#
# WHY --heap-size-hint DOES NOT HELP
#   --heap-size-hint only controls Julia's GC heap aggressiveness.  LLVM
#   (which does the actual native code generation) allocates memory directly
#   from the OS, completely outside Julia's GC.  The only way to reduce peak
#   LLVM RSS is to reduce the number of methods it has to compile (-O1 helps
#   because each method's IR is processed by fewer passes).
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
# not limit LLVM memory.  These thresholds reflect measured peak RSS:
#
#   Full build  (-O2, default): ~4–6 GiB LLVM
#   Lean build  (-O1):          ~2–3 GiB LLVM  (same warmup, fewer LLVM passes)
#   Skip:                       < 1.5 GiB free — nothing useful can be built
#
const FULL_MIN_GiB  = 3.5   # hard minimum for full (-O2) build
const FULL_WANT_GiB = 4.0   # auto-select: prefer full only if ≥ this
const LEAN_MIN_GiB  = 1.5   # hard minimum for lean (-O1) build; skip below
const LEAN_WANT_GiB = 2.5   # auto-select: use lean if ≥ this but < FULL_WANT

avail_gib   = available_memory_gib()
total_gib   = Sys.total_memory() / 2^30
force_full  = haskey(ENV, "JEXPRESSO_SYSIMAGE_FULL")
force_lean  = haskey(ENV, "JEXPRESSO_SYSIMAGE_LEAN")

println("Memory: $(round(avail_gib, digits=1)) GiB available / $(round(total_gib, digits=1)) GiB total")

# ── Decide: full / lean / skip ────────────────────────────────────────────────
if force_full && avail_gib < FULL_MIN_GiB
    # Hard abort: user explicitly asked for full but RAM is insufficient.
    println(stderr, """

    ERROR: Full sysimage build (-O2) requires at least $(FULL_MIN_GiB) GiB of free RAM.
           Available right now: $(round(avail_gib, digits=1)) GiB / $(round(total_gib, digits=1)) GiB total.

    Options:
      1. Close other applications to free ≥ $(FULL_MIN_GiB) GiB and retry:
             JEXPRESSO_SYSIMAGE_FULL=1 REBUILD_SYSIMAGE=1 ./run_coupled.sh 2 2

      2. Use lean mode (-O1, needs ≥ $(LEAN_MIN_GiB) GiB) — same warmup, slightly less
         optimised code, roughly half the LLVM memory:
             JEXPRESSO_SYSIMAGE_LEAN=1 REBUILD_SYSIMAGE=1 ./run_coupled.sh 2 2

      3. Let the script auto-pick (drop both override flags):
             REBUILD_SYSIMAGE=1 ./run_coupled.sh 2 2
    """)
    exit(1)

elseif avail_gib < LEAN_MIN_GiB
    # Not enough RAM for any useful build.  Exit gracefully (code 0).
    println(stderr, """

    WARNING: Not enough free RAM ($(round(avail_gib, digits=1)) GiB) to build any sysimage.
             Minimum for lean build: $(LEAN_MIN_GiB) GiB.
             Skipping sysimage — Julia will run without --sysimage (slower startup).

    To build a sysimage, close other applications to free ≥ $(LEAN_MIN_GiB) GiB, then:
        REBUILD_SYSIMAGE=1 ./run_coupled.sh 2 2
    """)
    exit(0)
end

# Auto-select mode: full if enough RAM, lean otherwise.
# Explicit override flags (JEXPRESSO_SYSIMAGE_FULL / _LEAN) take priority.
use_full_mode = (force_full && avail_gib >= FULL_MIN_GiB) ||
                (!force_lean && !force_full && avail_gib >= FULL_WANT_GiB)
use_lean_mode = !use_full_mode   # everything that is not full is lean

if use_full_mode
    println("Memory: $(round(avail_gib, digits=1)) GiB available — using FULL build mode (-O2).")
else
    println("""
    Memory: $(round(avail_gib, digits=1)) GiB available — using LEAN build mode (-O1).
      Both modes run the same warmup and eliminate first-timestep JIT.
      Lean uses roughly half the LLVM memory by skipping some optimisation passes.
      Runtime performance is slightly lower than full mode but still far faster
      than running without a sysimage.
      To force full mode (needs ≥ $(FULL_MIN_GiB) GiB free):
          JEXPRESSO_SYSIMAGE_FULL=1 REBUILD_SYSIMAGE=1 ./run_coupled.sh 2 2
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

@info "Building Jexpresso sysimage" equations case output=sysimage_path mode=(use_full_mode ? "full" : "lean") avail_gib=round(avail_gib, digits=1) cpu_target

# ── Warmup script (full and lean modes) ───────────────────────────────────────
warmup_script = nothing
if use_full_mode || use_lean_mode
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
end   # warmup script written for both full and lean modes

# ── Build the sysimage ────────────────────────────────────────────────────────
# Signal to Jexpresso.jl that it should skip __precompile__(false).
# This env var is inherited by every subprocess PackageCompiler spawns,
# including the --output-o compilation subprocess that loads `using Jexpresso`.
ENV["JEXPRESSO_BUILDING_SYSIMAGE"] = "1"

try
    # Lean mode passes -O1 to the LLVM compilation subprocess to cut its peak
    # RSS roughly in half compared to the default -O2.  Full mode omits -O1
    # so LLVM uses its default optimisation level for best runtime performance.
    build_args = use_lean_mode ? `-O1 --heap-size-hint=$(heap_gib)G` :
                                 `--heap-size-hint=$(heap_gib)G`

    kwargs = Dict{Symbol,Any}(
        :sysimage_path       => sysimage_path,
        :cpu_target          => cpu_target,
        :sysimage_build_args => build_args,
    )
    warmup_script !== nothing && (kwargs[:precompile_execution_file] = warmup_script)

    create_sysimage(["Jexpresso"]; kwargs...)

    mode_note = use_full_mode ?
        "Full sysimage (-O2) — package loading AND first-timestep JIT eliminated." :
        "Lean sysimage (-O1) — package loading AND first-timestep JIT eliminated.\n  │  Slightly less optimised than full mode; rebuild with JEXPRESSO_SYSIMAGE_FULL=1\n  │  on a machine with ≥ $(FULL_MIN_GiB) GiB free RAM for best performance."

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
