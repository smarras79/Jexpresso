# Install Jexpresso's optional AMR dependency (GridapP4est) into envs/amr/.
#
#   julia --project=. tools/setup_amr.jl
#
# Run it once after the root environment is instantiated. Afterwards launch
# AMR cases exactly like any other case; Jexpresso picks envs/amr/ up on its
# own when a case turns AMR on (src/kernel/Adaptivity/p4est_hooks.jl).
#
# After changing the MPI binding (INSTALL.md §5) run it again with
#   --rebuild   re-run the build steps of GridapP4est and its dependencies
#               (P4est_wrapper), which link against the MPI found at build time
#   --fresh     delete envs/amr/Manifest.toml first, so P4est_jll is re-resolved
#               to the variant of the new MPI (needed after a system <-> JLL
#               route change); implies --rebuild
#
# Not supported on Windows: P4est_wrapper's build step needs a p4est + MPI
# toolchain. Everything except AMR runs without this step.

using Pkg

const ROOT    = dirname(@__DIR__)
const AMR_ENV = joinpath(ROOT, "envs", "amr")
const FRESH   = "--fresh" in ARGS
const REBUILD = FRESH || "--rebuild" in ARGS

if Sys.iswindows()
    @warn "AMR (GridapP4est / P4est_wrapper) is not supported on Windows; " *
          "Jexpresso runs every non-AMR case without it."
end

# P4est_jll ships one build per MPI, picked when envs/amr's Manifest is
# resolved, from the MPIPreferences binding visible to that environment. Give
# it the root environment's binding so p4est and MPI.jl load the same libmpi.
root_prefs = joinpath(ROOT, "LocalPreferences.toml")
if isfile(root_prefs)
    cp(root_prefs, joinpath(AMR_ENV, "LocalPreferences.toml"); force = true)
    println(" # Copied the root MPI binding (LocalPreferences.toml) into envs/amr/")
else
    rm(joinpath(AMR_ENV, "LocalPreferences.toml"); force = true)
end

FRESH && rm(joinpath(AMR_ENV, "Manifest.toml"); force = true)

Pkg.activate(AMR_ENV)
withenv("JULIA_PKG_PRECOMPILE_AUTO" => "0") do
    Pkg.instantiate()
end
# Pkg.build builds GridapP4est and, depth first, every dependency with a build
# step -- P4est_wrapper among them.
REBUILD && Pkg.build("GridapP4est"; verbose = true)
Pkg.status("GridapP4est")

# Precompile GridapP4est and the extension the way a run loads them: from the
# root environment, with envs/amr appended to LOAD_PATH. Doing it here keeps
# the first AMR run (possibly on many MPI ranks at once) from compiling it.
println(" # Precompiling the AMR extension ...")
run(`$(Base.julia_cmd()) --project=$ROOT -e "using Jexpresso; Jexpresso._ensure_amr_loaded!(); println(\" # AMR extension loaded: \", Base.get_extension(Jexpresso, :JexpressoP4estExt))"`)
