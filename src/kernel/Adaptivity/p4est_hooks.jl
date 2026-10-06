# ─── AMR / p4est: optional extension hooks ──────────────────────────────
#
# GridapP4est (and through it P4est_wrapper/P4est_jll) is a WEAK dependency
# of Jexpresso: every call into it lives in ext/JexpressoP4estExt/impl.jl,
# which defines methods for the hook functions below. Without GridapP4est,
# Jexpresso loads and runs every non-AMR case — on Windows too, where
# P4est_wrapper does not build.
#
# The hooks are reached only by the AMR-style paths, all gated by
# `_adaptive_mesh_run(inputs)` (`:lamr`, `:ladapt`, `:lpreadapt`,
# `:linitial_refine`, `:lrestart_amr`):
#   * mesh.jl            :: mod_mesh_read_gmsh!  (octree model construction)
#   * sem_setup.jl       :: AMR restart branch   (checkpoint load, ad_lvl)
#   * TimeIntegrators.jl :: write_p4est_checkpoint
# Everything else the AMR path does to an octree model goes through Gridap /
# GridapDistributed generic functions (`Gridap.Adaptivity.adapt`,
# `GridapDistributed.redistribute`, `get_cell_gids`, ...) that GridapP4est
# extends, so it needs no hook.
#
# Loading. run.jl calls `_ensure_amr_loaded!()` as a TOP-LEVEL statement,
# before the driver, whenever the case is an AMR run. Doing it there and not
# lazily from inside the driver matters: methods defined by loading a package
# are invisible to code that is already running (world age), so GridapP4est's
# `adapt`/`redistribute` methods would not dispatch from a driver that started
# before the load. The in-driver call sites therefore only *check*
# (`_assert_amr_loaded()`), they never load.

const _GRIDAPP4EST_PKGID = Base.PkgId(
    Base.UUID("c2c8e14b-f5fd-423d-9666-1dd9ad120af9"), "GridapP4est")

# Optional environment that installs GridapP4est (see envs/amr/Project.toml).
const _AMR_ENV_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "envs", "amr"))

# Implementation shared by the package extension and the script-mode loader.
const _P4EST_IMPL_FILE = normpath(joinpath(@__DIR__, "..", "..", "..",
                                           "ext", "JexpressoP4estExt", "impl.jl"))

# Set by the extension's __init__ (package mode) or by
# _load_p4est_impl_script_mode! (script mode).
const _AMR_LOADED = Ref(false)

# GridapP4est's adaptivity flags (OctreeDistributedDiscreteModels.jl). They
# are plain Cint values handed to `Gridap.Adaptivity.adapt`; defined here so
# that user `initialize.jl` files and the flag bookkeeping in mesh.jl /
# Projection.jl work without GridapP4est loaded. Must stay equal to
# GridapP4est.{nothing,refine,coarsen}_flag — impl.jl asserts it on load.
const nothing_flag = Cint(0)
const refine_flag  = Cint(1)
const coarsen_flag = Cint(2)

"""
    amr_uniformly_refined_model(parts, coarse_model, nlevels)

`GridapP4est.UniformlyRefinedForestOfOctreesDiscreteModel(parts, coarse_model,
nlevels)`. Implemented in `ext/JexpressoP4estExt`.
"""
function amr_uniformly_refined_model end

"""
    amr_octree_model(parts, coarse_model)

`GridapP4est.OctreeDistributedDiscreteModel(parts, coarse_model)`.
Implemented in `ext/JexpressoP4estExt`.
"""
function amr_octree_model end

"""
    amr_copy_octree_model(model, dmodel)

A new `OctreeDistributedDiscreteModel` that shares `model`'s connectivity and
coarse model, owns a copy of its p4est forest (`pXest_copy`) and wraps the
distributed discrete model `dmodel`. Used by the AMR re-adapt path so
`Gridap.Adaptivity.adapt` can consume the copy while `model` stays intact.
Implemented in `ext/JexpressoP4estExt`.
"""
function amr_copy_octree_model end

"""
    write_p4est_checkpoint(output_dir, iter, partitioned_model)

Save the p4est forest topology to `output_dir/iter_N/iter_N.p4est`.
Called alongside each VTK write to enable AMR restarts.
Only called when `inputs[:lamr] == true`. Implemented in `ext/JexpressoP4estExt`.
"""
function write_p4est_checkpoint end

"""
    read_ad_lvl_from_p4est(pXest_type, ptr_pXest) -> Vector{TInt}

Walk the local trees of a p4est/p8est forest and return the level of every
local leaf quadrant, in p4est ordering.  This ordering matches the
Jexpresso mesh element ordering when the mesh is built directly from
the same forest (e.g. after an AMR restart via `load_p4est_checkpoint_model`).
Implemented in `ext/JexpressoP4estExt`.
"""
function read_ad_lvl_from_p4est end

"""
    load_p4est_checkpoint_model(base_model, forest_file) -> OctreeDistributedDiscreteModel

Load a p4est forest checkpoint saved by `write_p4est_checkpoint` and build a full
`OctreeDistributedDiscreteModel` from it.

`base_model` should be the coarse (or preadapted) `OctreeDistributedDiscreteModel`
built from the original .msh file — its `coarse_model` and `ptr_pXest_connectivity`
provide the geometric context for the loaded forest.  All MPI ranks call collectively.
Implemented in `ext/JexpressoP4estExt`.
"""
function load_p4est_checkpoint_model end

_amr_unavailable_msg() = """
    This case uses AMR (one of :lamr, :ladapt, :lpreadapt, :linitial_refine,
    :lrestart_amr is true), which needs the optional GridapP4est dependency.
    It is not installed in this environment. Install it once with

        julia --project=. tools/setup_amr.jl

    (not available on Windows: P4est_wrapper does not build there), or load
    GridapP4est yourself before running the case. See docs/amr_setup.md,
    "Installing the AMR extension"."""

"""
    Jexpresso._ensure_amr_loaded!()

Load GridapP4est and the AMR extension, if they are not loaded yet. Looks for
GridapP4est on the current `LOAD_PATH` first, then in `envs/amr/` (set up by
`tools/setup_amr.jl`). Call it at top level — not from inside a running
driver; see the world-age note in src/kernel/Adaptivity/p4est_hooks.jl.
run.jl does this for every AMR case.
"""
function _ensure_amr_loaded!()
    _AMR_LOADED[] && return nothing
    if Base.locate_package(_GRIDAPP4EST_PKGID) === nothing &&
       isfile(joinpath(_AMR_ENV_DIR, "Manifest.toml")) &&
       !(_AMR_ENV_DIR in LOAD_PATH)
        # Appended, not prepended: every package the active environment
        # already provides keeps coming from its own Manifest; envs/amr only
        # contributes what is missing there (GridapP4est, P4est_wrapper, ...).
        push!(LOAD_PATH, _AMR_ENV_DIR)
    end
    Base.locate_package(_GRIDAPP4EST_PKGID) === nothing &&
        error(_amr_unavailable_msg())
    gridapp4est = Base.require(_GRIDAPP4EST_PKGID)
    # Loading Jexpresso as a package: requiring GridapP4est has just run the
    # extension, whose __init__ set the flag. Running src/Jexpresso.jl as a
    # script: there is no package for an extension to attach to, so load the
    # same implementation into this module by hand.
    _AMR_LOADED[] || Base.invokelatest(_load_p4est_impl_script_mode!, gridapp4est)
    return nothing
end

function _load_p4est_impl_script_mode!(gridapp4est::Module)
    ext = Module(:JexpressoP4estExt)
    Core.eval(ext, :(const JX = $(@__MODULE__)))
    Core.eval(ext, :(const GridapP4est = $gridapp4est))
    Base.include(ext, _P4EST_IMPL_FILE)
    _AMR_LOADED[] = true
    return nothing
end

function _assert_amr_loaded()
    _AMR_LOADED[] && return nothing
    error(_amr_unavailable_msg() * """


        (GridapP4est must be loaded before the run starts: call
        Jexpresso._ensure_amr_loaded!() at top level first. run.jl does this
        for you.)""")
end
