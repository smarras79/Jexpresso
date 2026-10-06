# Implementation of Jexpresso's AMR hooks (src/kernel/Adaptivity/p4est_hooks.jl).
#
# Included into a module that binds
#   JX          => the Jexpresso module
#   GridapP4est => the GridapP4est module
# by ext/JexpressoP4estExt.jl (package mode) and by
# Jexpresso._load_p4est_impl_script_mode! (script mode). Use only those two
# names (and Base) here, fully qualified.

const P4est_wrapper = GridapP4est.P4est_wrapper

# Jexpresso hard-codes GridapP4est's flag values so it can run without it.
@assert JX.nothing_flag == GridapP4est.nothing_flag &&
        JX.refine_flag  == GridapP4est.refine_flag  &&
        JX.coarsen_flag == GridapP4est.coarsen_flag  "GridapP4est changed its adaptivity flag values; update src/kernel/Adaptivity/p4est_hooks.jl"

JX.amr_uniformly_refined_model(parts, coarse_model, nlevels) =
    GridapP4est.UniformlyRefinedForestOfOctreesDiscreteModel(parts, coarse_model, nlevels)

JX.amr_octree_model(parts, coarse_model) =
    GridapP4est.OctreeDistributedDiscreteModel(parts, coarse_model)

function JX.amr_copy_octree_model(model, dmodel)
    return GridapP4est.OctreeDistributedDiscreteModel(
        model.parts,
        dmodel,
        model.non_conforming_glue,
        model.coarse_model,
        model.ptr_pXest_connectivity,
        GridapP4est.pXest_copy(model.pXest_type, model.ptr_pXest),
        model.pXest_type,
        model.pXest_refinement_rule_type,
        model.owns_ptr_pXest_connectivity,
        model.gc_ref)
end

function JX.write_p4est_checkpoint(output_dir::String, iter::Int, partitioned_model)
    comm  = JX.MPI.COMM_WORLD
    rank  = JX.MPI.Comm_rank(comm)
    dir   = joinpath(output_dir, "iter_$(iter)")
    fname = joinpath(dir, "iter_$(iter).p4est")
    if rank == 0
        mkpath(dir)
    end
    JX.MPI.Barrier(comm)
    # save_data=0: no per-quadrant payload, forest topology only
    # Dispatch on 2D (p4est_save) vs 3D (p8est_save) to avoid passing the wrong struct type.
    if partitioned_model.pXest_type isa GridapP4est.P4estType
        JX.@outputrootonly P4est_wrapper.p4est_save(fname, partitioned_model.ptr_pXest, Cint(0))
    else
        JX.@outputrootonly P4est_wrapper.p8est_save(fname, partitioned_model.ptr_pXest, Cint(0))
    end
end

# Dispatches on 2D (`p4est_tree_t`/`p4est_quadrant_t`) vs 3D
# (`p8est_tree_t`/`p8est_quadrant_t`) — these have different memory layouts,
# so reading a 2D forest with the 3D struct types silently misreads garbage.
# Mirrors the same 2D/3D dispatch used by write_p4est_checkpoint /
# load_p4est_checkpoint_model.
function JX.read_ad_lvl_from_p4est(pXest_type, ptr_pXest)
    TreeT, QuadT = pXest_type isa GridapP4est.P4estType ?
        (P4est_wrapper.p4est_tree_t, P4est_wrapper.p4est_quadrant_t) :
        (P4est_wrapper.p8est_tree_t, P4est_wrapper.p8est_quadrant_t)

    TInt      = JX.TInt
    forest    = unsafe_load(ptr_pXest)
    trees_arr = unsafe_load(forest.trees)          # sc_array_t of {p4est,p8est}_tree_t
    levels    = TInt[]
    for t in forest.first_local_tree:forest.last_local_tree
        tree_ptr = Ptr{TreeT}(
            trees_arr.array + t * trees_arr.elem_size)
        tree   = unsafe_load(tree_ptr)
        n_quads = Int(tree.quadrants.elem_count)
        for q in 0:n_quads-1
            quad_ptr = Ptr{QuadT}(
                tree.quadrants.array + q * tree.quadrants.elem_size)
            quad = unsafe_load(quad_ptr)
            push!(levels, TInt(quad.level))
        end
    end
    return levels
end

function JX.load_p4est_checkpoint_model(base_model, forest_file::String)
    pXest_type = base_model.pXest_type
    parts      = base_model.parts

    # Load forest (MPI-collective). p4est_load/p8est_load also fill
    # *connectivity_ref with a freshly allocated connectivity that we leave
    # to be GCed — we use the Gridap-managed connectivity from base_model
    # throughout.
    # Dispatch on 2D (p4est_load) vs 3D (p8est_load) to avoid passing the
    # wrong struct type — mirrors write_p4est_checkpoint's save-side dispatch.
    if pXest_type isa GridapP4est.P4estType
        connectivity_ref = Ref{Ptr{P4est_wrapper.p4est_connectivity_t}}()
        loaded_ptr_pXest = P4est_wrapper.p4est_load(
            forest_file,
            JX.get_mpi_comm(),
            Csize_t(0),      # no per-quadrant data stored
            Cint(0),         # do not read payload
            C_NULL,
            connectivity_ref)
    else
        connectivity_ref = Ref{Ptr{P4est_wrapper.p8est_connectivity_t}}()
        loaded_ptr_pXest = P4est_wrapper.p8est_load(
            forest_file,
            JX.get_mpi_comm(),
            Csize_t(0),      # no per-quadrant data stored
            Cint(0),         # do not read payload
            C_NULL,
            connectivity_ref)
    end

    # Ghost layer and lnodes for a non-conforming (AMR) forest
    ptr_ghost  = GridapP4est.setup_pXest_ghost(pXest_type, loaded_ptr_pXest)
    ptr_lnodes = GridapP4est.setup_pXest_lnodes_nonconforming(pXest_type, loaded_ptr_pXest, ptr_ghost)

    # Build Gridap distributed mesh from the loaded forest.
    # base_model.coarse_model (GmshDiscreteModel) provides physical coordinates.
    fmodel, nc_glue = GridapP4est.setup_non_conforming_distributed_discrete_model(
        pXest_type,
        GridapP4est.PXestUniformRefinementRuleType(),
        parts,
        base_model.coarse_model,
        base_model.ptr_pXest_connectivity,
        loaded_ptr_pXest,
        ptr_ghost,
        ptr_lnodes)

    GridapP4est.pXest_ghost_destroy(pXest_type, ptr_ghost)
    GridapP4est.pXest_lnodes_destroy(pXest_type, ptr_lnodes)

    Dc = JX.num_cell_dims(base_model.dmodel)
    Dp = JX.num_point_dims(base_model.dmodel)

    return GridapP4est.OctreeDistributedDiscreteModel(Dc, Dp,
        parts,
        fmodel,
        nc_glue,
        base_model.coarse_model,
        base_model.ptr_pXest_connectivity,
        loaded_ptr_pXest,
        pXest_type,
        GridapP4est.PXestUniformRefinementRuleType(),
        false,         # does not own connectivity (base_model owns it)
        base_model)    # gc_ref: keep base_model alive
end
