# PERF: this file used to be wrapped in `module JeGeometry … end` with
# `__precompile__(false)` at the top, which disabled precompilation
# for the ENTIRE Jexpresso package (the toplevel `using Jexpresso`
# then fell through to "Skipping precompilation due to precompilable
# error" and re-loaded all 80+ source files at REPL launch). The
# wrapper bought us nothing — every name defined below is either a
# method on an external type (Gridap.Geometry.DiscreteModel,
# GridapGmsh.GmshDiscreteModel, Visualization.visualization_data) or
# a global const for a Gridap internal — none of them needed
# namespace isolation, and no caller outside this file referenced
# `JeGeometry.<x>` qualified. Flattening lets the rest of Jexpresso
# precompile cleanly while preserving the method definitions.
#
# Most of the symbols needed below are already in scope in the
# top-level Jexpresso module. The one extra pulled in here is the
# GridapDistributed internals (GenericDistributedDiscreteModel,
# compute_cell_graph) that aren't at the top of src/Jexpresso.jl.
using GridapDistributed: GenericDistributedDiscreteModel, compute_cell_graph

# Jexpresso's parallel-partition constructor for a
# `Gridap.Geometry.DiscreteModel`.  This used to be defined as a
# method on the imported `Gridap.Geometry.DiscreteModel` constructor
# itself (type piracy), and it had EXACTLY the same 3-arg positional
# signature as GridapDistributed's own
# `DiscreteModel(::AbstractArray, ::DiscreteModel, ::AbstractArray)`
# (RrOp1/src/Geometry.jl:377). Julia ≥ 1.10 refuses to precompile a
# module that overwrites an existing method, so the only way to keep
# the override was `__precompile__(false)` at the top of the (former)
# JeGeometry submodule — which killed precompilation of the entire
# Jexpresso package.
#
# Rename to a Jexpresso-local function (`je_DiscreteModel`) so it no
# longer collides with GridapDistributed's method table. The two call
# sites in src/kernel/mesh/mesh.jl that previously relied on dispatch
# now invoke `je_DiscreteModel(...)` explicitly.
function je_DiscreteModel(
    parts::AbstractArray,
    model::Geometry.DiscreteModel,
    cell_to_part::AbstractArray,
    cell_graph::SparseMatrixCSC = compute_cell_graph(model),
    new_param::Int = 1  # kept for caller-side compat; not used below
)
    ncells = num_cells(model)
    @assert length(cell_to_part) == ncells
    @assert size(cell_graph, 1) == ncells
    @assert size(cell_graph, 2) == ncells

    lcell_to_cell, lcell_to_part, gid_to_part = map(parts) do part
        cell_to_mask = fill(false, ncells)
        icell_to_jcells_ptrs = cell_graph.colptr
        icell_to_jcells_data = cell_graph.rowval
        for icell in 1:ncells
            if cell_to_part[icell] == part
                cell_to_mask[icell] = true
            end
        end
        lcell_to_cell = findall(cell_to_mask)
        lcell_to_part = zeros(Int32, length(lcell_to_cell))
        lcell_to_part .= cell_to_part[lcell_to_cell]
        lcell_to_cell, lcell_to_part, cell_to_part
    end |> tuple_of_arrays

    partition = map(parts, lcell_to_cell, lcell_to_part) do part, lcell_to_cell, lcell_to_part
        LocalIndices(ncells, part, lcell_to_cell, lcell_to_part)
    end

    assembly_neighbors(partition; symmetric=true)

    gids = PRange(partition)

    models = map(lcell_to_cell) do lcell_to_cell
        DiscreteModelPortion(model, lcell_to_cell)
    end

    # Incorporate new_param into the logic if needed
    # if new_param > 1
    #println("New parameter is greater than 1: ", new_param)
    # end

    GenericDistributedDiscreteModel(models, gids)
end

# function Geometry.DiscreteModelPortion(dmodel::DistributedDiscreteModel, cell_to_parent_cell::AbstractArray, parts::AbstractArray)

#     # gtopo = get_grid_topology(dmodel)
#     # glabels = get_face_labeling(dmodel)
#     # map(parts, dmodel, gtopo, glabels) do part, model, topo, lables

#     #     grid_p =  GridPortion(get_grid(model),cell_to_parent_cell)
#     #     topo_p, d_to_dface_to_parent_dface = _grid_topology_portion(topo,cell_to_parent_cell)
#     #     labels_p = _setup_labels_p(labels, d_to_dface_to_parent_dface)
#     #     model_p = DiscreteModel(grid_p,topo_p,labels_p)
#     #     DiscreteModelPortion(model_p,model,d_to_dface_to_parent_dface)
#     # end
# end


# ─── REMOVED: two more type-piracy overrides ──────────────────────────
#
# (1) GridapGmsh.GmshDiscreteModel(::Module; …) — overwrote
#     GridapGmsh/ZGvkt/src/GmshDiscreteModels.jl:30. No call site
#     inside Jexpresso ever invoked the `Module`-argument form (we
#     always pass either a filename string or
#     `(parts, filename)`), so the override was unreachable dead code
#     that just blocked precompilation. Dropped along with the six
#     `const _setup_* = GridapGmsh._setup_*` aliases that only the
#     override used.
#
# (2) Visualization.visualization_data(::GenericDistributedDiscreteModel{Dc},
#                                      ::AbstractString; labels=…) — overwrote
#     Gridap's distributed-model visualization hook to inject three
#     extra cell-data fields ("gid", "fgid", "part") into the VTU
#     output. The one reachable call site, `writevtk(partitioned_model,
#     "./refine/")` in src/kernel/Adaptivity/Projection.jl:799, is a
#     debug dump that only fires when `inputs[:linitial_refine] == true`
#     and falls back cleanly to GridapDistributed's default
#     `visualization_data` (refined-mesh VTU still written; missing
#     only the three debug fields).
#
# Both were guarded by the former `__precompile__(false)` at the top
# of the (now-flattened) `JeGeometry` submodule; with that opt-out
# gone, Julia ≥ 1.10 refuses to precompile the package.  If anyone
# later needs the extra cell-data fields back, the cleanest path is
# to define a Jexpresso-local `writevtk_with_part_info` that calls
# Gridap's default `visualization_data` and then mutates the returned
# tuple — same effect, no piracy.
# ─────────────────────────────────────────────────────────────────────

# Removed: `const DistributedVisualizationData = …` — the only
# reference was inside the now-removed `visualization_data` override.




# (no `end` here — the former `module JeGeometry … end` wrapper was
# removed; everything below now lives in the top-level Jexpresso
# module.)


function get_boundary_cells(model,nsd)
  facet_dim = 1
  if nsd == 3
    facet_dim = 2
  end
  bdy_facets = get_boundary_faces(model,nsd,facet_dim)
  topo = get_grid_topology(model)
  facet_cell_id = get_faces(topo, facet_dim, nsd)
  bdy_cells = facet_cell_id[bdy_facets]
  ret = unique(vcat(bdy_cells...))
  # @info ret``
  return ret
  # @info facet_cell_id
end

function get_boundary_faces(model,nsd,dim)
  if (dim == nsd)
    return get_boundary_cells(model,nsd)
  end
  labels = get_face_labeling(model)
  if nsd == 3
    Base.filter!(x -> !(x in ["internal", "hanging"]), labels.tag_to_name)
  elseif nsd == 2
    Base.filter!(x -> !(x in ["domain", "hanging"]), labels.tag_to_name)
  end
  facet_to_tag = get_face_tag_index(labels,labels.tag_to_name,dim)
  # @info facet_to_tag, labels.tag_to_name,  length(findall(x -> x>0, facet_to_tag))
  findall(x -> x>0, facet_to_tag)
end



