struct NoEmbedMeshFile
    path::String
end

function _gmsh_model_no_embed(gmsh; has_affine_map=nothing)
    Dc = GridapGmsh._setup_cell_dim(gmsh)
    Dp = Dc  # force Dp=Dc: avoids UnstructuredDiscreteModel{2,3} which breaks P4est
    node_to_coords = GridapGmsh._setup_node_coords(gmsh, Dp)
    vertex_to_node, node_to_vertex = GridapGmsh._setup_nodes_and_vertices(gmsh, node_to_coords)
    grid, cell_to_entity = GridapGmsh._setup_grid(gmsh, Dc, Dp, node_to_coords, node_to_vertex; has_affine_map)
    cell_to_vertices, vertex_to_node, node_to_vertex = GridapGmsh._setup_cell_to_vertices(grid, vertex_to_node, node_to_vertex)
    grid_topology = Gridap.Geometry.UnstructuredGridTopology(grid, cell_to_vertices, vertex_to_node)
    labeling = GridapGmsh._setup_labeling(gmsh, grid, grid_topology, cell_to_entity, vertex_to_node, node_to_vertex)
    Gridap.Geometry.UnstructuredDiscreteModel(grid, grid_topology, labeling)
end

# Not type piracy — NoEmbedMeshFile belongs to this module.
function GridapGmsh.GmshDiscreteModel(f::NoEmbedMeshFile; renumber=true, has_affine_map=nothing, kwargs...)
    GridapGmsh.GMSH_FOUND || error("GridapGmsh not loaded properly.")
    isfile(f.path) || error("Msh file not found: $(f.path)")
    gmsh = GridapGmsh.gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.option.setNumber("Mesh.MedImportGroupsOfNodes", 1)
    gmsh.open(f.path)
    renumber && gmsh.model.mesh.renumberNodes()
    renumber && gmsh.model.mesh.renumberElements()
    model = _gmsh_model_no_embed(gmsh; has_affine_map)
    gmsh.finalize()
    model
end

# Override of GridapGmsh.GmshDiscreteModel(gmsh::Module) restored from
# hw/giga_les. Forces Dp = Dc (embedding dim = cell dim).
#
# Without this override, GridapGmsh's default GmshDiscreteModel(gmsh::Module)
# reads the .msh's stored coordinate dimension as Dp - so a 2D mesh saved
# with z=0 columns becomes UnstructuredDiscreteModel{2,3}, not {2,2}. That
# silently propagates a 3D coordinate slot through node_to_coords, the
# grid, and downstream metric / DSS / boundary code paths. In serial the
# extra 3rd component is just zero and most computations still work; in
# parallel, GridapDistributed and PartitionedArrays make different
# decisions about ghost/halo construction and face-partition equality
# checks based on Dp, which manifests as accumulating DSS drift at
# interface DOFs (theta blowing up at t=100s with rho*theta negative).
#
# This override is invoked by the internal call chain underneath BOTH
# the serial GmshDiscreteModel(file) and the parallel
# GmshDiscreteModel(parts, file) constructors, so installing it here
# fixes both code paths.
function GridapGmsh.GmshDiscreteModel(gmsh::Module; has_affine_map=nothing, orient_if_simplex=nothing)
    Dc = GridapGmsh._setup_cell_dim(gmsh)
    Dp = Dc
    node_to_coords = GridapGmsh._setup_node_coords(gmsh, Dp)
    vertex_to_node, node_to_vertex = GridapGmsh._setup_nodes_and_vertices(gmsh, node_to_coords)
    grid, cell_to_entity = GridapGmsh._setup_grid(gmsh, Dc, Dp, node_to_coords, node_to_vertex; has_affine_map)
    cell_to_vertices, vertex_to_node, node_to_vertex = GridapGmsh._setup_cell_to_vertices(grid, vertex_to_node, node_to_vertex)
    grid_topology = Gridap.Geometry.UnstructuredGridTopology(grid, cell_to_vertices, vertex_to_node)
    labeling = GridapGmsh._setup_labeling(gmsh, grid, grid_topology, cell_to_entity, vertex_to_node, node_to_vertex)
    Gridap.Geometry.UnstructuredDiscreteModel(grid, grid_topology, labeling)
end

struct NoGhostParts{T<:AbstractVector} <: AbstractVector{eltype(T)}
    data::T
end
Base.size(p::NoGhostParts) = size(p.data)
Base.getindex(p::NoGhostParts, i::Integer) = p.data[i]
# Override map to use the underlying array's distributed map (avoids MPIArray scalar indexing)
Base.map(f::Function, p::NoGhostParts) = map(f, p.data)
Base.map(f::Function, p::NoGhostParts, args...) = map(f, p.data, args...)

# Override Gridap.Geometry.DiscreteModel(parts, smodel, cell_to_part) to
# build a partition WITHOUT ghost cells (drops Gridap's neighbor loop that
# adds ghost/halo cells to each part). Matches hw/giga_les behavior.
#
# Dispatching on parts::AbstractArray (not NoGhostParts) means this method
# also intercepts the internal Gridap call made by the parallel
# GmshDiscreteModel(parts, file) constructor in GridapGmsh — without this
# broader dispatch, that path silently uses Gridap's stock partitioner
# WITH ghost cells, which then makes every interface DOF receive an extra
# contribution at DSS time (owner's element + ghost element's element),
# accumulating drift that blows up runs in parallel only. See
# https://claude.ai/code for the diagnostic chain.
function Gridap.Geometry.DiscreteModel(
    parts::AbstractArray,
    model::Gridap.Geometry.DiscreteModel,
    cell_to_part::AbstractArray,
    cell_graph::SparseArrays.SparseMatrixCSC = GridapDistributed.compute_cell_graph(model)
)
    ncells = Gridap.Geometry.num_cells(model)
    @assert length(cell_to_part) == ncells
    @assert size(cell_graph, 1) == ncells
    @assert size(cell_graph, 2) == ncells

    lcell_to_cell, lcell_to_part, _ = map(parts) do part
        cell_to_mask = fill(false, ncells)
        for icell in 1:ncells
            if cell_to_part[icell] == part
                cell_to_mask[icell] = true
            end
        end
        lcell_to_cell = findall(cell_to_mask)
        lcell_to_part = zeros(Int32, length(lcell_to_cell))
        lcell_to_part .= cell_to_part[lcell_to_cell]
        lcell_to_cell, lcell_to_part, cell_to_part
    end |> PartitionedArrays.tuple_of_arrays

    partition = map(parts, lcell_to_cell, lcell_to_part) do part, lcell_to_cell, lcell_to_part
        PartitionedArrays.LocalIndices(ncells, part, lcell_to_cell, lcell_to_part)
    end

    PartitionedArrays.assembly_neighbors(partition; symmetric=true)
    gids = PartitionedArrays.PRange(partition)

    models = map(lcell_to_cell) do lcell_to_cell
        Gridap.Geometry.DiscreteModelPortion(model, lcell_to_cell)
    end

    GridapDistributed.GenericDistributedDiscreteModel(models, gids)
end


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


