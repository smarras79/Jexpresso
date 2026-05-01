# module JeGeometry
# using Gridap
# using Gridap.Arrays
# using Gridap.Arrays: Table
# using Gridap.Geometry
# using Gridap.Fields
# using Gridap.ReferenceFEs
# using Gridap.CellData
# using Gridap.Visualization
# using Gridap.Geometry: GridMock
# using GridapDistributed
# using GridapDistributed: GenericDistributedDiscreteModel, compute_cell_graph
# using PartitionedArrays
# using GridapGmsh
# using GridapP4est
# using P4est_wrapper
# using SparseArrays

# const DistributedVisualizationData = GridapDistributed.DistributedVisualizationData


# function Visualization.visualization_data(
#     model::GenericDistributedDiscreteModel{Dc},
#     filebase::AbstractString;
#     labels=get_face_labeling(model)) where Dc
  
#     cell_gids = get_cell_gids(model)
#     fact_gids = get_face_gids(model,Dc-1)
#     vd = map(local_views(model),partition(cell_gids), partition(fact_gids),labels.labels) do model,gids,fgids,labels
#       part = part_id(gids)
#       vd = visualization_data(model,filebase;labels=labels)
#       vd_cells = vd[end]
#       vd_facets = vd[Dc]
#       # @info part, size(vd)
#       push!(vd_cells.celldata, "gid" => local_to_global(gids))
#       push!(vd_facets.celldata, "fgid" => local_to_global(fgids))
#       push!(vd_cells.celldata, "part" => local_to_owner(gids))
#       # @info part, vd[end].celldata, vd_facets.celldata
#       vd
#     end
#     r = []
#     for i in 0:Dc
#       push!(r,DistributedVisualizationData(map(x->x[i+1],vd)))
#     end
#     r
# end

# end  # End of module JeGeometry


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

struct NoGhostParts{T<:AbstractVector} <: AbstractVector{eltype(T)}
    data::T
end
Base.size(p::NoGhostParts) = size(p.data)
Base.getindex(p::NoGhostParts, i::Integer) = p.data[i]
# Override map to use the underlying array's distributed map (avoids MPIArray scalar indexing)
Base.map(f::Function, p::NoGhostParts) = map(f, p.data)
Base.map(f::Function, p::NoGhostParts, args...) = map(f, p.data, args...)

# Not type piracy — dispatches on NoGhostParts which belongs to this module.
# Omits Gridap's neighbor loop that adds ghost/halo cells to each partition.
function Gridap.Geometry.DiscreteModel(
    parts::NoGhostParts,
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



