module MyGeometry

using MPI
using Gridap
using Gridap.Arrays
using Gridap.Arrays: Table
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Visualization
using Gridap.Geometry: GridMock
using GridapDistributed
using GridapDistributed: GenericDistributedDiscreteModel
using PartitionedArrays
using GridapGmsh
using GridapP4est
using P4est_wrapper
using Metis
using SparseArrays


# Define your custom version of DiscreteModel function
function Geometry.DiscreteModel(
    parts::AbstractArray,
    model::Geometry.DiscreteModel,
    cell_to_part::AbstractArray,
    cell_graph::SparseMatrixCSC = compute_cell_graph(model),
    new_param::Int = 1  # Add your new parameter here with a default value
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
                # pini = icell_to_jcells_ptrs[icell]
                # pend = icell_to_jcells_ptrs[icell + 1] - 1
                # for p in pini:pend
                #     jcell = icell_to_jcells_data[p]
                #     cell_to_mask[jcell] = true
                # end
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
        println("New parameter is greater than 1: ", new_param)
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


const _setup_cell_dim = GridapGmsh._setup_cell_dim
const _setup_node_coords = GridapGmsh._setup_node_coords
const _setup_nodes_and_vertices = GridapGmsh._setup_nodes_and_vertices
const _setup_cell_to_vertices = GridapGmsh._setup_cell_to_vertices
const _setup_grid =GridapGmsh._setup_grid
const _setup_labeling = GridapGmsh._setup_labeling

function GridapGmsh.GmshDiscreteModel(gmsh::Module)

    Dc = _setup_cell_dim(gmsh)
    Dp = Dc
    # @info Dp
    node_to_coords = _setup_node_coords(gmsh,Dp)
    nnodes = length(node_to_coords)
    vertex_to_node, node_to_vertex = _setup_nodes_and_vertices(gmsh,node_to_coords)
    grid, cell_to_entity = _setup_grid(gmsh,Dc,Dp,node_to_coords,node_to_vertex)
    cell_to_vertices, vertex_to_node, node_to_vertex = _setup_cell_to_vertices(grid,vertex_to_node,node_to_vertex)
    grid_topology = UnstructuredGridTopology(grid,cell_to_vertices,vertex_to_node)
    labeling = _setup_labeling(gmsh,grid,grid_topology,cell_to_entity,vertex_to_node,node_to_vertex)
    UnstructuredDiscreteModel(grid,grid_topology,labeling)
end

const DistributedVisualizationData = GridapDistributed.DistributedVisualizationData


function Visualization.visualization_data(
    model::GenericDistributedDiscreteModel{Dc},
    filebase::AbstractString;
    labels=get_face_labeling(model)) where Dc
  
    cell_gids = get_cell_gids(model)
    fact_gids = get_face_gids(model,Dc-1)
    vd = map(local_views(model),partition(cell_gids), partition(fact_gids),labels.labels) do model,gids,fgids,labels
      part = part_id(gids)
      vd = visualization_data(model,filebase;labels=labels)
      vd_cells = vd[end]
      vd_facets = vd[Dc]
      # @info part, size(vd)
      push!(vd_cells.celldata, "gid" => local_to_global(gids))
      push!(vd_facets.celldata, "fgid" => local_to_global(fgids))
      push!(vd_cells.celldata, "part" => local_to_owner(gids))
      # @info part, vd[end].celldata, vd_facets.celldata
      vd
    end
    r = []
    for i in 0:Dc
      push!(r,DistributedVisualizationData(map(x->x[i+1],vd)))
    end
    r
  end



const _unwrap_ghost_quadrants = GridapP4est._unwrap_ghost_quadrants
const unsafe_wrap = GridapP4est.unsafe_wrap
const _unwrap_global_first_quadrant = GridapP4est._unwrap_global_first_quadrant
const PXestType = GridapP4est.PXestType



# function GridapP4est.setup_cell_prange(pXest_type::PXestType,
#                            parts::AbstractVector{<:Integer},
#                            ptr_pXest,
#                            ptr_pXest_ghost)
#   comm = parts.comm
#   # @info "setup_cell_prange"

#   pXest_ghost = ptr_pXest_ghost[]
#   pXest       = ptr_pXest[]

#   # Obtain ghost quadrants
#   ptr_ghost_quadrants = _unwrap_ghost_quadrants(pXest_type,pXest_ghost)
#   proc_offsets = unsafe_wrap(Array, pXest_ghost.proc_offsets, pXest_ghost.mpisize+1)

#   global_first_quadrant = _unwrap_global_first_quadrant(pXest_type,pXest)

#   noids,firstgid,gho_to_glo,gho_to_own=map(parts) do part
#     gho_to_glo = Vector{Int}(undef, pXest_ghost.ghosts.elem_count)
#     gho_to_own = Vector{Int32}(undef, pXest_ghost.ghosts.elem_count)
#     k=1
#     for i=1:pXest_ghost.mpisize
#       for j=proc_offsets[i]:proc_offsets[i+1]-1
#         quadrant       = ptr_ghost_quadrants[j+1]
#         piggy3         = quadrant.p.piggy3
#         gho_to_glo[k]  = global_first_quadrant[i]+piggy3.local_num+1
#         gho_to_own[k] = Int32(i)
#         k=k+1
#       end
#     end
#     global_first_quadrant[part+1]-global_first_quadrant[part],global_first_quadrant[part]+1,gho_to_glo,gho_to_own
#   end |> tuple_of_arrays
#   ngids = global_first_quadrant[end]

#   partition=map(parts,noids,firstgid,gho_to_glo,gho_to_own) do part, noids, firstgid, gho_to_glo, gho_to_own
#     owner = part
#     own_indices=OwnIndices(ngids,owner,(collect(firstgid:firstgid+noids-1)))
#     ghost_indices=GhostIndices(ngids)
#     partition = OwnAndGhostIndices(own_indices,ghost_indices)
#   end
#   # This is required to provide the hint that the communication 
#   # pattern underlying partition is symmetric, so that we do not have 
#   # to execute the algorithm the reconstructs the reciprocal in the 
#   # communication graph
#   assembly_neighbors(partition;symmetric=true)
#   PRange(partition)

# end


# const fetch_vector_ghost_values_cache = GridapP4est.fetch_vector_ghost_values_cache
# const fetch_vector_ghost_values! = GridapP4est.fetch_vector_ghost_values!


# function GridapP4est.generate_cell_vertex_gids(ptr_pXest_lnodes, cell_prange)
#   pXest_lnodes = ptr_pXest_lnodes[]

#   nvertices = pXest_lnodes.vnodes
#   element_nodes = unsafe_wrap(Array,
#                               pXest_lnodes.element_nodes,
#                               pXest_lnodes.num_local_elements*nvertices)

#   nonlocal_nodes = unsafe_wrap(Array,
#                                pXest_lnodes.nonlocal_nodes,
#                                pXest_lnodes.num_local_nodes-pXest_lnodes.owned_count)

#   k = 1
#   cell_vertex_gids = map(partition(cell_prange)) do indices
#     n = length(local_to_own(indices))
#     ptrs = Vector{Int32}(undef,n+1)
#     ptrs[1] = 1
#     for i = 1:n
#       ptrs[i+1] = ptrs[i] + nvertices
#     end
#     k = 1
#     current = 1
#     data = Vector{Int}(undef, ptrs[n+1]-1)
#     for i = 1:pXest_lnodes.num_local_elements
#       for j = 1:nvertices
#         l = element_nodes[k+j-1]
#         if (l < pXest_lnodes.owned_count)
#           data[current] = pXest_lnodes.global_offset+l+1
#         else
#           data[current] = nonlocal_nodes[l-pXest_lnodes.owned_count+1]+1
#         end
#         current = current+1
#       end
#       k = k + nvertices
#     end
#     JaggedArray(data,ptrs)
#   end
#   fetch_cache = fetch_vector_ghost_values_cache(cell_vertex_gids,partition(cell_prange))
#   fetch_vector_ghost_values!(cell_vertex_gids,fetch_cache) |> wait
#   return cell_vertex_gids
# end


# # function setup_distributed_discrete_model(pXest_type::PXestType,
# #                                           parts,
# #                                           coarse_discrete_model,
# #                                           ptr_pXest_connectivity,
# #                                           ptr_pXest,
# #                                           ptr_pXest_ghost,
# #                                           ptr_pXest_lnodes)
  
# #    cell_prange = setup_cell_prange(pXest_type,parts,ptr_pXest,ptr_pXest_ghost)
# #    cell_vertex_gids=generate_cell_vertex_gids(ptr_pXest_lnodes,cell_prange)
# #    cell_corner_lids = generate_cell_corner_lids(cell_vertex_gids)
# #    cell_vertex_coordinates = generate_cell_vertex_coordinates(pXest_type,
# #                                               cell_corner_lids,
# #                                               ptr_pXest_connectivity,
# #                                               ptr_pXest,
# #                                               ptr_pXest_ghost)

# #    grid,topology = generate_grid_and_topology(pXest_type,
# #                                               cell_corner_lids,
# #                                               cell_vertex_coordinates)

# #    face_labeling = generate_face_labeling(pXest_type,
# #                                           parts,
# #                                           cell_prange,
# #                                           coarse_discrete_model,
# #                                           topology,
# #                                           ptr_pXest,
# #                                           ptr_pXest_ghost)

# #    local_models = map(grid,topology,face_labeling) do grid, topology, face_labeling
# #      Gridap.Geometry.UnstructuredDiscreteModel(grid,topology,face_labeling)
# #    end
# #    return GridapDistributed.DistributedDiscreteModel(local_models,cell_prange)
# # end


# const num_cell_dims = GridapP4est.num_cell_dims
# const pXest_tree_array_index = GridapP4est.pXest_tree_array_index
# const pXest_quadrant_array_index = GridapP4est.pXest_quadrant_array_index
# const pXest_num_quadrant_layers = GridapP4est.pXest_num_quadrant_layers
# const pXest_get_layer = GridapP4est.pXest_get_layer
# const pXest_cell_coords = GridapP4est.pXest_cell_coords
# const pXest_get_quadrant_and_layer_levels = GridapP4est.pXest_get_quadrant_and_layer_levels
# const pXest_get_quadrant_vertex_coordinates = GridapP4est.pXest_get_quadrant_vertex_coordinates



# """
#   Generate the geometrical cellwise vertex coordinates.
# """
# function GridapP4est.generate_cell_vertex_coordinates(pXest_type::PXestType,
#                                           cell_vertex_lids,
#                                           ptr_pXest_connectivity,
#                                           ptr_pXest,
#                                           ptr_pXest_ghost)

#   Dc = num_cell_dims(pXest_type)
#   # @info "generate_cell_vertex_coordinates"

#   PXEST_CORNERS = 2^Dc
#   pXest_ghost = ptr_pXest_ghost[]
#   pXest       = ptr_pXest[]

#   cell_vertex_coordinates = map(cell_vertex_lids) do cell_vertex_lids
#     data = Vector{Point{Dc,Float64}}(undef,length(cell_vertex_lids.data))
#     current = 1
#     vxy = Vector{Cdouble}(undef,Dc)
#     pvxy = pointer(vxy,1)
#     for itree = 1:pXest_ghost.num_trees
#       tree = pXest_tree_array_index(pXest_type, pXest, itree-1)[]
#       for cell=1:tree.quadrants.elem_count
#         quadrant=pXest_quadrant_array_index(pXest_type,tree,cell-1)[] 
#         # Loop over layers in the current column
#         for l=1:pXest_num_quadrant_layers(pXest_type,quadrant)
#           layer=pXest_get_layer(pXest_type, quadrant, pXest, l-1)
#           coords=pXest_cell_coords(pXest_type,quadrant,layer)
#           levels=pXest_get_quadrant_and_layer_levels(pXest_type,quadrant,layer)
          
#           for vertex=1:PXEST_CORNERS
#             pXest_get_quadrant_vertex_coordinates(pXest_type,
#                                                   ptr_pXest_connectivity,
#                                                   p4est_topidx_t(itree-1),
#                                                   coords,
#                                                   levels,
#                                                   Cint(vertex-1),
#                                                   pvxy)
#             data[current] = Point{Dc,Float64}(vxy...)
#             current=current+1
#           end
#         end 
#       end
#     end

#     ptrs = copy(cell_vertex_lids.ptrs)
#     return Gridap.Arrays.Table(data,ptrs)
#   end
#   return cell_vertex_coordinates
# end



# const P4P8estType = GridapP4est.P4P8estType
# const update_face_to_entity_with_ghost_data! = GridapP4est.update_face_to_entity_with_ghost_data!
# const P4EST_2_GRIDAP_FACET_2D = GridapP4est.P4EST_2_GRIDAP_FACET_2D
# const P4EST_2_GRIDAP_FACET_3D = GridapP4est.P4EST_2_GRIDAP_FACET_3D
# const cell_to_faces = GridapP4est.cell_to_faces


# const ITERATOR_RESTRICT_TO_BOUNDARY=Cint(100)
# const ITERATOR_RESTRICT_TO_INTERIOR=Cint(101)


# function GridapP4est.generate_face_labeling(pXest_type::P4P8estType,
#                                 parts,
#                                 cell_prange,
#                                 coarse_discrete_model::DiscreteModel{Dc,Dp},
#                                 topology,
#                                 ptr_pXest,
#                                 ptr_pXest_ghost) where {Dc,Dp}

#   pXest       = ptr_pXest[]
#   pXest_ghost = ptr_pXest_ghost[]

#   coarse_grid_topology = Gridap.Geometry.get_grid_topology(coarse_discrete_model)
#   coarse_grid_labeling = Gridap.Geometry.get_face_labeling(coarse_discrete_model)
#   coarse_cell_vertices = Gridap.Geometry.get_faces(coarse_grid_topology,Dc,0)
#   if (Dc==3)
#     coarse_cell_edgets = Gridap.Geometry.get_faces(coarse_grid_topology,Dc,1)
#   end
#   coarse_cell_facets   = Gridap.Geometry.get_faces(coarse_grid_topology,Dc,Dc-1)

#   owned_trees_offset=Vector{Int}(undef,pXest_ghost.num_trees+1)
#   owned_trees_offset[1]=0
#   for itree=1:pXest_ghost.num_trees
#     if Dc==2
#       tree = p4est_tree_array_index(pXest.trees, itree-1)[]
#     else
#       tree = p8est_tree_array_index(pXest.trees, itree-1)[]
#     end
#     owned_trees_offset[itree+1]=owned_trees_offset[itree]+tree.quadrants.elem_count
#   end
#   # @info "owned_trees_offset", owned_trees_offset
#  faces_to_entity=map(topology) do topology
#      num_cells = Gridap.Geometry.num_faces(topology,Dc)
#      # Iterate over corners
#      num_vertices=Gridap.Geometry.num_faces(topology,0)
#      vertex_to_entity=zeros(Int,num_vertices)
#      cell_vertices=Gridap.Geometry.get_faces(topology,Dc,0)

#      # Corner iterator callback
#      function jcorner_callback(pinfo     :: Ptr{p8est_iter_corner_info_t},
#                                user_data :: Ptr{Cvoid})
#         info=pinfo[]
#         if (Dc==2)
#           sides=Ptr{p4est_iter_corner_side_t}(info.sides.array)
#           CONNECT_CORNER=P4est_wrapper.P4EST_CONNECT_CORNER
#         else
#           sides=Ptr{p8est_iter_corner_side_t}(info.sides.array)
#           CONNECT_CORNER=P4est_wrapper.P8EST_CONNECT_CORNER
#         end
#         nsides=info.sides.elem_count
#         tree=sides[1].treeid+1
#         data=sides[1]
#         if data.is_ghost==1
#            ref_cell=pXest.local_num_quadrants+data.quadid+1
#         else
#            ref_cell=owned_trees_offset[tree]+data.quadid+1
#         end
#         if ref_cell > num_cells
#             return nothing
#         end
#         corner=sides[1].corner+1
#         ref_cornergid=cell_vertices[ref_cell][corner]
#         if (info.tree_boundary!=0 && info.tree_boundary==CONNECT_CORNER)
#               # The current corner is also a corner of the coarse mesh
#               coarse_cornergid=coarse_cell_vertices[tree][corner]
#               vertex_to_entity[ref_cornergid]=
#                  coarse_grid_labeling.d_to_dface_to_entity[1][coarse_cornergid]
#               @debug "[GLOBAL CORNER] vertex_to_entity[$(ref_cornergid)]=$(coarse_grid_labeling.d_to_dface_to_entity[1][coarse_cornergid])"
#         else
#           if vertex_to_entity[ref_cornergid]==0
#             # We are on the interior of a tree (if we did not touch it yet)
#             vertex_to_entity[ref_cornergid]=coarse_grid_labeling.d_to_dface_to_entity[Dc+1][tree]
#             @debug "[INTERIOR CORNER] vertex_to_entity[$(ref_cornergid)]=$(coarse_grid_labeling.d_to_dface_to_entity[Dc+1][tree])"
#           end
#         end
#         nothing
#      end

#      #  C-callable face callback
#      ccorner_callback = @cfunction($jcorner_callback,
#                                    Cvoid,
#                                    (Ptr{p8est_iter_corner_info_t},Ptr{Cvoid}))

#      cell_edgets=Gridap.Geometry.get_faces(topology,Dc,1)
#      num_edgets=Gridap.Geometry.num_faces(topology,1)
#      edget_to_entity=zeros(Int,num_edgets)
#     if (Dc==3)
#       # Edge iterator callback
#       function jedge_callback(pinfo     :: Ptr{p8est_iter_edge_info_t},
#                                user_data :: Ptr{Cvoid})
#         info=pinfo[]
#         sides=Ptr{p8est_iter_edge_side_t}(info.sides.array)
#         function process_edget(tree,edge,ref_cell,info)
#             if ref_cell > num_cells
#               return nothing
#             end
#            polytope=HEX
#            poly_faces=Gridap.ReferenceFEs.get_faces(polytope)
#            poly_edget_range=Gridap.ReferenceFEs.get_dimrange(polytope,1)
#            poly_first_edget=first(poly_edget_range)
#            poly_facet=poly_first_edget+edge-1
#            if (info.tree_boundary!=0 && info.tree_boundary==P4est_wrapper.P8EST_CONNECT_EDGE)
#              coarse_edgetgid=coarse_cell_edgets[tree][edge]
#              coarse_edgetgid_entity=coarse_grid_labeling.d_to_dface_to_entity[2][coarse_edgetgid]
#              # We are on the boundary of coarse mesh or inter-octree boundary
#              for poly_incident_face in poly_faces[poly_facet]
#                if poly_incident_face == poly_facet
#                  ref_edgetgid=cell_edgets[ref_cell][edge]
#                  edget_to_entity[ref_edgetgid]=coarse_edgetgid_entity
#                else
#                  ref_cornergid=cell_vertices[ref_cell][poly_incident_face]
#                  vertex_to_entity[ref_cornergid]=coarse_edgetgid_entity
#                end
#              end
#            else
#              # We are on the interior of the domain if we did not touch the edge yet
#              ref_edgetgid=cell_edgets[ref_cell][edge]
#              if (edget_to_entity[ref_edgetgid]==0)
#                edget_to_entity[ref_edgetgid]=coarse_grid_labeling.d_to_dface_to_entity[Dc+1][tree]
#              end
#            end
#         end 
  
#         nsides=info.sides.elem_count
#         for iside=1:nsides
#          edge=sides[iside].edge+1
#          tree=sides[iside].treeid+1
#          if (sides[iside].is_hanging == 0)
#            data=sides[iside].is.full
#            if data.is_ghost==1
#              ref_cell=pXest.local_num_quadrants+data.quadid+1
#            else
#              ref_cell=owned_trees_offset[tree]+data.quadid+1
#            end 
#            process_edget(tree,edge,ref_cell,info)
#          else 
#            for i=1:length(sides[iside].is.hanging.quadid)
#              quadid=sides[iside].is.hanging.quadid[i]
#              if (sides[iside].is.hanging.is_ghost[i]==1)
#                ref_cell=pXest.local_num_quadrants+quadid+1
#              else 
#                ref_cell=owned_trees_offset[tree]+quadid+1
#              end 
#              process_edget(tree,edge,ref_cell,info)
#            end 
#          end 
#         end
#         nothing
#       end

#       # C-callable edge callback
#       cedge_callback = @cfunction($jedge_callback,
#                                    Cvoid,
#                                    (Ptr{p8est_iter_edge_info_t},Ptr{Cvoid}))
#     end

#     # Iterate over faces
#     num_faces=Gridap.Geometry.num_faces(topology,Dc-1)
#     facet_to_entity=zeros(Int,num_faces)
#     cell_facets=Gridap.Geometry.get_faces(topology,Dc,Dc-1)

#      # Face iterator callback
#      function jface_callback(pinfo     :: Ptr{p8est_iter_face_info_t},
#                              user_data :: Ptr{Cvoid})
#         info=pinfo[]
#         if Dc==2
#           sides=Ptr{p4est_iter_face_side_t}(info.sides.array)
#         else
#           sides=Ptr{p8est_iter_face_side_t}(info.sides.array)
#         end
#         ptr_user_data=Ptr{Cint}(user_data)
#         iterator_mode=unsafe_wrap(Array, ptr_user_data, 1)
#         if (iterator_mode[1]==ITERATOR_RESTRICT_TO_BOUNDARY)
#           # If current face is NOT in the boundary
#           if (info.tree_boundary==0)
#             return nothing
#           end
#         else
#           # If current face is in the boundary 
#           if (info.tree_boundary!=0)
#             return nothing
#           end
#         end

#         function process_facet(tree,face,ref_cell,info)
#           if ref_cell > num_cells
#             return nothing
#           end
#           if Dc==2
#             gridap_facet=P4EST_2_GRIDAP_FACET_2D[face]
#           else
#             gridap_facet=P4EST_2_GRIDAP_FACET_3D[face]
#           end

#           polytope= Dc==2 ? QUAD : HEX
#           poly_faces=Gridap.ReferenceFEs.get_faces(polytope)
#           poly_facet_range=Gridap.ReferenceFEs.get_dimrange(polytope,Dc-1)
#           poly_first_facet=first(poly_facet_range)
#           poly_facet=poly_first_facet+gridap_facet-1

#           if (info.tree_boundary!=0)
#             # We are on the boundary of coarse mesh or inter-octree boundary
#             coarse_facetgid=coarse_cell_facets[tree][gridap_facet]
#             coarse_facetgid_entity=coarse_grid_labeling.d_to_dface_to_entity[Dc][coarse_facetgid]
#           else
#             coarse_facetgid_entity=coarse_grid_labeling.d_to_dface_to_entity[Dc+1][tree]
#           end 

#           for poly_incident_face in poly_faces[poly_facet]
#             if poly_incident_face == poly_facet
#               ref_facetgid=cell_facets[ref_cell][gridap_facet]
#               facet_to_entity[ref_facetgid]=coarse_facetgid_entity
#               @debug "[FACE] info.tree_boundary=$(info.tree_boundary) facet_to_entity[$(ref_facetgid)]=$(coarse_facetgid_entity)"
#             elseif (Dc==3 && poly_incident_face in Gridap.ReferenceFEs.get_dimrange(polytope,1))
#               poly_first_edget=first(Gridap.ReferenceFEs.get_dimrange(polytope,1))
#               edget=poly_incident_face-poly_first_edget+1
#               ref_edgetgid=cell_edgets[ref_cell][edget]
#               if (edget_to_entity[ref_edgetgid]==0)
#                 edget_to_entity[ref_edgetgid]=coarse_facetgid_entity
#                 @debug "[EDGE] info.tree_boundary=$(info.tree_boundary) edget_to_entity[$(ref_edgetgid)]=$(coarse_facetgid_entity)"
#               end
#             else
#               ref_cornergid=cell_vertices[ref_cell][poly_incident_face]
#               if (vertex_to_entity[ref_cornergid]==0)
#                 @debug "[CORNER ON FACE] info.tree_boundary=$(info.tree_boundary) vertex_to_entity[$(ref_cornergid)]=$(coarse_facetgid_entity)"
#                 vertex_to_entity[ref_cornergid]=coarse_facetgid_entity
#               end
#             end
#           end
#         end 

#         nsides=info.sides.elem_count
#         for iside=1:nsides
#           if (iside==2 && 
#             sides[iside].is_hanging == 0 &&
#             sides[1].is_hanging == 0)
#             break
#           end
#           face=sides[iside].face+1
#           tree=sides[iside].treeid+1
#           if (sides[iside].is_hanging == 0)
#             data=sides[iside].is.full
#             if data.is_ghost==1
#               ref_cell=pXest.local_num_quadrants+data.quadid+1
#             else
#               ref_cell=owned_trees_offset[tree]+data.quadid+1
#             end 
#             @debug "nsides=$(nsides) sides[$(iside)].is_hanging == $(sides[iside].is_hanging) process_facet(tree=$(tree),face=$(face),ref_cell=$(ref_cell))"
#             process_facet(tree,face,ref_cell,info)
#           else 
#             for i=1:length(sides[iside].is.hanging.quadid)
#               quadid=sides[iside].is.hanging.quadid[i]
#               if (sides[iside].is.hanging.is_ghost[i]==1)
#                 ref_cell=pXest.local_num_quadrants+quadid+1
#               else 
#                 ref_cell=owned_trees_offset[tree]+quadid+1
#               end
#               @debug "nsides=$(nsides) sides[$(iside)].is_hanging == $(sides[iside].is_hanging) process_facet(tree=$(tree),face=$(face),ref_cell=$(ref_cell))"
#               process_facet(tree,face,ref_cell,info)
#             end
#           end
#         end 
#         nothing
#     end

#     #  C-callable face callback
#     cface_callback = @cfunction($jface_callback,
#                                  Cvoid,
#                                  (Ptr{p8est_iter_face_info_t},Ptr{Cvoid}))

#     # Iterate over cells
#     # num_cells=Gridap.Geometry.num_faces(topology,Dc)
#     cell_to_entity=zeros(Int,num_cells)

#     # Cell iterator callback
#     function jcell_callback(pinfo     :: Ptr{p8est_iter_volume_info_t},
#                             user_data :: Ptr{Cvoid})
#       info=pinfo[]
#       tree=info.treeid+1
#       cell=owned_trees_offset[tree]+info.quadid+1
#       if cell > num_cells
#           return nothing
#       end
#       cell_to_entity[cell]=coarse_grid_labeling.d_to_dface_to_entity[Dc+1][tree]
#       nothing
#     end
#     ccell_callback = @cfunction($jcell_callback,
#                                  Cvoid,
#                                  (Ptr{p8est_iter_volume_info_t},Ptr{Cvoid}))

#     iterator_mode=Ref{Int}(ITERATOR_RESTRICT_TO_BOUNDARY)
#     if (Dc==2)
#        p4est_iterate(ptr_pXest,ptr_pXest_ghost,iterator_mode,C_NULL,cface_callback,C_NULL)
#        p4est_iterate(ptr_pXest,ptr_pXest_ghost,C_NULL,ccell_callback,C_NULL,ccorner_callback)
#        iterator_mode[]=ITERATOR_RESTRICT_TO_INTERIOR
#        p4est_iterate(ptr_pXest,ptr_pXest_ghost,iterator_mode,C_NULL,cface_callback,C_NULL)
#     else
#        p8est_iterate(ptr_pXest,ptr_pXest_ghost,iterator_mode,C_NULL,cface_callback,C_NULL,C_NULL)
#        p8est_iterate(ptr_pXest,ptr_pXest_ghost,iterator_mode,C_NULL,C_NULL,cedge_callback,C_NULL)
#        p8est_iterate(ptr_pXest,ptr_pXest_ghost,C_NULL,ccell_callback,C_NULL,C_NULL,ccorner_callback)
#        iterator_mode[]=ITERATOR_RESTRICT_TO_INTERIOR
#        p8est_iterate(ptr_pXest,ptr_pXest_ghost,iterator_mode,C_NULL,cface_callback,C_NULL,C_NULL)
#     end
#     if (Dc==2)
#       vertex_to_entity, facet_to_entity, cell_to_entity
#     else
#       vertex_to_entity, edget_to_entity, facet_to_entity, cell_to_entity
#     end
#   end

#   vertex_to_entity  = map(x->x[1]   , faces_to_entity)
#   if Dc==3
#     edget_to_entity = map(x->x[2]   , faces_to_entity)
#   end
#   facet_to_entity   = map(x->x[Dc]  , faces_to_entity)
#   cell_to_entity    = map(x->x[Dc+1], faces_to_entity)


#   polytope = Dc==2 ? QUAD : HEX
#   update_face_to_entity_with_ghost_data!(vertex_to_entity,
#                                           cell_prange,
#                                           num_faces(polytope,0),
#                                           cell_to_faces(topology,Dc,0))
#   if Dc==3
#     update_face_to_entity_with_ghost_data!(edget_to_entity,
#                                             cell_prange,
#                                             num_faces(polytope,1),
#                                             cell_to_faces(topology,Dc,1))   
#   end

#   update_face_to_entity_with_ghost_data!(facet_to_entity,
#                                          cell_prange,
#                                          num_faces(polytope,Dc-1),
#                                          cell_to_faces(topology,Dc,Dc-1))
 
#   update_face_to_entity_with_ghost_data!(cell_to_entity,
#                                         cell_prange,
#                                         num_faces(polytope,Dc),
#                                         cell_to_faces(topology,Dc,Dc))

# #  map(vertex_to_entity,facet_to_entity,cell_to_entity) do vertex_to_entity,facet_to_entity,cell_to_entity
# #    @assert all(vertex_to_entity .!= 0)
# #    @assert all(facet_to_entity  .!= 0)
# #    @assert all(cell_to_entity   .!= 0)
# #  end 

#   if (Dc==2)
#     faces_to_entity=[vertex_to_entity,facet_to_entity,cell_to_entity]
#   else
#     faces_to_entity=[vertex_to_entity,edget_to_entity,facet_to_entity,cell_to_entity]
#   end

#   face_labeling = map(faces_to_entity...) do faces_to_entity...
#     d_to_dface_to_entity       = Vector{Vector{Int}}(undef,Dc+1)
#     d_to_dface_to_entity[1]    = faces_to_entity[1]
#     if (Dc==3)
#       d_to_dface_to_entity[2]  = faces_to_entity[2]
#     end
#     d_to_dface_to_entity[Dc]   = faces_to_entity[Dc]
#     d_to_dface_to_entity[Dc+1] = faces_to_entity[Dc+1]
#     Gridap.Geometry.FaceLabeling(d_to_dface_to_entity,
#                                  copy(coarse_grid_labeling.tag_to_entities),
#                                  copy(coarse_grid_labeling.tag_to_name))
#   end
#   face_labeling
# end



# const num_cell_vertices = GridapP4est.num_cell_vertices
# const num_cell_edges = GridapP4est.num_cell_edges
# const num_cell_faces = GridapP4est.num_cell_faces
# const fetch_vector_ghost_values_cache = GridapP4est.fetch_vector_ghost_values_cache
# const hanging_lvertex_within_face_2d = GridapP4est.hanging_lvertex_within_face_2d
# const hanging_lvertex_within_face_3d = GridapP4est.hanging_lvertex_within_face_3d
# const p4est_face_corners = GridapP4est.p4est_face_corners
# const p8est_face_corners = GridapP4est.p8est_face_corners
# const _build_map_from_faces_to_cell_lface = GridapP4est._build_map_from_faces_to_cell_lface
# const _build_map_from_faces_edges_to_cell_lface_ledge = GridapP4est._build_map_from_faces_edges_to_cell_lface_ledge
# const pXest_2_gridap_vertex = GridapP4est.pXest_2_gridap_vertex
# const pXest_2_gridap_facet = GridapP4est.pXest_2_gridap_facet
# const p8est_2_gridap_edge = GridapP4est.p8est_2_gridap_edge
# const _faces_to_cell_element_nodes = GridapP4est._faces_to_cell_element_nodes
# const _edges_to_cell_element_nodes = GridapP4est._edges_to_cell_element_nodes
# const _vertices_to_cell_element_nodes = GridapP4est._vertices_to_cell_element_nodes
# const p4est_lnodes_decode = GridapP4est.p4est_lnodes_decode
# const pXest_lnodes_decode = GridapP4est.pXest_lnodes_decode
# const process_current_face! = GridapP4est.process_current_face!
# const hanging_lvertex_within_edge = GridapP4est.hanging_lvertex_within_edge
# const regular_lvertex_within_face = GridapP4est.regular_lvertex_within_face
# const subface_to_hanging_edges_within_subface = GridapP4est.subface_to_hanging_edges_within_subface
# const subface_to_hanging_edges_within_face = GridapP4est.subface_to_hanging_edges_within_face
# const regular_lvertex_within_edge = GridapP4est.regular_lvertex_within_edge
# const NonConformingGlue = GridapP4est.NonConformingGlue
# const hanging_edge_from_face_code = -3
# const hanging_vertex_code         = -2



# function GridapP4est.generate_cell_faces_and_non_conforming_glue(pXest_type::PXestType, 
#                                                      ptr_pXest_lnodes, 
#                                                      cell_prange) 
#   # @info "generate_cell_faces_and_non_conforming_glue"
#   Dc = num_cell_dims(pXest_type)
#   n_cell_vertices = num_cell_vertices(Val{Dc})
#   n_cell_edges    = num_cell_edges(Val{Dc})
#   n_cell_faces    = num_cell_faces(Val{Dc})
  
#   lnodes = ptr_pXest_lnodes[]
#   element_nodes = unsafe_wrap(Array, lnodes.element_nodes, lnodes.vnodes * lnodes.num_local_elements)
#   face_code = unsafe_wrap(Array, lnodes.face_code, lnodes.num_local_elements)
#   hanging_face = Vector{Cint}(undef, n_cell_faces)
#   face_code_with_ghosts = map(partition(cell_prange)) do indices
#       @assert length(face_code)==own_length(indices)
#       @assert own_length(indices)==lnodes.num_local_elements
#       face_code_with_ghosts=similar(face_code, local_length(indices))
#       face_code_with_ghosts[1:own_length(indices)] .= face_code
#       face_code_with_ghosts
#   end

# #   cache_face_code=fetch_vector_ghost_values_cache(face_code_with_ghosts, partition(cell_prange))
# #   fetch_vector_ghost_values!(face_code_with_ghosts, cache_face_code) |> wait

#   element_nodes_with_ghosts = map(partition(cell_prange)) do indices
#     nonlocal_nodes = unsafe_wrap(Array, lnodes.nonlocal_nodes, lnodes.num_local_nodes-lnodes.owned_count)
#     element_nodes_with_ghosts_data=similar(element_nodes, local_length(indices)*lnodes.vnodes)
#     for (i,node) in enumerate(element_nodes)
#       if (node<lnodes.owned_count)
#         element_nodes_with_ghosts_data[i] = lnodes.global_offset+node
#       else
#         element_nodes_with_ghosts_data[i] = nonlocal_nodes[node-lnodes.owned_count+1]
#       end 
#     end     
#     element_nodes_with_ghosts_ptrs = [i for i=1:lnodes.vnodes:length(element_nodes_with_ghosts_data)+1]
#     JaggedArray(element_nodes_with_ghosts_data,element_nodes_with_ghosts_ptrs)
#   end
# #   cache_element_nodes_with_ghosts=fetch_vector_ghost_values_cache(element_nodes_with_ghosts, partition(cell_prange))
# #   fetch_vector_ghost_values!(element_nodes_with_ghosts, cache_element_nodes_with_ghosts) |> wait

#   map(element_nodes_with_ghosts,face_code_with_ghosts,partition(cell_prange)) do element_nodes_with_ghosts, face_code_with_ghosts, indices
#     @debug "ENDES_WO_GHOSTS[$(part_id(indices))]: $(element_nodes)"
#     @debug "ENDES_WITH_GHOSTS[$(part_id(indices))]: $(element_nodes_with_ghosts.data)"
#     @debug "FCODS_WO_GHOSTS[$(part_id(indices))]: $(face_code)"
#     @debug "FCODS_WITH_GHOSTS[$(part_id(indices))]: $(face_code_with_ghosts)"
#   end

#   if (Dc==2)
#     hanging_lvertex_within_face=hanging_lvertex_within_face_2d
#     pXest_face_corners = p4est_face_corners 
#   else
#     hanging_lvertex_within_face=hanging_lvertex_within_face_3d
#     pXest_face_corners = p8est_face_corners
#   end 

#   if (Dc==3)
#     hanging_edge = Vector{Cint}(undef, n_cell_edges)
#   end

#   num_regular_faces,
#   num_hanging_faces,
#   gridap_cell_faces,
#   hanging_faces_glue = 
#       map(partition(cell_prange),
#           element_nodes_with_ghosts,
#           face_code_with_ghosts) do indices, 
#                                     element_nodes_with_ghosts, 
#                                     face_code_with_ghosts
#     @assert local_length(indices)==length(face_code_with_ghosts) 
 
#     num_local_elements = local_length(indices)
#     num_regular_faces = Vector{Int}(undef, Dc)
#     num_hanging_faces = Vector{Int}(undef, Dc)

#     # Count regular vertices
#     num_regular_faces[1] = 0
#     regular_vertices_p4est_to_gridap = Dict{Int,Int}()

#     num_regular_faces[Dc] = 0
#     regular_faces_p4est_to_gridap = Dict{Int,Int}()

#     if (Dc==2)
#       p4est_gface_to_gcell_p4est_lface = 
#          _build_map_from_faces_to_cell_lface(lnodes.vnodes, element_nodes_with_ghosts.data, face_code_with_ghosts)
#     else 
#       p4est_gface_to_gcell_p4est_lface, 
#          p4est_gedge_to_gcell_p4est_ledge = 
#            _build_map_from_faces_edges_to_cell_lface_ledge(pXest_type,
#                                                            lnodes.vnodes, 
#                                                            element_nodes_with_ghosts.data, 
#                                                            face_code_with_ghosts)          
#     end 

#     PXEST_2_GRIDAP_VERTEX = pXest_2_gridap_vertex(Val{Dc})
#     PXEST_2_GRIDAP_FACE   = pXest_2_gridap_facet(Val{Dc})
#     PXEST_2_GRIDAP_EDGE   = p8est_2_gridap_edge()

#     n = local_length(indices)
#     gridap_cell_vertices_ptrs = Vector{Int32}(undef,n+1)
#     gridap_cell_faces_ptrs = Vector{Int32}(undef,n+1)
#     gridap_cell_vertices_ptrs[1]=1
#     gridap_cell_faces_ptrs[1]=1

#     hanging_vertices_pairs_to_owner_face = Dict{Tuple{Int,Int},Int}()
#     hanging_faces_pairs_to_owner_face = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
#     for i=1:n
#       gridap_cell_vertices_ptrs[i+1]=gridap_cell_vertices_ptrs[i]+n_cell_vertices
#       gridap_cell_faces_ptrs[i+1]=gridap_cell_faces_ptrs[i]+n_cell_faces
#     end

#     gridap_cell_vertices_data = Vector{Int}(undef, num_local_elements * n_cell_vertices)
#     gridap_cell_vertices_data .= -1

#     gridap_cell_faces_data = Vector{Int}(undef, num_local_elements * n_cell_faces)
#     gridap_cell_faces_data .= -1

#     if (Dc==3)
#       num_regular_faces[2] = 0
#       gridap_cell_edges_ptrs = Vector{Int32}(undef,n+1)
#       gridap_cell_edges_ptrs[1]=1
#       for i=1:n
#         gridap_cell_edges_ptrs[i+1]=gridap_cell_edges_ptrs[i]+n_cell_edges
#       end
#       gridap_cell_edges_data = Vector{Int}(undef, num_local_elements * n_cell_edges)
#       gridap_cell_edges_data .= -1
#       hanging_edges_cell_ledge_to_owner_face_half = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
#       owner_edge_subedge_to_cell_ledge = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
#       hanging_vertices_pairs_to_owner_edge = Dict{Tuple{Int,Int},Int}()
#       regular_edges_p4est_to_gridap = Dict{Int,Int}()
#     end
    
#     faces_to_cell_element_nodes    = _faces_to_cell_element_nodes(pXest_type)
#     if (Dc==3)
#       edges_to_cell_element_nodes    = _edges_to_cell_element_nodes(pXest_type)
#     end
    
#     vertices_to_cell_element_nodes = _vertices_to_cell_element_nodes(pXest_type)

#     for cell = 1:num_local_elements
#       start_gridap_vertices = (cell - 1) * n_cell_vertices
#       start_gridap_faces    = (cell - 1) * n_cell_faces

#       s = (cell-1)*lnodes.vnodes + 1
#       e = cell*lnodes.vnodes
#       p4est_cell_nodes = view(element_nodes_with_ghosts.data, s:e)

#       p4est_cell_faces = view(p4est_cell_nodes,faces_to_cell_element_nodes)
#       p4est_cell_vertices = view(p4est_cell_nodes,vertices_to_cell_element_nodes)
    
#       gridap_cell_vertices = view(gridap_cell_vertices_data,
#         start_gridap_vertices+1:start_gridap_vertices+n_cell_vertices)
#       gridap_cell_faces = view(gridap_cell_faces_data,
#         start_gridap_faces+1:start_gridap_faces+n_cell_faces)

#       if (Dc==2)  
#         has_hanging = p4est_lnodes_decode(face_code_with_ghosts[cell], hanging_face)
#       else
#         has_hanging = pXest_lnodes_decode(pXest_type, face_code_with_ghosts[cell], hanging_face, hanging_edge)
#         start_gridap_edges = (cell-1)*n_cell_edges
#         gridap_cell_edges  = view(gridap_cell_edges_data, start_gridap_edges+1:start_gridap_edges+n_cell_edges)
#         p4est_cell_edges   = view(p4est_cell_nodes,edges_to_cell_element_nodes)
#       end
#       if has_hanging == 0
#         # All vertices/edges/faces of the current cell are regular 
#         # Process vertices
#         for (p4est_lvertex, p4est_gvertex) in enumerate(p4est_cell_vertices)
#           num_regular_faces[1] =
#             process_current_face!(gridap_cell_vertices,
#               regular_vertices_p4est_to_gridap,
#               num_regular_faces[1],
#               p4est_cell_vertices,
#               p4est_lvertex,
#               p4est_gvertex,
#               PXEST_2_GRIDAP_VERTEX)
#         end
        
#         if (Dc==3)
#           for (p4est_ledge, p4est_gedge) in enumerate(p4est_cell_edges)
#             num_regular_faces[2] =
#               process_current_face!(gridap_cell_edges,
#                 regular_edges_p4est_to_gridap,
#                 num_regular_faces[2],
#                 p4est_cell_edges,
#                 p4est_ledge,
#                 p4est_gedge,
#                 PXEST_2_GRIDAP_EDGE)
#           end
#         end 

#         # Process faces
#         for (p4est_lface, p4est_gface) in enumerate(p4est_cell_faces)
#           num_regular_faces[Dc] =
#             process_current_face!(gridap_cell_faces,
#               regular_faces_p4est_to_gridap,
#               num_regular_faces[Dc],
#               p4est_cell_faces,
#               p4est_lface,
#               p4est_gface,
#               PXEST_2_GRIDAP_FACE)
#         end
#       else
#         if (isa(pXest_type,P4P8estType))
#           # "Touch" hanging vertices before processing current cell
#           # This is required as we dont have any means to detect 
#           # a hanging vertex from a non-hanging face
#           for (p4est_lface, half) in enumerate(hanging_face)
#             if (half != -1)
#               hanging_vertex_lvertex_within_face = hanging_lvertex_within_face(half)
#               p4est_lvertex = pXest_face_corners[p4est_lface,
#                                                 hanging_vertex_lvertex_within_face+1]
#               gridap_cell_vertices[PXEST_2_GRIDAP_VERTEX[p4est_lvertex+1]] = hanging_vertex_code
#             end 
#           end
#         end

#         if (Dc==3)
#           for (p4est_ledge, half) in enumerate(hanging_edge)
#             if (half != -1 && half !=4)
#               hanging_vertex_lvertex_within_edge = hanging_lvertex_within_edge(half)
#               p4est_lvertex = p8est_edge_corners[p4est_ledge,
#                                                  hanging_vertex_lvertex_within_edge+1]
#               gridap_cell_vertices[PXEST_2_GRIDAP_VERTEX[p4est_lvertex+1]] = hanging_vertex_code

#               @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] cell=$(cell) hanging p4est_ledge=$(p4est_ledge) half=$(half) hanging p4est_lvertex=$(PXEST_2_GRIDAP_VERTEX[p4est_lvertex+1])"
#             end 
#           end
#         end 

#         # Current cell has at least one hanging face 
#         for (p4est_lface, half) in enumerate(hanging_face)
#           # Current face is NOT hanging
#           if (half == -1)
#             # Process vertices on the boundary of p4est_lface
#             for p4est_lvertex in pXest_face_corners[p4est_lface, :]
#               p4est_gvertex = p4est_cell_vertices[p4est_lvertex+1]
#               if (gridap_cell_vertices[p4est_lvertex+1] != hanging_vertex_code)
#                 num_regular_faces[1] =
#                   process_current_face!(gridap_cell_vertices,
#                     regular_vertices_p4est_to_gridap,
#                     num_regular_faces[1],
#                     p4est_cell_vertices,
#                     p4est_lvertex + 1,
#                     p4est_gvertex,
#                     PXEST_2_GRIDAP_VERTEX)
#               end
#             end
#             # Process non-hanging face
#             p4est_gface = p4est_cell_faces[p4est_lface]
#             num_regular_faces[Dc] =
#               process_current_face!(gridap_cell_faces,
#                 regular_faces_p4est_to_gridap,
#                 num_regular_faces[Dc],
#                 p4est_cell_faces,
#                 p4est_lface,
#                 p4est_gface,
#                 PXEST_2_GRIDAP_FACE)
#           else # Current face is hanging
#             owner_face = p4est_cell_faces[p4est_lface]

#             if half in 0:3
#               # Identify regular vertex and hanging vertex 
#               # Repeat code above for regular vertex 
#               # Special treatment for hanging vertex 
#               regular_vertex_lvertex_within_face = regular_lvertex_within_face(half)
#               hanging_vertex_lvertex_within_face = hanging_lvertex_within_face(half)

#               # Process regular vertex
#               p4est_regular_lvertex = pXest_face_corners[p4est_lface, regular_vertex_lvertex_within_face+1]
#               p4est_gvertex = p4est_cell_vertices[p4est_regular_lvertex+1]
#               num_regular_faces[1] =
#                 process_current_face!(gridap_cell_vertices,
#                   regular_vertices_p4est_to_gridap,
#                   num_regular_faces[1],
#                   p4est_cell_vertices,
#                   p4est_regular_lvertex + 1,
#                   p4est_gvertex,
#                   PXEST_2_GRIDAP_VERTEX)
              
#               # Process hanging vertex
#               p4est_hanging_lvertex = pXest_face_corners[p4est_lface, hanging_vertex_lvertex_within_face+1]
#               hanging_vertices_pairs_to_owner_face[(cell, PXEST_2_GRIDAP_VERTEX[p4est_hanging_lvertex+1])] = owner_face

#               # Process hanging face
#               hanging_faces_pairs_to_owner_face[(cell, PXEST_2_GRIDAP_FACE[p4est_lface])] = (owner_face,half+1)
#             else
#               # Anisotropic refinement
#               @assert half in 4:7
#               for regular_vertex_lvertex_within_face in p6est_half_to_regular_vertices[mod(half,4)+1,:]
#                 # Process regular vertex
#                 p4est_regular_lvertex = pXest_face_corners[p4est_lface, regular_vertex_lvertex_within_face+1]
#                 p4est_gvertex = p4est_cell_vertices[p4est_regular_lvertex+1]
#                 num_regular_faces[1] =
#                   process_current_face!(gridap_cell_vertices,
#                     regular_vertices_p4est_to_gridap,
#                     num_regular_faces[1],
#                     p4est_cell_vertices,
#                     p4est_regular_lvertex + 1,
#                     p4est_gvertex,
#                     PXEST_2_GRIDAP_VERTEX)
#               end
#               subface = mod(half,2)+1

#               # Process hanging face
#               hanging_faces_pairs_to_owner_face[(cell, PXEST_2_GRIDAP_FACE[p4est_lface])] = (owner_face,subface)
#             end  
           
#             if (Dc==3)
#               _subface_to_hanging_edges_within_subface = subface_to_hanging_edges_within_subface(pXest_type)
#               _subface_to_hanging_edges_within_face    = subface_to_hanging_edges_within_face(pXest_type)

#               for (i,ledge_within_face) in enumerate(_subface_to_hanging_edges_within_subface[mod(half,4)+1,:])
#                 p4est_ledge=p8est_face_edges[p4est_lface,ledge_within_face+1]

#                 @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  cell=$(cell) p4est_lface=$(p4est_lface) half=$(half) index=$(mod(half,4)+1) p4est_ledge=$(p4est_ledge) ledge_within_face=$(ledge_within_face) "

#                 gridap_ledge = PXEST_2_GRIDAP_EDGE[p4est_ledge+1]
#                 # Identify the two edges which are hanging within the face
#                 hanging_edges_cell_ledge_to_owner_face_half[(cell, gridap_ledge)] =
#                     (owner_face,-_subface_to_hanging_edges_within_face[mod(half,4)+1,i])
#                 gridap_cell_edges[gridap_ledge] = hanging_edge_from_face_code
#               end 
#             end 

#           end
#         end


#         if (Dc==3)
#           for (p4est_ledge, half) in enumerate(hanging_edge)
#             # Current edge is NOT hanging
#             if (half == -1)
#               # Process vertices on the boundary of p4est_ledge
#               for p4est_lvertex in p8est_edge_corners[p4est_ledge, :]
#                 p4est_gvertex = p4est_cell_vertices[p4est_lvertex+1]
#                 if (gridap_cell_vertices[p4est_lvertex+1] != hanging_vertex_code)
#                   num_regular_faces[1] =
#                     process_current_face!(gridap_cell_vertices,
#                       regular_vertices_p4est_to_gridap,
#                       num_regular_faces[1],
#                       p4est_cell_vertices,
#                       p4est_lvertex + 1,
#                       p4est_gvertex,
#                       PXEST_2_GRIDAP_VERTEX)
#                 end
#               end
#               # Process non-hanging edge
#               p4est_gedge = p4est_cell_edges[p4est_ledge]
#               num_regular_faces[2] =
#                 process_current_face!(gridap_cell_edges,
#                   regular_edges_p4est_to_gridap,
#                   num_regular_faces[2],
#                   p4est_cell_edges,
#                   p4est_ledge,
#                   p4est_gedge,
#                   PXEST_2_GRIDAP_EDGE)
#             else # Current edge is hanging
#               @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  cell=$(cell) hanging_edge=$(hanging_edge) p4est_ledge=$(p4est_ledge) gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge]]: $(gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge]]) half=$(half)"

#               if ( gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge]] != hanging_edge_from_face_code )
#                 # The present hanging edge cannot be within a coarser face
#                 @assert half != 4 

#                 # # Identify regular vertex and hanging vertex 
#                 # # Repeat code above for regular vertex 
#                 # # Special treatment for hanging vertex 
#                 regular_vertex_lvertex_within_edge = regular_lvertex_within_edge(half)
#                 hanging_vertex_lvertex_within_edge = hanging_lvertex_within_edge(half)

#                 # # Process regular vertex
#                 p4est_regular_lvertex = p8est_edge_corners[p4est_ledge, regular_vertex_lvertex_within_edge+1]
#                 p4est_gvertex = p4est_cell_vertices[p4est_regular_lvertex+1]
                
#                 num_regular_faces[1] =
#                   process_current_face!(gridap_cell_vertices,
#                     regular_vertices_p4est_to_gridap,
#                     num_regular_faces[1],
#                     p4est_cell_vertices,
#                     p4est_regular_lvertex + 1,
#                     p4est_gvertex,
#                     PXEST_2_GRIDAP_VERTEX)
                
#                 # Process hanging vertex
#                 p4est_hanging_lvertex = p8est_edge_corners[p4est_ledge, hanging_vertex_lvertex_within_edge+1]
#                 p4est_owner_edge = p4est_cell_edges[p4est_ledge]
#                 hanging_vertices_pairs_to_owner_edge[(cell, 
#                                                       PXEST_2_GRIDAP_VERTEX[p4est_hanging_lvertex+1])] = p4est_owner_edge

#                 # Process hanging edge
#                 subedge = regular_vertex_lvertex_within_edge+1
#                 owner_edge_subedge_pair=(p4est_owner_edge,subedge)

#                 @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  cell=$(cell) hanging_edge=$(hanging_edge) p4est_ledge=$(p4est_ledge) gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge]]: $(gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge]]) half=$(half) owner_edge_subedge_pair=$(owner_edge_subedge_pair)"

#                 gridap_ledge=PXEST_2_GRIDAP_EDGE[p4est_ledge]
#                 hanging_edges_cell_ledge_to_owner_face_half[(cell, gridap_ledge)] = owner_edge_subedge_pair
#                 if (!haskey(owner_edge_subedge_to_cell_ledge,owner_edge_subedge_pair))
#                   owner_edge_subedge_to_cell_ledge[owner_edge_subedge_pair] = (cell,gridap_ledge)
#                 end
#               end 
#             end 
#         end
#       end
#     end
#   end 

#     function is_ghost(cell)
#       cell>own_length(indices)
#     end

#     # Go over all touched hanging faces and start 
#     # assigning IDs from the last num_regular_faces ID
#     # For each hanging face, keep track of (owner_cell,lface)
#     # Go over all hanging faces 
#     # Detect if the owner face is in a ghost cell. 
#     # If not in a ghost cell or touched 
#     # Else 
#     #   The current face becomes a regular face 
#     # end 
#     hanging_faces_owner_cell_and_lface =
#       Vector{Tuple{Int,Int,Int}}(undef, length(keys(hanging_faces_pairs_to_owner_face)))
#     # @info hanging_faces_pairs_to_owner_face, regular_faces_p4est_to_gridap
#     num_hanging_faces[Dc] = 0
#     for key in keys(hanging_faces_pairs_to_owner_face)
#       (cell, lface) = key
#       (owner_p4est_gface, half) = hanging_faces_pairs_to_owner_face[key]
#       num_hanging_faces[Dc] += 1
#       start_gridap_faces = (cell - 1) * n_cell_faces
#       gridap_cell_faces_data[start_gridap_faces+lface] = num_regular_faces[Dc] + num_hanging_faces[Dc]
#       if (!(is_ghost(cell)) || haskey(regular_faces_p4est_to_gridap,owner_p4est_gface)) && (haskey(p4est_gface_to_gcell_p4est_lface,owner_p4est_gface))  
#         # @assert haskey(regular_faces_p4est_to_gridap,owner_p4est_gface)
#         (owner_cell, p4est_lface) = p4est_gface_to_gcell_p4est_lface[owner_p4est_gface]
#         hanging_faces_owner_cell_and_lface[num_hanging_faces[Dc]] =
#           (owner_cell, n_cell_vertices+n_cell_edges+PXEST_2_GRIDAP_FACE[p4est_lface], half)
#       else
#         # Glue info cannot be computed for this hanging face
#         hanging_faces_owner_cell_and_lface[num_hanging_faces[Dc]] = (-1,-1,-1)
#       end
#     end

#     @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  gridap_cell_faces_data: $(gridap_cell_faces_data)"


#     # Go over all touched hanging vertices and start 
#     # assigning IDs from the last num_regular_vertices ID
#     # For each hanging vertex, keep track of (owner_cell,lface)
#     num_hanging_faces[1] = 0
#     hanging_vertices_owner_cell_and_lface = Tuple{Int,Int,Int}[]
#     half=1
#     owner_p4est_gface_to_hanging_vertex = Dict{Int,Int}()
#     for key in keys(hanging_vertices_pairs_to_owner_face)
#       (cell, lvertex) = key
#       owner_p4est_gface = hanging_vertices_pairs_to_owner_face[key]
#       if !(haskey(owner_p4est_gface_to_hanging_vertex, owner_p4est_gface))
#         num_hanging_faces[1] += 1
#         owner_p4est_gface_to_hanging_vertex[owner_p4est_gface] = num_hanging_faces[1]
#         if (!is_ghost(cell) || (haskey(regular_faces_p4est_to_gridap,owner_p4est_gface))) && (haskey(p4est_gface_to_gcell_p4est_lface,owner_p4est_gface))
#           (owner_cell, p4est_lface) = p4est_gface_to_gcell_p4est_lface[owner_p4est_gface]
#           push!(hanging_vertices_owner_cell_and_lface,
#             (owner_cell, n_cell_vertices+n_cell_edges+PXEST_2_GRIDAP_FACE[p4est_lface],half))
#         else
#           push!(hanging_vertices_owner_cell_and_lface,(-1, -1,-1))
#         end
#       end 
#       start_gridap_vertices = (cell - 1) * n_cell_vertices
#       gridap_cell_vertices_data[start_gridap_vertices+lvertex] = num_regular_faces[1] +
#                                                           owner_p4est_gface_to_hanging_vertex[owner_p4est_gface] 
#     end

#     @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  gridap_cell_vertices_data: $(gridap_cell_vertices_data)"
      
#     if (Dc==3)
#       half=1
#       owner_p4est_gedge_to_hanging_vertex = Dict{Int,Int}()
#       for key in keys(hanging_vertices_pairs_to_owner_edge)
#         (cell, lvertex) = key
#         owner_p4est_gedge = hanging_vertices_pairs_to_owner_edge[key]
#         if !(haskey(owner_p4est_gedge_to_hanging_vertex, owner_p4est_gedge))
#           num_hanging_faces[1] += 1
#           owner_p4est_gedge_to_hanging_vertex[owner_p4est_gedge] = num_hanging_faces[1]
#           if (!is_ghost(cell) || (haskey(regular_edges_p4est_to_gridap,owner_p4est_gedge)))
#               (owner_cell, p4est_ledge) = p4est_gedge_to_gcell_p4est_ledge[owner_p4est_gedge]
#               push!(hanging_vertices_owner_cell_and_lface,
#                       (owner_cell, n_cell_vertices+PXEST_2_GRIDAP_EDGE[p4est_ledge],half))
#           else 
#             push!(hanging_vertices_owner_cell_and_lface,(-1, -1,-1))
#           end
#         end
#         start_gridap_vertices = (cell - 1) * n_cell_vertices
#         gridap_cell_vertices_data[start_gridap_vertices+lvertex] = num_regular_faces[1] +
#                                                 owner_p4est_gedge_to_hanging_vertex[owner_p4est_gedge]
#       end 

#       # Go over all touched hanging edges and start 
#       # assigning IDs from the last num_regular_edge ID
#       # For each hanging edge, keep track of (owner_cell,lface/ledge)
#       hanging_edges_owner_cell_and_lface = Tuple{Int,Int,Int}[]
#       owner_p4est_gface_half_to_hanging_edge = Dict{Tuple{Int,Int},Int}()
#       owner_p4est_gedge_subedge_to_hanging_edge = Dict{Tuple{Int,Int},Int}()
#       num_hanging_faces[2] = 0
#       ledge_to_cvertices = Gridap.ReferenceFEs.get_faces(HEX, 1, 0)
#       # The following loop needs gridap cell vertices to be already completed
#       for key in keys(hanging_edges_cell_ledge_to_owner_face_half)
#           (cell, ledge) = key
#           (owner_p4est_gface_or_gedge, half) = hanging_edges_cell_ledge_to_owner_face_half[key]
#           if (half<0) # hanging edge is within a coarser face 
#             owner_p4est_gface = owner_p4est_gface_or_gedge
#             if !(haskey(owner_p4est_gface_half_to_hanging_edge, (owner_p4est_gface,half)))
#               num_hanging_faces[2] += 1
#               owner_p4est_gface_half_to_hanging_edge[(owner_p4est_gface,half)] = num_hanging_faces[2]
#               if (!is_ghost(cell) || (haskey(regular_faces_p4est_to_gridap,owner_p4est_gface)))
#                 (owner_cell, p4est_lface) = p4est_gface_to_gcell_p4est_lface[owner_p4est_gface]
#                 push!(hanging_edges_owner_cell_and_lface,
#                       (owner_cell, n_cell_vertices+n_cell_edges+PXEST_2_GRIDAP_FACE[p4est_lface],half))
#               else
#                 push!(hanging_edges_owner_cell_and_lface,(-1, -1, -1))
#               end 
#             end
#             start_gridap_edges = (cell - 1) * n_cell_edges
#             gridap_cell_edges_data[start_gridap_edges+ledge] = num_regular_faces[2] + 
#                           owner_p4est_gface_half_to_hanging_edge[(owner_p4est_gface,half)]
#           else # hanging edge is within a coarser edge
#             @assert half==1 || half==2
#             owner_p4est_gedge = owner_p4est_gface_or_gedge
#             owner_gedge_pair = (owner_p4est_gedge,half)
#             if (haskey(owner_edge_subedge_to_cell_ledge,owner_gedge_pair))
#               (owner_cell, owner_cell_ledge) = owner_edge_subedge_to_cell_ledge[owner_gedge_pair]
#               if (owner_cell==cell)
#                 @assert owner_cell_ledge == ledge
#                 num_hanging_faces[2] += 1
#                 owner_p4est_gedge_subedge_to_hanging_edge[(owner_p4est_gedge,half)] = num_hanging_faces[2]
#                 if (!is_ghost(cell) || (haskey(regular_edges_p4est_to_gridap,owner_p4est_gedge)))
#                   (owner_cell, p4est_ledge) = p4est_gedge_to_gcell_p4est_ledge[owner_p4est_gedge]
#                   push!(hanging_edges_owner_cell_and_lface,
#                       (owner_cell, n_cell_vertices+PXEST_2_GRIDAP_EDGE[p4est_ledge],half))
#                 else
#                   push!(hanging_edges_owner_cell_and_lface,(-1, -1, -1))
#                 end
#                 start_gridap_edges = (cell - 1) * n_cell_edges
#                 gridap_cell_edges_data[start_gridap_edges+ledge] = num_regular_faces[2] + num_hanging_faces[2] 
#               end
#             end
#           end
#       end
#       @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] gridap_cell_edges_data: $(gridap_cell_edges_data)"
      
#       for key in keys(hanging_edges_cell_ledge_to_owner_face_half)
#         (cell, ledge) = key
#         (owner_p4est_gface_or_gedge, half) = hanging_edges_cell_ledge_to_owner_face_half[key]
#         @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] own_length=$(own_length(indices)) cell=$(cell) ledge=$(ledge) owner_p4est_gface_or_gedge=$(owner_p4est_gface_or_gedge) half=$(half)"
#         if (half>0) # hanging edge is within a coarser face 
#           @assert half==1 || half==2
#           owner_p4est_gedge = owner_p4est_gface_or_gedge
#           owner_gedge_pair = (owner_p4est_gedge,half)
#           if (haskey(owner_edge_subedge_to_cell_ledge,owner_gedge_pair))
#             (owner_cell, owner_cell_ledge) = owner_edge_subedge_to_cell_ledge[owner_gedge_pair]
#             if (owner_cell!=cell)
#               haskey_first_subedge  = haskey(owner_p4est_gedge_subedge_to_hanging_edge,(owner_p4est_gedge,1))
#               haskey_second_subedge = haskey(owner_p4est_gedge_subedge_to_hanging_edge,(owner_p4est_gedge,2))
#               if (!(is_ghost(cell)))
#                 @assert haskey_first_subedge && haskey_second_subedge
#               else
#                 @assert haskey_first_subedge || haskey_second_subedge
#               end 
#               if (haskey_first_subedge && haskey_second_subedge)
#                 # The following code is required as we may have edges 
#                 # with different orientations at the inter-octree boundaries
#                 match=true
#                 start_gridap_vertices_cell = (cell - 1) * n_cell_vertices
#                 start_gridap_vertices_cell_owner = (owner_cell - 1) * n_cell_vertices
#                 for lvertex_cell in ledge_to_cvertices[ledge]
#                   vertex_cell=gridap_cell_vertices_data[start_gridap_vertices_cell+lvertex_cell]
#                   found=false
#                   # Go over vertices of owner_cell_ledge in owner_cell 
#                   for lvertex_owner_cell in ledge_to_cvertices[owner_cell_ledge]
#                     vertex_owner_cell=gridap_cell_vertices_data[start_gridap_vertices_cell_owner+lvertex_owner_cell]
#                     if (vertex_owner_cell==vertex_cell)
#                       found=true
#                       break
#                     end
#                   end
#                   if (!found)
#                     match=false
#                     break
#                   end
#                 end
#                 if (match)
#                   owner_half=half
#                 else 
#                   owner_half=half==1 ? 2 : 1
#                 end 
#               elseif (haskey_first_subedge)
#                 owner_half=1
#               elseif (haskey_second_subedge)
#                 owner_half=2
#               end
#               @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] cell=$(cell) ledge=$(ledge) owner_p4est_gface_or_gedge=$(owner_p4est_gface_or_gedge) half=$(half) owner_half=$(owner_half)"
#               start_gridap_edges = (cell - 1) * n_cell_edges
#               gridap_cell_edges_data[start_gridap_edges+ledge] = 
#                  num_regular_faces[2] + owner_p4est_gedge_subedge_to_hanging_edge[(owner_p4est_gedge,owner_half)] 
#             end
#           end
#        end
#       end 
#       @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] gridap_cell_edges_data: $(gridap_cell_edges_data)"
#     end

#     gridap_cell_faces = Vector{JaggedArray}(undef,Dc)
#     gridap_cell_faces[1] = JaggedArray(gridap_cell_vertices_data,gridap_cell_vertices_ptrs)
#     if (Dc==3)
#       gridap_cell_faces[2] = JaggedArray(gridap_cell_edges_data,gridap_cell_edges_ptrs)
#     end   
#     gridap_cell_faces[Dc] = JaggedArray(gridap_cell_faces_data,gridap_cell_faces_ptrs)

#     hanging_faces_glue      = Vector{Vector{Tuple}}(undef,Dc)
#     hanging_faces_glue[1]   = hanging_vertices_owner_cell_and_lface
#     if (Dc==3)
#       hanging_faces_glue[2] = hanging_edges_owner_cell_and_lface
#     end 
#     hanging_faces_glue[Dc]  = hanging_faces_owner_cell_and_lface


#     return num_regular_faces, 
#            num_hanging_faces,
#            gridap_cell_faces,
#            hanging_faces_glue

#   end |> tuple_of_arrays

  
#   gridap_cell_faces_out  = Vector{MPIArray}(undef,Dc)
#   for i=1:Dc
#     gridap_cell_faces_out[i] = map(gridap_cell_faces) do  gridap_cell_faces
#       gridap_cell_faces[i]
#     end
#   end
#   non_conforming_glue=map(num_regular_faces,num_hanging_faces,hanging_faces_glue) do nrf, nhf, hfg
#     NonConformingGlue(nrf, nhf, hfg)
#   end 
#   gridap_cell_faces_out,non_conforming_glue
# end


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
  filter!(x -> x != "hanging", labels.tag_to_name)
  facet_to_tag = get_face_tag_index(labels,labels.tag_to_name,dim)
  # @info facet_to_tag, labels.tag_to_name,  length(findall(x -> x>0, facet_to_tag))
  findall(x -> x>0, facet_to_tag)
end




# function mod_mesh_adapt!(omesh::St_mesh, mesh::St_mesh, inputs::Dict, nparts, distribute)

#     # determine backend
#     backend = CPU()
#     comm = MPI.COMM_WORLD
#     rank = MPI.Comm_rank(comm)
    
#     #
#     # Read GMSH grid from file
#     #      
#     parts  = distribute(LinearIndices((nparts,)))
#     mesh.parts = distribute(LinearIndices((nparts,)))
#     mesh.nparts = nparts
#     ladpative = 1
#     if ladpative == 0
#         gmodel = GmshDiscreteModel(inputs[:gmsh_filename], renumber=true)
#         partitioned_model = GmshDiscreteModel(parts, inputs[:gmsh_filename], renumber=true)
#         # smodel = GmshDiscreteModel(inputs[:gmsh_filename], renumber=true)
#         # g = GridapDistributed.compute_cell_graph(smodel)
#         # cell_to_part = Metis.partition(g,4)
#         # @info cell_to_part

#         # coarse_discrete_model = GmshDiscreteModel(inputs[:gmsh_filename], renumber=true)
#         # num_uniform_refinements = 1
#         # partitioned_model = UniformlyRefinedForestOfOctreesDiscreteModel(parts,
#         #                                                    coarse_discrete_model,
#         #                                                    num_uniform_refinements)

#         # model = local_views(partitioned_model).item_ref[]
#     elseif ladpative == 1
#         gmodel = GmshDiscreteModel(inputs[:gmsh_filename], renumber=true)
#         # gmodel_type = typeof(gmodel)
#         # @info gmodel_type.parameters
#         partitioned_model_coarse = OctreeDistributedDiscreteModel(parts,gmodel)
#         ref_coarse_flags=map(parts,partition(get_cell_gids(partitioned_model_coarse.dmodel))) do rank,indices
#             flags=zeros(Cint,length(indices))
#             flags.=nothing_flag
#             # @info flags
#             if rank == 2
#                 # flags[10:5:end] .= refine_flag
#                 flags .= refine_flag
#                 # flags[1:5] .= refine_flag
#             end
#             flags
#         end
#         partitioned_model,glue_adapt=Gridap.Adaptivity.adapt(partitioned_model_coarse,ref_coarse_flags);
#         # ref_coarse_flags=map(parts,partition(get_cell_gids(partitioned_model.dmodel))) do rank,indices
#         #     flags=zeros(Cint,length(indices))
#         #     flags.=nothing_flag
#         #     # @info flags
#         #     if rank == 2
#         #         # flags[10:5:end] .= refine_flag
#         #         flags[1:end] .= refine_flag
#         #     end
#         #     flags
#         # end
#         # partitioned_model,glue=Gridap.Adaptivity.adapt(partitioned_model,ref_coarse_flags);
#         # partitioned_model, glue_redistribute = redistribute(partitioned_model)
#         map(parts, glue_adapt, partition(get_cell_gids(partitioned_model.dmodel)), partition(get_cell_gids(partitioned_model_coarse.dmodel))) do part, glue, indices_f, indices_c
#             # @info typeof(glue)
#             @info part, glue.n2o_cell_to_child_id
#             @info part, glue.n2o_faces_map[3], local_to_global(indices_f), local_to_global(indices_c), glue.is_refined, glue.o2n_faces_map
#             # @info rank, glue.hanging_faces_glue[1], "  abc  ", glue.hanging_faces_glue[2]
#         end
#         n2o_el_map = local_views(glue_adapt).item_ref[].n2o_faces_map[3]
#         # domain = Triangulation(partitioned_model.dmodel)
#         # cell_gids = get_cell_gids(partitioned_model)
#         # remove_ghost_cells(domain, cell_gids)
        
#         # num_cell = map(partition(cell_gids)) do gids
#         #     # @info global_to_own(gids)
#         #     num_cells = own_length(gids)
#         # end

#         # dmodel = map(parts, partitioned_model.dmodel.models, partition(cell_gids)) do part, model, gids
#         #     # cell_id = own_to_global(gids)
#         #     # @info part, cell_id
#         #     cell_id2 = global_to_own(gids)
#         #     dmodel = DiscreteModelPortion(model, filter(x -> x > 0, cell_id2))
#         # end
#         # part_model = GridapDistributed.GenericDistributedDiscreteModel(dmodel, cell_gids)
#     end


#     # MyGeometry.get_boundary_cells(model,mesh.nsd)

#     cell_gids = local_views(partition(get_cell_gids(partitioned_model))).item_ref[]
#     dmodel = local_views(partitioned_model.dmodel.models).item_ref[]
#     # @info rank, own_to_local(cell_gids), local_to_own(cell_gids), local_to_global(cell_gids)
#     model  = DiscreteModelPortion(dmodel, own_to_local(cell_gids))
#     topology      = get_grid_topology(model)
#     gtopology      = get_grid_topology(model)
#     mesh.nsd      = num_cell_dims(model)
    
#     POIN_flg = 0
#     EDGE_flg = 1
#     FACE_flg = 2
#     ELEM_flg = 3
    
#     if mesh.nsd == 3
#         mesh.SD = NSD_3D()

#         mesh.NNODES_EL  = 8
#         mesh.NEDGES_EL  = 12
#         mesh.NFACES_EL  = 6
#         mesh.EDGE_NODES = 2
#         mesh.FACE_NODES = 4
#     elseif mesh.nsd == 2

#         mesh.SD = NSD_2D()
#         ELEM_flg = FACE_flg
        
#         mesh.NNODES_EL  = 4
#         mesh.NEDGES_EL  = 4
#         mesh.NFACES_EL  = 1
#         mesh.EDGE_NODES = 2
#         mesh.FACE_NODES = 4
#     elseif mesh.nsd == 1
#         mesh.SD = NSD_1D()

#         mesh.NNODES_EL  = 2
#         mesh.NEDGES_EL  = 1
#         mesh.NFACES_EL  = 0
#         mesh.EDGE_NODES = 2
#         mesh.FACE_NODES = 0
#     else
#         error( " WRONG NSD: This is not theoretical physics: we only handle 1, 2, or 3 dimensions!")
#     end
    


#     # p2pp = get_faces(topology, 0, 0)
#     # @info rank, p2pp.data, p2pp.ptrs
#     # d_to_num_dfaces = [num_vertices(model), num_edges(model), local_views(num_cell).item_ref[]]
#     d_to_num_dfaces = [num_vertices(model), num_edges(model), num_cells(model)]
#     # @info d_to_num_dfaces
#     # @info edge2pedge

#     # Write the partitioned model to a VTK file
#     vtk_directory = "./coarse/" 
#     writevtk(partitioned_model_coarse, vtk_directory)
#     vtk_directory = "./refine/" 
#     writevtk(partitioned_model, vtk_directory)



#     #dump(topology)
#     #
#     # Mesh elements, nodes, faces, edges
#     #


#     p2pp = Geometry.get_face_to_parent_face(model,0)
#     # edge2pedge = Geometry.get_face_to_parent_face(model,1)
#     # face2pface = Geometry.get_face_to_parent_face(model,2)
#     # @info point2ppoint, edge2pedge, face2pface

#     pgids = local_views(partition(get_face_gids(partitioned_model,POIN_flg))).item_ref[]
#     edgids = local_views(partition(get_face_gids(partitioned_model,EDGE_flg))).item_ref[]
#     fgids = local_views(partition(get_face_gids(partitioned_model,FACE_flg))).item_ref[]
#     elgids = local_views(partition(get_face_gids(partitioned_model,ELEM_flg))).item_ref[]
#     # @info rank, own_to_local(pgids), "a ", local_to_own(pgids),  "a ", local_to_global(pgids), "a ", p2pp
#     point2ppoint = local_to_global(pgids)[p2pp]
#     edge2pedge = local_to_global(edgids)
#     face2pface = local_to_global(fgids)
#     elm2pelm = local_to_global(elgids)
#     # @info rank, p2pp, point2ppoint


#     mesh.gnpoin_linear = num_faces(partitioned_model.dmodel,POIN_flg)    
#     mesh.gnpoin        = mesh.gnpoin_linear         #This will be updated for the high order grid
#     mesh.gnedges       = num_faces(partitioned_model,EDGE_flg)
#     mesh.gnfaces       = num_faces(partitioned_model,FACE_flg)   
#     mesh.gnelem        = num_faces(partitioned_model,ELEM_flg)


#     mesh.npoin_linear = num_faces(model,POIN_flg)    
#     mesh.npoin        = mesh.npoin_linear         #This will be updated for the high order grid
#     mesh.nedges       = num_faces(model,EDGE_flg)
#     mesh.nfaces       = num_faces(model,FACE_flg)   
#     mesh.nelem        = num_faces(model,ELEM_flg)


    
#     if (ladpative == 1)
#         mesh.nelem_bdy    = length(MyGeometry.get_boundary_cells(model,mesh.nsd))
#         mesh.nfaces_bdy   = length(MyGeometry.get_boundary_faces(model,mesh.nsd,FACE_flg))
#         mesh.nedges_bdy   = length(MyGeometry.get_boundary_faces(model,mesh.nsd,EDGE_flg))
#     else 
#         mesh.nelem_bdy    = count(get_isboundary_face(topology,mesh.nsd))
#         mesh.nfaces_bdy   = count(get_isboundary_face(topology,mesh.nsd-1))
#         mesh.nedges_bdy   = count(get_isboundary_face(topology,mesh.nsd-2))
#     end
#     mesh.nelem_int    = mesh.nelem - mesh.nelem_bdy
#     mesh.nfaces_int   = mesh.nfaces - mesh.nfaces_bdy
#     mesh.nedges_int   = mesh.nedges - mesh.nedges_bdy

#     #get_isboundary_face(topology,mesh.nsd-1)
#     println(" # GMSH LINEAR GRID PROPERTIES")
#     println(" # N. Global points         : ", mesh.gnpoin_linear)
#     println(" # N. Global elements       : ", mesh.gnelem)
#     println(" # N. Global edges          : ", mesh.gnedges)
#     println(" # N. Global faces          : ", mesh.gnfaces)    
#     println(" # N. points         : ", mesh.npoin_linear)
#     println(" # N. elements       : ", mesh.nelem)
#     println(" # N. edges          : ", mesh.nedges)
#     println(" # N. faces          : ", mesh.nfaces)    
#     println(" # N. internal elem  : ", mesh.nelem_int)
#     println(" # N. internal edges : ", mesh.nedges_int) 
#     println(" # N. internal faces : ", mesh.nfaces_int)    
#     println(" # N. boundary elem  : ", mesh.nelem_bdy)
#     println(" # N. boundary edges : ", mesh.nedges_bdy)
#     println(" # N. boundary faces : ", mesh.nfaces_bdy)
#     println(" # GMSH LINEAR GRID PROPERTIES ...................... END")
    
#     ngl                     = mesh.nop + 1
#     tot_linear_poin         = mesh.npoin_linear
    
#     tot_edges_internal_nodes = mesh.nedges*(ngl-2)
#     tot_faces_internal_nodes = mesh.nfaces*(ngl-2)*(ngl-2)
#     tot_vol_internal_nodes   = mesh.nelem*(ngl-2)^(mesh.nsd)
    
#     el_edges_internal_nodes = mesh.NEDGES_EL*(ngl-2)
#     el_faces_internal_nodes = mesh.NFACES_EL*(ngl-2)*(ngl-2)
#     el_vol_internal_nodes   = (ngl-2)^(mesh.nsd)
    
#     #Update number of grid points from linear count to total high-order points
#     mesh.npoin = tot_linear_poin + tot_edges_internal_nodes + tot_faces_internal_nodes + (mesh.nsd - 2)*tot_vol_internal_nodes
    
#     if (mesh.nop > 1)
#         println(" # GMSH HIGH-ORDER GRID PROPERTIES")
#         println(" # N. edges internal points   : ", tot_edges_internal_nodes)
#         println(" # N. faces internal points   : ", tot_faces_internal_nodes)
#         println(" # N. volumes internal points : ", tot_vol_internal_nodes)
#         println(" # N. total high order points : ", mesh.npoin)
#         println(" # GMSH HIGH-ORDER GRID PROPERTIES ...................... END")
#     end
    
#     #
#     # Resize as needed
#     #
#     mesh.x = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
#     mesh.y = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
#     mesh.z = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))

#     mesh.ip2gip = KernelAbstractions.zeros(backend, TInt, Int64(mesh.npoin))
#     mesh.gip2owner = KernelAbstractions.ones(backend, TInt, Int64(mesh.npoin))*local_views(parts).item_ref[]
    
    
#     mesh.conn_edge_el = KernelAbstractions.zeros(backend, TInt, 2, Int64(mesh.NEDGES_EL), Int64(mesh.nelem))    
#     mesh.conn_face_el = KernelAbstractions.zeros(backend, TInt,  4, Int64(mesh.NFACES_EL), Int64(mesh.nelem))  
#     mesh.bdy_edge_in_elem = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy))  
#     mesh.poin_in_edge = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges), Int64(mesh.ngl))
#     mesh.poin_in_bdy_edge = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy), Int64(mesh.ngl))
    
#     mesh.poin_in_face     = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nfaces), Int64(mesh.ngl), Int64(mesh.ngl))
#     mesh.edge_type        = Array{Union{Nothing, String}}(nothing, Int64(mesh.nedges))
#     mesh.bdy_edge_type    = Array{Union{Nothing, String}}(nothing, Int64(mesh.nedges_bdy))
#     mesh.bdy_edge_type_id = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy))  
    
#     if mesh.nsd > 2
#         mesh.poin_in_bdy_face = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nfaces_bdy), Int64(mesh.ngl), Int64(mesh.ngl))
#         mesh.face_type = Array{Union{Nothing, String}}(nothing, Int64(mesh.nfaces))
#         mesh.bdy_face_type = Array{Union{Nothing, String}}(nothing, Int64(mesh.nfaces_bdy))
#         mesh.bdy_face_in_elem = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nfaces_bdy))
#     end
#     mesh.npoin_el         = mesh.NNODES_EL + el_edges_internal_nodes + el_faces_internal_nodes + (mesh.nsd - 2)*el_vol_internal_nodes
#     mesh.conn = KernelAbstractions.zeros(backend,TInt, Int64(mesh.nelem), Int64(mesh.npoin_el))
    
#     #
#     # Connectivity matrices
#     #
#     mesh.cell_node_ids     = get_cell_node_ids(get_grid(model))
#     mesh.conn_unique_faces = get_face_nodes(model, FACE_flg) #faces --> 4 nodes
#     mesh.conn_unique_edges = get_face_nodes(model, EDGE_flg) #edges --> 2 nodes

#     mesh.cell_edge_ids     = get_faces(topology, mesh.nsd, 1) #edge map from local to global numbering i.e. iedge_g = cell_edge_ids[1:NELEM][1:NEDGES_EL]
#     mesh.cell_face_ids     = get_faces(topology, mesh.nsd, mesh.nsd-1) #face map from local to global numbering i.e. iface_g = cell_face_ids[1:NELEM][1:NFACE_EL]
#     mesh.face_edge_ids     = get_faces(topology,mesh.nsd-1, 1)
#     mesh.edge_g_color::Array{Int64, 1} = zeros(Int64, mesh.nedges)
#     # @info rank,node_global



#     if (mesh.nsd == 1)
#         nothing
#     elseif (mesh.nsd == 2)
    
#         mesh.connijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl),1)
    
#         for iel = 1:mesh.nelem
#             oiel = n2o_el_map[iel]
#             mesh.conn[iel, 1] = omesh.cell_node_ids[oiel][1]
#             mesh.conn[iel, 2] = omesh.cell_node_ids[oiel][2]
#             mesh.conn[iel, 3] = omesh.cell_node_ids[oiel][4]
#             mesh.conn[iel, 4] = omesh.cell_node_ids[oiel][3]

#             #
#             # 3-----4
#             # |     |
#             # |     |
#             # 1-----2
#             #
#             mesh.connijk[iel, 1,      1] = omesh.cell_node_ids[oiel][2]
#             mesh.connijk[iel, 1,    ngl] = omesh.cell_node_ids[oiel][1]
#             mesh.connijk[iel, ngl,  ngl] = omesh.cell_node_ids[oiel][3]
#             mesh.connijk[iel, ngl,    1] = omesh.cell_node_ids[oiel][4]
            
#             # @printf(" [1,1] [ngl, 1] [1, ngl] [ngl, ngl] %d %d %d %d\n", mesh.connijk[iel, 1, 1], mesh.connijk[iel, ngl, 1] , mesh.connijk[iel, 1,ngl], mesh.connijk[iel, ngl, ngl] )
#         end
#         #
#         # Fill in elements dictionary needed by NodeOrdering.jl
#         #
#         elements = Dict(
#             kk => mesh.conn[kk, 1:4]
#             for kk = 1:mesh.nelem)
#         element_types = Dict(
#             kk => :Quad4
#             for kk = 1:mesh.nelem)
        
#         #
#         # Rewrite coordinates in RCM order:
#         #
#         # filename = "./COORDS_LO_" + rank + ".dat" 
#         open("./COORDS_LO_$rank.dat", "w") do f
#             for ip = 1:mesh.npoin_linear
                
#                 mesh.x[ip] = get_node_coordinates(get_grid(model))[ip][1]
#                 mesh.y[ip] = get_node_coordinates(get_grid(model))[ip][2]
                
#                 mesh.ip2gip[ip] = point2ppoint[ip]
#                 # mesh.gip2owner[ip] = 1
#                 @printf(f, " %.6f %.6f 0.000000 %d %d\n", mesh.x[ip],  mesh.y[ip], ip, point2ppoint[ip])
#             end
#         end #f

#     elseif (mesh.nsd == 3)
        
#         mesh.connijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl), Int64(mesh.ngl))
#         mesh.conn_edgesijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.NEDGES_EL))

#         for iel = 1:mesh.nelem
#             #CGNS numbering: OK ref: HEXA...
#             mesh.conn[iel, 1] = mesh.cell_node_ids[iel][1]#9
#             mesh.conn[iel, 2] = mesh.cell_node_ids[iel][5]#11
#             mesh.conn[iel, 3] = mesh.cell_node_ids[iel][6]#6
#             mesh.conn[iel, 4] = mesh.cell_node_ids[iel][2]#1
#             mesh.conn[iel, 5] = mesh.cell_node_ids[iel][3]#10
#             mesh.conn[iel, 6] = mesh.cell_node_ids[iel][7]#12
#             mesh.conn[iel, 7] = mesh.cell_node_ids[iel][8]#5
#             mesh.conn[iel, 8] = mesh.cell_node_ids[iel][4]#4

#             #OK
#             mesh.connijk[iel, 1, 1, 1]       = mesh.cell_node_ids[iel][2]
#             mesh.connijk[iel, ngl, 1, 1]     = mesh.cell_node_ids[iel][1]
#             mesh.connijk[iel, ngl, ngl, 1]   = mesh.cell_node_ids[iel][5]
#             mesh.connijk[iel, 1, ngl, 1]     = mesh.cell_node_ids[iel][6]
#             mesh.connijk[iel, 1, 1, ngl]     = mesh.cell_node_ids[iel][4]
#             mesh.connijk[iel, ngl, 1, ngl]   = mesh.cell_node_ids[iel][3]
#             mesh.connijk[iel, ngl, ngl, ngl] = mesh.cell_node_ids[iel][7]
#             mesh.connijk[iel, 1, ngl, ngl]   = mesh.cell_node_ids[iel][8]
            
#         end
        
#         #
#         # Fill in elements dictionary needed by NodeOrdering.jl
#         #
#         elements = Dict(
#             kk => mesh.conn[kk, 1:8]
#             for kk = 1:mesh.nelem)
#         element_types = Dict(
#             kk => :Hexa8
#             for kk = 1:mesh.nelem)
        
#         #
#         #Use NodeNumbering.jl
#         #
#         #adjacency = create_adjacency_graph(elements, element_types)
#         #degrees = node_degrees(adjacency)
#         #neworder = RCM(adjacency, degrees, tot_linear_poin, tot_linear_poin)
#         #finalorder = renumbering(neworder)
#         #RCM_adjacency = create_RCM_adjacency(adjacency, finalorder)
#         #newmatrix = adjacency_visualization(RCM_adjacency)
#         #display(UnicodePlots.heatmap(newmatrix))
        
        
#         #
#         # Rewrite coordinates in RCM order:
#         #
#         open("./COORDS_LO_$rank.dat", "w") do f
#             #open("./COORDS_LO.dat", "w") do f
#             for ip = 1:mesh.npoin_linear
#                 mesh.x[ip] = get_node_coordinates(get_grid(model))[ip][1]
#                 mesh.y[ip] = get_node_coordinates(get_grid(model))[ip][2]
#                 mesh.z[ip] = get_node_coordinates(get_grid(model))[ip][3]
#                 mesh.ip2gip[ip] = point2ppoint[ip]
#                 # mesh.gip2owner[ip] = 1
#                 @printf(f, " %.6f %.6f %.6f %d %d\n", mesh.x[ip],  mesh.y[ip], mesh.z[ip], ip, point2ppoint[ip])
#             end
#         end #f
#     end


#     #
#     # Add high-order points to edges, faces, and elements (volumes)
#     #
#     # initialize LGL struct and buyild Gauss-Lobatto-xxx points
#     lgl = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.nop, backend)

#     println(" # POPULATE GRID with SPECTRAL NODES ............................ ")
#     #
#     # Edges
#     #
#     populate_conn_edge_el!(mesh, mesh.SD)
#     add_high_order_nodes_edges!(mesh, lgl, mesh.SD, backend, edge2pedge)

#     #
#     # Faces
#     #
#     populate_conn_face_el!(mesh, mesh.SD)
#     add_high_order_nodes_faces!(mesh, lgl, mesh.SD, face2pface)

#     #
#     # Volume
#     #
#     # NOTICE: in 2D we consider only edges. faces are the elements.
#     #         
#     add_high_order_nodes_volumes!(mesh, lgl, mesh.SD, elm2pelm)

#     mesh.gnpoin = MPI.Allreduce(maximum(mesh.ip2gip), MPI.MAX, comm)
#     # @info mesh.gnpoin
#     mesh.gip2owner = find_gip_owner(mesh.ip2gip)
    
#     for ip = mesh.npoin_linear+1:mesh.npoin
#         mesh.x[ip] = mesh.x_ho[ip]
#         mesh.y[ip] = mesh.y_ho[ip]
#         mesh.z[ip] = 0.0
#         if (mesh.nsd > 2)
#             mesh.z[ip] = mesh.z_ho[ip]
#         end
#     end
    
#     mesh.xmax = maximum(mesh.x)
#     mesh.xmin = minimum(mesh.x)
#     mesh.ymax = maximum(mesh.y)
#     mesh.ymin = minimum(mesh.y)
#     if (mesh.nsd > 2)
#         mesh.zmax = maximum(mesh.z)
#         mesh.zmin = minimum(mesh.z)
#     end

#     for ip = 1: mesh.npoin
#         if mesh.gip2owner[ip] != rank+1
#             # @info mesh.x[ip], mesh.y[ip], mesh.gip2owner[ip], rank+1
#         end
#     end

#     #----------------------------------------------------------------------
#     # Extract boundary edges and faces nodes:
#     #----------------------------------------------------------------------
#     #
#     # Bdy edges
#     #
#     if mesh.nsd == 2
#         # isboundary_edge = compute_isboundary_face(topology, EDGE_flg)
#         isboundary_edge = fill(false, mesh.nedges)  
        
#         # @info isboundary_edge
#         #
#         # Get labels contained in the current GMSH grid:
#         #
#         n_semi_inf = 0
#         labels = get_face_labeling(model)
#         # @info rank, labels.tag_to_name
#         for ilabel in labels.tag_to_name
#             edges_to_tag  = get_face_tag_index(labels,ilabel,EDGE_flg)
#             idx_edges_inflow = findall( x -> x == 1, edges_to_tag)
#             #    
#             # Tag the boundary edge with its type as defined in the user-provided GMSH file:
#             #
#             for idx in idx_edges_inflow
#                 mesh.edge_type[idx] = ilabel
#                 isboundary_edge[idx] = true
#             end
#             # @info mesh.edge_type
#         end
#         # @info rank, isboundary_edge
#         iedge_bdy = 1
#         for iedge = 1:mesh.nedges #total nedges
#             if isboundary_edge[iedge] == true
#                 # if rank == 1
#                 #     @info mesh.x[mesh.poin_in_edge[iedge, 1]], mesh.y[mesh.poin_in_edge[iedge, 1]]
#                 # end
#                 for igl = 1:mesh.ngl
#                     mesh.poin_in_bdy_edge[iedge_bdy, igl] = mesh.poin_in_edge[iedge, igl]
#                     mesh.bdy_edge_type[iedge_bdy] = mesh.edge_type[iedge]
#                 end
#                 if (mesh.bdy_edge_type[iedge_bdy] == "Laguerre")
#                     n_semi_inf += 1
#                 end
#                 iedge_bdy += 1
#             end
#         end
#         for iel = 1:mesh.nelem
#             for iedge_bdy = 1:mesh.nedges_bdy
#                 if issubset(mesh.poin_in_bdy_edge[iedge_bdy, :], mesh.connijk[iel, :, :])
#                     mesh.bdy_edge_in_elem[iedge_bdy] = iel
#                 end
#             end
#         end
#         # build mesh data structs for Laguerre semi-infinite elements
#         if ("Laguerre" in mesh.bdy_edge_type)
#             gr = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta],backend) 
#             factorx = inputs[:xfac_laguerre]#0.1
#             factory = inputs[:yfac_laguerre]#0.025
#             mesh.connijk_lag = KernelAbstractions.zeros(backend, TInt, Int64(n_semi_inf), Int64(mesh.ngl), Int64(mesh.ngr),1)
#             bdy_normals = zeros(n_semi_inf, 2)
#             bdy_tangents = zeros(n_semi_inf, 2)
#             e_iter = 1
#             iter = mesh.npoin + 1
#             x_new = KernelAbstractions.zeros(backend, TFloat, mesh.npoin + n_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)+mesh.ngr-1)
#             y_new = KernelAbstractions.zeros(backend, TFloat, mesh.npoin + n_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)+mesh.ngr-1)
#             x_new[1:mesh.npoin] .= mesh.x[:]
#             y_new[1:mesh.npoin] .= mesh.y[:]
#             for iedge = 1:size(mesh.bdy_edge_type,1)
#                 if (mesh.bdy_edge_type[iedge] == "Laguerre") 
#                     iel = mesh.bdy_edge_in_elem[iedge]
#                     #find tangent and normal vectors to the boundary
#                     ip = mesh.poin_in_bdy_edge[iedge,1]
#                     ip1 = mesh.poin_in_bdy_edge[iedge,2]
#                     #tangent vector 
#                     x = mesh.x[ip]
#                     x1 = mesh.x[ip1]
#                     y = mesh.y[ip]
#                     y1 = mesh.y[ip1]
#                     tan = [x-x1, y-y1]
#                     # deduce normal vector components
#                     if (tan[2] > 1e-7)
#                         x2 = 1.0
#                         y2 = -x2*tan[1]/tan[2]
#                     else
#                         y2 = 1.0
#                         x2 = -y2*tan[2]/tan[1]
#                     end
#                     nor = [x2,y2]
#                     # generate unit versions of tangent and normal vectors
#                     modu = sqrt(tan[1]^2+tan[2]^2)
#                     tan = tan * (1/modu)
#                     modu = sqrt(nor[1]^2+nor[2]^2)
#                     nor = nor * (1/modu)
#                     #make sure normal is outward facing
#                     l = 1
#                     m = 1
#                     l1 = 1
#                     m1 = 1
#                     for ii=1:mesh.ngl
#                         for jj=1:mesh.ngl
#                             if (mesh.connijk[iel,ii,jj] == ip)
#                                 l=ii
#                                 m=jj
#                             end
#                             if (mesh.connijk[iel,ii,jj] == ip1)
#                                 l1 = ii
#                                 m1 = jj
#                             end
#                         end
#                     end
#                     if (l == l1)
#                         ip2 = mesh.connijk[iel,3,m]
#                     else
#                         ip2 = mesh.connijk[iel,l,3]
#                     end
#                     v = [mesh.x[ip2]-x, mesh.y[ip2]-y]
#                     if (dot(v,nor) > 0.0)
#                         nor .= -nor
#                     end
#                     bdy_normals[e_iter,:] .= nor
#                     bdy_tangents[e_iter,:] .= tan
#                     for i=1:mesh.ngl
#                         ip = mesh.poin_in_bdy_edge[iedge,i]
#                         mesh.connijk_lag[e_iter,i,1] = ip
#                         for j=2:mesh.ngr
# 			                if (inputs[:xscale]==1.0)
#                                 x_temp = mesh.x[ip] + nor[1]*gr.ξ[j]*factorx
#                             else
#                                 x_temp = mesh.x[ip] + nor[1]*gr.ξ[j]*factorx/(inputs[:xscale] * 0.5)
# 			                end
#                             if (inputs[:yscale] == 1.0)
# 			                    y_temp = mesh.y[ip] + nor[2]*gr.ξ[j]*factory
# 			                else 
#                                 y_temp = mesh.y[ip] + nor[2]*gr.ξ[j]*factory/(inputs[:yscale] * 0.5)
#                             end
# 			                matched = 0
#                             if (i == mesh.ngl || i == 1)
#                                 iter_end = 0
#                                 while (matched == 0 && iter_end == 0)
#                                     for e_check = 1:n_semi_inf
#                                         for i_check =1:mesh.ngl
#                                             for j_check =1:mesh.ngr
#                                                 ip_check = mesh.connijk_lag[e_check,i_check,j_check]
#                                                 if (ip_check != 0.0 && e_check != e_iter)
#                                                     if (AlmostEqual(x_temp,x_new[ip_check]) && AlmostEqual(y_temp,y_new[ip_check]))
#                                                         mesh.connijk_lag[e_iter,i,j] = ip_check
#                                                         matched = 1
#                                                     end
#                                                 end
#                                             end
#                                         end 
#                                     end
#                                     iter_end = 1
#                                 end    
#                             else
#                                 x_new[iter] = x_temp#mesh.x[ip] + nor[1]*gr.ξ[j]*factorx
#                                 y_new[iter] = y_temp#mesh.y[ip] + nor[2]*gr.ξ[j]*factory
#                                 mesh.connijk_lag[e_iter,i,j] = iter
#                                 iter += 1
#                                 matched = 1
#                             end
#                             if (matched == 0)
#                                 x_new[iter] = x_temp#mesh.x[ip] + nor[1]*gr.ξ[j]*factorx
#                                 y_new[iter] = y_temp#mesh.y[ip] + nor[2]*gr.ξ[j]*factory
#                                 mesh.connijk_lag[e_iter,i,j] = iter
#                                 iter += 1 
#                             end
                            
#                             #@info nor[1],nor[2],x_new[iter],y_new[iter], mesh.x[ip],mesh.y[ip]
#                         end
#                     end
#                     e_iter += 1
#                 end
#             end
#             #@info mesh.npoin, iter - 1, mesh.ngr, n_semi_inf, e_iter - 1
#             mesh.npoin_original = mesh.npoin
#             mesh.npoin = iter -1
#             mesh.x = x_new
#             mesh.y = y_new
#             mesh.z = KernelAbstractions.zeros(backend, TFloat, mesh.npoin)
#             mesh.nelem_semi_inf = n_semi_inf
#         end

#     elseif mesh.nsd > 2
#         # isboundary_face = compute_isboundary_face(topology, FACE_flg)
#         isboundary_face = fill(false, mesh.nfaces)  
#         #
#         # Get labels contained in the current GMSH grid:
#         #
#         labels = get_face_labeling(model)
#         for ilabel in labels.tag_to_name
#             faces_to_tag  = get_face_tag_index(labels,ilabel,FACE_flg)
#             idx_faces_inflow = findall( x -> x == 1, faces_to_tag)
#             #    
#             # Tag the boundary edge with its type as defined in the user-provided GMSH file:
#             #
#             for idx in idx_faces_inflow
#                 mesh.face_type[idx] = ilabel
#                 isboundary_face[idx] = true
#             end
#         end
#         get_bdy_poin_in_face_on_edges!(mesh, @view(isboundary_face[:]), mesh.SD)
#         iface_bdy = 1
#         for iface in findall(x -> x == true, isboundary_face) #total nedges
#             # if isboundary_face[iface] == true
#             for igl = 1:mesh.ngl
#                 for jgl = 1:mesh.ngl
#                     mesh.poin_in_bdy_face[iface_bdy, igl,jgl] = mesh.poin_in_face[iface, igl,jgl]
#                     mesh.bdy_face_type[iface_bdy] = mesh.face_type[iface]
#                     #@info "face point number", mesh.poin_in_face[iface,igl,jgl],iface,igl,jgl
#                 end
#             end
#             iface_bdy += 1
#             # end
#         end
#         for iel = 1:mesh.nelem
#             for iface_bdy = 1:mesh.nfaces_bdy
#                 if issubset(mesh.poin_in_bdy_face[iface_bdy, :,:], mesh.connijk[iel, :, :, :])
#                     mesh.bdy_face_in_elem[iface_bdy] = iel
#                 end
#             end
#         end
#         #=for iface =1:mesh.nfaces_bdy
#             for i=1:mesh.ngl
#                 for j=1:mesh.ngl
#                     ip = mesh.poin_in_bdy_face[iface,i,j]
#                     @info "bdy points coords", mesh.x[ip],mesh.y[ip],mesh.z[ip]
#                 end
#             end
#         end=#
#     end


end  # End of module MyGeometry