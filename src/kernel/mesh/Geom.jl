module JeGeometry
__precompile__(false)
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
using SparseArrays

# Define your custom version of DiscreteModel function
function Gridap.Geometry.DiscreteModel(
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


const _setup_cell_dim = GridapGmsh._setup_cell_dim
const _setup_node_coords = GridapGmsh._setup_node_coords
const _setup_nodes_and_vertices = GridapGmsh._setup_nodes_and_vertices
const _setup_cell_to_vertices = GridapGmsh._setup_cell_to_vertices
const _setup_grid =GridapGmsh._setup_grid
const _setup_labeling = GridapGmsh._setup_labeling

function GridapGmsh.GmshDiscreteModel(gmsh::Module; has_affine_map=nothing, orient_if_simplex=nothing)

    Dc = _setup_cell_dim(gmsh)
    Dp = Dc
    # @info Dp
    node_to_coords = _setup_node_coords(gmsh,Dp)
    nnodes = length(node_to_coords)
    vertex_to_node, node_to_vertex = _setup_nodes_and_vertices(gmsh,node_to_coords)
    grid, cell_to_entity = _setup_grid(gmsh,Dc,Dp,node_to_coords,node_to_vertex;has_affine_map)
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




end  # End of module JeGeometry


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



mutable struct AssemblerCache
    global_max_index::Int
    index_a::Vector{Int}
    owner_a::Vector{Int}

    # Index communication buffers
    recv_idx_buffers::Vector{Vector{Int}}
    # combined_recv_idx::Vector{Int}

    # Send-back buffers
    recvback_idx_buffers::Vector{Vector{Int}}
    # combined_recv_back_idx::Vector{Int}

    # Compact global index mapping
    global_to_compact::Dict{Int, Int}     # Maps global index -> compact index
    local_to_compact::Dict{Int, Int}      # Maps local index -> compact index (precomputed)

    sum_array_1D::Vector{Float64}
    sum_array_2D::Matrix{Float64}

    # auxiliary
    send_i::Vector{Vector{Int}} 
    send_data_buffers::Vector{Vector{Float64}}
    recv_data_buffers::Vector{Vector{Float64}}
    send_data_sizes::Vector{Int}
    recv_data_sizes::Vector{Int}
    # i_local::Dict{Int, Vector{Int}}

end

function setup_assembler(SD, a, index_a, owner_a)

    if SD == NSD_1D() return nothing end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)

    global_max_index = maximum(index_a)

    m = size(a, 2)

    # i_local = Dict{Int,  Vector{Int}}()
    # for (i, idx) in enumerate(index_a)
    #     owner = owner_a[i]
    #     if owner == rank
    #         idx_i = get!(i_local, idx, Int[])
    #         push!(idx_i, i)
    #     end
    # end
    # filter!(p -> !isempty(p.second), i_local)

    send_idx = Dict(i => Int[] for i in 0:rank_sz-1)
    send_i = [Int[] for i in 0:rank_sz-1]
    for (i, idx) in enumerate(index_a)
        owner = owner_a[i]
        if owner != rank
            buf_idx = get!(send_idx, owner, Int[])
            push!(buf_idx, idx)
            push!(send_i[owner+1], i)
        end
    end


    send_idx_sizes = [length(send_idx[i]) for i in 0:rank_sz-1]
    recv_idx_sizes = MPI.Alltoall(MPI.UBuffer(send_idx_sizes, 1), comm)
    MPI.Barrier(comm)

    # Prepare buffers for sending and receiving data
    send_idx_buffers = [send_idx[i] for i in 0:rank_sz-1]
    recv_idx_buffers = [Vector{Int}(undef, recv_idx_sizes[i+1]) for i in 0:rank_sz-1]

    # Communicate data
    requests = MPI.Request[]
    for i in 0:rank_sz-1
        if send_idx_sizes[i+1] > 0
            push!(requests, MPI.Isend(send_idx_buffers[i+1], i, 1, comm))
        end
        if recv_idx_sizes[i+1] > 0
            push!(requests, MPI.Irecv!(recv_idx_buffers[i+1], i, 1, comm))
        end
    end

    # Wait for all communication to complete
    MPI.Waitall!(requests)

    # Combine received data into a single vector
    # combined_recv_idx = Int[]
    # for i in 0:rank_sz-1
    #     if recv_idx_sizes[i+1] > 0
    #         append!(combined_recv_idx, recv_idx_buffers[i+1])
    #     end
    # end

    # send data back to original ranks
    sendback_idx = Dict(i => Int[] for i in 0:rank_sz-1)
    for i in 0:rank_sz-1
        if recv_idx_sizes[i+1] > 0
            buf_idx = get!(sendback_idx, i, Int[])
            for idx in recv_idx_buffers[i+1]
                push!(buf_idx, idx)
            end
        end
    end
    sendback_idx_sizes = recv_idx_sizes
    recvback_idx_sizes = send_idx_sizes
    MPI.Barrier(comm)



    # Prepare buffers for sending and receiving back data
    sendback_idx_buffers = [sendback_idx[i] for i in 0:rank_sz-1]
    recvback_idx_buffers = [Vector{Int}(undef, length(send_idx[i])) for i in 0:rank_sz-1]


    # Communicate back data
    requests_back = MPI.Request[]
    for i in 0:rank_sz-1
        if sendback_idx_sizes[i+1] > 0
            push!(requests_back, MPI.Isend(sendback_idx_buffers[i+1], i, 3, comm))
        end
        if recvback_idx_sizes[i+1] > 0
            push!(requests_back, MPI.Irecv!(recvback_idx_buffers[i+1], i, 3, comm))
        end
    end


    # Wait for all communication to complete
    MPI.Waitall!(requests_back)



    # Combine received data into a single vector
    # combined_recv_back_idx = Int[]
    # for i in 0:rank_sz-1
    #     if recvback_idx_sizes[i+1] > 0
    #         append!(combined_recv_back_idx, recvback_idx_buffers[i+1])
    #     end
    # end

    needed_indices = Set{Int}()
    
    # 1. Add repeated indices from local index_a (owned by this rank)
    local_index_counts = Dict{Int, Int}()
    for (i, idx) in enumerate(index_a)
        if owner_a[i] == rank
            local_index_counts[idx] = get(local_index_counts, idx, 0) + 1
        end
    end
    for (idx, count) in local_index_counts
        if count > 1
            push!(needed_indices, idx)
        end
    end
    
    # 2. Add all indices from recv_idx_buffers
    for i in 0:rank_sz-1
        if recv_idx_sizes[i+1] > 0
            for idx in recv_idx_buffers[i+1]
                push!(needed_indices, idx)
            end
        end
    end
    
    # 3. Add all indices from recvback_idx_buffers
    for i in 0:rank_sz-1
        if recvback_idx_sizes[i+1] > 0
            for idx in recvback_idx_buffers[i+1]
                push!(needed_indices, idx)
            end
        end
    end
    
    # Create compact mapping
    compact_global_indices = sort(collect(needed_indices))
    global_to_compact = Dict{Int, Int}()
    for (compact_idx, global_idx) in enumerate(compact_global_indices)
        global_to_compact[global_idx] = compact_idx
    end
    
    compact_size = length(compact_global_indices)
    
    # Allocate compact sum arrays
    sum_array_1D = zeros(Float64, compact_size)
    sum_array_2D = zeros(Float64, compact_size, m)

    # create local to compact Dict
    local_to_compact = Dict{Int, Int}()

    for (local_idx, global_idx) in enumerate(index_a)
        compact_idx = get(global_to_compact, global_idx, 0)
        if compact_idx > 0
            local_to_compact[local_idx] = compact_idx
        end
    end


    # sum_array_1D = zeros(Float64, global_max_index)
    # sum_array_2D = zeros(Float64, global_max_index, m)

    
    send_data_sizes = [send_idx_sizes[i+1] * m for i in 0:rank_sz-1]
    recv_data_sizes = [recv_idx_sizes[i+1] * m for i in 0:rank_sz-1]

    send_data_buffers = [zeros(Float64, send_idx_sizes[i+1] * m) for i in 0:rank_sz-1]
    recv_data_buffers = [zeros(Float64, recv_idx_sizes[i+1] * m) for i in 0:rank_sz-1]



    cache = AssemblerCache(global_max_index, index_a, owner_a,
            recv_idx_buffers, recvback_idx_buffers,
            global_to_compact, local_to_compact,
            sum_array_1D, sum_array_2D,
            send_i,send_data_buffers,recv_data_buffers, send_data_sizes, recv_data_sizes)
    return cache
end


function assemble_mpi!(a, cache::AssemblerCache)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    T = eltype(a)

    is1D = ndims(a) == 1
    n = size(a, 1)
    m = is1D ? 1 : size(a, 2)

    if is1D
        fill!(cache.sum_array_1D, zero(T))
    else
        fill!(cache.sum_array_2D, zero(T))

    end




    @inbounds for (i, idx) in cache.local_to_compact
        owner = cache.owner_a[i]
        if owner == rank
            if is1D
                cache.sum_array_1D[idx] += a[i]
            else
                for j = 1:m
                    cache.sum_array_2D[idx, j] += a[i, j]
                end
            end
        end
    end


    for i in 0:rank_sz-1
        fill!(cache.send_data_buffers[i+1], zero(T))
    end

    @inbounds for owner = 0:rank_sz-1
            buf_data = cache.send_data_buffers[owner+1]
            send_i_local = cache.send_i[owner+1]
            if is1D
                for (i,idx) in enumerate(send_i_local)
                    buf_data[i] = a[idx]
                end
            else

                for (i,idx) in enumerate(send_i_local)
                    for j = 1:m
                    buf_data[(i-1)*m + j] = a[idx,j]
                end
            end
        end
    end
    # send_data_sizes = [length(cache.send_data_buffers[i+1]) for i in 0:rank_sz-1]
    # recv_data_sizes = MPI.Alltoall(MPI.UBuffer(send_data_sizes, 1), comm)
    MPI.Barrier(comm)

    # Prepare buffers for sending and receiving data
    # recv_data_buffers = [Vector{T}(undef, recv_data_sizes[i+1]) for i in 0:rank_sz-1]

    for i in 0:rank_sz-1
        fill!(cache.recv_data_buffers[i+1], zero(T))
    end


    # Communicate data
    requests = MPI.Request[]
    for i in 0:rank_sz-1
        if cache.send_data_sizes[i+1] > 0
            push!(requests, MPI.Isend(cache.send_data_buffers[i+1], i, 0, comm))
        end
        if cache.recv_data_sizes[i+1] > 0
            push!(requests, MPI.Irecv!(cache.recv_data_buffers[i+1], i, 0, comm))
        end
    end

    # Wait for all communication to complete
    MPI.Waitall!(requests)

    # Combine received data into a single vector
    # combined_recv_data = vcat(cache.recv_data_buffers...)
    # combined_recv_data = T[]
    # for i in 0:rank_sz-1
    #     if recv_data_sizes[i+1] > 0
    #         append!(combined_recv_data, cache.recv_data_buffers[i+1])
    #     end
    # end
    @inbounds for rk in 0:rank_sz-1
        if cache.recv_data_sizes[rk+1] > 0
            buffer = cache.recv_data_buffers[rk+1]
            for (i, idx) in enumerate(cache.recv_idx_buffers[rk+1])
                compact_idx = cache.global_to_compact[idx]
                if is1D
                    cache.sum_array_1D[compact_idx] += buffer[i]
                else
                    for j = 1:m
                        cache.sum_array_2D[compact_idx, j] += buffer[(i-1)*m+j]
                    end
                end
            end
        end
    end

    # send data back to original ranks
    sendback_data_buffers = cache.recv_data_buffers
    @inbounds for i in 0:rank_sz-1
        if cache.recv_data_sizes[i+1] > 0
            buf_data = sendback_data_buffers[i+1]
            for (j, idx) in enumerate(cache.recv_idx_buffers[i+1])
                compact_idx = cache.global_to_compact[idx]
                if is1D
                    buf_data[j] = cache.sum_array_1D[compact_idx]
                else
                    for k = 1:m
                        buf_data[(j-1)*m+k] = cache.sum_array_2D[compact_idx, k]
                    end
                end
            end
        end
    end
    sendback_data_sizes = cache.recv_data_sizes
    recvback_data_sizes = cache.send_data_sizes
    MPI.Barrier(comm)



    # Prepare buffers for sending and receiving back data
    recvback_data_buffers = cache.send_data_buffers


    # Communicate back data
    requests_back = MPI.Request[]
    @inbounds for i in 0:rank_sz-1
        if sendback_data_sizes[i+1] > 0
            push!(requests_back, MPI.Isend(sendback_data_buffers[i+1], i, 2, comm))
        end
        if recvback_data_sizes[i+1] > 0
            push!(requests_back, MPI.Irecv!(recvback_data_buffers[i+1], i, 2, comm))
        end
    end


    # Wait for all communication to complete
    MPI.Waitall!(requests_back)


    # Combine received data into a single vector
    # combined_recv_back_data = vcat(recvback_data_buffers...)
#     @time begin 
#     combined_recv_back_data = T[]
#     for i in 0:rank_sz-1
#         if recvback_data_sizes[i+1] > 0
#             append!(combined_recv_back_data, recvback_data_buffers[i+1])
#         end
#     end
# end
    @inbounds for rk = 0:rank_sz-1
        if recvback_data_sizes[rk+1] > 0
            buffer = recvback_data_buffers[rk+1]
            for (i, idx) in enumerate(cache.recvback_idx_buffers[rk+1])
                compact_idx = cache.global_to_compact[idx]
                if is1D
                    cache.sum_array_1D[compact_idx] = buffer[i]
                else
                    for j = 1:m
                        cache.sum_array_2D[compact_idx, j] = buffer[(i-1)*m+j]
                    end
                end
            end
        end
    end


    @inbounds for (i, idx) in cache.local_to_compact
        if is1D
            a[i] = cache.sum_array_1D[idx]
        else
            for j = 1:m
                a[i, j] = cache.sum_array_2D[idx, j]
            end
        end
    end
end
