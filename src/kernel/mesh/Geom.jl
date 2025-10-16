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




mutable struct CyclingReverseDict
    mapping::Dict{Int, Vector{Int}}
    counters::Dict{Int, Int}
    repeated_keys::Vector{Int}
end

function CyclingReverseDict(a::Vector)
    mapping = Dict{Int, Vector{Int}}()
    repeated_keys = Int[]
    for (i, val) in enumerate(a)
        if haskey(mapping, val)
            if length(mapping[val]) == 1
                push!(repeated_keys, val)
            end
            push!(mapping[val], i)
        else
            mapping[val] = [i]
        end
    end
    CyclingReverseDict(mapping, Dict{Int, Int}(), repeated_keys)
end

function get_repeated_keys(crd::CyclingReverseDict)
    return crd.repeated_keys
end

function get_vals(crd::CyclingReverseDict, key::Int; all = false, first = false, order = true)
    indices = crd.mapping[key]
    if all == true
        return indices
    end
    if first == true
        return indices[1]
    end
    if order == true
        # Get current counter for this key (default to 0)
        counter = get(crd.counters, key, 0)
        # Get the index to return
        result = indices[counter % length(indices) + 1]
        # Increment counter
        crd.counters[key] = counter + 1
        return result
    end
end

# for non-periodic only
mutable struct AssemblerCache
    # Index communication buffers
    recv_idx_buffers::Vector{Vector{Int}}
    # combined_recv_idx::Vector{Int}

    # Send-back buffers
    recvback_idx_buffers::Vector{Vector{Int}}

    # auxiliary
    send_i::Vector{Vector{Int}} 
    send_data_buffers::Vector{Vector{Float64}}
    recv_data_buffers::Vector{Vector{Float64}}
    send_data_sizes::Vector{Int}
    recv_data_sizes::Vector{Int}
    # i_local::Dict{Int, Vector{Int}}

    # Preallocated requests
    requests::MPI.MultiRequest
    requests_back::MPI.MultiRequest
end

function setup_assembler(SD, a, index_a, owner_a)

    if SD == NSD_1D() return nothing end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)

    global_max_index = maximum(index_a)

    m = size(a, 2)


    send_idx = Dict(i => Int[] for i in 0:rank_sz-1)
    send_i = [Int[] for i in 0:rank_sz-1]
    # send list remote
    for (i, idx) in enumerate(index_a)
        owner = owner_a[i]
        if owner != rank
            buf_idx = get!(send_idx, owner, Int[])
            push!(buf_idx, idx)
            push!(send_i[owner+1], i)
        end
    end
    # send list local (for periodic)
    a_g2l_idx      = CyclingReverseDict(index_a)
    a_g2l_repeated = get_repeated_keys(a_g2l_idx)
    for idx in a_g2l_repeated
        local_idx = get_vals(a_g2l_idx, idx; all = true)[2:end]
        for i in local_idx
            owner = owner_a[i]
            if owner == rank
                buf_idx = get!(send_idx, owner, Int[])
                push!(buf_idx, idx)
                push!(send_i[owner+1], i)
            end
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
    MPI.Waitall(requests)


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
    MPI.Waitall(requests_back)



    # needed_indices = Set{Int}()
    
    # # 2. Add all indices from recv_idx_buffers
    # for i in 0:rank_sz-1
    #     if recv_idx_sizes[i+1] > 0
    #         for idx in recv_idx_buffers[i+1]
    #             push!(needed_indices, idx)
    #         end
    #     end
    # end
    
    # # 3. Add all indices from recvback_idx_buffers
    # for i in 0:rank_sz-1
    #     if recvback_idx_sizes[i+1] > 0
    #         for idx in recvback_idx_buffers[i+1]
    #             push!(needed_indices, idx)
    #         end
    #     end
    # end


    # change recv_idx_buffers and recvback_idx_buffers to local_idx
    for rk in 1:rank_sz
        recv_idx_buffers_rk     = recv_idx_buffers[rk]
        recvback_idx_buffers_rk = recvback_idx_buffers[rk]
        for (i,idx) in enumerate(recv_idx_buffers_rk)
            recv_idx_buffers_rk[i] = get_vals(a_g2l_idx, idx; first = true)
        end
        if rk-1 == rank
            for (i,idx) in enumerate(recvback_idx_buffers_rk)
                recvback_idx_buffers_rk[i] = send_i[rk][i]
            end
        else
            for (i,idx) in enumerate(recvback_idx_buffers_rk)
                recvback_idx_buffers_rk[i] = get_vals(a_g2l_idx, idx)
            end
        end
    end

    send_data_sizes = [send_idx_sizes[i+1] * m for i in 0:rank_sz-1]
    recv_data_sizes = [recv_idx_sizes[i+1] * m for i in 0:rank_sz-1]

    send_data_buffers = [zeros(Float64, send_idx_sizes[i+1] * m) for i in 0:rank_sz-1]
    recv_data_buffers = [zeros(Float64, recv_idx_sizes[i+1] * m) for i in 0:rank_sz-1]

    # Preallocate requests for assemble function
    n_req = sum(send_data_sizes .> 0) + sum(recv_data_sizes .> 0)
    n_req_back = sum(recv_data_sizes .> 0) + sum(send_data_sizes .> 0)

    return AssemblerCache( recv_idx_buffers, recvback_idx_buffers,
            send_i,send_data_buffers,recv_data_buffers, send_data_sizes, recv_data_sizes,
            MPI.MultiRequest(n_req),
            MPI.MultiRequest(n_req_back))
end


function assemble_mpi!(a, cache::AssemblerCache)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    T = eltype(a)

    is1D = ndims(a) == 1
    n = size(a, 1)
    m = is1D ? 1 : size(a, 2)

    for i in 0:rank_sz-1
        fill!(cache.send_data_buffers[i+1], zero(T))
    end

    @inbounds for owner = 0:rank_sz-1
        buf_data = cache.send_data_buffers[owner+1]
        send_i_local = cache.send_i[owner+1]

        for (i,idx) in enumerate(send_i_local)
            for j = 1:m
                buf_data[(i-1)*m + j] = a[idx,j]
            end
        end
    end
    # MPI.Barrier(comm)


    for i in 0:rank_sz-1
        fill!(cache.recv_data_buffers[i+1], zero(T))
    end


    # Communicate data
    req_idx = 1
    @inbounds for i in 0:rank_sz-1
        if cache.send_data_sizes[i+1] > 0
            MPI.Isend(cache.send_data_buffers[i+1], comm, cache.requests[req_idx]; dest=i, tag=0)
            req_idx += 1
        end
        if cache.recv_data_sizes[i+1] > 0
            MPI.Irecv!(cache.recv_data_buffers[i+1], comm, cache.requests[req_idx]; source=i, tag=0)
            req_idx += 1
        end
    end

    # Wait for all communication to complete
    MPI.Waitall(cache.requests)

    @inbounds for rk in 0:rank_sz-1
        if cache.recv_data_sizes[rk+1] > 0
            buffer = cache.recv_data_buffers[rk+1]
            for (i, local_idx) in enumerate(cache.recv_idx_buffers[rk+1])
                for j = 1:m
                    a[local_idx, j] += buffer[(i-1)*m+j]
                end
            end
        end
    end

    # send data back to original ranks
    sendback_data_buffers = cache.recv_data_buffers
    @inbounds for rk in 0:rank_sz-1
        if cache.recv_data_sizes[rk+1] > 0
            buf_data = sendback_data_buffers[rk+1]
            for (i, local_idx) in enumerate(cache.recv_idx_buffers[rk+1])
                for j = 1:m
                    buf_data[(i-1)*m+j] = a[local_idx, j]
                end
            end
        end
    end
    sendback_data_sizes = cache.recv_data_sizes
    recvback_data_sizes = cache.send_data_sizes
    # MPI.Barrier(comm)



    # Prepare buffers for sending and receiving back data
    recvback_data_buffers = cache.send_data_buffers


    # Communicate back data
    req_idx = 1
    @inbounds for i in 0:rank_sz-1
        if sendback_data_sizes[i+1] > 0
            MPI.Isend(sendback_data_buffers[i+1], comm, cache.requests_back[req_idx]; dest=i, tag=2)
            req_idx += 1
        end
        if recvback_data_sizes[i+1] > 0
            MPI.Irecv!(recvback_data_buffers[i+1], comm, cache.requests_back[req_idx]; source=i, tag=2)
            req_idx += 1
        end
    end


    # Wait for all communication to complete
    MPI.Waitall(cache.requests_back)


    @inbounds for rk = 0:rank_sz-1
        if recvback_data_sizes[rk+1] > 0
            buffer = recvback_data_buffers[rk+1]
            for (i, local_idx) in enumerate(cache.recvback_idx_buffers[rk+1])
                for j = 1:m
                    a[local_idx, j] = buffer[(i-1)*m+j]
                end
            end
        end
    end
end