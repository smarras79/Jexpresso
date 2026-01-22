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

function setup_global_numbering_extra_dim(ip2gip, gip2owner, npoin, npoin_ang, npoin_total)

    comm = MPI.COMM_WORLD

    ip2gip_extra = KernelAbstractions.zeros(CPU(),Int64,npoin_total)
    for ip = 1:npoin
        gip = ip2gip[ip]
        for ip_ext = 1:npoin_ang
            idx_ip = (ip-1)*(npoin_ang) + ip_ext
            idx_gip = (gip-1)*(npoin_ang) + ip_ext
            #=if (ip != gip)
                @info ip, gip
            end=#
            ip2gip_extra[idx_ip] = idx_gip
        end
    end
    
    gnpoin    = MPI.Allreduce(maximum(ip2gip_extra), MPI.MAX, comm)
    gip2owner_extra = find_gip_owner(ip2gip_extra)
    gip2ip    = KernelAbstractions.zeros(CPU(), TInt, gnpoin)
    @info gnpoin, npoin_total
    for (ip, gip) in enumerate(ip2gip_extra)
        gip2ip[gip] = ip
    end

    return ip2gip_extra, gip2owner_extra, gnpoin
end

function setup_global_numbering_adaptive_angular_scalable(
    ip2gip, gip2owner, mesh, connijk_spa,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    n_spa, n_non_global_nodes
)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nproc = MPI.Comm_size(comm)
    
    nelem = mesh.nelem
    ngl = mesh.ngl

    # =========================================================================
    # PHASE 1: Build local unique signatures
    # =========================================================================
    
    local_points = Dict{NTuple{3,Float64}, Int}()
    signature_list = NTuple{3,Float64}[]
    
    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel, i, j, k]
            gip = ip2gip[ip]
            
            for e_ext = 1:extra_meshes_extra_nelems[iel]
                for jθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                    for iθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                        ip_ext = extra_meshes_connijk[iel][e_ext, iθ, jθ]
                        θ = extra_meshes_coords[iel][1, ip_ext]
                        ϕ = extra_meshes_coords[iel][2, ip_ext]
                        ip_spa = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                        
                        sig = (Float64(gip), round(θ, digits=12), round(ϕ, digits=12))
                        
                        if !haskey(local_points, sig)
                            local_points[sig] = ip_spa
                            push!(signature_list, sig)
                        end
                    end
                end
            end
        end
    end
    sort!(signature_list)
    n_local = length(signature_list)
    
    @info "[Rank $rank] Local unique points: $n_local"

    # =========================================================================
    # PHASE 2: Identify processor boundary spatial-angular points
    # =========================================================================
    
    # A spatial-angular point is on a processor boundary if:
    # 1. Its spatial node (gip) is NOT owned by this processor, OR
    # 2. Its spatial node (gip) IS owned but exists on other processors too
    
    # Find spatial nodes on processor boundaries
    processor_boundary_spatial = Set{Int}()
    
    for ip = 1:mesh.npoin
        gip = ip2gip[ip]
        owner = gip2owner[ip]
        
        # If this processor doesn't own this spatial node, it's on a boundary
        if owner != rank
            push!(processor_boundary_spatial, gip)
        end
    end

    # Additionally, find spatial nodes this processor owns but are shared
    # (These are nodes owned by this rank but present on other ranks too)
    owned_spatial = Set{Int}()
    for ip = 1:mesh.npoin
        gip = ip2gip[ip]
        if gip2owner[ip] == rank
            push!(owned_spatial, gip)
        end
    end

    local_spatial = collect(processor_boundary_spatial)
    n_local = Int32(length(local_spatial))

    # Gather the counts from all processors
    counts = MPI.Allgather([n_local], comm)

    # Prepare the receive buffer
    total_count = sum(counts)

    # Handle empty case - need a typed array even if empty
    if isempty(local_spatial)
        local_spatial = Int[]  # Empty array with correct type
    end

    # Use Allgatherv to gather variable-length data
    if total_count > 0
        all_owned_spatial = MPI.Allgatherv(local_spatial, counts, comm)
    else
        # All ranks are empty
        all_owned_spatial = Int[]
    end
    
    for other_rank = 0:(nproc-1)
        if other_rank == 0
            offset = 0
        else
            offset = sum(counts[1:other_rank])
        end

        if other_rank == rank
            continue
        end
        
        other_not_owned = Set(all_owned_spatial[offset+1:counts[other_rank+1]])
        
        # Find intersection: spatial nodes owned by us but also present on other rank
        shared = intersect(owned_spatial, other_not_owned)
        
        union!(processor_boundary_spatial, shared)
        
    end

    @info "[Rank $rank] Found $(length(processor_boundary_spatial)) spatial nodes on processor boundaries"
    
    # =========================================================================
    # PHASE 3: Extract processor boundary spatial-angular signatures
    # =========================================================================
    
    processor_boundary_sigs = Dict{NTuple{3,Float64}, Int}()
    
    for sig in signature_list
        gip_spatial = Int(sig[1])  # First element is the global spatial ID
        
        if gip_spatial in processor_boundary_spatial
            idx = findfirst(x -> x == sig, signature_list)
            processor_boundary_sigs[sig] = idx  # Store local index for now
        end
    end
    
    @info "[Rank $rank] Found $(length(processor_boundary_sigs)) spatial-angular points on processor boundaries"
    
    # =========================================================================
    # PHASE 4: Exchange boundary signatures and resolve duplicates
    # =========================================================================
    
    # Gather all processor boundary signatures from all ranks
    local_keys = collect(keys(processor_boundary_sigs))
    n_local = length(local_keys)

    # Flatten tuples to a matrix (3 x n_local)
    local_data = hcat([collect(k) for k in local_keys]...)  # or stack them appropriately
   
    # Gather counts from all processors
    counts = MPI.Allgather([n_local], comm)

    # Prepare receive buffer
    total_count = sum(counts)
    
    all_data = MPI.Allgatherv(local_data[:], Int32.(counts .* 3), comm)

    # Reshape back to tuples
    all_boundary_sigs = [Tuple(all_data[i:i+2]) for i in 1:3:length(all_data)]
    #all_boundary_sigs = MPI.Allgather(collect(keys(processor_boundary_sigs)), comm)
    
    # Find globally shared signatures
    global_boundary_sigs = Set{NTuple{3,Float64}}()
    sig_counts = Dict{NTuple{3,Float64}, Int}()

    for sig in all_boundary_sigs
        sig_counts[sig] = get(sig_counts, sig, 0) + 1
    end
    
    # A signature is shared if it appears on multiple processors
    for (sig, count) in sig_counts
        if count > 1
            push!(global_boundary_sigs, sig)
        end
    end
    
    @info "[Rank $rank] Found $(length(global_boundary_sigs)) globally shared spatial-angular points"
    
    # =========================================================================
    # PHASE 5: Assign tentative global IDs using parallel prefix sum
    # =========================================================================
    
    # Each processor gets a contiguous range of global IDs
    n_local = length(signature_list)
    local_count = n_local
    offset = MPI.Exscan(local_count, MPI.SUM, comm)
    
    if rank == 0
        offset = 0
    end
    
    # Build tentative mapping
    sig_to_tentative_gid = Dict{NTuple{3,Float64}, Int}()
    for (idx, sig) in enumerate(signature_list)
        sig_to_tentative_gid[sig] = offset + idx
    end
    
    @info "[Rank $rank] Tentative global ID range: [$(offset+1), $(offset+n_local)]"

    # =========================================================================
    # PHASE 6: Resolve shared signatures (take minimum GID)
    # =========================================================================
    
    # For shared signatures, gather all tentative GIDs and take minimum
    shared_sig_resolution = Dict{NTuple{3,Float64}, Int}()
    
    # Each processor broadcasts its tentative GIDs for shared signatures
    local_shared_tentative = Dict{NTuple{3,Float64}, Int}()
    for sig in global_boundary_sigs
        if haskey(sig_to_tentative_gid, sig)
            local_shared_tentative[sig] = sig_to_tentative_gid[sig]
        end
    end
    
    # Serialize the dictionary
    local_buffer = IOBuffer()
    serialize(local_buffer, local_shared_tentative)
    local_data = take!(local_buffer)

    # Gather the sizes
    n_local = Int32(length(local_data))
    counts = MPI.Allgather([n_local], comm)

    # Gather the serialized data
    total_count = sum(counts)
    if total_count > 0
        all_data = MPI.Allgatherv(local_data, counts, comm)
    
        # Deserialize on each rank to get all dictionaries
        all_shared_tentative = Dict{NTuple{3,Float64}, Int}[]
        offset = 1
        for count in counts
            if count > 0
                chunk = all_data[offset:offset+count-1]
                push!(all_shared_tentative, deserialize(IOBuffer(chunk)))
            else
                push!(all_shared_tentative, Dict{NTuple{3,Float64}, Int}())
            end
            offset += count
        end
    else
        all_shared_tentative = [Dict{NTuple{3,Float64}, Int}() for _ in 1:length(counts)]
    end

    
    #all_shared_tentative = MPI.Allgather(local_shared_tentative, comm)
    
    # Resolve: take minimum tentative GID for each shared signature
    for proc_shared in all_shared_tentative
        for (sig, tentative_gid) in proc_shared
            if haskey(shared_sig_resolution, sig)
                shared_sig_resolution[sig] = min(shared_sig_resolution[sig], tentative_gid)
            else
                shared_sig_resolution[sig] = tentative_gid
            end
        end
    end
    
    @info "[Rank $rank] Resolved $(length(shared_sig_resolution)) shared signatures"

    # =========================================================================
    # PHASE 7: Build final global IDs (with duplicates for shared points)
    # =========================================================================
    
    # For shared points, use resolved GID; for interior points, use tentative GID
    sig_to_final_gid = Dict{NTuple{3,Float64}, Int}()
    
    for sig in signature_list
        if haskey(shared_sig_resolution, sig)
            sig_to_final_gid[sig] = shared_sig_resolution[sig]
        else
            sig_to_final_gid[sig] = sig_to_tentative_gid[sig]
        end
    end

    # =========================================================================
    # PHASE 8: Compact numbering to remove gaps
    # =========================================================================
    
    # Collect all used global IDs across all processors
    #local_used_gids = Set(values(sig_to_final_gid))
    local_used_gids = collect(values(sig_to_final_gid))
    n_local = Int32(length(local_used_gids))

    # Gather the counts from all processors
    counts = MPI.Allgather([n_local], comm)

    # Prepare the receive buffer
    total_count = sum(counts)

    # Handle empty case - need a typed array even if empty
    if isempty(local_used_gids)
        local_used_gids = Int[]  # Empty array with correct type
    end

    # Use Allgatherv to gather variable-length data
    if total_count > 0
        all_used_gids_list = MPI.Allgatherv(local_used_gids, counts, comm)
    else
        # All ranks are empty
        all_used_gids_list = Int[]
    end
    
    #all_used_gids_list = MPI.Allgather(collect(local_used_gids), comm)
    
    global_used_gids = Set{Int}(all_used_gids_list)
    
    # Create compaction map: old_gid -> new_gid (1:gnpoin with no gaps)
    sorted_gids = sort(collect(global_used_gids))
    old_to_new = Dict{Int, Int}()
    for (new_id, old_id) in enumerate(sorted_gids)
        old_to_new[old_id] = new_id
    end
    
    gnpoin = length(sorted_gids)
    
    @info "[Rank $rank] Compacted to $gnpoin global spatial-angular points (no gaps)"

    # =========================================================================
    # PHASE 9: Apply global numbering to local spatial-angular points
    # =========================================================================
    
    ip2gip_spa = zeros(Int64, n_spa)
    
    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel, i, j, k]
            gip = ip2gip[ip]
            
            for e_ext = 1:extra_meshes_extra_nelems[iel]
                for jθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                    for iθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                        ip_ext = extra_meshes_connijk[iel][e_ext, iθ, jθ]
                        θ = extra_meshes_coords[iel][1, ip_ext]
                        ϕ = extra_meshes_coords[iel][2, ip_ext]
                        ip_spa = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                        
                        sig = (Float64(gip), round(θ, digits=12), round(ϕ, digits=12))
                        
                        # Get final compacted global ID
                        old_gid = sig_to_final_gid[sig]
                        new_gid = old_to_new[old_gid]
                        
                        ip2gip_spa[ip_spa] = new_gid
                    end
                end
            end
        end
    end
    
    # Verify compactness
    verify_compact_numbering(ip2gip_spa, gnpoin, rank, comm)
    
    gip2owner_spa = find_gip_owner_spa(ip2gip_spa, n_spa, gnpoin, comm)
    
    gip2ip = zeros(Int, gnpoin)
    for (ip, gid) in enumerate(ip2gip_spa)
        if gid > 0 && gip2owner_spa[gid] == rank
            gip2ip[gid] = ip
        end
    end

    return ip2gip_spa, gip2ip, gip2owner_spa, gnpoin
end

function find_gip_owner_spa(ip2gip_spa, n_spa, gnpoin, comm)
    """
    Determine which processor owns each global spatial-angular point
    
    Ownership rule: The processor with the minimum rank that has this point owns it
    This ensures deterministic, consistent ownership across all processors
    """
    rank = MPI.Comm_rank(comm)
    
    # Each processor claims ownership of its points
    local_ownership = fill(typemax(Int), gnpoin)
    
    for ip = 1:n_spa
        gid = ip2gip_spa[ip]
        if gid > 0 && gid <= gnpoin
            local_ownership[gid] = rank
        end
    end
    
    # Global reduction: take minimum rank (lowest rank wins)
    gip2owner_spa = similar(local_ownership)
    MPI.Allreduce!(local_ownership, gip2owner_spa, MPI.MIN, comm)
    
    # Verify: all points should have an owner
    unowned = findall(x -> x == typemax(Int), gip2owner_spa)
    if !isempty(unowned) && rank == 0
        @warn "Found $(length(unowned)) unowned spatial-angular points!"
    end
    
    return gip2owner_spa
end

function verify_compact_numbering(ip2gip_spa, gnpoin, rank, comm)
    """
    Verify that numbering is compact (1:gnpoin with no gaps)
    """
    local_present = falses(gnpoin)
    for gid in ip2gip_spa
        if gid > 0 && gid <= gnpoin
            local_present[gid] = true
        end
    end
    
    
    local_present_bool = convert(Vector{Bool}, local_present)
    global_present_bool = MPI.Allreduce(local_present_bool, MPI.LOR, comm)
    global_present = BitVector(global_present_bool)
    
    has_gaps = false
    missing_ids = Int[]
    for i = 1:gnpoin
        if !global_present[i]
            push!(missing_ids, i)
            has_gaps = true
        end
    end
    
    if has_gaps
        if rank == 0
            @warn "Gaps found in global numbering at IDs: $(missing_ids[1:min(10, length(missing_ids))])..."
        end
    else
        if rank == 0
            @info "✓ Global numbering is compact: 1:$gnpoin with no gaps"
        end
    end
    
    return !has_gaps
end

function setup_assembler(SD, a, index_a, owner_a)

    if SD == NSD_1D() return nothing end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)

    if rank_sz == 1 return nothing end

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


    sum_array_1D = zeros(Float64, global_max_index)
    sum_array_2D = zeros(Float64, global_max_index, m)

    send_data_sizes = [send_idx_sizes[i+1] * m for i in 0:rank_sz-1]
    recv_data_sizes = [recv_idx_sizes[i+1] * m for i in 0:rank_sz-1]
    send_data_buffers = [zeros(Float64, send_idx_sizes[i+1] * m) for i in 0:rank_sz-1]
    recv_data_buffers = [zeros(Float64, recv_idx_sizes[i+1] * m) for i in 0:rank_sz-1]

    cache = AssemblerCache(global_max_index, index_a, owner_a, recv_idx_buffers,
            recvback_idx_buffers, sum_array_1D, sum_array_2D,
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




    @inbounds for (i, idx) in enumerate(cache.index_a)
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
                if is1D
                    cache.sum_array_1D[idx] += buffer[i]
                else
                    for j = 1:m
                        cache.sum_array_2D[idx, j] += buffer[(i-1)*m+j]
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
                if is1D
                    buf_data[j] = cache.sum_array_1D[idx]
                else
                    for k = 1:m
                        buf_data[(j-1)*m+k] = cache.sum_array_2D[idx, k]
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
                if is1D
                    cache.sum_array_1D[idx] = buffer[i]
                else
                    for j = 1:m
                        cache.sum_array_2D[idx, j] = buffer[(i-1)*m+j]
                    end
                end
            end
        end
    end


    @inbounds for (i, idx) in enumerate(cache.index_a)
        if is1D
            a[i] = cache.sum_array_1D[idx]
        else
            for j = 1:m
                a[i, j] = cache.sum_array_2D[idx, j]
            end
        end
    end
end

function assemble_mpi_global_vec!(a, cache::AssemblerCache)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    T = eltype(a)

    n = size(a, 1)

    fill!(cache.sum_array_1D, zero(T))




    @inbounds for (i, idx) in enumerate(cache.index_a)
        owner = cache.owner_a[i]
        if owner == rank
            cache.sum_array_1D[idx] += a[i]
        end
    end
    
    for i in 0:rank_sz-1
        fill!(cache.send_data_buffers[i+1], zero(T))
    end

    @inbounds for owner = 0:rank_sz-1
        buf_data = cache.send_data_buffers[owner+1]
        send_i_local = cache.send_i[owner+1]
        for (i,idx) in enumerate(send_i_local)
            buf_data[i] = a[idx]
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
                cache.sum_array_1D[idx] += buffer[i]
            end
        end
    end

    # send data back to original ranks
    sendback_data_buffers = cache.recv_data_buffers
    @inbounds for i in 0:rank_sz-1
        if cache.recv_data_sizes[i+1] > 0
            buf_data = sendback_data_buffers[i+1]
            for (j, idx) in enumerate(cache.recv_idx_buffers[i+1])
                buf_data[j] = cache.sum_array_1D[idx]
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
                cache.sum_array_1D[idx] = buffer[i]
            end
        end
    end
    
    @inbounds for (i, idx) in enumerate(cache.index_a)
        a[i] = cache.sum_array_1D[idx]
    end
end

