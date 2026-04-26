module JeGeometry
#__precompile__(false)
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
    n_spa, n_non_global_nodes, nc_non_global_nodes,
)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nproc = MPI.Comm_size(comm)
    nelem = mesh.nelem
    ngl   = mesh.ngl

    # ── Phase 1: Build local unique signatures separated by type ─────────────
    hanging_node_set = Set{Int}(nc_non_global_nodes)

    local_free_sigs    = Dict{NTuple{3,Float64}, Int}()  # sig → local ip_spa
    local_hanging_sigs = Dict{NTuple{3,Float64}, Int}()

    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip  = mesh.connijk[iel, i, j, k]
            gip = ip2gip[ip]
            for e_ext = 1:extra_meshes_extra_nelems[iel]
                nop = extra_meshes_extra_nops[iel][e_ext]
                for jθ = 1:nop+1, iθ = 1:nop+1
                    ip_ext = extra_meshes_connijk[iel][e_ext, iθ, jθ]
                    θ      = extra_meshes_coords[iel][1, ip_ext]
                    ϕ      = extra_meshes_coords[iel][2, ip_ext]
                    ip_spa = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                    sig    = (Float64(gip), round(θ, digits=12), round(ϕ, digits=12))

                    if ip_spa in hanging_node_set
                        get!(local_hanging_sigs, sig, ip_spa)
                    else
                        get!(local_free_sigs, sig, ip_spa)
                    end
                end
            end
        end
    end

    # Sort for deterministic prefix-sum ordering
    free_sig_list    = sort!(collect(keys(local_free_sigs)))
    hanging_sig_list = sort!(collect(keys(local_hanging_sigs)))
    n_local_free    = length(free_sig_list)
    n_local_hanging = length(hanging_sig_list)

    # ── Phase 2: Identify processor-boundary spatial nodes ───────────────────
    # Build owned and non-owned sets from local mesh data only — no communication
    owned_spatial    = Set{Int}()
    nonowned_spatial = Set{Int}()
    for ip = 1:mesh.npoin
        gip = ip2gip[ip]
        if gip2owner[ip] == rank
            push!(owned_spatial, gip)
        else
            push!(nonowned_spatial, gip)
        end
    end

    # Exchange non-owned sets to find what each rank shares with others.
    # Send our non-owned list (= nodes we have but others own).
    # After Allgatherv, find which of our owned nodes appear as non-owned elsewhere.
    local_nonowned     = collect(nonowned_spatial)
    n_local_no         = Int32(length(local_nonowned))
    no_counts          = MPI.Allgather([n_local_no], comm)
    all_nonowned       = MPI.Allgatherv(local_nonowned, no_counts, comm)

    processor_boundary_spatial = copy(nonowned_spatial)
    all_nonowned_set   = Set(all_nonowned)
    for gip in owned_spatial
        if gip in all_nonowned_set
            push!(processor_boundary_spatial, gip)
        end
    end

    # ── Phase 3: Boundary free signatures — O(1) lookup via Dict ─────────────
    # Replaced findfirst O(n) scan with direct Dict membership test
    boundary_free_sigs = Set{NTuple{3,Float64}}()
    for sig in free_sig_list
        gip_spatial = Int(sig[1])
        if gip_spatial in processor_boundary_spatial
            push!(boundary_free_sigs, sig)
        end
    end

    # ── Phase 4: Exchange boundary signatures — flat Float64, no serialization
    local_bsig_data = Vector{Float64}(undef, 3 * length(boundary_free_sigs))
    for (k, sig) in enumerate(boundary_free_sigs)
        local_bsig_data[3k-2] = sig[1]
        local_bsig_data[3k-1] = sig[2]
        local_bsig_data[3k]   = sig[3]
    end
    n_local_bsig  = Int32(length(boundary_free_sigs))
    bsig_counts   = MPI.Allgather([n_local_bsig], comm)
    all_bsig_data = MPI.Allgatherv(local_bsig_data, Int32.(bsig_counts .* 3), comm)

    # Count occurrences to find truly shared signatures
    sig_count = Dict{NTuple{3,Float64}, Int}()
    for k = 1:3:length(all_bsig_data)
        sig = (all_bsig_data[k], all_bsig_data[k+1], all_bsig_data[k+2])
        sig_count[sig] = get(sig_count, sig, 0) + 1
    end
    globally_shared_sigs = Set(sig for (sig, c) in sig_count if c > 1)

    # ── Phase 5: Assign tentative GIDs via prefix sum — free nodes ───────────
    offset_free = MPI.Exscan(n_local_free, MPI.SUM, comm)
    offset_free = (rank == 0) ? 0 : offset_free

    sig_to_tentative = Dict{NTuple{3,Float64}, Int}()
    sizehint!(sig_to_tentative, n_local_free)
    for (idx, sig) in enumerate(free_sig_list)
        sig_to_tentative[sig] = offset_free + idx
    end

    # ── Phase 6: Resolve shared signatures — flat exchange, no serialization ──
    # Send (sig_f1, sig_f2, sig_f3, tentative_gid) for each shared sig we own
    local_shared = NTuple{3,Float64}[]
    local_gids   = Int[]
    for sig in globally_shared_sigs
        if haskey(sig_to_tentative, sig)
            push!(local_shared, sig)
            push!(local_gids, sig_to_tentative[sig])
        end
    end

    n_shared = Int32(length(local_shared))
    # Pack as Float64: [f1, f2, f3, gid_as_float64] per entry — 4 values each
    packed = Vector{Float64}(undef, 4 * n_shared)
    for (k, (sig, gid)) in enumerate(zip(local_shared, local_gids))
        packed[4k-3] = sig[1]
        packed[4k-2] = sig[2]
        packed[4k-1] = sig[3]
        packed[4k]   = Float64(gid)
    end
    shared_counts  = MPI.Allgather([n_shared], comm)
    all_packed     = MPI.Allgatherv(packed, Int32.(shared_counts .* 4), comm)

    # Resolve: minimum tentative GID wins
    shared_resolution = Dict{NTuple{3,Float64}, Int}()
    for k = 1:4:length(all_packed)
        sig = (all_packed[k], all_packed[k+1], all_packed[k+2])
        gid = Int(all_packed[k+3])
        if haskey(shared_resolution, sig)
            shared_resolution[sig] = min(shared_resolution[sig], gid)
        else
            shared_resolution[sig] = gid
        end
    end

    # ── Phase 7: Final GIDs for free nodes ───────────────────────────────────
    sig_to_final_free = Dict{NTuple{3,Float64}, Int}()
    sizehint!(sig_to_final_free, n_local_free)
    for sig in free_sig_list
        sig_to_final_free[sig] = get(shared_resolution, sig, sig_to_tentative[sig])
    end

    # ── Phase 8: Compact free numbering without global gather ────────────────
    # Key insight: instead of gathering all GIDs globally, use a second prefix sum.
    # Count how many of our free signatures are the canonical owner
    # (i.e., our tentative GID == resolved GID, meaning we assigned the minimum).
    # This tells us how many unique free DOFs we contribute to the global count.
    #
    # For shared nodes: only the rank whose tentative == resolved GID is canonical.
    # For interior nodes: always canonical.

    n_canonical_free = 0
    for sig in free_sig_list
        tentative = sig_to_tentative[sig]
        resolved  = sig_to_final_free[sig]
        if tentative == resolved
            n_canonical_free += 1
        end
    end

    # Prefix sum gives each rank its offset into the canonical free numbering
    canonical_offset = MPI.Exscan(n_canonical_free, MPI.SUM, comm)
    canonical_offset = (rank == 0) ? 0 : canonical_offset

    # Assign compact GIDs to canonical nodes; share resolved GIDs for shared nodes
    # via a second exchange so non-canonical ranks get the final compact GID
    canonical_sig_to_compact = Dict{NTuple{3,Float64}, Int}()
    canon_idx = canonical_offset
    for sig in free_sig_list
        tentative = sig_to_tentative[sig]
        resolved  = sig_to_final_free[sig]
        if tentative == resolved
            canon_idx += 1
            canonical_sig_to_compact[sig] = canon_idx
        end
    end

    gnpoin_free = MPI.Allreduce(n_canonical_free, MPI.SUM, comm)

    # Exchange compact GIDs for shared signatures
    # Pack: [f1, f2, f3, compact_gid] for each sig we are canonical for
    n_canon = Int32(length(canonical_sig_to_compact))
    canon_packed = Vector{Float64}(undef, 4 * n_canon)
    for (k, (sig, cgid)) in enumerate(canonical_sig_to_compact)
        canon_packed[4k-3] = sig[1]
        canon_packed[4k-2] = sig[2]
        canon_packed[4k-1] = sig[3]
        canon_packed[4k]   = Float64(cgid)
    end
    canon_counts    = MPI.Allgather([n_canon], comm)
    all_canon_packed = MPI.Allgatherv(canon_packed, Int32.(canon_counts .* 4), comm)

    # Build final compact GID lookup for all free signatures
    sig_to_compact_gid = Dict{NTuple{3,Float64}, Int}()
    sizehint!(sig_to_compact_gid, n_local_free)
    # First add our own canonical entries
    merge!(sig_to_compact_gid, canonical_sig_to_compact)
    # Then fill in compact GIDs received from other ranks for shared non-canonical sigs
    for k = 1:4:length(all_canon_packed)
        sig  = (all_canon_packed[k], all_canon_packed[k+1], all_canon_packed[k+2])
        cgid = Int(all_canon_packed[k+3])
        if !haskey(sig_to_compact_gid, sig) && haskey(sig_to_final_free, sig)
            sig_to_compact_gid[sig] = cgid
        end
    end

    # ── Phase 9: Hanging node GIDs via prefix sum — no global gather ──────────
    offset_hanging = MPI.Exscan(n_local_hanging, MPI.SUM, comm)
    offset_hanging = (rank == 0) ? 0 : offset_hanging

    sig_to_hanging_gid = Dict{NTuple{3,Float64}, Int}()
    sizehint!(sig_to_hanging_gid, n_local_hanging)
    for (idx, sig) in enumerate(hanging_sig_list)
        sig_to_hanging_gid[sig] = gnpoin_free + offset_hanging + idx
    end

    gnpoin_hanging = MPI.Allreduce(n_local_hanging, MPI.SUM, comm)
    gnpoin = gnpoin_free + gnpoin_hanging

    # ── Phase 10: Apply global numbering ─────────────────────────────────────
    ip2gip_spa = zeros(Int64, n_spa)

    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip  = mesh.connijk[iel, i, j, k]
            gip = ip2gip[ip]
            for e_ext = 1:extra_meshes_extra_nelems[iel]
                nop = extra_meshes_extra_nops[iel][e_ext]
                for jθ = 1:nop+1, iθ = 1:nop+1
                    ip_ext = extra_meshes_connijk[iel][e_ext, iθ, jθ]
                    θ      = extra_meshes_coords[iel][1, ip_ext]
                    ϕ      = extra_meshes_coords[iel][2, ip_ext]
                    ip_spa = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                    sig    = (Float64(gip), round(θ, digits=12), round(ϕ, digits=12))

                    if ip_spa in hanging_node_set
                        ip2gip_spa[ip_spa] = sig_to_hanging_gid[sig]
                    else
                        ip2gip_spa[ip_spa] = sig_to_compact_gid[sig]
                    end
                end
            end
        end
    end

    # ── Phase 11: Owner map — local computation only ─────────────────────────
    # gip2owner_spa: for each global ID, which rank owns it.
    # A rank owns a free DOF if it was canonical (tentative == resolved).
    # A rank owns a hanging DOF always (hanging nodes are local by construction).
    # Build by exchanging (gid, owner_rank) pairs for all DOFs this rank owns.

    n_local_dofs = Int32(n_spa)
    dof_counts   = MPI.Allgather([n_local_dofs], comm)

    # Pack (gid, rank) for all local DOFs
    local_claims = Vector{Int}(undef, 2 * n_spa)
    for ip = 1:n_spa
        local_claims[2ip-1] = ip2gip_spa[ip]
        local_claims[2ip]   = rank
    end

    all_claims = MPI.Allgatherv(local_claims, Int32.(dof_counts .* 2), comm)

    # Resolve: minimum rank wins for each gid
    gip2owner_spa = fill(typemax(Int), gnpoin)
    for k = 1:2:length(all_claims)
        gid  = all_claims[k]
        ownr = all_claims[k+1]
        if ownr < gip2owner_spa[gid]
            gip2owner_spa[gid] = ownr
        end
    end

# Verify no unowned points in debug mode
@assert all(x -> x != typemax(Int), gip2owner_spa) "Unowned global DOFs detected"

gip2ip = zeros(Int, gnpoin)
for (ip, gid) in enumerate(ip2gip_spa)
    if gid > 0 && gip2owner_spa[gid] == rank
        gip2ip[gid] = ip
    end
end

    return ip2gip_spa, gip2ip, gip2owner_spa, gnpoin
end

function verify_hanging_node_numbering_fast(ip2gip_spa, n_spa, gnpoin_free, gnpoin,
                                             hanging_node_set, rank)
    # Local checks only — no MPI communication
    free_gids    = Set{Int}()
    hanging_gids = Set{Int}()
    ok = true

    for ip = 1:n_spa
        gid = ip2gip_spa[ip]
        if ip in hanging_node_set
            if gid <= gnpoin_free || gid > gnpoin
                @warn "[Rank $rank] Hanging node $ip has GID $gid outside [$gnpoin_free+1, $gnpoin]"
                ok = false
            end
            if gid in hanging_gids
                @warn "[Rank $rank] Duplicate hanging GID $gid"
                ok = false
            end
            push!(hanging_gids, gid)
        else
            if gid < 1 || gid > gnpoin_free
                @warn "[Rank $rank] Free node $ip has GID $gid outside [1, $gnpoin_free]"
                ok = false
            end
            if gid in free_gids
                @warn "[Rank $rank] Duplicate free GID $gid on rank $rank"
                ok = false
            end
            push!(free_gids, gid)
        end
    end

    if ok && rank == 0
        @info "Numbering verification passed (local checks)"
    end
end

