using .JeGeometry
using Gridap
using Gridap.Arrays
using Gridap.Arrays: Table
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Geometry: GridMock
using GridapDistributed
using PartitionedArrays
using GridapGmsh
using GridapP4est
using Metis
using Test

export St_mesh
export mod_mesh_mesh_driver
export mod_mesh_build_mesh!
export mod_mesh_read_gmsh!

const VERTEX_NODES = UInt64(1)
const EDGE_NODES   = UInt64(2)
const FACE_NODES   = UInt64(4)

Base.@kwdef mutable struct St_mesh{TInt, TFloat, backend}

    x = KernelAbstractions.zeros(backend, TFloat, 2)
    y = KernelAbstractions.zeros(backend, TFloat, 2)
    z = KernelAbstractions.zeros(backend, TFloat, 2)

    x_ho = KernelAbstractions.zeros(backend, TFloat, 2)
    y_ho = KernelAbstractions.zeros(backend, TFloat, 2)
    z_ho = KernelAbstractions.zeros(backend, TFloat, 2)

    Δx = KernelAbstractions.zeros(backend, TFloat, 2)
    Δy = KernelAbstractions.zeros(backend, TFloat, 2)
    Δz = KernelAbstractions.zeros(backend, TFloat, 2)
    
    xmin::Union{TFloat, Missing} = -1.0;
    xmax::Union{TFloat, Missing} = +1.0;
    
    ymin::Union{TFloat, Missing} = -1.0;
    ymax::Union{TFloat, Missing} = +1.0;

    zmin::Union{TFloat, Missing} = -1.0;
    zmax::Union{TFloat, Missing} = +1.0;
    
    npx::Union{TInt, Missing} = 1
    npy::Union{TInt, Missing} = 1
    npz::Union{TInt, Missing} = 1
    
    nelem::Union{TInt, Missing} = 1
    nelem_semi_inf::Union{TInt, Missing} = 0# Semi infinite elements for Laguerre BC
    nelem_int::Union{TInt, Missing} = 1    # internal elements
    npoin::Union{TInt, Missing} = 1        # This is updated after populating with high-order nodes
    npoin_original::Union{TInt, Missing} =1# Storage for original npoin if modified for Laguerre semi_inf
    npoin_linear::Union{TInt, Missing} = 1 # This is always the original number of the first-order grid
    nelem_bdy::Union{TInt, Missing} = 1    # bdy elements
    
    nedges::Union{TInt, Missing} = 1       # total number of edges
    nedges_bdy::Union{TInt, Missing} = 1   # bdy edges
    nedges_int::Union{TInt, Missing} = 1   # internal edges
    
    nfaces::Union{TInt, Missing} = 1       # total number of faces
    nfaces_bdy::Union{TInt, Missing} = 1   # bdy faces
    nfaces_int::Union{TInt, Missing} = 1   # internal faces

    # global for MPI
    gnelem::Union{TInt, Missing} = 1
    gnpoin::Union{TInt, Missing} = 1        # This is updated after populating with high-order nodes
    gnpoin_linear::Union{TInt, Missing} = 1 # This is always the original number of the first-order grid
    gnedges::Union{TInt, Missing} = 1       # total number of edges
    gnfaces::Union{TInt, Missing} = 1       # total number of faces
    
    nsd::Union{TInt, Missing} = 1              # number of space dim
    nop::Union{TInt, Missing} = 4              # poly order
    ngl::Union{TInt, Missing} = nop + 1        # number of quad point 
    ngr::Union{TInt, Missing} = 0              # nop_gr
    lLaguerre::Union{Bool, Missing} = false # whether or not Laguerre boundaries are in the mesh
    npoin_el::Union{TInt, Missing} = 1         # Total number of points in the reference element
    

    NNODES_EL::Union{TInt, Missing}  =  2^nsd
    NEDGES_EL::Union{TInt, Missing}  = 12
    NFACES_EL::Union{TInt, Missing}  =  6
    EDGE_NODES::Union{TInt, Missing} =  2
    FACE_NODES::Union{TInt, Missing} =  4
    
    #low and high order connectivity tables
    cell_node_ids::Table{Int64,Vector{Int64},Vector{Int64}}    = Gridap.Arrays.Table(zeros(nelem), zeros(1))
    cell_node_ids_ho::Table{Int64,Vector{Int64},Vector{Int64}} = Gridap.Arrays.Table(zeros(nelem), zeros(1))
    cell_edge_ids::Table{Int64,Vector{Int64},Vector{Int64}}    = Gridap.Arrays.Table(zeros(nelem), zeros(1))    
    cell_face_ids::Table{Int64,Vector{Int64},Vector{Int64}}    = Gridap.Arrays.Table(zeros(nelem), zeros(1))
    face_edge_ids::Table{Int64,Vector{Int64},Vector{Int64}}    = Gridap.Arrays.Table(zeros(nfaces), zeros(1))
    facet_cell_ids::Table{Int64,Vector{Int64},Vector{Int64}}   = Gridap.Arrays.Table(zeros(nfaces), zeros(1))

    connijk_lag = KernelAbstractions.zeros(backend,TInt, 0, 0, 0, 0)
    connijk =  KernelAbstractions.zeros(backend,TInt, 0, 0, 0, 0)
    conn_edgesijk::Array{Int64,2} = KernelAbstractions.zeros(backend, TInt, 0, 0)    # edge analogue of connijk
    conn_facesijk::Array{Int64,2} = KernelAbstractions.zeros(backend, TInt, 0, 0)    # face analogue of connijk

    el_min = KernelAbstractions.zeros(backend,TFloat, 0, 0)
    el_max = KernelAbstractions.zeros(backend,TFloat, 0, 0)

    conn::Array{TInt,2}  = KernelAbstractions.zeros(backend, TInt, 0, 0)
    conn_unique_edges    = Array{TInt}(undef,  1, 2)
    conn_unique_edges1   = Array{Int64}(undef,  1, 2)
    conn_unique_faces    = Array{TInt}(undef,  1, 4)
    poin_in_edge         = Array{TInt}(undef, 0, 0)
    conn_edge_el         = Array{TInt}(undef, 0, 0, 0)
    poin_in_face         = Array{TInt}(undef, 0, 0, 0)
    conn_face_el         = Array{TInt}(undef, 0, 0, 0)
    face_in_elem         = Array{TInt}(undef, 0, 0, 0)

    edge_g_color::Array{Int64, 1} = zeros(Int64, 1)
    
    #MPI variables
    rank                = 0
    ip2gip              = KernelAbstractions.zeros(backend, TInt, 0)
    gip2owner           = KernelAbstractions.zeros(backend, TInt, 0)
    gip2ip              = KernelAbstractions.zeros(backend, TInt, 0)
    parts               = 1
    nparts              = 1
    
    #Auxiliary arrays for boundary conditions
    
    bdy_edge_in_elem  =  KernelAbstractions.zeros(backend, TInt, 0)
    poin_in_bdy_edge  =  KernelAbstractions.zeros(backend, TInt, 0, 0)
    bdy_face_in_elem  =  KernelAbstractions.zeros(backend, TInt, 0)
    poin_in_bdy_face  =  KernelAbstractions.zeros(backend, TInt, 0, 0, 0)
    edge_type     = Array{Union{Nothing, String}}(nothing, 1)
    face_type     = Array{Union{Nothing, String}}(nothing, 1)
    bdy_edge_type = Array{Union{Nothing, String}}(nothing, 1)
    bdy_face_type = Array{Union{Nothing, String}}(nothing, 1)
    bdy_edge_type_id  =  KernelAbstractions.zeros(backend, TInt, 0)

    Δelem        = KernelAbstractions.zeros(backend, TInt, 0)
    Δelem_s      = 0.0
    Δelem_l      = 0.0
    Δeffective_s = 0.0
    Δeffective_l = 0.0
        
    SD::AbstractSpaceDimensions

    # for AMR
    ad_lvl = KernelAbstractions.zeros(backend, TInt, 0)
    num_hanging_facets::Union{TInt, Missing} = 0
    non_conforming_facets                = Vector{Vector{TInt}}(undef,num_hanging_facets)
    non_conforming_facets_parents_ghost  = Vector{Vector{TInt}}(undef,num_hanging_facets)
    non_conforming_facets_children_ghost = Vector{Vector{TInt}}(undef,num_hanging_facets)

    pgip_ghost = KernelAbstractions.zeros(backend, TInt, 0)
    pgip_owner = KernelAbstractions.zeros(backend, TInt, 0)
    cgip_ghost = KernelAbstractions.zeros(backend, TInt, 0)
    cgip_owner = KernelAbstractions.zeros(backend, TInt, 0)

    msg_suppress::Bool = false

end

const get_d_to_face_to_parent_face = Gridap.Adaptivity.get_d_to_face_to_parent_face


function mod_mesh_read_gmsh!(mesh::St_mesh, inputs::Dict, nparts, distribute, adapt_flags = nothing, partitioned_model_coarse = nothing, omesh = nothing)

    # determine backend
    backend = CPU()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    mpi_size = MPI.Comm_size(comm)
    
    #
    # Read GMSH grid from file
    #      
    parts           = distribute(LinearIndices((nparts,)))
    mesh.parts      = distribute(LinearIndices((nparts,)))
    mesh.nparts     = nparts
    mesh.rank       = rank
    ladaptive       = inputs[:ladapt]
    linitial_refine = inputs[:linitial_refine]
    lamr_mesh       = !isnothing(adapt_flags)
    if isnothing(adapt_flags)
    
        if ladaptive == false && linitial_refine == false
            partitioned_model = GmshDiscreteModel(parts, inputs[:gmsh_filename], renumber=true)
            model = local_views(partitioned_model).item_ref[]
        elseif linitial_refine == true
            gmodel = GmshDiscreteModel(inputs[:gmsh_filename], renumber=true)
            partitioned_model = UniformlyRefinedForestOfOctreesDiscreteModel(parts, gmodel, inputs[:init_refine_lvl])
            cell_gids = local_views(partition(get_cell_gids(partitioned_model))).item_ref[]
            dmodel = local_views(partitioned_model).item_ref[]
            model  = DiscreteModelPortion(dmodel, own_to_local(cell_gids))
        elseif ladaptive == true && linitial_refine == false
            gmodel = GmshDiscreteModel(inputs[:gmsh_filename], renumber=true)
            partitioned_model_coarse = OctreeDistributedDiscreteModel(parts,gmodel)

            ref_coarse_flags = map(parts,partition(get_cell_gids(partitioned_model_coarse.dmodel))) do rank,indices
                flags = zeros(Cint,length(indices))
                flags.=nothing_flag
                # @info flags
                # flags[250] = refine_flag
                # if rank == 2
                    # flags[1:3:end] .= refine_flag
                    # flags[1] = refine_flag
                # end
                flags
            end
            partitioned_model, glue_adapt=Gridap.Adaptivity.adapt(partitioned_model_coarse,ref_coarse_flags);
            # partitioned_model, glue_redistribute = redistribute(partitioned_model)
            # glue_adapt = get_adaptivity_glue(partitioned_model.dmodel)
            cell_gids = local_views(partition(get_cell_gids(partitioned_model))).item_ref[]
            dmodel = local_views(partitioned_model.dmodel.models).item_ref[]
            cmodel = local_views(partitioned_model_coarse.dmodel.models).item_ref[]
            cell_gids_c = local_views(partition(get_cell_gids(partitioned_model_coarse))).item_ref[]
            # @info rank, own_to_local(cell_gids), local_to_own(cell_gids), local_to_global(cell_gids)
            model  = DiscreteModelPortion(dmodel, own_to_local(cell_gids))
            model_coarse  = DiscreteModelPortion(cmodel, own_to_local(cell_gids_c))
            dtopology      = get_grid_topology(dmodel)
        end
    else
        ref_coarse_flags = map(parts,partition(get_cell_gids(partitioned_model_coarse.dmodel))) do rank,indices
            flags  = zeros(Cint,length(indices))
            flags .= nothing_flag
            flags[1:length(adapt_flags)] = adapt_flags
            flags
        end

        partitioned_model,glue_adapt=Gridap.Adaptivity.adapt(partitioned_model_coarse,ref_coarse_flags);
        cell_gids = local_views(partition(get_cell_gids(partitioned_model))).item_ref[]
        dmodel = local_views(partitioned_model.dmodel.models).item_ref[]
        cmodel = local_views(partitioned_model_coarse.dmodel.models).item_ref[]
        cell_gids_c = local_views(partition(get_cell_gids(partitioned_model_coarse))).item_ref[]
        # @info rank, own_to_local(cell_gids), local_to_own(cell_gids), local_to_global(cell_gids)
        model  = DiscreteModelPortion(dmodel, own_to_local(cell_gids))
        model_coarse  = DiscreteModelPortion(cmodel, own_to_local(cell_gids_c))
        dtopology      = get_grid_topology(dmodel)

        mesh.msg_suppress = true
    end

    topology      = get_grid_topology(model)
    mesh.nsd      = num_cell_dims(model)
    
    POIN_flg = 0
    EDGE_flg = 1
    FACE_flg = 2
    ELEM_flg = 3
    
    if mesh.nsd == 3
        mesh.SD = NSD_3D()

        mesh.NNODES_EL  = 8
        mesh.NEDGES_EL  = 12
        mesh.NFACES_EL  = 6
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 4
    elseif mesh.nsd == 2

        mesh.SD = NSD_2D()
        ELEM_flg = FACE_flg
        
        mesh.NNODES_EL  = 4
        mesh.NEDGES_EL  = 4
        mesh.NFACES_EL  = 1
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 4
    elseif mesh.nsd == 1
        mesh.SD = NSD_1D()

        mesh.NNODES_EL  = 2
        mesh.NEDGES_EL  = 1
        mesh.NFACES_EL  = 0
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 0
    else
        error( " WRONG NSD: This is not theoretical physics: we only handle 1, 2, or 3 dimensions!")
    end
    


    d_to_num_dfaces = [num_vertices(model), num_edges(model), num_cells(model)]

    # Write the partitioned model to a VTK file
    # vtk_directory = "./coarse/" 
    # writevtk(partitioned_model_coarse, vtk_directory)
    vtk_directory = "./refine/"
    if ladaptive == true
        writevtk(partitioned_model.dmodel, vtk_directory)
    end



    #dump(topology)
    #
    # Mesh elements, nodes, faces, edges
    #
    if (ladaptive == true) || (linitial_refine == true) 
        p2pp   = Geometry.get_face_to_parent_face(model,POIN_flg)
        eg2peg = Geometry.get_face_to_parent_face(model,EDGE_flg)
        f2pf   = Geometry.get_face_to_parent_face(model,FACE_flg)
        e2pe   = Geometry.get_face_to_parent_face(model,ELEM_flg)
        # @info point2ppoint, edge2pedge, face2pface

        pgids  = local_views(partition(get_face_gids(partitioned_model,POIN_flg))).item_ref[]
        edgids = local_views(partition(get_face_gids(partitioned_model,EDGE_flg))).item_ref[]
        fgids  = local_views(partition(get_face_gids(partitioned_model,FACE_flg))).item_ref[]
        elgids = local_views(partition(get_face_gids(partitioned_model,ELEM_flg))).item_ref[]
        # @info rank, own_to_local(pgids), "a ", local_to_own(pgids),  "a ", local_to_global(pgids), "a ", p2pp
        point2ppoint = local_to_global(pgids)[p2pp]
        edge2pedge   = local_to_global(edgids)[eg2peg]
        face2pface   = local_to_global(fgids)[f2pf]
        elm2pelm     = local_to_global(elgids)[e2pe]
    else
        point2ppoint = Geometry.get_face_to_parent_face(model,POIN_flg)
        edge2pedge   = Geometry.get_face_to_parent_face(model,EDGE_flg)
        face2pface   = Geometry.get_face_to_parent_face(model,FACE_flg)
        elm2pelm    = Geometry.get_face_to_parent_face(model,ELEM_flg)
    end
    # @info rank, p2pp, point2ppoint
    if ladaptive == true
        hanging_vert_glue  = local_views(partitioned_model.non_conforming_glue).item_ref[].hanging_faces_glue[1]
        hanging_facet_glue = local_views(partitioned_model.non_conforming_glue).item_ref[].hanging_faces_glue[mesh.nsd]
        num_regular_facets = local_views(partitioned_model.non_conforming_glue).item_ref[].num_regular_faces[mesh.nsd]
        num_hanging_facets = local_views(partitioned_model.non_conforming_glue).item_ref[].num_hanging_faces[mesh.nsd]
    else
        num_hanging_facets = 0
    end
    # @info rank, hanging_vert_glue, local_to_global(elgids), "a", hanging_facet_glue,num_hanging_facets


    mesh.gnpoin_linear = num_faces(partitioned_model,POIN_flg)    
    mesh.gnpoin        = mesh.gnpoin_linear         #This will be updated for the high order grid
    mesh.gnedges       = num_faces(partitioned_model,EDGE_flg)
    mesh.gnfaces       = num_faces(partitioned_model,FACE_flg)   
    mesh.gnelem        = num_faces(partitioned_model,ELEM_flg)


    mesh.npoin_linear = num_faces(model,POIN_flg)    
    mesh.npoin        = mesh.npoin_linear         #This will be updated for the high order grid
    mesh.nedges       = num_faces(model,EDGE_flg)
    mesh.nfaces       = num_faces(model,FACE_flg)   
    mesh.nelem        = num_faces(model,ELEM_flg)


    
    # if (ladaptive == 1)
        mesh.nelem_bdy    = length(get_boundary_cells(model,mesh.nsd))
        mesh.nfaces_bdy   = length(get_boundary_faces(model,mesh.nsd,FACE_flg))
        mesh.nedges_bdy   = length(get_boundary_faces(model,mesh.nsd,EDGE_flg))
    # else 
        # mesh.nelem_bdy    = count(get_isboundary_face(topology,mesh.nsd))
        # mesh.nfaces_bdy   = count(get_isboundary_face(topology,FACE_flg))
        # mesh.nedges_bdy   = count(get_isboundary_face(topology,EDGE_flg))
    # end
    mesh.nelem_int    = mesh.nelem - mesh.nelem_bdy
    mesh.nfaces_int   = mesh.nfaces - mesh.nfaces_bdy
    mesh.nedges_int   = mesh.nedges - mesh.nedges_bdy

    #get_isboundary_face(topology,mesh.nsd-1)
    if !lamr_mesh
        println_rank(" # GMSH LINEAR GRID PROPERTIES"; msg_rank = rank, suppress = mesh.msg_suppress)
        println_rank(" # N. Global points         : ", mesh.gnpoin_linear; msg_rank = rank, suppress = mesh.msg_suppress)
        println_rank(" # N. Global elements       : ", mesh.gnelem; msg_rank = rank, suppress = mesh.msg_suppress)
        println_rank(" # N. Global edges          : ", mesh.gnedges; msg_rank = rank, suppress = mesh.msg_suppress)
        println_rank(" # N. Global faces          : ", mesh.gnfaces; msg_rank = rank, suppress = mesh.msg_suppress)
        MPI.Barrier(comm)
        if mesh.msg_suppress == false
            for i = 0 : mpi_size
                if i == rank
                    println("   # Rank                       : ", rank)
                    println("     # N. points                : ", mesh.npoin_linear)
                    println("     # N. elements              : ", mesh.nelem)
                    println("     # N. edges                 : ", mesh.nedges)
                    println("     # N. faces                 : ", mesh.nfaces)    
                    println("     # N. internal elem         : ", mesh.nelem_int)
                    println("     # N. internal edges        : ", mesh.nedges_int) 
                    println("     # N. internal faces        : ", mesh.nfaces_int)    
                    println("     # N. boundary elem         : ", mesh.nelem_bdy)
                    println("     # N. boundary edges        : ", mesh.nedges_bdy)
                    println("     # N. boundary faces        : ", mesh.nfaces_bdy)
                end
                MPI.Barrier(comm)
            end
        end
        println_rank(" # GMSH LINEAR GRID PROPERTIES ...................... END"; msg_rank = rank, suppress = mesh.msg_suppress)
    end

    ngl                     = mesh.nop + 1
    tot_linear_poin         = mesh.npoin_linear
    
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    tot_faces_internal_nodes = mesh.nfaces*(ngl-2)*(ngl-2)
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)^(mesh.nsd)
    
    el_edges_internal_nodes = mesh.NEDGES_EL*(ngl-2)
    el_faces_internal_nodes = mesh.NFACES_EL*(ngl-2)*(ngl-2)
    el_vol_internal_nodes   = (ngl-2)^(mesh.nsd)
    
    #Update number of grid points from linear count to total high-order points
    mesh.npoin = tot_linear_poin + tot_edges_internal_nodes + tot_faces_internal_nodes + (mesh.nsd - 2)*tot_vol_internal_nodes
    
    if (mesh.nop > 1) && (!lamr_mesh)
        println_rank(" # GMSH HIGH-ORDER GRID PROPERTIES"; msg_rank = rank, suppress = mesh.msg_suppress)
        MPI.Barrier(comm)
        if mesh.msg_suppress == false
            for i = 0 : mpi_size
                if i == rank
                    println("   # Rank                         : ", rank)
                    println("     # N. edges internal points   : ", tot_edges_internal_nodes)
                    println("     # N. faces internal points   : ", tot_faces_internal_nodes)
                    println("     # N. volumes internal points : ", tot_vol_internal_nodes)
                    println("     # N. total high order points : ", mesh.npoin)
                end
                MPI.Barrier(comm)
            end
        end
        println_rank(" # GMSH HIGH-ORDER GRID PROPERTIES ...................... END"; msg_rank = rank, suppress = mesh.msg_suppress)
    end
    
    #
    # Resize as needed
    #
    mesh.x = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
    mesh.y = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
    mesh.z = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))

    mesh.ip2gip    = KernelAbstractions.zeros(backend, TInt, Int64(mesh.npoin))
    mesh.gip2owner = KernelAbstractions.ones(backend, TInt, Int64(mesh.npoin))*local_views(parts).item_ref[]
    
    
    mesh.conn_edge_el     = KernelAbstractions.zeros(backend, TInt, 2, Int64(mesh.NEDGES_EL), Int64(mesh.nelem))    
    mesh.conn_face_el     = KernelAbstractions.zeros(backend, TInt,  4, Int64(mesh.NFACES_EL), Int64(mesh.nelem))  
    mesh.bdy_edge_in_elem = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy))  
    mesh.poin_in_edge     = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges), Int64(mesh.ngl))
    mesh.poin_in_bdy_edge = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy), Int64(mesh.ngl))
    
    mesh.poin_in_face     = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nfaces), Int64(mesh.ngl), Int64(mesh.ngl))
    mesh.edge_type        = Array{Union{Nothing, String}}(nothing, Int64(mesh.nedges))
    mesh.bdy_edge_type    = Array{Union{Nothing, String}}(nothing, Int64(mesh.nedges_bdy))
    mesh.bdy_edge_type_id = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy))  
    
    if mesh.nsd > 2
        mesh.poin_in_bdy_face = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nfaces_bdy), Int64(mesh.ngl), Int64(mesh.ngl))
        mesh.face_type        = Array{Union{Nothing, String}}(nothing, Int64(mesh.nfaces))
        mesh.bdy_face_type    = Array{Union{Nothing, String}}(nothing, Int64(mesh.nfaces_bdy))
        mesh.bdy_face_in_elem = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nfaces_bdy))
    end
    mesh.npoin_el         = mesh.NNODES_EL + el_edges_internal_nodes + el_faces_internal_nodes + (mesh.nsd - 2)*el_vol_internal_nodes
    mesh.conn = KernelAbstractions.zeros(backend,TInt, Int64(mesh.nelem), Int64(mesh.npoin_el))
    mesh.el_max = KernelAbstractions.zeros(backend,TFloat, Int64(mesh.nelem), 3)
    mesh.el_min = KernelAbstractions.zeros(backend,TFloat, Int64(mesh.nelem), 3)
    #
    # Connectivity matrices
    #
    mesh.cell_node_ids     = get_cell_node_ids(get_grid(model))
    mesh.conn_unique_faces = get_face_nodes(model, FACE_flg) #faces --> 4 nodes
    mesh.conn_unique_edges = get_face_nodes(model, EDGE_flg) #edges --> 2 nodes

    mesh.cell_edge_ids     = get_faces(topology, mesh.nsd, 1) #edge map from local to global numbering i.e. iedge_g = cell_edge_ids[1:NELEM][1:NEDGES_EL]
    mesh.cell_face_ids     = get_faces(topology, mesh.nsd, mesh.nsd-1) #face map from local to global numbering i.e. iface_g = cell_face_ids[1:NELEM][1:NFACE_EL]
    mesh.face_edge_ids     = get_faces(topology,mesh.nsd-1, 1)
    mesh.facet_cell_ids    = get_faces(topology,mesh.nsd-1, mesh.nsd)
    # @info mesh.facet_cell_ids
    mesh.edge_g_color::Array{Int64, 1} = zeros(Int64, mesh.nedges)

    

    #
    # element refinement level
    #
    mesh.ad_lvl = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem))
    if ladaptive == true
        nelem_c        = num_faces(model_coarse,ELEM_flg)

        glue = local_views(glue_adapt).item_ref[]
        r_c_flags = local_views(ref_coarse_flags).item_ref[]
        if isnothing(adapt_flags)
            for i in 1:nelem_c
                elem_idx = glue.o2n_faces_map[i]
                for j in elem_idx
                    mesh.ad_lvl[j] = 0
                    if r_c_flags[i] == 1
                        mesh.ad_lvl[j] += 1
                    elseif r_c_flags[i] == 2
                        mesh.ad_lvl[j] -= 1
                    end
                end
            end
        else
            for i = 1:omesh.nelem
                elem_idx = glue.o2n_faces_map[i]
                for j in elem_idx
                    mesh.ad_lvl[j] = omesh.ad_lvl[i]
                    if r_c_flags[i] == 1
                        mesh.ad_lvl[j] += 1
                    elseif r_c_flags[i] == 2
                        mesh.ad_lvl[j] = mesh.ad_lvl[j] > 0 ? mesh.ad_lvl[j] - 1 : 0
                    end
                end
            end
        end
    end
    # @info rank, mesh.ad_lvl


    if (mesh.nsd == 1)
        nothing
    elseif (mesh.nsd == 2)
    
        mesh.connijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl),1)
    
        for iel = 1:mesh.nelem
            mesh.conn[iel, 1] = mesh.cell_node_ids[iel][1]
            mesh.conn[iel, 2] = mesh.cell_node_ids[iel][2]
            mesh.conn[iel, 3] = mesh.cell_node_ids[iel][4]
            mesh.conn[iel, 4] = mesh.cell_node_ids[iel][3]

            #
            # 1-----3
            # |     |
            # |     |
            # 2-----4
            #
            mesh.connijk[iel, 1,      1] = mesh.cell_node_ids[iel][2]
            mesh.connijk[iel, 1,    ngl] = mesh.cell_node_ids[iel][1]
            mesh.connijk[iel, ngl,  ngl] = mesh.cell_node_ids[iel][3]
            mesh.connijk[iel, ngl,    1] = mesh.cell_node_ids[iel][4]
            
            # @printf(" [1,1] [ngl, 1] [1, ngl] [ngl, ngl] %d %d %d %d\n", mesh.connijk[iel, 1, 1], mesh.connijk[iel, ngl, 1] , mesh.connijk[iel, 1,ngl], mesh.connijk[iel, ngl, ngl] )
        end
        #
        # Fill in elements dictionary needed by NodeOrdering.jl
        #
        elements = Dict(
            kk => mesh.conn[kk, 1:4]
            for kk = 1:mesh.nelem)
        element_types = Dict(
            kk => :Quad4
            for kk = 1:mesh.nelem)
        
        #
        # Rewrite coordinates in RCM order:
        #
        # filename = "./COORDS_LO_" + rank + ".dat" 
        # open("./COORDS_LO_$rank.dat", "w") do f
            for ip = 1:mesh.npoin_linear
                
                mesh.x[ip] = get_node_coordinates(get_grid(model))[ip][1]
                mesh.y[ip] = get_node_coordinates(get_grid(model))[ip][2]
                
                mesh.ip2gip[ip] = point2ppoint[ip]
                # mesh.gip2owner[ip] = 1
                # @printf(f, " %.6f %.6f 0.000000 %d %d\n", mesh.x[ip],  mesh.y[ip], ip, point2ppoint[ip])
            end
        # end #f

    elseif (mesh.nsd == 3)
        
        mesh.connijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl), Int64(mesh.ngl))
        mesh.conn_edgesijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.NEDGES_EL))

        for iel = 1:mesh.nelem
            #CGNS numbering: OK ref: HEXA...
            mesh.conn[iel, 1] = mesh.cell_node_ids[iel][1]#9
            mesh.conn[iel, 2] = mesh.cell_node_ids[iel][5]#11
            mesh.conn[iel, 3] = mesh.cell_node_ids[iel][6]#6
            mesh.conn[iel, 4] = mesh.cell_node_ids[iel][2]#1
            mesh.conn[iel, 5] = mesh.cell_node_ids[iel][3]#10
            mesh.conn[iel, 6] = mesh.cell_node_ids[iel][7]#12
            mesh.conn[iel, 7] = mesh.cell_node_ids[iel][8]#5
            mesh.conn[iel, 8] = mesh.cell_node_ids[iel][4]#4

            #OK
            mesh.connijk[iel, 1, 1, 1]       = mesh.cell_node_ids[iel][2]
            mesh.connijk[iel, ngl, 1, 1]     = mesh.cell_node_ids[iel][1]
            mesh.connijk[iel, ngl, ngl, 1]   = mesh.cell_node_ids[iel][5]
            mesh.connijk[iel, 1, ngl, 1]     = mesh.cell_node_ids[iel][6]
            mesh.connijk[iel, 1, 1, ngl]     = mesh.cell_node_ids[iel][4]
            mesh.connijk[iel, ngl, 1, ngl]   = mesh.cell_node_ids[iel][3]
            mesh.connijk[iel, ngl, ngl, ngl] = mesh.cell_node_ids[iel][7]
            mesh.connijk[iel, 1, ngl, ngl]   = mesh.cell_node_ids[iel][8]
            
        end
        
        #
        # Fill in elements dictionary needed by NodeOrdering.jl
        #
        elements = Dict(
            kk => mesh.conn[kk, 1:8]
            for kk = 1:mesh.nelem)
        element_types = Dict(
            kk => :Hexa8
            for kk = 1:mesh.nelem)
        
        #
        #Use NodeNumbering.jl
        #
        #adjacency = create_adjacency_graph(elements, element_types)
        #degrees = node_degrees(adjacency)
        #neworder = RCM(adjacency, degrees, tot_linear_poin, tot_linear_poin)
        #finalorder = renumbering(neworder)
        #RCM_adjacency = create_RCM_adjacency(adjacency, finalorder)
        #newmatrix = adjacency_visualization(RCM_adjacency)
        #display(UnicodePlots.heatmap(newmatrix))
        
        
        #
        # Rewrite coordinates in RCM order:
        #
        # open("./COORDS_LO_$rank.dat", "w") do f
            #open("./COORDS_LO.dat", "w") do f
            for ip = 1:mesh.npoin_linear
                mesh.x[ip] = get_node_coordinates(get_grid(model))[ip][1]
                mesh.y[ip] = get_node_coordinates(get_grid(model))[ip][2]
                mesh.z[ip] = get_node_coordinates(get_grid(model))[ip][3]
                mesh.ip2gip[ip] = point2ppoint[ip]
                # mesh.gip2owner[ip] = 1
                # @printf(f, " %.6f %.6f %.6f %d %d\n", mesh.x[ip],  mesh.y[ip], mesh.z[ip], ip, point2ppoint[ip])
            end
        # end #f
    end


    #
    # Add high-order points to edges, faces, and elements (volumes)
    #
    # initialize LGL struct and buyild Gauss-Lobatto-xxx points
    lgl = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.nop, backend)

    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ "; msg_rank = rank, suppress = mesh.msg_suppress)
    #
    # Edges
    #
    populate_conn_edge_el!(mesh, mesh.SD)
    add_high_order_nodes_edges!(mesh, lgl, mesh.SD, backend, edge2pedge)

    #
    # Faces
    #
    populate_conn_face_el!(mesh, mesh.SD)
    add_high_order_nodes_faces!(mesh, lgl, mesh.SD, face2pface)

    #
    # Volume
    #
    # NOTICE: in 2D we consider only edges. faces are the elements.
    #         
    add_high_order_nodes_volumes!(mesh, lgl, mesh.SD, elm2pelm)

    
    for ip = mesh.npoin_linear+1:mesh.npoin
        mesh.x[ip] = mesh.x_ho[ip]
        mesh.y[ip] = mesh.y_ho[ip]
        mesh.z[ip] = 0.0
        if (mesh.nsd > 2)
            mesh.z[ip] = mesh.z_ho[ip]
        end
    end

    mesh.xmax = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
    mesh.xmin = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)
    mesh.ymax = MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm)
    mesh.ymin = MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm)
    if (mesh.nsd > 2)
        mesh.zmax = MPI.Allreduce(maximum(mesh.z), MPI.MAX, comm)
        mesh.zmin = MPI.Allreduce(minimum(mesh.z), MPI.MIN, comm)
    end

    # for ip = 1: mesh.npoin
    #     if mesh.gip2owner[ip] != rank+1
    #         # @info mesh.x[ip], mesh.y[ip], mesh.gip2owner[ip], rank+1
    #     end
    # end



    gpelm_ghost = KernelAbstractions.zeros(backend, TInt, 0)
    gpfacets_ghost = KernelAbstractions.zeros(backend, TInt, 0)
    gpfacets_owner = KernelAbstractions.zeros(backend, TInt, 0)
    gcelm_ghost = KernelAbstractions.zeros(backend, TInt, 0)
    gcfacets_ghost = KernelAbstractions.zeros(backend, TInt, 0)
    gcfacets_owner = KernelAbstractions.zeros(backend, TInt, 0)
    mesh.num_hanging_facets = num_hanging_facets
    pedge2edge = Dict(val => idx for (idx, val) in enumerate(edge2pedge))
    pface2face = Dict(val => idx for (idx, val) in enumerate(face2pface))
    pelm2elm = Dict(val => idx for (idx, val) in enumerate(elm2pelm))
    # mesh.non_conforming_facets = [KernelAbstractions.zeros(backend, TInt, 0, 0, 0, 0) for _ in 1:num_hanging_facets]
    if ladaptive == true
     
        cell_fecet_pids = get_faces(dtopology, mesh.nsd, mesh.nsd-1) #edge map from local to global numbering i.e. iedge_g = cell_edge_ids[1:NELEM][1:NEDGES_EL]
        facet_cell_pids = get_faces(dtopology, mesh.nsd-1, mesh.nsd) #edge map from local to global numbering i.e. iedge_g = cell_edge_ids[1:NELEM][1:NEDGES_EL]
        if mesh.nsd == 2
            offset = 4
            for idx in 1: num_hanging_facets
                cfacet = idx+num_regular_facets
                cid = facet_cell_pids[cfacet][1]
                pid, lfacetid, half = hanging_facet_glue[idx]
                pfacet = cell_fecet_pids[pid][lfacetid-offset]
                gfacet_p = local_to_global(edgids)[pfacet]
                gfacet_c = local_to_global(edgids)[cfacet]
                if (lfacetid == 7) || (lfacetid == 8)
                    half = 3-half
                end
                # own child facet, ghost parent facet 
                if (cfacet ∈ eg2peg) && (pfacet ∉ eg2peg)
                    # add ghost ip
                    gpid = local_to_global(elgids)[pid]
                    # push!(gpfacets_ghost, gfacet_p)
                    push!(gpelm_ghost, gpid)
                    push!(gpfacets_ghost, lfacetid - offset)
                    push!(gpfacets_owner, local_to_owner(edgids)[pfacet]-1)
                    pfacet = -pfacet
                    push!(mesh.non_conforming_facets_parents_ghost, [cid, lfacetid - offset, half])
                    continue
                end
                # ghost child facet, own parent facet
                if (cfacet ∉ eg2peg) && (pfacet ∈ eg2peg)
                    gcid = local_to_global(elgids)[cid]
                    push!(gcelm_ghost, gcid)
                    push!(gcfacets_ghost, lfacetid - offset)
                    # push!(gcfacets_ghost, gfacet_c)
                    push!(gcfacets_owner, local_to_owner(edgids)[cfacet]-1)
                    cfacet = -cfacet
                    push!(mesh.non_conforming_facets_children_ghost, [pid, lfacetid - offset, half])
                    continue
                end
                # ghost child facet, ghost parent facet
                if (cfacet ∉ eg2peg) && (pfacet ∉ eg2peg)
                    pfacet = -pfacet
                    cfacet = -cfacet
                    continue
                end
                # @info rank, cfacet, pfacet, [pid, lfacetid, half], gfacet_c, gfacet_p
                # mesh.non_conforming_facets[idx] = [cfacet, pfacet, lfacetid-offset, half]
                push!(mesh.non_conforming_facets, [cfacet, cid, pfacet, pid, lfacetid - offset, half])
            end
            # reorder mesh.non_conforming_facets_parents_ghost and mesh.non_conforming_facets_children_ghost in rank orders\
            sorted_idx = sortperm(gpfacets_owner)
            sort!(gpfacets_owner)
            gpelm_ghost    .= gpelm_ghost[sorted_idx]
            gpfacets_ghost .= gpfacets_ghost[sorted_idx]
            mesh.non_conforming_facets_parents_ghost .= mesh.non_conforming_facets_parents_ghost[sorted_idx]

            sorted_idx = sortperm(gcfacets_owner)
            sort!(gcfacets_owner)
            gcelm_ghost    .= gcelm_ghost[sorted_idx]
            gcfacets_ghost .= gcfacets_ghost[sorted_idx]
            mesh.non_conforming_facets_children_ghost .= mesh.non_conforming_facets_children_ghost[sorted_idx]
            # @info "edge2pedge", rank, edge2pedge
            ghost_p_or_c = 1
            mesh.pgip_ghost, mesh.pgip_owner = get_ghost_ips(gpelm_ghost, gpfacets_ghost, gpfacets_owner, mesh.connijk, pelm2elm, mesh.ip2gip, ngl, ghost_p_or_c, comm)
            ghost_p_or_c = 2
            mesh.cgip_ghost, mesh.cgip_owner = get_ghost_ips(gcelm_ghost, gcfacets_ghost, gcfacets_owner, mesh.connijk, pelm2elm, mesh.ip2gip, ngl, ghost_p_or_c, comm)

        elseif mesh.nsd == 3
            offset = 20
            for idx in 1: num_hanging_facets
                half_1 = 0
                half_2 = 0
                cfacet     = idx+num_regular_facets
                cfacet = idx+num_regular_facets
                cid = facet_cell_pids[cfacet][1]
                pid, lfacetid, half = hanging_facet_glue[idx]
                
                pfacet     = cell_fecet_pids[pid][lfacetid-offset]
                gfacet_p   = local_to_global(fgids)[pfacet]
                gfacet_c   = local_to_global(fgids)[cfacet]
                if (cfacet ∈ f2pf) && (pfacet ∉ f2pf)
                    pfacet = -pfacet
                end
                if (cfacet ∉ f2pf) && (pfacet ∈ f2pf)
                    cfacet = -cfacet
                end
                if (cfacet ∉ f2pf) && (pfacet ∉ f2pf)
                    pfacet = -pfacet
                    cfacet = -cfacet
                end
                if (lfacetid-offset == 1) || (lfacetid-offset == 2) 
                    if half == 1
                        half_1 = 1
                        half_2 = 2
                    elseif half == 2
                        half_1 = 2
                        half_2 = 2
                    elseif half == 3
                        half_1 = 1
                        half_2 = 1
                    elseif half == 4
                        half_1 = 2
                        half_2 = 1
                    end
                elseif (lfacetid-offset == 3) || (lfacetid-offset == 4) 
                    if half == 1
                        half_1 = 1
                        half_2 = 2
                    elseif half == 2
                        half_1 = 2
                        half_2 = 2
                    elseif half == 3
                        half_1 = 1
                        half_2 = 1
                    elseif half == 4
                        half_1 = 2
                        half_2 = 1
                    end
                elseif (lfacetid-offset == 5) || (lfacetid-offset == 6)
                    if half == 1
                        half_1 = 2
                        half_2 = 2
                    elseif half == 2
                        half_1 = 2
                        half_2 = 1
                    elseif half == 3
                        half_1 = 1
                        half_2 = 2
                    elseif half == 4
                        half_1 = 1
                        half_2 = 1
                    end
                end
                push!(mesh.non_conforming_facets, [cfacet, cid, pfacet, pid, lfacetid - offset, half_1, half_2])
                # @info rank, cfacet, pfacet, facet_glue
                # comm_ip = intersect(mesh.conn[pid,:], mesh.conn[cid,:])
                # @info "coords: ",  mesh.x[mesh.conn[pid,1:8]], mesh.y[mesh.conn[pid,1:8]], mesh.z[mesh.conn[pid,1:8]]
                # @info "comm_coord", mesh.x[comm_ip], mesh.y[comm_ip], mesh.z[comm_ip] 
                for k = 1:ngl
                    for j = 1:ngl
                        for i = 1:ngl
                            pip = mesh.connijk[pid, k, j, i]
                            for n = 1:ngl
                                for m = 1:ngl
                                    for l = 1:ngl
                                        cip = mesh.connijk[cid, n, m, l]
                                        if pip == cip
                                            if lfacetid-offset == 1
                                                # @info "front, ", half, (k, i), half_1, half_2
                                                x1 = k
                                                x2 = i
                                            elseif lfacetid-offset == 2
                                                # @info "back, ", half, (k, i), half_1, half_2
                                                x1 = k
                                                x2 = i
                                            elseif lfacetid-offset == 3
                                                # @info "bottom, ", half, (k, j), half_1, half_2
                                                x1 = k
                                                x2 = j
                                            elseif lfacetid-offset == 4
                                                # @info "top, ", half, (k, j), half_1, half_2
                                                x1 = k
                                                x2 = j
                                            elseif lfacetid-offset == 5
                                                # @info "right, ", half, (j, i), half_1, half_2
                                                x1 = j
                                                x2 = i
                                            elseif lfacetid-offset == 6
                                                # @info "left, ", half, (j, i), half_1, half_2
                                                x1 = j
                                                x2 = i
                                            end
                                            if x1==1
                                                half1 = 2
                                            else
                                                half1 = 1
                                            end
                                            if x2 == 1
                                                half2 = 2
                                            else
                                                half2 = 1
                                            end
                                            @test (half1 == half_1) && (half2 == half_2) 
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    mesh.gnpoin    = MPI.Allreduce(maximum(mesh.ip2gip), MPI.MAX, comm)
    mesh.gip2owner = find_gip_owner(mesh.ip2gip)
    mesh.gip2ip    = KernelAbstractions.zeros(backend, TInt, mesh.gnpoin)

    for (ip, gip) in enumerate(mesh.ip2gip)
        mesh.gip2ip[gip] = ip
    end
    #----------------------------------------------------------------------
    # Extract boundary edges and faces nodes:
    #----------------------------------------------------------------------
    #
    # Bdy edges
    #
    if mesh.nsd == 2
        # isboundary_edge = compute_isboundary_face(topology, EDGE_flg)
        isboundary_edge = fill(false, mesh.nedges)  
        # boundary_edges  = findall(x -> size(x,1) == 1, mesh.facet_cell_ids)
        # isboundary_edge[boundary_edges] .= true
        
        # @info isboundary_edge
        #
        # Get labels contained in the current GMSH grid:
        #
        n_semi_inf = 0
        labels = get_face_labeling(model)
        # @info rank, labels.tag_to_name
        for ilabel in labels.tag_to_name
            edges_to_tag  = get_face_tag_index(labels,ilabel,EDGE_flg)
            idx_edges_inflow = findall( x -> x == 1, edges_to_tag)
            #    
            # Tag the boundary edge with its type as defined in the user-provided GMSH file:
            #
            for idx in idx_edges_inflow
                mesh.edge_type[idx] = ilabel
                isboundary_edge[idx] = true
            end
            # @info mesh.edge_type
        end
        iedge_bdy = 1
        for iedge = 1:mesh.nedges #total nedges
            if isboundary_edge[iedge] == true
                # if rank == 1
                #     @info mesh.x[mesh.poin_in_edge[iedge, 1]], mesh.y[mesh.poin_in_edge[iedge, 1]]
                # end
                for igl = 1:mesh.ngl
                    mesh.poin_in_bdy_edge[iedge_bdy, igl] = mesh.poin_in_edge[iedge, igl]
                    mesh.bdy_edge_type[iedge_bdy] = mesh.edge_type[iedge]
                    mesh.bdy_edge_in_elem[iedge_bdy] = mesh.facet_cell_ids[iedge][1]
                    # if (size(mesh.facet_cell_ids[iedge],1) ≠ 1)
                    #     s = """
                    #     Check boundary elements! size(mesh.facet_cell_ids[iedge],1) ≠ 1 
                    #         """
                
                    #     @error s
                    # end
                end
                if (mesh.bdy_edge_type[iedge_bdy] == "Laguerre")
                    n_semi_inf += 1
                end
                iedge_bdy += 1
            end
        end
        n_semi_infg = MPI.Allreduce(n_semi_inf, MPI.SUM, comm)
        if (n_semi_infg > 0) 
            mesh.lLaguerre = true
        end
        # build mesh data structs for Laguerre semi-infinite elements
        if (mesh.lLaguerre && n_semi_inf > 0)
            mesh.gip2ip = KernelAbstractions.zeros(backend, TInt, mesh.gnpoin + n_semi_infg*(mesh.ngl-1)*(mesh.ngr-1)+mesh.ngr-1)
            ip2gip_new = KernelAbstractions.zeros(backend, TInt, mesh.npoin + n_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)+mesh.ngr-1)
            ip2gip_new[1:mesh.npoin] .= mesh.ip2gip[:]
            gr = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta],backend) 
            factorx = inputs[:xfac_laguerre]#0.1
            factory = inputs[:yfac_laguerre]#0.025
            mesh.connijk_lag = KernelAbstractions.zeros(backend, TInt, Int64(n_semi_inf), Int64(mesh.ngl), Int64(mesh.ngr),1)
            bdy_normals = zeros(n_semi_inf, 2)
            bdy_tangents = zeros(n_semi_inf, 2)
            e_iter = 1
            iter = mesh.npoin + 1
            x_new = KernelAbstractions.zeros(backend, TFloat, mesh.npoin + n_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)+mesh.ngr-1)
            y_new = KernelAbstractions.zeros(backend, TFloat, mesh.npoin + n_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)+mesh.ngr-1)
            x_new[1:mesh.npoin] .= mesh.x[:]
            y_new[1:mesh.npoin] .= mesh.y[:]
            for iedge = 1:size(mesh.bdy_edge_type,1)
                if (mesh.bdy_edge_type[iedge] == "Laguerre") 
                    iel = mesh.bdy_edge_in_elem[iedge]
                    #find tangent and normal vectors to the boundary
                    ip = mesh.poin_in_bdy_edge[iedge,1]
                    ip1 = mesh.poin_in_bdy_edge[iedge,2]
                    ipend = mesh.poin_in_bdy_edge[iedge,mesh.ngl]
                    #tangent vector 
                    x = mesh.x[ip]
                    x1 = mesh.x[ip1]
                    y = mesh.y[ip]
                    y1 = mesh.y[ip1]
                    xend = mesh.x[ipend]
                    yend = mesh.y[ipend]
                    tan = [x-x1, y-y1]
                    # deduce normal vector components
                    if (tan[2] > 1e-7)
                        x2 = 1.0
                        y2 = -x2*tan[1]/tan[2]
                    else
                        y2 = 1.0
                        x2 = -y2*tan[2]/tan[1]
                    end
                    nor = [x2,y2]
                    # generate unit versions of tangent and normal vectors
                    modu = sqrt(tan[1]^2+tan[2]^2)
                    tan = tan * (1/modu)
                    modu = sqrt(nor[1]^2+nor[2]^2)
                    nor = nor * (1/modu)
                    #make sure normal is outward facing
                    l = 1
                    m = 1
                    l1 = 1
                    m1 = 1
                    for ii=1:mesh.ngl
                        for jj=1:mesh.ngl
                            if (mesh.connijk[iel,ii,jj] == ip)
                                l=ii
                                m=jj
                            end
                            if (mesh.connijk[iel,ii,jj] == ip1)
                                l1 = ii
                                m1 = jj
                            end
                        end
                    end
                    if (l == l1)
                        ip2 = mesh.connijk[iel,3,m]
                    else
                        ip2 = mesh.connijk[iel,l,3]
                    end
                    v = [mesh.x[ip2]-x, mesh.y[ip2]-y]
                    if (dot(v,nor) > 0.0)
                        nor .= -nor
                    end
                    
                    bdy_normals[e_iter,:] .= nor
                    bdy_tangents[e_iter,:] .= tan
                    if (dot(tan,[1,0]) == 0)
                        ord = Int64(round((max(y,yend)-mesh.ymin)/abs(y-yend)))
                    else
                        ord = Int64(round((max(x,xend)-mesh.xmin)/abs(x-xend)))
                    end
                    
                    for i=1:mesh.ngl
                        ip = mesh.poin_in_bdy_edge[iedge,i]
                        mesh.connijk_lag[e_iter,i,1] = ip
                        for j=2:mesh.ngr
                            gip = Int64(mesh.gnpoin + (ord-1)*(mesh.ngl-1)*(mesh.ngr-1) + (mesh.ngr-1)*(i-1) + j-1)
			    #@info ord, gip, mesh.gnpoin, x, xend, mesh.xmax, mesh.xmin
                            if (inputs[:xscale]==1.0)
                                x_temp = mesh.x[ip] + nor[1]*gr.ξ[j]*factorx
                            else
                                x_temp = mesh.x[ip] + nor[1]*gr.ξ[j]*factorx/(inputs[:xscale] * 0.5)
			    end
                            if (inputs[:yscale] == 1.0)
			                    y_temp = mesh.y[ip] + nor[2]*gr.ξ[j]*factory
			                else 
                                y_temp = mesh.y[ip] + nor[2]*gr.ξ[j]*factory/(inputs[:yscale] * 0.5)
                            end
			                matched = 0
                            if (i == mesh.ngl || i == 1)
                                iter_end = 0
                                while (matched == 0 && iter_end == 0)
                                    for e_check = 1:n_semi_inf
                                        for i_check =1:mesh.ngl
                                            for j_check =1:mesh.ngr
                                                ip_check = mesh.connijk_lag[e_check,i_check,j_check]
                                                if (ip_check != 0.0 && e_check != e_iter)
                                                    if (AlmostEqual(x_temp,x_new[ip_check]) && AlmostEqual(y_temp,y_new[ip_check]))
                                                        mesh.connijk_lag[e_iter,i,j] = ip_check
                                                        matched = 1
                                                    end
                                                end
                                            end
                                        end 
                                    end
                                    iter_end = 1
                                end    
                            else
                                x_new[iter] = x_temp#mesh.x[ip] + nor[1]*gr.ξ[j]*factorx
                                y_new[iter] = y_temp#mesh.y[ip] + nor[2]*gr.ξ[j]*factory
                                ip2gip_new[iter] = gip
                                mesh.connijk_lag[e_iter,i,j] = iter
                                mesh.gip2ip[gip] = iter
                                @info gip, iter, rank, x_temp, y_temp, rank
                                iter += 1
                                matched = 1
                                
                            end
                            if (matched == 0)
                                x_new[iter] = x_temp#mesh.x[ip] + nor[1]*gr.ξ[j]*factorx
                                y_new[iter] = y_temp#mesh.y[ip] + nor[2]*gr.ξ[j]*factory
                                ip2gip_new[iter] = gip
                                mesh.connijk_lag[e_iter,i,j] = iter
                                mesh.gip2ip[gip] = iter
                                @info gip, iter, x_temp, y_temp,rank
                                iter += 1 
                            end
                            
                            #@info nor[1],nor[2],x_new[iter],y_new[iter], mesh.x[ip],mesh.y[ip]
                        end
                    end
                    e_iter += 1
                end
            end
            #@info mesh.npoin, iter - 1, mesh.ngr, n_semi_inf, e_iter - 1
            mesh.npoin_original = mesh.npoin
            mesh.npoin = iter -1
            mesh.gnpoin += n_semi_infg*(mesh.ngl-1)*(mesh.ngr-1)
            mesh.x = x_new
            mesh.y = y_new
            mesh.ip2gip = ip2gip_new
            mesh.z = KernelAbstractions.zeros(backend, TFloat, mesh.npoin)
           # mesh.nelem_semi_inf = n_semi_inf
        end
        if (mesh.lLaguerre)
            mesh.nelem_semi_inf = n_semi_inf
            mesh.gip2owner = find_gip_owner(mesh.ip2gip)
            for (ip, gip) in enumerate(mesh.ip2gip)
                mesh.gip2ip[gip] = ip
            end
            mesh.gnpoin    = MPI.Allreduce(maximum(mesh.ip2gip), MPI.MAX, comm)
            
        end
    elseif mesh.nsd > 2
        # isboundary_face = compute_isboundary_face(topology, FACE_flg)
        isboundary_face = fill(false, mesh.nfaces)
        # boundary_faces  = findall(x -> size(x,1) == 1, mesh.facet_cell_ids)
        # isboundary_face[boundary_faces] .= true
        #
        # Get labels contained in the current GMSH grid:
        #
        labels = get_face_labeling(model)
        # @info labels.tag_to_name, labels.tag_to_entities, labels.d_to_dface_to_entity
        for ilabel in labels.tag_to_name
            if ilabel == "internal"
                continue
            end
            faces_to_tag  = get_face_tag_index(labels,ilabel,FACE_flg)
            # @info "faces_to_tag", faces_to_tag
            idx_faces_inflow = findall( x -> x == 1, faces_to_tag)
            #    
            # Tag the boundary edge with its type as defined in the user-provided GMSH file:
            #
            for idx in idx_faces_inflow
                mesh.face_type[idx]  = ilabel
                isboundary_face[idx] = true
            end
        end
        get_bdy_poin_in_face_on_edges!(mesh, @view(isboundary_face[:]), mesh.SD)
        # @info isboundary_face
        iface_bdy = 1
        for iface in findall(x -> x == true, isboundary_face) #total nedges
            # if isboundary_face[iface] == true
                for igl = 1:mesh.ngl
                    for jgl = 1:mesh.ngl
                        mesh.poin_in_bdy_face[iface_bdy, igl,jgl] = mesh.poin_in_face[iface, igl,jgl]
                        mesh.bdy_face_type[iface_bdy] = mesh.face_type[iface]
                        mesh.bdy_face_in_elem[iface_bdy] = mesh.facet_cell_ids[iface][1]
                        if (size(mesh.facet_cell_ids[iface],1) ≠ 1)
                            s = """
                            Check boundary elements! size(mesh.facet_cell_ids[iface],1) ≠ 1 
                                """
                            @info iface_bdy, iface, mesh.x[mesh.poin_in_face[iface, igl,jgl]], mesh.y[mesh.poin_in_face[iface, igl,jgl]], mesh.z[mesh.poin_in_face[iface, igl,jgl]]
                            @error s
                        end
                        # @info "face point number", mesh.poin_in_face[iface,igl,jgl],iface,igl,jgl
                    end
                end
                iface_bdy += 1
            # end
        end
        #=for iface =1:mesh.nfaces_bdy
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    ip = mesh.poin_in_bdy_face[iface,i,j]
                    @info "bdy points coords", mesh.x[ip],mesh.y[ip],mesh.z[ip]
                end
            end
        end=#
    end

    #----------------------------------------------------------------------
    # END Extract boundary edges and faces nodes
    #----------------------------------------------------------------------


    #----------------------------------------------------------------------
    # periodicity_restructure for MPI
    #----------------------------------------------------------------------
    if mesh.nsd > 2

        nor1 = [1.0, 0.0, 0.0]
        nor2 = [0.0, 1.0, 0.0]
        nor3 = [0.0, 0.0, 1.0]
        if ("periodicx" in mesh.bdy_face_type)
            finder = false
            iface_bdy = 1
            while (finder == false)
                if (mesh.bdy_face_type[iface_bdy] == "periodicx")
                    ip = mesh.poin_in_bdy_face[iface_bdy,1,1]
                    ip1 = mesh.poin_in_bdy_face[iface_bdy,1,2]
                    ip2 = mesh.poin_in_bdy_face[iface_bdy,2,1]
                    t1 = [mesh.x[ip] - mesh.x[ip1],mesh.y[ip] - mesh.y[ip1], mesh.z[ip] - mesh.z[ip1]]
                    t2 = [mesh.x[ip] - mesh.x[ip2],mesh.y[ip] - mesh.y[ip2], mesh.z[ip] - mesh.z[ip2]]
                    s1 = t1[2]*t2[3] - t1[3]*t2[2]
                    s2 = t1[3]*t2[1] - t1[1]*t2[3]
                    s3 = t1[1]*t2[2] - t1[2]*t2[1]
                    mag = sqrt(s1^2 + s2^2 + s3^2)
                    nor1 .= [s1/mag, s2/mag, s3/mag]
                    finder = true
                else
                    iface_bdy +=1
                end
            end
        end
        if ("periodicz" in mesh.bdy_face_type)
            finder = false
            iface_bdy = 1
            while (finder == false)
                if (mesh.bdy_face_type[iface_bdy] == "periodicz")
                    ip = mesh.poin_in_bdy_face[iface_bdy,1,1]
                    ip1 = mesh.poin_in_bdy_face[iface_bdy,1,2]
                    ip2 = mesh.poin_in_bdy_face[iface_bdy,2,1]
                    t1 = [mesh.x[ip] - mesh.x[ip1],mesh.y[ip] - mesh.y[ip1], mesh.z[ip] - mesh.z[ip1]]
                    t2 = [mesh.x[ip] - mesh.x[ip2],mesh.y[ip] - mesh.y[ip2], mesh.z[ip] - mesh.z[ip2]]
                    s1 = t1[2]*t2[3] - t1[3]*t2[2]
                    s2 = t1[3]*t2[1] - t1[1]*t2[3]
                    s3 = t1[1]*t2[2] - t1[2]*t2[1]
                    mag = sqrt(s1^2 + s2^2 + s3^2)
                    nor2 = [s1/mag, s2/mag, s3/mag] 
                    finder = true
                else
                    iface_bdy +=1
                end
            end
        end
        if ("periodicy" in mesh.bdy_face_type)
            finder = false
            iface_bdy = 1
            while (finder == false)
                if (mesh.bdy_face_type[iface_bdy] == "periodicy")
                    ip = mesh.poin_in_bdy_face[iface_bdy,1,1]
                    ip1 = mesh.poin_in_bdy_face[iface_bdy,1,2]
                    ip2 = mesh.poin_in_bdy_face[iface_bdy,2,1]
                    t1 = [mesh.x[ip] - mesh.x[ip1],mesh.y[ip] - mesh.y[ip1], mesh.z[ip] - mesh.z[ip1]]
                    t2 = [mesh.x[ip] - mesh.x[ip2],mesh.y[ip] - mesh.y[ip2], mesh.z[ip] - mesh.z[ip2]]
                    s1 = t1[2]*t2[3] - t1[3]*t2[2]
                    s2 = t1[3]*t2[1] - t1[1]*t2[3]
                    s3 = t1[1]*t2[2] - t1[2]*t2[1]
                    mag = sqrt(s1^2 + s2^2 + s3^2)
                    nor3 = [s1/mag, s2/mag, s3/mag]
                    finder = true
                else
                    iface_bdy +=1
                end
            end
        end
        restructure4periodicity_3D(mesh, nor1, "periodicx")
        # @info mesh.ip2gip
        restructure4periodicity_3D(mesh, nor2, "periodicz")
        restructure4periodicity_3D(mesh, nor3, "periodicy")
        # @info mesh.ip2gip

    end
    #----------------------------------------------------------------------
    # END periodicity_restructure for MPI
    #----------------------------------------------------------------------


    #
    #
    # Free memory of obsolete arrays
    #
    mesh.x_ho = zeros(1)
    mesh.y_ho = zeros(1)
    mesh.z_ho = zeros(1)
    GC.gc()
    #
    # END Free memory of obsolete arrays
    #
    
    #=open("./COORDS_GLOBAL.dat", "w") do f
        for ip = 1:mesh.npoin
            #@printf(" %.6f %.6f %.6f %d\n", mesh.x[ip],  mesh.y[ip], mesh.z[ip], ip)
            @printf(f, " %.6f %.6f %.6f %d\n", mesh.x[ip],  mesh.y[ip], mesh.z[ip], ip)
        end
    end #f


    # =#
   
    # write_vtk_grid_only(mesh.SD, mesh, "VTK_grid", "./", parts, nparts)

    #
    # gridapDistributed test on gmsh
    #
    # @mystop("my stop at mesh.jl L135")
    # end gridapDistributed test on gmsh


    #show(stdout, "text/plain", mesh.conn')
    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ DONE"; msg_rank = rank, suppress = mesh.msg_suppress)
    if isnothing(adapt_flags)
        return partitioned_model
    else
        return partitioned_model, glue.n2o_faces_map[mesh.nsd+1]
    end

    #writevtk(model,"gmsh_grid")
end


function restructure4periodicity_3D(mesh, norm, periodic_direction)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    per_ip = Int[]
    ngl = mesh.ngl
    for iface_bdy =1:size(mesh.bdy_face_type,1)
        for k=1:ngl
            for l=1:ngl
                ip = mesh.poin_in_bdy_face[iface_bdy,k,l]
                if (mesh.bdy_face_type[iface_bdy] == periodic_direction)
                    per_ip = [per_ip; ip]
                end
            end
        end
    end
    ### remove duplicates
    unique!(per_ip)
    x_local  = mesh.x[per_ip]
    y_local  = mesh.y[per_ip]
    z_local  = mesh.z[per_ip]
    per_gip  = mesh.ip2gip[per_ip]
    ip_owner = mesh.gip2owner[per_ip]
    # @info  mesh.x[per_ip]

    # Gather arrays onto the root processor (rank 0)
    root = 0

    # Gather per_gip
    buffer_sz::Int32    = size(per_ip, 1)
    # @info rank, buffer_sz
    recv_counts  = MPI.Gather(buffer_sz, 0, comm)
    # total_counts = sum(recv_counts)
    # if total_counts == 0
        # return
    # end
    # else
    # if total_counts == 0
    #     return
    # end
    
    x_gather     = MPI.gather(x_local, comm)
    y_gather     = MPI.gather(y_local, comm)
    z_gather     = MPI.gather(z_local, comm)
    gathered_per = MPI.gather(per_gip, comm)
    owner_gather = MPI.gather(ip_owner, comm)
    if mesh.rank == root
    
    # On the root processor, combine and remove duplicates
        # Concatenate gathered arrays
        x              = vcat(x_gather...)
        y              = vcat(y_gather...)
        z              = vcat(z_gather...)
        global_per_gip = vcat(gathered_per...)
        owner          = vcat(owner_gather...)

        sz = size(global_per_gip,1)
        for i = 1:sz
            i1 = i+1
            for i1 = (i+1):sz
                if global_per_gip[i] == global_per_gip[i1]
                    continue
                end
                vec = [x[i] - x[i1], y[i] - y[i1], z[i] - z[i1]]
                # @info vec, norm
                if (determine_colinearity(vec, norm))
                    xt = x[i1]
                    yt = y[i1]
                    zt = z[i1]
                    xi = x[i]
                    yi = y[i]
                    zi = z[i]
                    if (yi == 0 && yt == 0 && zi == 0 && zt == 0)
                        comp1 = xi < xt
                    elseif (yi == 0 && yt == 0)
                        comp1 = xi*abs(zi) < xt*abs(zt)
                    elseif (zi == 0 && zt == 0)
                        comp1 = xi*abs(zi) < xt*abs(yt)
                    else
                        comp1 = xi*abs(yi*zi) < xt*abs(yt*zt)
                    end
                    if (xi ==0 && xt == 0 && zi == 0 && zt ==0)
                        comp2 = yi < yt
                    elseif (xi == 0 && xt == 0)
                        comp2 = yi*abs(zi) < yt*abs(zt)
                    elseif (zi == 0 && zt == 0)
                        comp2 = yi*abs(xi) < yt*abs(xt)
                    else
                        comp2 = yi*abs(xi*zi) < yt*abs(xt*zt)
                    end
                    if (xi == 0 && xt == 0 && yi == 0 && yt ==0)
                        comp3 = zi < zt
                    elseif (xi == 0 && xt == 0)
                        comp3 = zi*abs(yi) < zt*abs(yt)
                    elseif (yi == 0 && yt == 0)
                        comp3 = zi*abs(xi) < zt*abs(xt)
                    else
                        comp3 = zi*abs(xi*yi) < zt*abs(xt*yt)
                    end
                    # @info "found", global_per_gip[i], global_per_gip[i1]
                    if (comp1 || comp2 || comp3)    
                        global_per_gip[i1] = global_per_gip[i]
                        if owner[i1] != owner[i]
                            owner[i1] = owner[i]
                        end
                    else
                        global_per_gip[i] = global_per_gip[i1]
                        if owner[i1] != owner[i]
                            owner[i] = owner[i1]
                        end
                    end
                    # break
                else
                    continue
                end
            end
        end
        # do something for global_per_gip
        s_gip_vbuf   = VBuffer(global_per_gip, recv_counts)
        s_owner_vbuf = VBuffer(owner, recv_counts)
    else
        s_gip_vbuf   = VBuffer(nothing)
        s_owner_vbuf = VBuffer(nothing)
    end
    MPI.Barrier(comm)
    per_ip_updated = MPI.Scatterv!(s_gip_vbuf,zeros(eltype(per_gip), buffer_sz), 0, comm)
    owner_updated  = MPI.Scatterv!(s_owner_vbuf,zeros(eltype(ip_owner), buffer_sz), 0, comm)
    # per_ip_updated = MPI.Scatterv!(global_per_gip,buffer_sz, 0, comm)
        
    mesh.ip2gip[per_ip]    .= per_ip_updated
    mesh.gip2owner[per_ip] .= owner_updated
end

function find_gip_owner(a)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    
    # Gather all elements from all ranks to rank 0
    all_elements = MPI.gather(a, comm)
    
    if rank == 0
        # Flatten the gathered list
        flat_elements = vcat(all_elements...)
        # all_owners = [i for i in 1:size for _ in 1:length(all_elements[i])]
        all_owners = [i for i in 0:size-1 for _ in 1:length(all_elements[i+1])]

        # Create a dictionary to store the smallest rank for each element
        element_owner_map = Dict{Int, Int}()
        for i in 1:length(flat_elements)
            el = flat_elements[i]
            owner = all_owners[i]
            if !haskey(element_owner_map, el) || owner < element_owner_map[el]
                element_owner_map[el] = owner
            end
        end

        # Map the owners back to the original elements in each rank's vector
        all_owners_result = [element_owner_map[el] for el in flat_elements]
        
        # Split the ownership result according to the original vectors
        chunked_owners = [all_owners_result[sum(length.(all_elements)[1:i-1])+1:sum(length.(all_elements)[1:i])] for i in 1:size]
    else
        chunked_owners = nothing
    end

    # Scatter the ownership chunks back to all ranks
    element_owners = MPI.scatter(chunked_owners, comm)

    return element_owners
end

function determine_colinearity(vec1, vec2)
    if AlmostEqual(norm(vec1), 0.0)
        return false
    end
    # For 2D vectors
    if size(vec1, 1) < 3
        if AlmostEqual(vec1[1], 0.0) && AlmostEqual(vec2[1], 0.0)
            return abs(vec1[2]) > 1e-7 && abs(vec2[2]) > 1e-7
        elseif AlmostEqual(vec1[2], 0.0) && AlmostEqual(vec2[2], 0.0)
            return abs(vec1[1]) > 1e-7 && abs(vec2[1]) > 1e-7
        elseif abs(vec1[2]) > 1e-7 && abs(vec2[2]) > 1e-7 && abs(vec1[1]) > 1e-7 && abs(vec2[1]) > 1e-7
            return AlmostEqual((vec1[1] + 1e-16) / (vec2[1] + 1e-16), (vec1[2] + 1e-16) / (vec2[2] + 1e-16))
        else
            return false
        end

    # For 3D vectors
    else
        return abs(vec1[2] * vec2[3] - vec1[3] * vec2[2]) < 1e-7 &&
               abs(vec1[3] * vec2[1] - vec1[1] * vec2[3]) < 1e-7 &&
               abs(vec1[1] * vec2[2] - vec1[2] * vec2[1]) < 1e-7
    end
end

function populate_conn_edge_el!(mesh::St_mesh, SD::NSD_2D)
    
    for iel = 1:mesh.nelem
        #
        # CGNS numbering
        #
        ip1 = mesh.cell_node_ids[iel][4]
        ip2 = mesh.cell_node_ids[iel][3]
        ip3 = mesh.cell_node_ids[iel][1]
        ip4 = mesh.cell_node_ids[iel][2]
        
	    # Edges bottom face:
	    iedg_el = 1
        mesh.conn_edge_el[1, iedg_el, iel] = ip1
	    mesh.conn_edge_el[2, iedg_el, iel] = ip2
        iedg_el = 2
        mesh.conn_edge_el[1, iedg_el, iel] = ip2
        mesh.conn_edge_el[2, iedg_el, iel] = ip3
        iedg_el = 3
        mesh.conn_edge_el[1, iedg_el, iel] = ip3
        mesh.conn_edge_el[2, iedg_el, iel] = ip4
        iedg_el = 4
        mesh.conn_edge_el[1, iedg_el, iel] = ip4
        mesh.conn_edge_el[2, iedg_el, iel] = ip1
    end
    
end #populate_edge_el!

function populate_conn_edge_el!(mesh::St_mesh, SD::NSD_3D)

    mesh.conn_unique_edges1 = Array{Int64}(undef, mesh.nedges, 2)
    
    cache_edge_ids = array_cache(mesh.cell_edge_ids) # allocation here
    for iel = 1:mesh.nelem
        #
        # CGNS numbering
        #
        ip1 = mesh.connijk[iel, 1, 1, 1]
        ip2 = mesh.connijk[iel, mesh.ngl, 1, 1]
        ip3 = mesh.connijk[iel, mesh.ngl, mesh.ngl, 1]
        ip4 = mesh.connijk[iel, 1, mesh.ngl, 1]
        ip5 = mesh.connijk[iel, 1, 1, mesh.ngl]
        ip6 = mesh.connijk[iel, mesh.ngl, 1, mesh.ngl]
        ip7 = mesh.connijk[iel, mesh.ngl, mesh.ngl, mesh.ngl]
        ip8 = mesh.connijk[iel, 1, mesh.ngl, mesh.ngl]
        
        edge_ids = getindex!(cache_edge_ids, mesh.cell_edge_ids, iel)
        
    end

end #populate_edge_el!

function populate_conn_face_el!(mesh::St_mesh, SD::NSD_2D)
    
    for iel = 1:mesh.nelem
        #
        # CGNS numbering
        #
        ip1 = mesh.cell_node_ids[iel][4]
        ip2 = mesh.cell_node_ids[iel][3]
        ip3 = mesh.cell_node_ids[iel][1]
        ip4 = mesh.cell_node_ids[iel][2]
        
        #
        # Local faces node connectivity:
        # i.e. what nodes belong to a given local face in iel:
        #
        face_el = 1
        mesh.conn_face_el[1, face_el, iel] = ip1
        mesh.conn_face_el[2, face_el, iel] = ip2
        mesh.conn_face_el[3, face_el, iel] = ip3
        mesh.conn_face_el[4, face_el, iel] = ip4
    end
    
end #populate_face_el

function populate_conn_face_el!(mesh::St_mesh, SD::NSD_3D)
    
    for iel = 1:mesh.nelem
        #
        # CGNS numbering
        #
        ip1 = mesh.connijk[iel, mesh.ngl, 1, 1]
        ip2 = mesh.connijk[iel, mesh.ngl, mesh.ngl, 1]
        ip3 = mesh.connijk[iel, 1, mesh.ngl, 1]
        ip4 = mesh.connijk[iel, 1, 1, 1]
        ip5 = mesh.connijk[iel, mesh.ngl, 1, mesh.ngl]
        ip6 = mesh.connijk[iel, mesh.ngl, mesh.ngl, mesh.ngl]
        ip7 = mesh.connijk[iel, 1, mesh.ngl, mesh.ngl]
        ip8 = mesh.connijk[iel, 1, 1, mesh.ngl]
        
        #
        # Local faces node connectivity:
        # i.e. what nodes belong to a given local face in iel:
        #
        face_el = 1
        mesh.conn_face_el[1, face_el, iel] = ip1
        mesh.conn_face_el[2, face_el, iel] = ip4
        mesh.conn_face_el[3, face_el, iel] = ip3
        mesh.conn_face_el[4, face_el, iel] = ip2

        face_el = 2
        mesh.conn_face_el[1, face_el, iel] = ip1
        mesh.conn_face_el[2, face_el, iel] = ip2
        mesh.conn_face_el[3, face_el, iel] = ip6
        mesh.conn_face_el[4, face_el, iel] = ip5
        
        face_el = 3
        mesh.conn_face_el[1, face_el, iel] = ip2
        mesh.conn_face_el[2, face_el, iel] = ip3
        mesh.conn_face_el[3, face_el, iel] = ip7
        mesh.conn_face_el[4, face_el, iel] = ip6
        
        face_el = 4
        mesh.conn_face_el[1, face_el, iel] = ip3
        mesh.conn_face_el[2, face_el, iel] = ip4
        mesh.conn_face_el[3, face_el, iel] = ip8
        mesh.conn_face_el[4, face_el, iel] = ip7
        
        face_el = 5
        mesh.conn_face_el[1, face_el, iel] = ip1
        mesh.conn_face_el[2, face_el, iel] = ip5
        mesh.conn_face_el[3, face_el, iel] = ip8
        mesh.conn_face_el[4, face_el, iel] = ip4
        
        face_el = 6
        mesh.conn_face_el[1, face_el, iel] = ip5
        mesh.conn_face_el[2, face_el, iel] = ip6
        mesh.conn_face_el[3, face_el, iel] = ip7
        mesh.conn_face_el[4, face_el, iel] = ip8
    end
    
end #populate_face_el


function  add_high_order_nodes!(mesh::St_mesh) end

function  add_high_order_nodes_1D_native_mesh!(mesh::St_mesh, interpolation_nodes, backend)
    
    if (mesh.nop < 2) return end
    
    println(" # POPULATE 1D GRID with SPECTRAL NODES ............................ ")
    println(" # ...")
    
    lgl = basis_structs_ξ_ω!(interpolation_nodes, mesh.nop, backend)
    
    x1, x2 = TFloat(0.0), TFloat(0.0)    
    ξ::typeof(lgl.ξ[1]) = 0.0

    ngl                      = mesh.nop + 1
    tot_linear_poin          = mesh.npoin_linear
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)
    el_internal_nodes        = ngl - 2
    el_nodes                 = ngl
    
    #Increase number of grid points from linear count to total high-order points
    mesh.npoin = mesh.npoin_linear + tot_vol_internal_nodes
    resize!(mesh.x, (mesh.npoin))
    
    mesh.connijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.ngl), 1, 1)

    #
    # First pass: build coordinates and store IP into poin_in_edge[iedge_g, l]
    #
    ip = tot_linear_poin + 1
    for iel_g = 1:mesh.nelem

        ip1 = iel_g
        ip2 = iel_g + 1
        
        mesh.conn[iel_g, 1], mesh.conn[iel_g, ngl] = ip1, ip2
        mesh.connijk[iel_g, 1, 1, 1], mesh.connijk[iel_g, ngl, 1, 1] = ip1, ip2
        x1, x2 = mesh.x[ip1], mesh.x[ip2]
        
        iconn = 1
        for l=2:ngl-1
            ξ = lgl.ξ[l];
            
            mesh.x[ip] = x1*(1.0 - ξ)*0.5 + x2*(1.0 + ξ)*0.5;
            
            mesh.conn[iel_g, l] = ip #OK
            mesh.connijk[iel_g, l, 1, 1] = ip #OK
            iconn = iconn + 1
            
            ip = ip + 1
        end
    end
    
    println(" # POPULATE 1D GRID with SPECTRAL NODES ............................ DONE")
    return 
end


function  add_high_order_nodes_edges!(mesh::St_mesh, lgl, SD::NSD_2D, backend, edge2pedge)
    
    if (mesh.nop < 2) return end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ EDGES"; msg_rank = rank, suppress = mesh.msg_suppress)
    println_rank(" # ..."; msg_rank = rank, suppress = mesh.msg_suppress)
    
    x1, y1 = TFloat(0.0), TFloat(0.0)
    x2, y2 = TFloat(0.0), TFloat(0.0)
    
    ξ::typeof(lgl.ξ[1]) = 0.0

    ngl                      = mesh.nop + 1
    tot_linear_poin          = mesh.npoin_linear
    gtot_linear_poin          = mesh.gnpoin_linear
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)*(ngl-2)
    el_edges_internal_nodes  = mesh.NEDGES_EL*(ngl-2)
    
    #Increase number of grid points from linear count to total high-order points
    mesh.npoin = mesh.npoin_linear + tot_edges_internal_nodes + tot_vol_internal_nodes

    if length(mesh.x_ho) < mesh.npoin
        #resize!(mesh.x_ho, (mesh.npoin))
        mesh.x_ho = KernelAbstractions.allocate(backend, TFloat, mesh.npoin)
    end
    if length(mesh.y_ho) < mesh.npoin        
        #resize!(mesh.y_ho, (mesh.npoin))
       mesh.y_ho = KernelAbstractions.allocate(backend, TFloat, mesh.npoin)
    end
    
    #poin_in_edge::Array{TInt, 2}  = zeros(mesh.nedges, mesh.ngl)
    # open("./COORDS_HO_edges_$rank.dat", "w") do f
        #
        # First pass: build coordinates and store IP into poin_in_edge[iedge_g, l]
        #
        ip = tot_linear_poin + 1
        gip::Int64 = 0
        for iedge_g = 1:mesh.nedges
            
            ip1 = mesh.conn_unique_edges[iedge_g][1]
            ip2 = mesh.conn_unique_edges[iedge_g][2]
            gip = gtot_linear_poin + 1 + (edge2pedge[iedge_g] - 1) * (ngl - 2)
            
            mesh.poin_in_edge[iedge_g,        1] = ip1
            mesh.poin_in_edge[iedge_g, mesh.ngl] = ip2
            
            x1, y1 = mesh.x[ip1], mesh.y[ip1]
            x2, y2 = mesh.x[ip2], mesh.y[ip2]
            
            #@printf(" %d: (ip1, ip2) = (%d %d) ", iedge_g, ip1, ip2)
            for l=2:ngl-1
                ξ = lgl.ξ[l];
                
                mesh.x_ho[ip] = x1*(1.0 - ξ)*0.5 + x2*(1.0 + ξ)*0.5;
	            mesh.y_ho[ip] = y1*(1.0 - ξ)*0.5 + y2*(1.0 + ξ)*0.5;
                
                mesh.poin_in_edge[iedge_g, l] = ip
                mesh.ip2gip[ip] = gip
                # mesh.gip2owner[ip] = 1
                
                #@printf(" lgl %d: %d %d ", l, iedge_g, mesh.poin_in_edge[iedge_g, l])
            #    @printf(f, " %.6f %.6f 0.000000 %d %d %d\n", mesh.x_ho[ip],  mesh.y_ho[ip], ip, gip, edge2pedge[iedge_g])
                ip  = ip + 1
                gip = gip + 1
            end
        end
    # end #do f
    #show(stdout, "text/plain", poin_in_edge)
    #@info "-----2D edges"
    
    #
    # Second pass: populate mesh.conn[∀ elem, 1:4+el_edges_internal_nodes]\n")
    #
    cache_edge_ids = array_cache(mesh.cell_edge_ids) # allocation here  
    for iel = 1:mesh.nelem
        edge_ids = getindex!(cache_edge_ids, mesh.cell_edge_ids, iel)
        #show(stdout, "text/plain",edge_ids)
        iconn = 1

        iedge_el = 1
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        # @info ip1, ip2, iedge_el, iedge_g, iel
        if (mesh.conn[iel,1] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel, 2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
        end

        iedge_el = 4
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        # @info ip1, ip2, iedge_el, iedge_g, iel
        if (mesh.conn[iel,2] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end 
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
        end

        iedge_el = 2
        iedge_g = edge_ids[iedge_el] 
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        # @info ip1, ip2, iedge_el, iedge_g, iel
        if (mesh.conn[iel,3] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
        end
        iedge_el = 3
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        # @info ip1, ip2, iedge_el, iedge_g, iel
        if (mesh.conn[iel,4] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end 
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
        end
    end
    #show(stdout, "text/plain", mesh.conn')
    
    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ EDGES DONE"; msg_rank = rank, suppress = mesh.msg_suppress)
    
    return 
end


function  add_high_order_nodes_edges!(mesh::St_mesh, lgl, SD::NSD_3D, backend, edge2pedge)
    
    if (mesh.nop < 2) return end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ EDGES"; msg_rank = rank, suppress = mesh.msg_suppress)
    println_rank(" # ..."; msg_rank = rank, suppress = mesh.msg_suppress)
    
    x1, y1, z1 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x2, y2, z2 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    
    ξ::typeof(lgl.ξ[1]) = 0.0

    ngl                      = mesh.nop + 1
    tot_linear_poin          = mesh.npoin_linear
    gtot_linear_poin          = mesh.gnpoin_linear
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    tot_faces_internal_nodes = mesh.nfaces*(ngl-2)*(ngl-2)
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)*(ngl-2)*(ngl-2)
    el_edges_internal_nodes  = mesh.NEDGES_EL*(ngl-2)
    
    #Increase number of grid points from linear count to total high-order points
    mesh.npoin = mesh.npoin_linear + tot_edges_internal_nodes + tot_faces_internal_nodes + tot_vol_internal_nodes

    if length(mesh.x_ho) < mesh.npoin
        resize!(mesh.x_ho, (mesh.npoin))
    end
    if length(mesh.y_ho) < mesh.npoin        
        resize!(mesh.y_ho, (mesh.npoin))
    end
    if length(mesh.z_ho) < mesh.npoin        
        resize!(mesh.z_ho, (mesh.npoin))
    end
    #
    # CGNS numbering
    #
    edge_g_color::Array{Int64, 1} = zeros(Int64, mesh.nedges)
    #poin_in_edge::Array{Int64, 2}  = zeros(mesh.nedges, mesh.ngl)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    open("./COORDS_HO_edges_$rank.dat", "w") do f
    # open("./COORDS_HO_edges.dat", "w") do f
        #
        # First pass: build coordinates and store IP into poin_in_edge[iedge_g, l]
        #
        ip = tot_linear_poin + 1
        gip::Int64 = 0
        for iedge_g = 1:mesh.nedges
            #
            # First pass: build coordinates and store IP into poin_in_edge[iedge_g, l]
            #
            ip1 = mesh.conn_unique_edges[iedge_g][1]
            ip2 = mesh.conn_unique_edges[iedge_g][2]
            gip = gtot_linear_poin + 1 + (edge2pedge[iedge_g] - 1) * (ngl - 2)
            #ip1 = mesh.conn_edge_el[1, iedge_el, iel]
            #ip2 = mesh.conn_edge_el[2, iedge_el, iel]
            mesh.poin_in_edge[iedge_g,        1] = ip1
            mesh.poin_in_edge[iedge_g, mesh.ngl] = ip2
            
            x1, y1, z1 = mesh.x[ip1], mesh.y[ip1], mesh.z[ip1]
            x2, y2, z2 = mesh.x[ip2], mesh.y[ip2], mesh.z[ip2]
            
            #@printf(" iedge_g=%d -> (ip1, ip2) = (%d %d) \n", iedge_g, ip1, ip2)
            for l=2:ngl-1
                ξ = lgl.ξ[l];
                
                mesh.x_ho[ip] = x1*(1.0 - ξ)*0.5 + x2*(1.0 + ξ)*0.5;
                mesh.y_ho[ip] = y1*(1.0 - ξ)*0.5 + y2*(1.0 + ξ)*0.5;
                mesh.z_ho[ip] = z1*(1.0 - ξ)*0.5 + z2*(1.0 + ξ)*0.5;
                
                mesh.poin_in_edge[iedge_g, l] = ip
                mesh.ip2gip[ip] = gip
                
                #@printf(" lgl %d: %d %d ", l, iedge_g, mesh.poin_in_edge[iedge_g, l])
                @printf(f, " %.6f %.6f %.6f %d %d\n", mesh.x_ho[ip],  mesh.y_ho[ip], mesh.z_ho[ip], ip, gip)
                #@printf( " %.6f %.6f %.6f %d\n", mesh.x_ho[ip],  mesh.y_ho[ip], mesh.z_ho[ip], ip)
                ip  = ip + 1
                gip = gip + 1
            end
        end
    end #end f
    #show(stdout, "text/plain", mesh.poin_in_edge)
    #@info "-----3D edges"
        
    #
    # Second pass: populate mesh.conn[1:8+el_edges_internal_nodes, ∀ elem]\n")
    #
    cache_edge_ids = array_cache(mesh.cell_edge_ids) # allocation here
    for iel = 1:mesh.nelem
        edge_ids = getindex!(cache_edge_ids, mesh.cell_edge_ids, iel)
        iconn = 1
        iedge_el = 10
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,1] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,1,l,1] = ip #OK
        end
        iedge_el = 8
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,2] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel, 2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,1,ngl,l] = ip
        end
        iedge_el = 12
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,3] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,1,l,ngl] = ip
        end
        iedge_el = 6
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,4] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,1,1,l] = ip
        end
        iedge_el = 1
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,1] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,ngl-l+1,1,1] = ip
        end
        iedge_el = 3
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,2] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,ngl-l+1,ngl,1] = ip
        end
        iedge_el = 4
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,3] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,ngl-l+1,ngl,ngl] = ip
        end
        iedge_el = 2
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,4] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,ngl-l+1,1,ngl] = ip
        end
        iedge_el = 9
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,5] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,ngl,l,1] = ip
        end
        iedge_el = 7
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,6] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,ngl,ngl,l] = ip
        end
        iedge_el = 11
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,7] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,ngl,l,ngl] = ip
        end
        iedge_el = 5
        iedge_g = edge_ids[iedge_el]
        ip1 = mesh.conn_unique_edges[iedge_g][1]
        ip2 = mesh.conn_unique_edges[iedge_g][2]
        if (mesh.conn[iel,4] == ip1)
            starter = 2
            ender = ngl-1
            stepper =1
        else
            starter = ngl-1
            ender = 2
            stepper =-1
        end
        for l = starter:stepper:ender
            ip = mesh.poin_in_edge[iedge_g, l]
            mesh.conn[iel,2^mesh.nsd + iconn] = ip #OK
            iconn = iconn + 1
            mesh.connijk[iel,ngl,1,l] = ip
        end
    end
    #show(stdout, "text/plain", mesh.conn')


    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ EDGES DONE"; msg_rank = rank, suppress = mesh.msg_suppress)
    return 
end


function  add_high_order_nodes_faces!(mesh::St_mesh, lgl, SD::NSD_2D, face2pface)

    if (mesh.nop < 2) return end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ FACES"; msg_rank = rank, suppress = mesh.msg_suppress)
    println_rank(" # ..."; msg_rank = rank, suppress = mesh.msg_suppress)
    
    x1, y1 = TFloat(0.0), TFloat(0.0)
    x2, y2 = TFloat(0.0), TFloat(0.0)
    x3, y3 = TFloat(0.0), TFloat(0.0)
    x4, y4 = TFloat(0.0), TFloat(0.0)
    
    ξ::typeof(lgl.ξ[1]) = 0.0
    ζ::typeof(lgl.ξ[1]) = 0.0
    
    ngl                      = mesh.nop + 1
    tot_linear_poin          = mesh.npoin_linear
    gtot_linear_poin          = mesh.gnpoin_linear
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    gtot_edges_internal_nodes = mesh.gnedges*(ngl-2)
    tot_faces_internal_nodes = mesh.nfaces*(ngl-2)*(ngl-2)
    
    el_edges_internal_nodes  = mesh.NEDGES_EL*(ngl-2)
    el_faces_internal_nodes  = mesh.NFACES_EL*(ngl-2)*(ngl-2)
    
    #Increase number of grid points from linear count to total high-order points
    mesh.npoin = mesh.npoin_linear + tot_edges_internal_nodes + tot_faces_internal_nodes
    
    if length(mesh.x_ho) < mesh.npoin
        setize!(mesh.x_ho, (mesh.npoin))
    end
    if length(mesh.y_ho) < mesh.npoin
        resize!(mesh.y_ho, (mesh.npoin))
    end

    #open("./COORDS_HO_faces.dat", "w") do f
        #
        # First pass:
        #
        ip  = tot_linear_poin + tot_edges_internal_nodes + 1
        gip::Int64  = 0
        for iface_g = 1:mesh.nelem #NOTICE: in 2D the faces are the elements themselves
            iel = iface_g
            gip = gtot_linear_poin + gtot_edges_internal_nodes + 1 + (face2pface[iel] - 1) * (ngl-2) * (ngl - 2)
            #GGNS numbering
            ip1 = mesh.cell_node_ids[iel][1]
            ip2 = mesh.cell_node_ids[iel][2]
            ip3 = mesh.cell_node_ids[iel][4]
            ip4 = mesh.cell_node_ids[iel][3]

            mesh.poin_in_face[iface_g, 1, 1]     = ip1
            mesh.poin_in_face[iface_g, ngl, 1]   = ip2
            mesh.poin_in_face[iface_g, ngl, ngl] = ip4
            mesh.poin_in_face[iface_g, 1, ngl]   = ip3
            
            x1, y1 = mesh.x[ip1], mesh.y[ip1]
            x2, y2 = mesh.x[ip2], mesh.y[ip2]
            x3, y3 = mesh.x[ip3], mesh.y[ip3]
            x4, y4 = mesh.x[ip4], mesh.y[ip4]
            
            for l=2:ngl-1
                ξ = lgl.ξ[l];
                
                for m=2:ngl-1
                    ζ = lgl.ξ[m];
                    
                    mesh.x_ho[ip] = (x1*(1 - ξ)*(1 - ζ)*0.25
                                        + x2*(1 + ξ)*(1 - ζ)*0.25
                                + x3*(1 + ξ)*(1 + ζ)*0.25			
                                + x4*(1 - ξ)*(1 + ζ)*0.25)
                    
                    mesh.y_ho[ip] =  (y1*(1 - ξ)*(1 - ζ)*0.25
		                      + y2*(1 + ξ)*(1 - ζ)*0.25
		                      + y3*(1 + ξ)*(1 + ζ)*0.25
		                      + y4*(1 - ξ)*(1 + ζ)*0.25)

                    mesh.poin_in_face[iface_g, l, m] = ip
                    #NEW ORDERING
                    mesh.connijk[iel, m, ngl-l+1] = ip

                    mesh.ip2gip[ip] = gip
                    #OLD ORDERING
                    #mesh.connijk[iel, m, l] = ip
      #              @printf(f, " %.6f %.6f 0.000000 %d\n", mesh.x_ho[ip],  mesh.y_ho[ip], ip)
                    
	                ip = ip + 1
	                gip = gip + 1
                end
            end
        end
    #end #do f

    #
    # Second pass: populate mesh.conn[1:8+el_edges_internal_nodes+el_faces_internal_nodes, ∀ elem]\n")
    #
    cache_face_ids = array_cache(mesh.cell_face_ids) # allocation here    
    for iel = 1:mesh.nelem
        iface_g = iel
        iconn = 1
        iterate = 1
        starter = 2
        ender = ngl-1
        while (iconn <= (ngl-2)^2)
            m=starter
            for l=starter:ender
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn] = ip
                iconn = iconn + 1
            end
            if (iconn > (ngl-2)^2) break end 
            l=ender
            for m=starter+1:ender
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn] = ip
                iconn = iconn + 1
            end
            if (iconn > (ngl-2)^2) break end
            m=ender
            for l=ender-1:-1:starter
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn] = ip
                iconn = iconn + 1
            end
            if (iconn > (ngl-2)^2) break end
            l=starter
            for m=ender-1:-1:starter+1
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn] = ip
                iconn = iconn + 1
            end
            starter = starter+1
            ender = ender-1 
        end 
        
    end
    
    for iel = 1:mesh.nelem
        iter =1
        for m=ngl-1:-1:2
            mesh.connijk[iel,1,m] = mesh.conn[iel,4+iter]
            iter = iter+1
        end
        for m=2:ngl-1
            mesh.connijk[iel,m,1] = mesh.conn[iel,4+iter]
            iter = iter+1
        end
        for m=2:ngl-1
            mesh.connijk[iel,ngl,m] = mesh.conn[iel,4+iter]
            iter = iter+1
        end
        for m=ngl-1:-1:2
            mesh.connijk[iel,m,ngl] = mesh.conn[iel,4+iter]
            iter = iter+1
        end
    end
    for iel =1:mesh.nelem
        #      show(stdout, "text/plain", mesh.connijk[iel,:,:]')
    end

    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ FACES DONE"; msg_rank = rank, suppress = mesh.msg_suppress)

end

function  add_high_order_nodes_faces!(mesh::St_mesh, lgl, SD::NSD_3D, face2pface)

    if (mesh.nop < 2) return end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ FACES"; msg_rank = rank, suppress = mesh.msg_suppress)
    
    x1, y1, z1 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x2, y2, z2 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x3, y3, z3 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x4, y4, z4 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    
    ξ::typeof(lgl.ξ[1]) = 0.0
    ζ::typeof(lgl.ξ[1]) = 0.0
    
    ngl                      = mesh.nop + 1
    tot_linear_poin          = mesh.npoin_linear
    gtot_linear_poin          = mesh.gnpoin_linear
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    gtot_edges_internal_nodes = mesh.gnedges*(ngl-2)
    tot_faces_internal_nodes = mesh.nfaces*(ngl-2)*(ngl-2)
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)*(ngl-2)*(ngl-2)

    el_edges_internal_nodes  = mesh.NEDGES_EL*(ngl-2)
    el_faces_internal_nodes  = mesh.NFACES_EL*(ngl-2)*(ngl-2)
    
    #Increase number of grid points from linear count to total high-order points
    mesh.npoin = mesh.npoin_linear + tot_edges_internal_nodes + tot_faces_internal_nodes + tot_vol_internal_nodes
    

    if length(mesh.x_ho) < mesh.npoin
        resize!(mesh.x_ho, (mesh.npoin))
    end
    if length(mesh.y_ho) < mesh.npoin
        resize!(mesh.y_ho, (mesh.npoin))
    end
    if length(mesh.z_ho) < mesh.npoin
        resize!(mesh.z_ho, (mesh.npoin))
    end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    open("./COORDS_HO_faces_$rank.dat", "w") do f
    #open("./COORDS_HO_faces.dat", "w") do f
        #
        # First pass:
        #
        ip  = tot_linear_poin + tot_edges_internal_nodes + 1
        gip::Int64  = 0
        for iface_g = 1:mesh.nfaces
            gip = gtot_linear_poin + gtot_edges_internal_nodes + 1 + (face2pface[iface_g] - 1) * (ngl-2) * (ngl - 2)
            
            #GGNS numbering
            ip1 = mesh.conn_unique_faces[iface_g][1]
            ip2 = mesh.conn_unique_faces[iface_g][2]
            ip3 = mesh.conn_unique_faces[iface_g][4]
            ip4 = mesh.conn_unique_faces[iface_g][3]

            mesh.poin_in_face[iface_g, 1, 1]     = ip1
            mesh.poin_in_face[iface_g, ngl, 1]   = ip2
            mesh.poin_in_face[iface_g, ngl, ngl] = ip3#ip4
            mesh.poin_in_face[iface_g, 1, ngl]   = ip4#ip3
            
            x1, y1, z1 = mesh.x[ip1], mesh.y[ip1], mesh.z[ip1]
            x2, y2, z2 = mesh.x[ip2], mesh.y[ip2], mesh.z[ip2]
            x3, y3, z3 = mesh.x[ip3], mesh.y[ip3], mesh.z[ip3]
            x4, y4, z4 = mesh.x[ip4], mesh.y[ip4], mesh.z[ip4]
            

            for l=2:ngl-1
                ξ = lgl.ξ[l];
                
                for m=2:ngl-1
                    ζ = lgl.ξ[m];
                    
	            mesh.x_ho[ip] =  (x1*(1 - ξ)*(1 - ζ)*0.25
                                      + x2*(1 + ξ)*(1 - ζ)*0.25
		                      + x3*(1 + ξ)*(1 + ζ)*0.25			
		                      + x4*(1 - ξ)*(1 + ζ)*0.25)
                    
                    mesh.y_ho[ip] =  (y1*(1 - ξ)*(1 - ζ)*0.25
                                    + y2*(1 + ξ)*(1 - ζ)*0.25
                                    + y3*(1 + ξ)*(1 + ζ)*0.25
                                    + y4*(1 - ξ)*(1 + ζ)*0.25)
                    
                    mesh.z_ho[ip] =  (z1*(1 - ξ)*(1 - ζ)*0.25
                                    + z2*(1 + ξ)*(1 - ζ)*0.25
                                    + z3*(1 + ξ)*(1 + ζ)*0.25
                                    + z4*(1 - ξ)*(1 + ζ)*0.25)

                    mesh.poin_in_face[iface_g, l, m] = ip
                    mesh.ip2gip[ip] = gip
                    
                    @printf(f, " %.6f %.6f %.6f %d %d\n", mesh.x_ho[ip],  mesh.y_ho[ip], mesh.z_ho[ip], ip, gip)

                    # mesh.connijk[iel, m, ngl-l+1] = ip #<==== need to build this 
                    
	                ip = ip + 1
                    gip = gip + 1
                end
            end
        end
    end #do f

    #
    # Second pass: populate mesh.conn[1:8+el_edges_internal_nodes+el_faces_internal_nodes, ∀ elem]\n")
    #
    cache_face_ids = array_cache(mesh.cell_face_ids) # allocation here    
    for iel = 1:mesh.nelem
        iconn = 0
        face_ids = getindex!(cache_face_ids, mesh.cell_face_ids, iel)
        iterate = 1
        starter = 2
        ender = ngl-1
        iface_el = 6 # left
        iface_g = face_ids[iface_el]
        iconn_face =1
        while (iconn_face <= (ngl-2)^2)
            m=starter
            for l=starter:ender
                ip = mesh.poin_in_face[iface_g, m, l]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn_face] = ip
                iconn_face=iconn_face+1
                mesh.connijk[iel,1,l,m] = ip # OK for nop=3. need to check nop>3 
            end
            if (iconn_face > (ngl-2)^2) break end
            l=ender
            for m=starter+1:ender
                ip = mesh.poin_in_face[iface_g, m, l]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,1,l,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            m=ender
            for l=ender-1:-1:starter
                ip = mesh.poin_in_face[iface_g, m, l]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,1,l,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            l=starter
            for m=ender-1:-1:starter+1
                ip = mesh.poin_in_face[iface_g, m, l]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,1,l,m] = ip
            end
            starter = starter+1
            ender = ender-1
        end
        iconn = iconn+iconn_face-1
        iconn_face = 1 
        iterate = 1
        starter = 2
        ender = ngl-1
        iface_el = 3 # bottom
        iface_g = face_ids[iface_el]
        while (iconn_face <= (ngl-2)^2)
            l=ender
            for m=starter:ender
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,m,1] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            m=ender
            for l=ender-1:-1:starter
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,m,1] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            l=starter
            for m=ender-1:-1:starter
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn+iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,m,1] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            m=starter
            for l=starter+1:ender-1
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,m,1] = ip
            end
            starter = starter+1
            ender = ender-1
        end
        iconn = iconn+iconn_face-1
        iconn_face = 1 
        iterate = 1
        starter = 2
        ender = ngl-1
        iface_el = 2 # back
        iface_g = face_ids[iface_el]
        while (iconn_face <= (ngl-2)^2)
            l=ender
            for m=starter:ender
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,ngl,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            m=ender
            for l=ender-1:-1:starter
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,ngl,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            l=starter
            for m=ender-1:-1:starter
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn+iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,ngl,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            m=starter
            for l=starter+1:ender-1
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,ngl,m] = ip
            end
            starter = starter+1
            ender = ender-1
        end
        iconn = iconn+iconn_face-1
        iconn_face = 1 
        iterate = 1
        starter = 2
        ender = ngl-1
        iface_el = 4 # top
        iface_g = face_ids[iface_el]
        while (iconn_face <= (ngl-2)^2)
            l=ender
            for m=ender:-1:starter
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,m,ngl] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            m=starter
            for l=ender-1:-1:starter
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,m,ngl] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            l=starter
            for m=starter+1:ender
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,m,ngl] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            m=ender
            for l=starter+1:ender-1
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn+iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,m,ngl] = ip
            end
            starter = starter+1
            ender = ender-1
        end
        iconn = iconn+iconn_face-1
        iconn_face = 1 
        iterate = 1
        starter = 2
        ender = ngl-1
        iface_el = 1 # front
        iface_g = face_ids[iface_el]
        while (iconn_face <= (ngl-2)^2)
            l=ender
            for m=ender:-1:starter
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,1,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            m=starter
            for l=ender-1:-1:starter
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,1,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            l=starter
            for m=starter+1:ender
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,1,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            m=ender
            for l=starter+1:ender-1
                ip = mesh.poin_in_face[iface_g, l, m]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn+iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl-l+1,1,m] = ip
            end
            starter = starter+1
            ender = ender-1
        end
        iconn = iconn+iconn_face-1
        iconn_face = 1 
        iterate = 1
        starter = 2
        ender = ngl-1
        iface_el = 5 # right
        iface_g = face_ids[iface_el]
        while (iconn_face <= (ngl-2)^2)
            m=starter
            for l=starter:ender
                ip = mesh.poin_in_face[iface_g, m, l]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl,l,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            l=ender
            for m=starter+1:ender
                ip = mesh.poin_in_face[iface_g, m, l]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl,l,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            m=ender
            for l=ender-1:-1:starter
                ip = mesh.poin_in_face[iface_g, m, l]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn + iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl,l,m] = ip
            end
            if (iconn_face > (ngl-2)^2) break end
            l=starter
            for m=ender-1:-1:starter+1
                ip = mesh.poin_in_face[iface_g, m, l]
                mesh.conn[iel, 2^mesh.nsd + el_edges_internal_nodes + iconn+iconn_face] = ip
                iconn_face = iconn_face + 1
                mesh.connijk[iel,ngl,l,m] = ip
            end
            starter = starter+1
            ender = ender-1
        end
        #=istart=2^mesh.nsd + el_edges_internal_nodes
        iter =1 
        iconn_face =1
        starter = 2
        ender = ngl-1 =#
        #=for iface_el = 1:length(face_ids)
        iface_g = face_ids[iface_el]
        for l = 2:ngl-1
        for m = 2:ngl-1
        ip = mesh.poin_in_face[iface_g, l, m]
        mesh.conn[ie, 8 + el_edges_internal_nodes + iconn] = ip
        iconn = iconn + 1
        end
        end
        end=#
    end
    #show(stdout, "text/plain", mesh.conn')
    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ FACES DONE"; msg_rank = rank, suppress = mesh.msg_suppress)

end

function  add_high_order_nodes_volumes!(mesh::St_mesh, lgl, SD::NSD_2D, elm2pelm)
    nothing
end

function  add_high_order_nodes_volumes!(mesh::St_mesh, lgl, SD::NSD_3D, elm2pelm)

    if (mesh.nop < 2) return end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ VOLUMES"; msg_rank = rank, suppress = mesh.msg_suppress)
    println_rank(" # ..."; msg_rank = rank, suppress = mesh.msg_suppress)
    
    x1, y1, z1 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x2, y2, z2 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x3, y3, z3 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x4, y4, z4 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x5, y5, z5 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x6, y6, z6 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x7, y7, z7 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    x8, y8, z8 = TFloat(0.0), TFloat(0.0), TFloat(0.0)
    
    ξ::typeof(lgl.ξ[1]) = 0.0
    η::typeof(lgl.ξ[1]) = 0.0
    ζ::typeof(lgl.ξ[1]) = 0.0

    ngl                      = mesh.nop + 1
    tot_linear_poin          = mesh.npoin_linear
    gtot_linear_poin          = mesh.gnpoin_linear
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    gtot_edges_internal_nodes = mesh.gnedges*(ngl-2)
    tot_faces_internal_nodes = mesh.nfaces*(ngl-2)*(ngl-2)
    gtot_faces_internal_nodes = mesh.gnfaces*(ngl-2)*(ngl-2)
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)*(ngl-2)*(ngl-2)

    el_edges_internal_nodes  = mesh.NEDGES_EL*(ngl-2)
    el_faces_internal_nodes  = mesh.NFACES_EL*(ngl-2)*(ngl-2)
    el_vol_internal_nodes    = (ngl-2)*(ngl-2)*(ngl-2)
    conn_vol_poin::Array{TInt, 4}  = zeros(mesh.ngl, mesh.ngl, mesh.ngl,mesh.nelem)
    #Increase number of grid points from linear count to titak high-order points
    mesh.npoin = tot_linear_poin + tot_edges_internal_nodes + tot_faces_internal_nodes + tot_vol_internal_nodes
    
    if length(mesh.x_ho) < mesh.npoin
        resize!(mesh.x_ho, (mesh.npoin))
    end
    if length(mesh.y_ho) < mesh.npoin
        resize!(mesh.y_ho, (mesh.npoin))
    end
    if length(mesh.z_ho) < mesh.npoin
        resize!(mesh.z_ho, (mesh.npoin))
    end

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    gip::Int64  = 0
    open("./COORDS_HO_faces_$rank.dat", "w") do f
    #open("./COORDS_HO_vol.dat", "w") do f
        ip  = tot_linear_poin + tot_edges_internal_nodes + tot_faces_internal_nodes + 1
        for iel = 1:mesh.nelem
            gip = gtot_linear_poin + gtot_edges_internal_nodes + gtot_faces_internal_nodes + 1 + (elm2pelm[iel] - 1) * (ngl-2) * (ngl - 2) * (ngl - 2)
            iconn = 1
            
            #
            # CGNS numbering
            #
            ip1 = mesh.cell_node_ids[iel][2]
            ip2 = mesh.cell_node_ids[iel][6]
            ip3 = mesh.cell_node_ids[iel][8]
            ip4 = mesh.cell_node_ids[iel][4]
            ip5 = mesh.cell_node_ids[iel][1]
            ip6 = mesh.cell_node_ids[iel][5]
            ip7 = mesh.cell_node_ids[iel][7]
            ip8 = mesh.cell_node_ids[iel][3]
            
            x1, y1, z1 = mesh.x[ip1], mesh.y[ip1], mesh.z[ip1]
            x2, y2, z2 = mesh.x[ip2], mesh.y[ip2], mesh.z[ip2]
            x3, y3, z3 = mesh.x[ip3], mesh.y[ip3], mesh.z[ip3]
            x4, y4, z4 = mesh.x[ip4], mesh.y[ip4], mesh.z[ip4]     
            x5, y5, z5 = mesh.x[ip5], mesh.y[ip5], mesh.z[ip5]
            x6, y6, z6 = mesh.x[ip6], mesh.y[ip6], mesh.z[ip6]
            x7, y7, z7 = mesh.x[ip7], mesh.y[ip7], mesh.z[ip7]
            x8, y8, z8 = mesh.x[ip8], mesh.y[ip8], mesh.z[ip8]
            
            for l=2:ngl-1
                ξ = lgl.ξ[l];
                
                for m=2:ngl-1
                    η = lgl.ξ[m];
                    
                    for n=2:ngl-1
                        ζ = lgl.ξ[n];
                        
	                mesh.x_ho[ip] = (x1*(1 - ξ)*(1 - η)*(1 - ζ)*0.125
			                 + x2*(1 + ξ)*(1 - η)*(1 - ζ)*0.125
			                 + x3*(1 + ξ)*(1 + η)*(1 - ζ)*0.125
			                 + x4*(1 - ξ)*(1 + η)*(1 - ζ)*0.125
			                 + x5*(1 - ξ)*(1 - η)*(1 + ζ)*0.125
			                 + x6*(1 + ξ)*(1 - η)*(1 + ζ)*0.125
			                 + x7*(1 + ξ)*(1 + η)*(1 + ζ)*0.125
			                 + x8*(1 - ξ)*(1 + η)*(1 + ζ)*0.125)
                        
	                mesh.y_ho[ip] = (y1*(1 - ξ)*(1 - η)*(1 - ζ)*0.125
			                 + y2*(1 + ξ)*(1 - η)*(1 - ζ)*0.125
			                 + y3*(1 + ξ)*(1 + η)*(1 - ζ)*0.125
			                 + y4*(1 - ξ)*(1 + η)*(1 - ζ)*0.125
			                 + y5*(1 - ξ)*(1 - η)*(1 + ζ)*0.125
			                 + y6*(1 + ξ)*(1 - η)*(1 + ζ)*0.125
			                 + y7*(1 + ξ)*(1 + η)*(1 + ζ)*0.125
			                 + y8*(1 - ξ)*(1 + η)*(1 + ζ)*0.125)
                        
	                mesh.z_ho[ip] = (z1*(1 - ξ)*(1 - η)*(1 - ζ)*0.125
			                 + z2*(1 + ξ)*(1 - η)*(1 - ζ)*0.125
			                 + z3*(1 + ξ)*(1 + η)*(1 - ζ)*0.125
			                 + z4*(1 - ξ)*(1 + η)*(1 - ζ)*0.125
			                 + z5*(1 - ξ)*(1 - η)*(1 + ζ)*0.125
			                 + z6*(1 + ξ)*(1 - η)*(1 + ζ)*0.125
			                 + z7*(1 + ξ)*(1 + η)*(1 + ζ)*0.125
			                 + z8*(1 - ξ)*(1 + η)*(1 + ζ)*0.125)

                        #mesh.conn[iel, 8 + el_edges_internal_nodes + el_faces_internal_nodes + iconn] = ip
                        conn_vol_poin[l,m,n,iel] = ip
                        mesh.connijk[iel,n,l,m] = ip
                        mesh.ip2gip[ip] = gip

                        @printf(f, " %.6f %.6f %.6f %d %d\n", mesh.x_ho[ip],  mesh.y_ho[ip], mesh.z_ho[ip], ip, gip)

                        ip = ip + 1
                        gip = gip + 1
                        iconn = iconn + 1
                    end
                end
            end
        end
    end # do f 
    #open("./CONNIJK.dat", "w") do f
        for iel =1:mesh.nelem
            iconn =1
            for n=2:ngl-1
                starter=2
                ender = ngl-1
                iconn_level=1
                while (iconn_level <= (ngl-2)^2)
                    m=starter
                    for l=starter:ender
                        ip = conn_vol_poin[l,m,n,iel] 
                        mesh.conn[iel, 8 + el_edges_internal_nodes + el_faces_internal_nodes + iconn] = ip
                        iconn=iconn+1
                        iconn_level=iconn_level+1
                    end 
                    if (iconn_level > (ngl-2)^2) break end
                    l=ender
                    for m=starter+1:ender
                        ip = conn_vol_poin[l,m,n,iel]
                        mesh.conn[iel, 8 + el_edges_internal_nodes + el_faces_internal_nodes + iconn] = ip
                        iconn=iconn+1
                        iconn_level=iconn_level+1
                    end 
                    if (iconn_level > (ngl-2)^2) break end
                    m=ender
                    for l=ender-1:-1:starter
                        ip = conn_vol_poin[l,m,n,iel]
                        mesh.conn[iel, 8 + el_edges_internal_nodes + el_faces_internal_nodes + iconn] = ip
                        iconn=iconn+1
                        iconn_level=iconn_level+1
                    end 
                    if (iconn_level > (ngl-2)^2) break end
                    l=starter
                    for m=ender-1:starter+1
                        ip = conn_vol_poin[l,m,n,iel]
                        mesh.conn[iel, 8 + el_edges_internal_nodes + el_faces_internal_nodes + iconn] = ip
                        iconn=iconn+1
                        iconn_level=iconn_level+1
                    end
                    starter=starter+1
                    ender = ender-1
                end
            end
            #for i = 1:1:ngl
            #    for j = 1:1:ngl
            #        for k = 1:1:ngl
            #            ip = mesh.connijk[iel,i,j,k]
            #            @printf(f, " %d %d %d %d %d\n", iel, i, j, k, ip)
            #        end
            #    end
            #end
        end
    #end # do f


    #show(stdout, "text/plain", mesh.conn')
    #for iel = 1:mesh.nelem
    #    show(stdout, "text/plain", mesh.connijk[iel,:,:,:])
    #end

    println_rank(" # POPULATE GRID with SPECTRAL NODES ............................ VOLUMES DONE"; msg_rank = rank, suppress = mesh.msg_suppress)

end

function mod_mesh_build_mesh!(mesh::St_mesh, interpolation_nodes, backend)

    if (mesh.nsd > 1)
        @error(" USE GMSH to build a higher-dimensional grid!")
    end
    
    println(" # BUILD LINEAR CARTESIAN GRID ............................")
    
    mesh.npoin_linear = mesh.npx
    mesh.npoin        = mesh.npoin_linear #This will be updated for high order grids
    mesh.nelem        = mesh.npx - 1
    
    Δx::TFloat=0.0
    resize!(mesh.Δx, mesh.nelem)
    
    Δx = abs(mesh.xmax - mesh.xmin)/(mesh.nelem)
    mesh.npoin = mesh.npx

    mesh.x[1] = mesh.xmin
    for i = 2:mesh.npx
        mesh.x[i] = mesh.x[i-1] + Δx
        mesh.Δx[i-1] = Δx #Constant for the sake of simplicity in 1D problems. This may change later
    end
    mesh.NNODES_EL  = 2
    
    println(" # 1D NATIVE LINEAR GRID PROPERTIES")
    println(" # N. elements       : ", mesh.nelem)
    println(" # N. points         : ", mesh.npoin_linear)
    println(" # 1D NATIVE LINEAR GRID PROPERTIES ...................... END")
    
    ngl                     = mesh.nop + 1
    tot_linear_poin         = mesh.npoin_linear    
    tot_vol_internal_nodes  = mesh.nelem*(ngl-2)  
    el_vol_internal_nodes   = (ngl-2)
    
    #Update number of grid points from linear count to total high-order points
    mesh.npoin = tot_linear_poin + tot_vol_internal_nodes
    
    if (mesh.nop > 1)
        println(" # 1D NATIVE HIGH-ORDER GRID PROPERTIES")
        println(" # N. volumes internal points : ", tot_vol_internal_nodes)
        println(" # N. total high order points : ", mesh.npoin)
        println(" # 1D NATIVE HIGH-ORDER GRID PROPERTIES ...................... END")
    end
    
    
    # Resize (using resize! from ElasticArrays) as needed
    resize!(mesh.x, (mesh.npoin))
    mesh.npoin_el = ngl
    #allocate mesh.conn and reshape it
    mesh.conn = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.npoin_el))
    mesh.connijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.ngl), 1, 1)
    for iel = 1:mesh.nelem
        mesh.conn[iel, 1] = iel
        mesh.conn[iel, 2] = iel + 1

        mesh.connijk[iel, 1] = iel
        mesh.connijk[iel, mesh.ngl] = iel + 1
    end
    
    #Add high-order nodes
    add_high_order_nodes_1D_native_mesh!(mesh, interpolation_nodes, backend)

    mesh.nelem_semi_inf = 0
    if (inputs[:llaguerre_1d_right]) mesh.nelem_semi_inf +=1 end 
    if (inputs[:llaguerre_1d_left]) mesh.nelem_semi_inf +=1 end
    if (mesh.nelem_semi_inf == 0) mesh.nelem_semi_inf = 1 end  
    mesh.connijk_lag = KernelAbstractions.zeros(backend,TInt, Int64(mesh.nelem_semi_inf), Int64(mesh.ngr), 1, 1)
    mesh.npoin_original = mesh.npoin
    if (inputs[:llaguerre_1d_right])
        x = KernelAbstractions.zeros(backend, TFloat, mesh.npoin+mesh.ngr-1)      
        x[1:mesh.npoin] .= mesh.x[1:mesh.npoin] 
        gr = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta],backend)
        mesh.connijk_lag[1,1,1] = mesh.npoin_linear 
        for i=2:mesh.ngr
            ip = mesh.npoin+i-1
            mesh.connijk_lag[1,i,1] = ip
            x[ip] = mesh.xmax + inputs[:yfac_laguerre]*gr.ξ[i]
        end    
        mesh.npoin = mesh.npoin + mesh.ngr-1
        mesh.x = x
    end
    if (inputs[:llaguerre_1d_left])
        e = min(2,mesh.nelem_semi_inf)
        x = zeros(Float64,mesh.npoin+mesh.ngr-1)
        x[1:mesh.npoin] .= mesh.x[1:mesh.npoin]
        gr = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta],backend)
        mesh.connijk_lag[e,1,1] = 1
        for i=2:mesh.ngr
            ip = mesh.npoin+i-1
            mesh.connijk_lag[e,i,1] = ip
            x[ip] = mesh.xmin - inputs[:yfac_laguerre]*gr.ξ[i]
        end
        mesh.npoin = mesh.npoin + mesh.ngr-1
        mesh.x = x
    end 
    #plot_1d_grid(mesh)
    resize!(mesh.y, (mesh.npoin))
    println(" # BUILD LINEAR CARTESIAN GRID ............................ DONE")
    
end


function mod_mesh_mesh_driver(inputs::Dict, nparts, distribute, adapt_flags = nothing, partitioned_model_coarse = nothing, omesh = nothing)
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if (haskey(inputs, :lread_gmsh) && inputs[:lread_gmsh]==true)
        
        println_rank(" # Read gmsh grid and populate with high-order points "; msg_rank = rank, suppress = omesh == !isnothing)
        
        # Initialize mesh struct: the arrays length will be increased in mod_mesh_read_gmsh
        mesh = St_mesh{TInt,TFloat, CPU()}(nsd=TInt(inputs[:nsd]),
                                    nop=TInt(inputs[:nop]),
                                    ngr=TInt(inputs[:nop_laguerre]+1),
                                    SD=NSD_1D())
        
        # Read gmsh grid using the GridapGmsh reader
        n2o_ele_map = nothing
        if isnothing(adapt_flags)
            partitioned_model = mod_mesh_read_gmsh!(mesh, inputs, nparts, distribute)
        else
            partitioned_model, n2o_ele_map = mod_mesh_read_gmsh!(mesh, inputs, nparts, distribute, adapt_flags, partitioned_model_coarse, omesh)
        end

        
        println_rank(" # Read gmsh grid and populate with high-order points ........................ DONE"; msg_rank = rank, suppress = mesh.msg_suppress)
        
    else
        
        println(" # Build native grid")
        # Initialize mesh struct for native structured grid:
        if (haskey(inputs, :nsd))
            
            if (inputs[:nsd]==1)
                println(" # ... build 1D grid ")
                mesh = St_mesh{TInt,TFloat, CPU()}(x = KernelAbstractions.zeros(CPU(),TFloat,Int64(inputs[:npx])),
                                            npx  = TInt(inputs[:npx]),
                                            xmin = TFloat(inputs[:xmin]), xmax = TFloat(inputs[:xmax]),
                                            nop=TInt(inputs[:nop]),
                                            connijk = KernelAbstractions.zeros(CPU(), TInt,  Int64(inputs[:nelx]), Int64(inputs[:nop]+1), 1, 1),
                                            ngr=TInt(inputs[:nop_laguerre]+1),
                                            SD=NSD_1D())
                
            elseif (inputs[:nsd]==2)
                println(" # ... build 2D grid ")
                mesh = St_mesh{TInt,TFloat, CPU()}(x =  KernelAbstractions.zeros(CPU(),TFloat,Int64(inputs[:npx])),
                                            z = zeros(TInt(inputs[:npz])),
                                            npx  = TInt(inputs[:npx]),
                                            npz  = TInt(inputs[:npz]), 
                                            xmin = TFloat(inputs[:xmin]), xmax = TFloat(inputs[:xmax]),
                                            zmin = TFloat(inputs[:zmin]), zmax = TFloat(inputs[:zmax]),
                                            nop=TInt(inputs[:nop]),
                                            SD=NSD_2D())
                
            elseif (inputs[:nsd]==3)
                println(" # ... build 3D grid ")
                mesh = St_mesh{TInt,TFloat, CPU()}(x = KernelAbstractions.zeros(CPU(),TFloat, Int64(inputs[:npx])),
                                            y = zeros(TInt(inputs[:npy])),
                                            z = zeros(TInt(inputs[:npz])),
                                            npx  = TInt(inputs[:npx]),
                                            npy  = TInt(inputs[:npy]),
                                            npz  = TInt(inputs[:npz]), 
                                            xmin = TFloat(inputs[:xmin]), xmax = TFloat(inputs[:xmax]),
                                            ymin = TFloat(inputs[:ymin]), ymax = TFloat(inputs[:ymax]),
                                            zmin = TFloat(inputs[:zmin]), zmax = TFloat(inputs[:zmax]),
                                            nop=TInt(inputs[:nop]),
                                            SD=NSD_3D())
            else
                @error( " INPUT ERROR: nsd must be an integer in [1, 2, 3] ")
            end
            
        else
            
            #
            # Default grid is 1D if `nsd` is not defined in user_input.jl
            #
            println(" # ... build DEFAULT 1D grid")
            println(" # ...... DEFINE NSD in your input dictionary if you want a different grid!")
            mesh = St_mesh{TInt,TFloat, CPU()}(x = KernelAbstractions.zeros(CPU(),TFloat,Int64(inputs[:npx])),
                                        npx  = Int64(inputs[:npx]),
                                        xmin = TFloat(inputs[:xmin]), xmax = TFloat(inputs[:xmax]),
                                        nop=Int64(inputs[:nop]),
                                        ngr=Int64(inputs[:nop_laguerre]+1),
                                        SD=NSD_1D())
        end
        mod_mesh_build_mesh!(mesh,  inputs[:interpolation_nodes], CPU())
        
        #Write structured grid to VTK
        #vtkfile = vtk_grid("mySTRUCTURED_GRID", mesh.x, mesh.y, mesh.z) # 3-D
        #outfiles = vtk_save(vtkfile)
        
        println(" # Build native grid ........................ DONE")
    end
    
    if (mesh.nsd == 1)
        mesh.SD = NSD_1D()
    elseif (mesh.nsd == 2)
        mesh.SD = NSD_2D()
    elseif (mesh.nsd == 3)
        mesh.SD = NSD_3D()
    else
        error(" Drivers.jl: Number of space dimnnsions unknow! CHECK Your grid!")
    end

    #
    # Element size:
    # 
    #   mesh.Δelem[1:nelem] --> smallest distance between two corner points inside each element
    #   mesh.Δelem_smallest --> size of the smallest element inside a grid (==minimum(mesh.Δelem))
    #   mesh.Δeffective     --> mesh.Δelem_smallest/mesh.nop
    #
    compute_element_size_driver(mesh, mesh.SD, Float64, CPU())
    
    if isnothing(adapt_flags)
        return mesh, partitioned_model
    else
        return mesh, partitioned_model, n2o_ele_map
    end
        
end

#-----------------------------------------------------------------------
# This subroutine computes size of the element hexa 
# and the mean spacing !between the lgl points inside the element.
#
# The following local node numbering is assumed
#
#  3----------4
#  |          |
#  |          |
#  |          |
#  |          |
#  1----------2
#
# The following local node numbering is assumed
#
#      7---------8
#     /| Top    /|
#    / |       / |
#   5---------6  |
#   |  |      |  |
#   |  3------|--4
#   | /  Bott | /
#   |/        |/
#   1---------2
#
# WARNING: this only gives an estimate if the grid is not cartesian
#
#----------------------------------------------------------------------
function compute_element_size_driver(mesh::St_mesh, SD, T, backend)
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    mesh.Δelem = KernelAbstractions.zeros(backend, T, mesh.nelem)
    for ie = 1:mesh.nelem
         compute_element_size!(SD, ie, mesh::St_mesh, T)
    end
    mesh.Δelem_s      = minimum(mesh.Δelem)    
    mesh.Δelem_l      = maximum(mesh.Δelem)
    mesh.Δeffective_s = TFloat(mesh.Δelem_s/mesh.nop)
    mesh.Δeffective_l = TFloat(mesh.Δelem_l/mesh.nop)

    println_rank(" # "; msg_rank = rank, suppress = mesh.msg_suppress)
    println_rank(" # ELEMENT SIZES:"; msg_rank = rank, suppress = mesh.msg_suppress)
    println_rank(" #   The smallest element has size: ", mesh.Δelem_s, " and effective resolution ", mesh.Δeffective_s; msg_rank = rank, suppress = mesh.msg_suppress)
    println_rank(" #   The biggest  element has size: ", mesh.Δelem_l, " and effective resolution ", mesh.Δeffective_l; msg_rank = rank, suppress = mesh.msg_suppress)
    println_rank(" # "; msg_rank = rank, suppress = mesh.msg_suppress)
end

#------------------------------------------------------------------------------------
#Computes element size assuming flow in a straght-sided cube
#------------------------------------------------------------------------------------
function compute_element_size!(SD::NSD_1D, ie, mesh::St_mesh, T) nothing end
    
function compute_element_size!(SD::NSD_2D, ie, mesh::St_mesh, T)
    
    #local arrays
    ngl   = mesh.ngl
    
    x     = zeros(T, 4)
    y     = zeros(T, 4)
    inode = zeros(TInt, 4)
    
    inode[1] = mesh.connijk[ie, 1,   ngl]
    inode[2] = mesh.connijk[ie, 1,     1]
    inode[3] = mesh.connijk[ie, ngl, ngl]
    inode[4] = mesh.connijk[ie, ngl,   1]
    
    #Store Coordinates
    for m = 1:4
        x[m] = mesh.x[inode[m]]
        y[m] = mesh.y[inode[m]]
        #@info m, x[m], y[m]
    end
    
    # Nodes distances as if it were linear
    dx = maximum(x) - minimum(x)
    dy = maximum(y) - minimum(y)
    
    mesh.Δelem[ie] = min(dx, dy)           #shortest distance of two points corner within a given element
    #mesh.Δelem_largest[ie]  = max(dx, dy) #longest distance of two points corner within a given element
    
end

function compute_element_size!(SD::NSD_3D, ie, mesh::St_mesh, T)
    
    #local arrays
    ngl   = mesh.ngl
    
    x     = zeros(T, 8)
    y     = zeros(T, 8)
    z     = zeros(T, 8)
    inode = zeros(TInt, 8)
    
    inode[1] = mesh.connijk[ie, 1,     1, ngl]
    inode[2] = mesh.connijk[ie, 1,     1,   1]
    inode[3] = mesh.connijk[ie, ngl,   1, ngl]
    inode[4] = mesh.connijk[ie, ngl,   1,   1]
    inode[5] = mesh.connijk[ie, 1,   ngl, ngl]
    inode[6] = mesh.connijk[ie, 1,   ngl,   1]
    inode[7] = mesh.connijk[ie, ngl, ngl, ngl]
    inode[8] = mesh.connijk[ie, ngl, ngl,   1]
    
    #Store Coordinates
    for m = 1:8
        x[m] = mesh.x[inode[m]]
        y[m] = mesh.y[inode[m]]
        z[m] = mesh.z[inode[m]]
        # @info m, x[m], y[m], z[m]
    end

    #Element sizes (as if it were linear)
    dx = maximum(x) - minimum(x)
    dy = maximum(y) - minimum(y)
    dz = maximum(z) - minimum(z)
    
    mesh.Δelem[ie] = min(dx, dy, dz)            #shortest distance of two points corner within a given element
    #mesh.Δelem_largest[ie]  = max(dx, dy, dz)  #longest distance of two points corner within a given element
    
end

function get_bdy_poin_in_face_on_edges!(mesh::St_mesh, isboundary_face, SD::NSD_3D)
    
    for iface_g in  findall(x -> x == true, isboundary_face)
        
        #GGNS numbering
        ip1 = mesh.conn_unique_faces[iface_g][1]
        ip2 = mesh.conn_unique_faces[iface_g][2]
        ip3 = mesh.conn_unique_faces[iface_g][4]
        ip4 = mesh.conn_unique_faces[iface_g][3]
            
        ###find edges belonging to this face and populate edges on poin_in_face array
        for iedge in mesh.face_edge_ids[iface_g]
            ipe1 = mesh.conn_unique_edges[iedge][1]
            ipe2 = mesh.conn_unique_edges[iedge][2]
            if (ipe1 == ip1 && ipe2 == ip2) || (ipe1 == ip2 && ipe2 == ip1)
                for i=2:mesh.ngl-1
                    mesh.poin_in_face[iface_g,i,1] = mesh.poin_in_edge[iedge,i]
                end
            elseif (ipe1 == ip2 && ipe2 == ip3) || (ipe1 == ip3 && ipe2 == ip2)
                for i=2:mesh.ngl-1
                    mesh.poin_in_face[iface_g,mesh.ngl,i] = mesh.poin_in_edge[iedge,i]
                end
            elseif (ipe1 == ip3 && ipe2 == ip4) || (ipe1 == ip4 && ipe2 == ip3)
                for i=2:mesh.ngl-1
                    mesh.poin_in_face[iface_g,i,mesh.ngl] = mesh.poin_in_edge[iedge,i]
                end
            elseif (ipe1 == ip4 && ipe2 == ip1) || (ipe1 == ip1 && ipe2 == ip4)
                for i=2:mesh.ngl-1
                    mesh.poin_in_face[iface_g,1,i] = mesh.poin_in_edge[iedge,i]
                end
            end
        end
    end
end

function send_and_receive(data2send, send_targets, comm)
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    T = eltype(data2send)

    # Validate inputs
    if length(data2send) != length(send_targets)
        @info rank, data2send, send_targets
        error("data2send and send_targets must have the same size")
    end

    # Organize data to send
    send_data = Dict(i => T[] for i in 0:size-1)
    send_data_order = Dict(i => T[] for i in 0:size-1)
    for (data, owner) in zip(data2send, send_targets)
        if owner < 0 || owner >= size
            error("send_targets contains invalid rank")
        end
        push!(send_data[owner], data)
    end

    # Prepare send sizes and receive sizes
    send_sizes = [length(send_data[i]) for i in 0:size-1]
    recv_sizes = MPI.Alltoall(MPI.UBuffer(send_sizes, 1), comm)

    # Prepare buffers for sending and receiving data
    send_buffers = [send_data[i] for i in 0:size-1]
    recv_buffers = [Vector{T}(undef, recv_sizes[i+1]) for i in 0:size-1]

    # Communicate data
    requests = MPI.Request[]
    for i in 0:size-1
        if send_sizes[i+1] > 0
            push!(requests, MPI.Isend(send_buffers[i+1], i, 0, comm))
        end
        if recv_sizes[i+1] > 0
            push!(requests, MPI.Irecv!(recv_buffers[i+1], i, 0, comm))
        end
    end

    # Wait for all communication to complete
    MPI.Waitall!(requests)

    # Combine received data into a single vector
    combined_recv_data = T[]
    original_senders   = Int[]
    for i in 0:size-1
        if recv_sizes[i+1] > 0
            append!(combined_recv_data, recv_buffers[i+1])
            append!(original_senders, fill(i, recv_sizes[i+1]))
        end
    end

    # # Output the combined data and source list
    # println("Processor $rank received combined data: ", combined_recv_data)
    # println("Processor $rank data sources: ", original_senders)

    return combined_recv_data, original_senders
end

function get_ghost_ips(gelm_ghost, gfacets_ghost, gfacets_owner, conn, pelm2elm, ip2gip, ngl, ghost_p_or_c, comm)
    rank = MPI.Comm_rank(comm)
    # @info "gfacets_ghost, gfacets_owner", rank, gfacets_ghost, gfacets_owner
    gelm_recv, original_senders = send_and_receive(gelm_ghost, gfacets_owner, comm)
    gfacets_recv = send_and_receive(gfacets_ghost, gfacets_owner, comm)[1]
    # @info "gfacets_recv, original_senders", rank,  gfacets_recv, original_senders
    # println("Rank $rank gfacets_recv: $gfacets_recv")
    lcells = [pelm2elm[x] for x in gelm_recv]
    # @info "lcells", rank, lcells 
    lfacets = gfacets_recv
    # @info "lfacets", rank, lfacets 
    # lfacets = global_to_local(facetsids)[vcat(gfacets_recv...)]
    nlfacets    = size(lfacets,1)
    ips_send    = KernelAbstractions.zeros(CPU(), TInt, nlfacets * ngl)
    ips_targets = KernelAbstractions.zeros(CPU(), TInt, nlfacets * ngl)
    IP          = KernelAbstractions.zeros(CPU(), TInt, ngl)
    cnt = 1
    for (i, (lfacet, lcell)) in enumerate(zip(lfacets, lcells))

        if (lfacet == 1)
            m = ngl
            n = 1:ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, 1, n]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, ngl, n]
            end
        elseif (lfacet == 2)
            m = 1
            n = 1:ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, ngl, n]
            elseif ghost_p_or_c == 2
                # @info rank, IP, conn, lcell, m 
                IP .= conn[lcell, 1, n]
            end
        elseif (lfacet == 3)
            m = 1:ngl
            n = 1
            if ghost_p_or_c == 1
                IP .= conn[lcell, m, ngl]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, m, 1]
            end
        elseif (lfacet == 4)
            m = 1:ngl
            n = ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, m, 1]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, m, ngl]
            end
        end
        for ip in IP
            # @info rank, ip, cnt,i, lfacet, lcell
            ips_send[cnt] = ip2gip[ip]
            ips_targets[cnt] = original_senders[i]
            cnt += 1
        end
    end
    # @info "ips_send, ips_targets", rank, ips_send, ips_targets
    ips_recv, ips_owner  = send_and_receive(ips_send, ips_targets, comm)
    # @info "ips_recv, ips_owner", rank, ips_recv, ips_owner

    return ips_recv, ips_owner
end
