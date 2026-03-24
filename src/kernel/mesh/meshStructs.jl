export St_mesh

Base.@kwdef mutable struct St_extra_mesh{TInt, TFloat, NSD, dims1, dims2, dims3, dims4, dims5, nelem, npoin, backend}

    extra_coords  = KernelAbstractions.zeros(backend,TFloat, dims1)
    extra_coords_cart = KernelAbstractions.zeros(backend,TFloat, dims5)
    extra_connijk = KernelAbstractions.zeros(backend,TInt, dims2)
    extra_nelem::Union{TInt, Missing} = nelem
    extra_npoin::Union{TInt, Missing} = npoin
    extra_nop = KernelAbstractions.zeros(backend,TInt, dims3)
    extra_metrics = allocate_metrics(NSD, dims4[1], dims4[2], dims4[3], TFloat, backend)
    Minv = KernelAbstractions.zeros(backend,TFloat, npoin)
    ωθ = KernelAbstractions.zeros(backend,TInt, 5)
    ωϕ = KernelAbstractions.zeros(backend,TInt, 5)
    ψ = KernelAbstractions.zeros(backend,TInt, 5,5)
    dψ = KernelAbstractions.zeros(backend,TInt, 5,5)
    ref_level = KernelAbstractions.zeros(backend,TInt,nelem)
end

Base.@kwdef mutable struct St_mesh{TInt, TFloat, backend}

    x      = KernelAbstractions.zeros(backend, TFloat, 2)
    y      = KernelAbstractions.zeros(backend, TFloat, 2)
    z      = KernelAbstractions.zeros(backend, TFloat, 2)
    coords = KernelAbstractions.zeros(backend, TFloat, 2, 1)
    
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

    conn::Array{TInt,2}   = KernelAbstractions.zeros(backend, TInt, 0, 0)
    conn_unique_edges     = Array{TInt}(undef,  1, 2)
    conn_unique_edges1    = Array{Int64}(undef,  1, 2)
    conn_unique_faces     = Array{TInt}(undef,  1, 4)
    poin_in_edge          = Array{TInt}(undef, 0, 0)
    internal_poin_in_edge = Array{TInt}(undef, 0, 0)
    conn_edge_el          = Array{TInt}(undef, 0, 0, 0)
    poin_in_face          = Array{TInt}(undef, 0, 0, 0)
    conn_face_el          = Array{TInt}(undef, 0, 0, 0)
    face_in_elem          = Array{TInt}(undef, 0, 0, 0)

    # Skeleton arrays needed by "element learning"
    lengthΓ  = 0  #non-repeated bdy points from mesh.poin_in_bdy_edge    
    lengthO  = 0  #all internal, including edges, but without domain's bdy
    length∂τ = 0
    lengthIo = 0
    length∂O = 0
    Γ::Array{TInt, 1}  = KernelAbstractions.zeros(backend, TInt, lengthΓ)
    O::Array{TInt, 1}  = KernelAbstractions.zeros(backend, TInt, lengthO)
    ∂τ::Array{TInt, 1} = KernelAbstractions.zeros(backend, TInt, length∂τ)
    Io::Array{TInt, 1} = KernelAbstractions.zeros(backend, TInt, lengthIo)
    ∂O::Array{TInt, 1} = KernelAbstractions.zeros(backend, TInt, length∂O)
        
    edge_g_color::Array{Int64, 1} = zeros(Int64, 1)
    
    #MPI variables
    rank                = 0
    ip2gip              = KernelAbstractions.zeros(backend, TInt, 0)
    gip2owner           = KernelAbstractions.zeros(backend, TInt, 0)
    gip2ip              = KernelAbstractions.zeros(backend, TInt, 0)
    parts               = 1
    nparts              = 1
    
    #Auxiliary arrays for boundary conditions
    edge_in_elem              = KernelAbstractions.zeros(backend, TInt, 0)
    bdy_edge_in_elem          = KernelAbstractions.zeros(backend, TInt, 0)
    poin_in_bdy_edge          = KernelAbstractions.zeros(backend, TInt, 0, 0)
    internal_poin_in_bdy_edge = KernelAbstractions.zeros(backend, TInt, 0, 0)
    internal_poin_in_elem     = KernelAbstractions.zeros(backend, TInt, 0, 0)
    bdy_face_in_elem          = KernelAbstractions.zeros(backend, TInt, 0)
    poin_in_bdy_face          = KernelAbstractions.zeros(backend, TInt, 0, 0, 0)
    elem_to_face              = KernelAbstractions.zeros(backend, TInt, 0, 0, 0, 0, 0)
    elem_to_edge              = KernelAbstractions.zeros(backend, TInt, 0, 0, 0, 0)
    edge_type                 = Array{Union{Nothing, String}}(nothing, 1)
    face_type                 = Array{Union{Nothing, String}}(nothing, 1)
    bdy_edge_type             = Array{Union{Nothing, String}}(nothing, 1)
    bdy_face_type             = Array{Union{Nothing, String}}(nothing, 1)
    bdy_edge_type_id          = KernelAbstractions.zeros(backend, TInt, 0)

    Δelem        = KernelAbstractions.zeros(backend, TInt, 0)
    Δelem_s      = 0.0
    Δelem_l      = 0.0
    Δeffective_s = 0.0
    Δeffective_l = 0.0
    extra_mesh = Array{St_extra_mesh}(undef, 0, 0)
    
    SD::AbstractSpaceDimensions

    # for AMR
    ad_lvl = KernelAbstractions.zeros(backend, TInt, 0)
    num_hanging_facets::Union{TInt, Missing} = 0
    non_conforming_facets                    = Vector{Vector{TInt}}(undef,num_hanging_facets)
    non_conforming_facets_parents_ghost      = Vector{Vector{TInt}}(undef,num_hanging_facets)
    non_conforming_facets_children_ghost     = Vector{Vector{TInt}}(undef,num_hanging_facets)

    pgip_ghost = KernelAbstractions.zeros(backend, TInt, 0)
    pgip_owner = KernelAbstractions.zeros(backend, TInt, 0)
    pgip_local = KernelAbstractions.zeros(backend, TInt, 0)
    cgip_ghost = KernelAbstractions.zeros(backend, TInt, 0)
    cgip_owner = KernelAbstractions.zeros(backend, TInt, 0)
    cgip_local = KernelAbstractions.zeros(backend, TInt, 0)

    q_local_c  = KernelAbstractions.zeros(backend, TFloat, 0)


    # non_conforming_facets arrays
    cip         = KernelAbstractions.zeros(backend, TInt, 0)
    pip         = KernelAbstractions.zeros(backend, TInt, 0)
    lfid        = KernelAbstractions.zeros(backend, TInt, 0)
    half1       = KernelAbstractions.zeros(backend, TInt, 0)
    half2       = KernelAbstractions.zeros(backend, TInt, 0)
    IPc_list    = KernelAbstractions.zeros(backend, TInt, 0, 0)
    IPp_list    = KernelAbstractions.zeros(backend, TInt, 0, 0)
    IPc_list_pg = KernelAbstractions.zeros(backend, TInt, 0, 0)
    IPp_list_cg = KernelAbstractions.zeros(backend, TInt, 0, 0)
    # non_conforming_facets arrays for own child facet, ghost parent facet 
    cip_pg   = KernelAbstractions.zeros(backend, TInt, 0)
    lfid_pg  = KernelAbstractions.zeros(backend, TInt, 0)
    half1_pg = KernelAbstractions.zeros(backend, TInt, 0)
    half2_pg = KernelAbstractions.zeros(backend, TInt, 0)
    IPpg_list= KernelAbstractions.zeros(backend, TInt, 0, 0)
    # non_conforming_facets arrays for ghost child facet, own parent facet
    pip_cg   = KernelAbstractions.zeros(backend, TInt, 0)
    lfid_cg  = KernelAbstractions.zeros(backend, TInt, 0)
    half1_cg = KernelAbstractions.zeros(backend, TInt, 0)
    half2_cg = KernelAbstractions.zeros(backend, TInt, 0)
    IPcg_list= KernelAbstractions.zeros(backend, TInt, 0, 0)

    num_ncf::Union{TInt, Missing}    = 0
    num_ncf_pg::Union{TInt, Missing} = 0
    num_ncf_cg::Union{TInt, Missing} = 0

    lneed_redistribute::Bool = false

    msg_suppress::Bool = false

end
