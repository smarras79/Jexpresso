using Gridap
using GridapGmsh
using Gridap.Arrays
using Gridap.Arrays: Table
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Geometry: GridMock



export mod_mesh_mesh_driver_extra
export mod_mesh_read_gmsh_extra!

const VERTEX_NODES = UInt64(1)
const EDGE_NODES   = UInt64(2)
const FACE_NODES   = UInt64(4)


function mod_mesh_read_gmsh_extra!(mesh::St_mesh, inputs::Dict)

    # determine backend
    backend = CPU()
    
    #
    # Read GMSH grid from file
    #
    model         = GmshDiscreteModel(inputs[:gmsh_filename_extra], renumber=true)
    topology      = get_grid_topology(model)
    mesh.nsd      = num_cell_dims(model)
    d_to_num_dfaces = [num_vertices(model), num_edges(model), num_cells(model)]
    
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
    
    #dump(topology)
    #
    # Mesh elements, nodes, faces, edges
    #
    mesh.npoin_linear = num_faces(model,POIN_flg)    
    mesh.npoin        = mesh.npoin_linear         #This will be updated for the high order grid
    mesh.nedges       = num_faces(model,EDGE_flg)
    mesh.nfaces       = num_faces(model,FACE_flg)   
    mesh.nelem        = num_faces(model,ELEM_flg)
    
    mesh.nelem_bdy    = count(get_isboundary_face(topology,mesh.nsd))
    mesh.nelem_int    = mesh.nelem - mesh.nelem_bdy
    mesh.nfaces_bdy   = count(get_isboundary_face(topology,mesh.nsd-1))
    mesh.nfaces_int   = mesh.nfaces - mesh.nfaces_bdy
    mesh.nedges_bdy   = count(get_isboundary_face(topology,mesh.nsd-2))
    mesh.nedges_int   = mesh.nedges - mesh.nedges_bdy
    
    #get_isboundary_face(topology,mesh.nsd-1)
    
    println(" # GMSH LINEAR GRID PROPERTIES")
    println(" # N. points         : ", mesh.npoin_linear)
    println(" # N. elements       : ", mesh.nelem)
    println(" # N. edges          : ", mesh.nedges)
    println(" # N. faces          : ", mesh.nfaces)    
    println(" # N. internal elem  : ", mesh.nelem_int)
    println(" # N. internal edges : ", mesh.nedges_int) 
    println(" # N. internal faces : ", mesh.nfaces_int)    
    println(" # N. boundary elem  : ", mesh.nelem_bdy)
    println(" # N. boundary edges : ", mesh.nedges_bdy)
    println(" # N. boundary faces : ", mesh.nfaces_bdy)
    println(" # GMSH LINEAR GRID PROPERTIES ...................... END")
    
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
    
    if (mesh.nop > 1)
        println(" # GMSH HIGH-ORDER GRID PROPERTIES")
        println(" # N. edges internal points   : ", tot_edges_internal_nodes)
        println(" # N. faces internal points   : ", tot_faces_internal_nodes)
        println(" # N. volumes internal points : ", tot_vol_internal_nodes)
        println(" # N. total high order points : ", mesh.npoin)
        println(" # GMSH HIGH-ORDER GRID PROPERTIES ...................... END")
    end
    
    #
    # Resize as needed
    #
    mesh.x = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
    mesh.y = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
    mesh.z = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
    
    mesh.conn_edge_el = KernelAbstractions.zeros(backend, TInt, 2, Int64(mesh.NEDGES_EL), Int64(mesh.nelem))    
    mesh.conn_face_el = KernelAbstractions.zeros(backend, TInt,  4, Int64(mesh.NFACES_EL), Int64(mesh.nelem))  
    mesh.bdy_edge_in_elem = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy))  
    mesh.poin_in_edge = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges), Int64(mesh.ngl))
    mesh.poin_in_bdy_edge = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy), Int64(mesh.ngl))
    
    mesh.poin_in_face = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nfaces), Int64(mesh.ngl), Int64(mesh.ngl))
    mesh.edge_type     = Array{Union{Nothing, String}}(nothing, Int64(mesh.nedges))
    mesh.bdy_edge_type                    = Array{Union{Nothing, String}}(nothing, Int64(mesh.nedges_bdy))
    mesh.bdy_edge_type_id = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy))  
    
    if mesh.nsd > 2
        mesh.poin_in_bdy_face = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nfaces_bdy), Int64(mesh.ngl), Int64(mesh.ngl))
        mesh.face_type = Array{Union{Nothing, String}}(nothing, Int64(mesh.nfaces))
        mesh.bdy_face_type = Array{Union{Nothing, String}}(nothing, Int64(mesh.nfaces_bdy))
        mesh.bdy_face_in_elem = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nfaces_bdy))
    end
    mesh.npoin_el         = mesh.NNODES_EL + el_edges_internal_nodes + el_faces_internal_nodes + (mesh.nsd - 2)*el_vol_internal_nodes
    mesh.conn = KernelAbstractions.zeros(backend,TInt, Int64(mesh.nelem), Int64(mesh.npoin_el))
    
    #
    # Connectivity matrices
    #
    mesh.cell_node_ids     = model.grid.cell_node_ids
    mesh.conn_unique_faces = get_face_nodes(model, FACE_flg) #faces --> 4 nodes
    mesh.conn_unique_edges = get_face_nodes(model, EDGE_flg) #edges --> 2 nodes

    mesh.cell_edge_ids     = get_faces(topology, mesh.nsd, 1) #edge map from local to global numbering i.e. iedge_g = cell_edge_ids[1:NELEM][1:NEDGES_EL]
    mesh.cell_face_ids     = get_faces(topology, mesh.nsd, mesh.nsd-1) #face map from local to global numbering i.e. iface_g = cell_face_ids[1:NELEM][1:NFACE_EL]
    mesh.face_edge_ids     = get_faces(topology,mesh.nsd-1, 1)
    mesh.facet_cell_ids     = get_faces(topology,mesh.nsd-1, mesh.nsd)
    # @info mesh.facet_cell_ids
    mesh.edge_g_color::Array{Int64, 1} = zeros(Int64, mesh.nedges)


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
            # 3-----4
            # |     |
            # |     |
            # 1-----2
            #
            mesh.connijk[iel, 1,      1] = mesh.cell_node_ids[iel][2]
            mesh.connijk[iel, 1,    ngl] = mesh.cell_node_ids[iel][1]
            mesh.connijk[iel, ngl,  ngl] = mesh.cell_node_ids[iel][3]
            mesh.connijk[iel, ngl,    1] = mesh.cell_node_ids[iel][4]
            
            #@printf(" [1,1] [ngl, 1] [1, ngl] [ngl, ngl] %d %d %d %d\n", mesh.connijk[iel, 1, 1], mesh.connijk[iel, ngl, 1] , mesh.connijk[iel, 1,ngl], mesh.connijk[iel, ngl, ngl] )
            
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
        #open("./COORDS_LO.dat", "w") do f
            for ip = 1:mesh.npoin_linear
                
                mesh.x[ip] = model.grid.node_coordinates[ip][1]
                mesh.y[ip] = model.grid.node_coordinates[ip][2]
                
        #        @printf(f, " %.6f %.6f 0.000000 %d\n", mesh.x[ip],  mesh.y[ip], ip)
            end
        #end #f

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
        #open("./COORDS_LO.dat", "w") do f
            for ip = 1:mesh.npoin_linear
                mesh.x[ip] = model.grid.node_coordinates[ip][1]
                mesh.y[ip] = model.grid.node_coordinates[ip][2]
                mesh.z[ip] = model.grid.node_coordinates[ip][3]
        #        @printf(f, " %.6f %.6f %.6f %d\n", mesh.x[ip],  mesh.y[ip], mesh.z[ip], ip)
            end
        #end #f
    end


    #
    # Add high-order points to edges, faces, and elements (volumes)
    #
    # initialize LGL struct and buyild Gauss-Lobatto-xxx points
    lgl = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.nop, backend)

    println(" # POPULATE GRID with SPECTRAL NODES ............................ ")
    #
    # Edges
    #
    populate_conn_edge_el!(mesh, mesh.SD)
    add_high_order_nodes_edges!(mesh, lgl, mesh.SD, backend)

    #
    # Faces
    #
    populate_conn_face_el!(mesh, mesh.SD)
    add_high_order_nodes_faces!(mesh, lgl, mesh.SD)

    #
    # Volume
    #
    # NOTICE: in 2D we consider only edges. faces are the elements.
    #         
    add_high_order_nodes_volumes!(mesh, lgl, mesh.SD)
    
    for ip = mesh.npoin_linear+1:mesh.npoin
        mesh.x[ip] = mesh.x_ho[ip]
        mesh.y[ip] = mesh.y_ho[ip]
        mesh.z[ip] = 0.0
        if (mesh.nsd > 2)
            mesh.z[ip] = mesh.z_ho[ip]
        end
    end
    
    mesh.xmax = maximum(mesh.x)
    mesh.xmin = minimum(mesh.x)
    mesh.ymax = maximum(mesh.y)
    mesh.ymin = minimum(mesh.y)
    if (mesh.nsd > 2)
        mesh.zmax = maximum(mesh.z)
        mesh.zmin = minimum(mesh.z)
    end
    
    #----------------------------------------------------------------------
    # Extract boundary edges and faces nodes:
    #----------------------------------------------------------------------
    #
    # Bdy edges
    #
    if mesh.nsd == 2
        isboundary_edge = compute_isboundary_face(topology, EDGE_flg)
        #
        # Get labels contained in the current GMSH grid:
        #
        n_semi_inf = 0
        labels = get_face_labeling(model)
        for ilabel in labels.tag_to_name
            edges_to_tag  = get_face_tag_index(labels,ilabel,EDGE_flg)
            idx_edges_inflow = findall( x -> x == 1, edges_to_tag)
            #    
            # Tag the boundary edge with its type as defined in the user-provided GMSH file:
            #
            for idx in idx_edges_inflow
                mesh.edge_type[idx] = ilabel
            end
        end
        iedge_bdy = 1
        for iedge = 1:mesh.nedges #total nedges
            if isboundary_edge[iedge] == true
                for igl = 1:mesh.ngl
                    mesh.poin_in_bdy_edge[iedge_bdy, igl] = mesh.poin_in_edge[iedge, igl]
                    mesh.bdy_edge_type[iedge_bdy] = mesh.edge_type[iedge]
                    mesh.bdy_edge_in_elem[iedge_bdy] = mesh.facet_cell_ids[iedge][1]
                    if (size(mesh.facet_cell_ids[iedge],1) ≠ 1)
                        s = """
                        Check boundary elements! size(mesh.facet_cell_ids[iedge],1) ≠ 1 
                            """
                
                        @error s
                    end
                end
                if (mesh.bdy_edge_type[iedge_bdy] == "Laguerre")
                    n_semi_inf += 1
                end
                iedge_bdy += 1
            end
        end
        # build mesh data structs for Laguerre semi-infinite elements
        if ("Laguerre" in mesh.bdy_edge_type)
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
                    #tangent vector 
                    x = mesh.x[ip]
                    x1 = mesh.x[ip1]
                    y = mesh.y[ip]
                    y1 = mesh.y[ip1]
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
                    for i=1:mesh.ngl
                        ip = mesh.poin_in_bdy_edge[iedge,i]
                        mesh.connijk_lag[e_iter,i,1] = ip
                        for j=2:mesh.ngr
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
                                mesh.connijk_lag[e_iter,i,j] = iter
                                iter += 1
                                matched = 1
                            end
                            if (matched == 0)
                                x_new[iter] = x_temp#mesh.x[ip] + nor[1]*gr.ξ[j]*factorx
                                y_new[iter] = y_temp#mesh.y[ip] + nor[2]*gr.ξ[j]*factory
                                mesh.connijk_lag[e_iter,i,j] = iter
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
            mesh.x = x_new
            mesh.y = y_new
            mesh.z = KernelAbstractions.zeros(backend, TFloat, mesh.npoin)
            mesh.nelem_semi_inf = n_semi_inf
        end

    elseif mesh.nsd > 2
        isboundary_face = compute_isboundary_face(topology, FACE_flg)
        #
        # Get labels contained in the current GMSH grid:
        #
        labels = get_face_labeling(model)
        for ilabel in labels.tag_to_name
            faces_to_tag  = get_face_tag_index(labels,ilabel,FACE_flg)
            idx_faces_inflow = findall( x -> x == 1, faces_to_tag)
            #    
            # Tag the boundary edge with its type as defined in the user-provided GMSH file:
            #
            for idx in idx_faces_inflow
                mesh.face_type[idx] = ilabel
            end
        end
        get_bdy_poin_in_face_on_edges!(mesh, @view(isboundary_face[:]), mesh.SD)
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
                    
                            @error s
                        end
                        #@info "face point number", mesh.poin_in_face[iface,igl,jgl],iface,igl,jgl
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
#
#
# Free memory of obsolete arrays
#
mesh.x_ho = zeros(1)
mesh.y_ho = zeros(1)
mesh.z_ho = zeros(1)
#resize!(mesh.x_ho, 1)
#resize!(mesh.y_ho, 1)
#resize!(mesh.z_ho, 1)
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
=#

#show(stdout, "text/plain", mesh.conn')
println(" # POPULATE GRID with SPECTRAL NODES ............................ DONE")

#writevtk(model,"gmsh_grid")
end

function mod_mesh_mesh_driver_extra(inputs::Dict)
    if (haskey(inputs, :lread_gmsh_extra) && inputs[:lread_gmsh_extra]==true)
        
        println(" # Read gmsh grid and populate with high-order points ")
        
        # Initialize mesh struct: the arrays length will be increased in mod_mesh_read_gmsh
        mesh = St_mesh{TInt,TFloat, CPU()}(nsd=TInt(inputs[:nsd]),
                                    nop=TInt(inputs[:nop]),
                                    ngr=TInt(inputs[:nop_laguerre]+1),
                                    SD=NSD_1D())
        
        # Read gmsh grid using the GridapGmsh reader
        mod_mesh_read_gmsh_extra!(mesh, inputs)
        
        println(" # Read gmsh grid and populate with high-order points ........................ DONE")
        
    else
        
        println(" # Build native grid")
        # Initialize mesh struct for native structured grid:
        if (haskey(inputs, :nsd_extra))
            
            if (inputs[:nsd_extra]==1)
                println(" # ... build 1D grid ")
                mesh = St_mesh{TInt,TFloat, CPU()}(x = KernelAbstractions.zeros(CPU(),TFloat,Int64(inputs[:npx_extra])),
                                            npx  = TInt(inputs[:npx_extra]),
                                            xmin = TFloat(inputs[:xmin_extra]), xmax = TFloat(inputs[:xmax_extra]),
                                            nop=TInt(inputs[:nop_extra]),
                                            connijk = KernelAbstractions.zeros(CPU(), TInt,  Int64(inputs[:nelx_extra]), Int64(inputs[:nop_extra]+1), 1, 1),
                                            ngr=TInt(inputs[:nop_laguerre_extra]+1),
                                            SD=NSD_1D())
                
            elseif (inputs[:nsd_extra]==2)
                println(" # ... build 2D grid ")
                mesh = St_mesh{TInt,TFloat, CPU()}(x =  KernelAbstractions.zeros(CPU(),TFloat,Int64(inputs[:npx_extra])),
                                            z = zeros(TInt(inputs[:npz_extra])),
                                            npx  = TInt(inputs[:npx_extra]),
                                            npz  = TInt(inputs[:npz_extra]), 
                                            xmin = TFloat(inputs[:xmin_extra]), xmax = TFloat(inputs[:xmax_extra]),
                                            zmin = TFloat(inputs[:zmin_extra]), zmax = TFloat(inputs[:zmax_extra]),
                                            nop=TInt(inputs[:nop_extra]),
                                            SD=NSD_2D())
                
            elseif (inputs[:nsd_extra]==3)
                println(" # ... build 3D grid ")
                mesh = St_mesh{TInt,TFloat, CPU()}(x = KernelAbstractions.zeros(CPU(),TFloat, Int64(inputs[:npx_extra])),
                                            y = zeros(TInt(inputs[:npy_extra])),
                                            z = zeros(TInt(inputs[:npz_extra])),
                                            npx  = TInt(inputs[:npx_extra]),
                                            npy  = TInt(inputs[:npy_extra]),
                                            npz  = TInt(inputs[:npz_extra]), 
                                            xmin = TFloat(inputs[:xmin_extra]), xmax = TFloat(inputs[:xmax_extra]),
                                            ymin = TFloat(inputs[:ymin_extra]), ymax = TFloat(inputs[:ymax_extra]),
                                            zmin = TFloat(inputs[:zmin_extra]), zmax = TFloat(inputs[:zmax_extra]),
                                            nop=TInt(inputs[:nop_extra]),
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
            mesh = St_mesh{TInt,TFloat, CPU()}(x = KernelAbstractions.zeros(CPU(),TFloat,Int64(inputs[:npx_extra])),
                                        npx  = Int64(inputs[:npx_extra]),
                                        xmin = TFloat(inputs[:xmin_extra]), xmax = TFloat(inputs[:xmax_extra]),
                                        nop=Int64(inputs[:nop_extra]),
                                        ngr=Int64(inputs[:nop_laguerre_extra]+1),
                                        SD=NSD_1D())
        end
        mod_mesh_build_mesh!(mesh,  inputs[:interpolation_nodes_extra], CPU())
        
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
    
    return mesh
    
end
