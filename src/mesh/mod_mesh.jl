using Test
using DelimitedFiles
using Gridap
using Gridap.Arrays
using Gridap.Arrays: Table
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Geometry: GridMock
using GridapGmsh
using LinearAlgebra
using Printf
using Revise
using StaticArrays

export St_mesh

export mod_mesh_build_mesh!
export mod_mesh_read_gmsh!

const POIN         = UInt8(0)
const EDGE         = UInt8(1)
const FACE         = UInt8(2)
const ELEM         = UInt8(3)

const VERTEX_NODES = UInt8(1)
const EDGE_NODES   = UInt8(2)
const FACE_NODES   = UInt8(4)


abstract type At_geo_entity end


include("../basis/basis_structs.jl")

Base.@kwdef mutable struct St_mesh{TInt, TFloat}
    
    x::Union{Array{TFloat}, Missing} = zeros(2)
    y::Union{Array{TFloat}, Missing} = zeros(2)
    z::Union{Array{TFloat}, Missing} = zeros(2)
    
    x_ho::Union{Array{TFloat}, Missing} = zeros(2)
    y_ho::Union{Array{TFloat}, Missing} = zeros(2)
    z_ho::Union{Array{TFloat}, Missing} = zeros(2)
    
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
    npoin::Union{TInt, Missing} = 1        #This is updated after populating with high-order nodes
    npoin_linear::Union{TInt, Missing} = 1 #This is always the original number of the first-order grid
    
    nedges::Union{TInt, Missing} = 1     # total number of edges
    nedges_bdy::Union{TInt, Missing} = 1 # bdy edges
    nedges_int::Union{TInt, Missing} = 1 # internal edges

    nfaces::Union{TInt, Missing} = 1     # total number of faces
    nfaces_bdy::Union{TInt, Missing} = 1 # bdy faces
    nfaces_int::Union{TInt, Missing} = 1 # internal faces
    
    nsd::Union{TInt, Missing} = 1
    nop::Union{TInt, Missing} = 4
    ngl::Union{TInt, Missing} = nop + 1

    NEDGES_EL::Union{TInt, Missing}  = 12
    NFACES_EL::Union{TInt, Missing}  =  6
    EDGE_NODES::Union{TInt, Missing} =  2
    FACE_NODES::Union{TInt, Missing} =  4
    
    #low and high order connectivity tables
    cell_node_ids::Table{Int64,Vector{Int64},Vector{Int64}} = Gridap.Arrays.Table(zeros(nelem), zeros(8))
    cell_node_ids_ho::Table{Int64,Vector{Int64},Vector{Int64}} = Gridap.Arrays.Table(zeros(nelem), zeros(8))
    
    conn_ho           = Array{Int64, 1}(undef, nelem*ngl)
    conn_unique_edges = Array{Int64, 2}(undef,  1, 2)
    conn_unique_faces = Array{Int64, 2}(undef,  1, 4)
    conn_edge_el      = Array{Int64, 3}(undef,  nelem, 12, 2)
    conn_face_el      = Array{Int64, 3}(undef,  nelem,  6, 4)
    face_in_elem      = Array{Int64, 3}(undef,  nelem,  6, 2)    
end

function mod_mesh_build_mesh!(mesh::St_mesh)

    Δx::TFloat=0.0
    Δy::TFloat=0.0
    Δz::TFloat=0.0   
    
    if (mesh.nsd == 1)
        mesh.npy = mesh.npz = 1
        Δx = abs(mesh.xmax - mesh.xmin)/(mesh.npx - 1)
        Δy = 0.0
        Δz = 0.0
    elseif (mesh.nsd == 2)
        mesh.npz = 1
        Δx = abs(mesh.xmax - mesh.xmin)/(mesh.npx - 1)
        Δy = abs(mesh.ymax - mesh.ymin)/(mesh.npy - 1)
        Δz = 0.0
    else
        Δx = abs(mesh.xmax - mesh.xmin)/(mesh.npx - 1)
        Δy = abs(mesh.ymax - mesh.ymin)/(mesh.npy - 1)
        Δz = abs(mesh.zmax - mesh.zmin)/(mesh.npz - 1)
    end
    mesh.npoin = mesh.npx*mesh.npy*mesh.npz
    
    for i = 1:mesh.npx
        mesh.x[i] = (i - 1)*Δx
    end
    for j = 1:mesh.npy
        mesh.y[j] = (j - 1)*Δy
    end
    for k = 1:mesh.npz
        mesh.z[k] = (k - 1)*Δz
    end
    
end

function mod_mesh_read_gmsh!(mesh::St_mesh, gmsh_filename::String)

    #
    # Read GMSH grid from file
    #
    model    = GmshDiscreteModel(gmsh_filename, renumber=true)
    topology = get_grid_topology(model)
    mesh.nsd = num_cell_dims(model)
    
    if mesh.nsd == 3
        mesh.NEDGES_EL  = 12
        mesh.NFACES_EL  = 6
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 4
    elseif mesh.nsd == 2
        mesh.NEDGES_EL  = 4
        mesh.NFACES_EL  = 1
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 4
    elseif mesh.nsd == 1
        mesh.NEDGES_EL  = 1
        mesh.NFACES_EL  = 0
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 0
    else
        error( " WRONG NSD: This is not theoretical physics: we only handle 1, 2, or 3 dimensions!")
    end
    
    
    
    #@info topology.vertex_coordinates
    
    #dump(topology)
    #
    # Mesh elements, nodes, faces, edges
    #
    mesh.npoin_linear = num_faces(model,POIN)
    mesh.npoin        = mesh.npoin_linear     #This will be updated for the high order grid
    mesh.nedges       = num_faces(model,EDGE)
    mesh.nfaces       = num_faces(model,FACE)
    mesh.nelem        = num_faces(model,ELEM)
    mesh.nfaces_bdy   = count(get_isboundary_face(topology,mesh.nsd-1))
    mesh.nfaces_int   = mesh.nfaces - mesh.nfaces_bdy
    mesh.nedges_bdy   = count(get_isboundary_face(topology,mesh.nsd-2))
    mesh.nedges_int   = mesh.nedges - mesh.nedges_bdy
    
    println(" # GMSH LINEAR GRID PROPERTIES")
    println(" # N. elements       : ", mesh.nelem)
    println(" # N. points         : ", mesh.npoin_linear)
    println(" # N. edges          : ", mesh.nedges)
    println(" # N. internal edges : ", mesh.nedges_int)
    println(" # N. boundary edges : ", mesh.nedges_bdy)
    println(" # N. faces          : ", mesh.nfaces) 
    println(" # N. internal faces : ", mesh.nfaces_int)
    println(" # N. boundary faces : ", mesh.nfaces_bdy)
    println(" # GMSH LINEAR GRID PROPERTIES ...................... END")

    #
    # Coordinates arrays
    # 
    resize!(mesh.x, mesh.npoin_linear)
    resize!(mesh.y, mesh.npoin_linear)
    resize!(mesh.z, mesh.npoin_linear)
    
    #
    # Connectivity matrices
    #
    mesh.cell_node_ids     = model.grid.cell_node_ids
    mesh.conn_unique_faces = get_face_nodes(model, FACE) #faces --> 4 nodes
    mesh.conn_unique_edges = get_face_nodes(model, EDGE) #edges --> 2 nodes
    ngl = mesh.ngl
    
    resize!(mesh.conn_ho, mesh.nelem*(mesh.ngl)^(mesh.nsd))
    if (mesh.nsd == 1)
        mesh.conn_ho = reshape(mesh.conn_ho, mesh.nelem, mesh.ngl)
    elseif (mesh.nsd == 2)
        mesh.conn_ho = reshape(mesh.conn_ho, mesh.nelem, mesh.ngl, mesh.ngl)
    elseif (mesh.nsd == 3)
        mesh.conn_ho = reshape(mesh.conn_ho, mesh.nelem, mesh.ngl, mesh.ngl, mesh.ngl)

        for iel = 1:mesh.nelem
            mesh.conn_ho[iel,1,1,1]                      = mesh.cell_node_ids[iel][3]
            mesh.conn_ho[iel,mesh.ngl,1,1]               = mesh.cell_node_ids[iel][2]
            mesh.conn_ho[iel,mesh.ngl,1,mesh.ngl]        = mesh.cell_node_ids[iel][1]
            mesh.conn_ho[iel,1,1,mesh.ngl]               = mesh.cell_node_ids[iel][4]

            mesh.conn_ho[iel,1,mesh.ngl,1]               = mesh.cell_node_ids[iel][7]
            mesh.conn_ho[iel,mesh.ngl,mesh.ngl,1]        = mesh.cell_node_ids[iel][6]
            mesh.conn_ho[iel,mesh.ngl,mesh.ngl,mesh.ngl] = mesh.cell_node_ids[iel][5]
            mesh.conn_ho[iel,1,mesh.ngl,mesh.ngl]        = mesh.cell_node_ids[iel][8]
        end
    #=@info size(get_isboundary_face(topology,mesh.nsd-1))
    for i=1:length(get_isboundary_face(topology,mesh.nsd-1))
        #Get nodes of each element's face
        if get_isboundary_face(topology,mesh.nsd-1)[i] == true
    #        @info get_face_nodes(model,EDGE) #edges
        end
    end=#
    end
    
    open("./COORDS_LO.dat", "w") do f

        for ip = 1:mesh.npoin_linear
            mesh.x[ip] = model.grid.node_coordinates[ip][1]
            mesh.y[ip] = model.grid.node_coordinates[ip][2]
            mesh.z[ip] = model.grid.node_coordinates[ip][3]
            
            @printf(f, " %.6f %.6f %.6f %d\n", mesh.x[ip],  mesh.y[ip], mesh.z[ip], ip)
        end
    end

    #
    # Add high-order points to edges, faces, and elements (volumes)
    #
    
    # initialize LGL struct and buyild Gauss-Lobatto-xxx points
    Legendre = St_legendre{Float64}(0.0, 0.0, 0.0, 0.0)
    lgl      = St_lgl{Float64}(zeros(mesh.nop+1),
                               zeros(mesh.nop+1))
    build_lgl!(Legendre, lgl, mesh.nop)
    
    add_high_order_nodes_edges!(mesh, lgl)
    add_high_order_nodes_faces!(mesh, lgl)
    add_high_order_nodes_volumes!(mesh, lgl)
    
    #writevtk(model,"gmsh_grid")
end

function populate_conn_edge_el!(mesh::St_mesh)
    
    for iel = 1:mesh.nelem
        
	# Edges bottom face:
	iedg_el = 1
        mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][2]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][6]       
	iedg_el = 2
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][6]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][8]
	iedg_el = 3
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][8]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][4]
	iedg_el = 4
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][4]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][2]

        #Edges top face
	iedg_el = 5
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][1]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][5]
	iedg_el = 6
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][5]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][7]
	iedg_el = 7
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][7]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][3]
	iedg_el = 8
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][3]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][1]
	
	#Vertical edges
	iedg_el = 9
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][3]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][4]
	iedg_el = 10
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][2]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][1]
	iedg_el = 11
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][7]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][8]
	iedg_el = 12
	mesh.conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][6]
	mesh.conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][5]
        
    end
    
end #populate_edge_el!


function populate_conn_face_el!(mesh::St_mesh)

    for iel = 1:mesh.nelem
         
        #
        # Local faces node connectivity:
        # i.e. what nodes belong to a given local face in iel:
        #
        mesh.conn_face_el[iel, 1, 1] = mesh.cell_node_ids[iel][1]
        mesh.conn_face_el[iel, 1, 2] = mesh.cell_node_ids[iel][3]
        mesh.conn_face_el[iel, 1, 3] = mesh.cell_node_ids[iel][4]
        mesh.conn_face_el[iel, 1, 4] = mesh.cell_node_ids[iel][2]

        mesh.conn_face_el[iel, 3, 1] = mesh.cell_node_ids[iel][5]
        mesh.conn_face_el[iel, 3, 2] = mesh.cell_node_ids[iel][6]
        mesh.conn_face_el[iel, 3, 3] = mesh.cell_node_ids[iel][8]
        mesh.conn_face_el[iel, 3, 4] = mesh.cell_node_ids[iel][7]
        
        mesh.conn_face_el[iel, 2, 1] = mesh.cell_node_ids[iel][1]
        mesh.conn_face_el[iel, 2, 2] = mesh.cell_node_ids[iel][2]
        mesh.conn_face_el[iel, 2, 3] = mesh.cell_node_ids[iel][6]
        mesh.conn_face_el[iel, 2, 4] = mesh.cell_node_ids[iel][5]
        
        mesh.conn_face_el[iel, 4, 1] = mesh.cell_node_ids[iel][3]
        mesh.conn_face_el[iel, 4, 2] = mesh.cell_node_ids[iel][7]
        mesh.conn_face_el[iel, 4, 3] = mesh.cell_node_ids[iel][8]
        mesh.conn_face_el[iel, 4, 4] = mesh.cell_node_ids[iel][4]
        
        mesh.conn_face_el[iel, 5, 1] = mesh.cell_node_ids[iel][2]
        mesh.conn_face_el[iel, 5, 2] = mesh.cell_node_ids[iel][4]
        mesh.conn_face_el[iel, 5, 3] = mesh.cell_node_ids[iel][8]
        mesh.conn_face_el[iel, 5, 4] = mesh.cell_node_ids[iel][6]
        
        mesh.conn_face_el[iel, 6, 1] = mesh.cell_node_ids[iel][1]
        mesh.conn_face_el[iel, 6, 2] = mesh.cell_node_ids[iel][5]
        mesh.conn_face_el[iel, 6, 3] = mesh.cell_node_ids[iel][7]
        mesh.conn_face_el[iel, 6, 4] = mesh.cell_node_ids[iel][3]
        
    end
    
end #populate_face_el


function populate_face_in_elem!(face_in_elem::Array{Int64, 3}, nelem, NFACES_EL, conn_face_el_sort::Array{Int64, 3})

    iface = Int64(0)
    for ifac = 1:NFACES_EL
        for iel = 1:nelem	    
	    for jfac = 1:NFACES_EL
                for jel = iel:nelem
		    
		    if(     conn_face_el_sort[iel,ifac,1] === conn_face_el_sort[jel,jfac,1] && 
			    conn_face_el_sort[iel,ifac,2] === conn_face_el_sort[jel,jfac,2] && 
			    conn_face_el_sort[iel,ifac,3] === conn_face_el_sort[jel,jfac,3] && 
			    conn_face_el_sort[iel,ifac,4] === conn_face_el_sort[jel,jfac,4] && 
			    iel ≠ jel) 
			
			face_in_elem[iel,ifac,1] = iel
			face_in_elem[iel,ifac,2] = jel
			
			face_in_elem[jel,jfac,1] = jel
			face_in_elem[jel,jfac,2] = iel
			
			#@info " SHARED FACE:  face %d of ELEMENT %d is shared with face %d of ELEMENT %d - (%d %d %d %d) = (%d %d %d %d)\n", ifac+1,iel+1, jfac+1, jel+1, conn_face_el_sort[iel,ifac,1], conn_face_el_sort[iel,ifac,2], conn_face_el_sort[iel,ifac,3], conn_face_el_sort[iel,ifac,4],  conn_face_el_sort[jel,jfac,1], conn_face_el_sort[jel,jfac,2], conn_face_el_sort[jel,jfac,3], conn_face_el_sort[jel,jfac,4]
			iface = iface + 1;
		    end
		end
	    end
	end
    end     
    nfaces_int = Int64(iface)
    
end

function  add_high_order_nodes!(mesh::St_mesh) end

function  add_high_order_nodes_edges!(mesh::St_mesh, lgl::St_lgl)

    if (mesh.nop < 2) return end
    
    x1, y1, z1 = Float64(0.0), Float64(0.0), Float64(0.0)
    x2, y2, z2 = Float64(0.0), Float64(0.0), Float64(0.0)
    
    ξ::typeof(lgl.ξ[1]) = 0.0

    ngl                      = mesh.nop + 1
    tot_linear_poin          = mesh.npoin_linear
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    tot_faces_internal_nodes = mesh.nfaces*(ngl-2)*(ngl-2)
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)*(ngl-2)*(ngl-2)

    #Increase number of grid points from linear count to titak high-order points
    mesh.npoin = mesh.npoin_linear + tot_edges_internal_nodes + tot_faces_internal_nodes + tot_vol_internal_nodes

    if length(mesh.x_ho) < mesh.npoin
        resize!(mesh.x_ho, mesh.npoin)
    end
    if length(mesh.y_ho) < mesh.npoin        
        resize!(mesh.y_ho, mesh.npoin)
    end
    if length(mesh.z_ho) < mesh.npoin        
        resize!(mesh.z_ho, mesh.npoin)
    end
     
    open("./COORDS_HO_edges.dat", "w") do f

        ip  = tot_linear_poin + 1
        for iedge_g = 1:mesh.nedges
            
            ip1 = min(mesh.conn_unique_edges[iedge_g][1], mesh.conn_unique_edges[iedge_g][2])
            ip2 = max(mesh.conn_unique_edges[iedge_g][1], mesh.conn_unique_edges[iedge_g][2])

            x1, y1, z1 = mesh.x[ip1], mesh.y[ip1], mesh.z[ip1]
            x2, y2, z2 = mesh.x[ip2], mesh.y[ip2], mesh.z[ip2]

            for l=2:ngl-1
                ξ = lgl.ξ[l];
                
                mesh.x_ho[ip] = x1*(1.0 - ξ)*0.5 + x2*(1.0 + ξ)*0.5;
	        mesh.y_ho[ip] = y1*(1.0 - ξ)*0.5 + y2*(1.0 + ξ)*0.5;
	        mesh.z_ho[ip] = z1*(1.0 - ξ)*0.5 + z2*(1.0 + ξ)*0.5;
                
                @printf(f, " %.6f %.6f %.6f %d\n", mesh.x_ho[ip],  mesh.y_ho[ip], mesh.z_ho[ip], ip)
                
                ip = ip + 1
            end
        end
        
    end #do f
    
    println(" # POPULATE GRID with SPECTRAL NODES ............................ EDGES DONE")
    
end

function  add_high_order_nodes_faces!(mesh::St_mesh, lgl::St_lgl)

    if (mesh.nop < 2) return end
    
    
    x1, y1, z1 = Float64(0.0), Float64(0.0), Float64(0.0)
    x2, y2, z2 = Float64(0.0), Float64(0.0), Float64(0.0)
    x3, y3, z3 = Float64(0.0), Float64(0.0), Float64(0.0)
    x4, y4, z4 = Float64(0.0), Float64(0.0), Float64(0.0)
    
    ξ::typeof(lgl.ξ[1]) = 0.0
    ζ::typeof(lgl.ξ[1]) = 0.0
    
    ngl                      = mesh.nop + 1
    tot_linear_poin          = mesh.npoin_linear
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    tot_faces_internal_nodes = mesh.nfaces*(ngl-2)*(ngl-2)
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)*(ngl-2)*(ngl-2)

    #Increase number of grid points from linear count to titak high-order points
    mesh.npoin = tot_linear_poin + tot_edges_internal_nodes + tot_faces_internal_nodes + tot_vol_internal_nodes

    if length(mesh.x_ho) < mesh.npoin
        resize!(mesh.x_ho, mesh.npoin)
    end
    if length(mesh.y_ho) < mesh.npoin
        resize!(mesh.y_ho, mesh.npoin)
    end
    if length(mesh.z_ho) < mesh.npoin
        resize!(mesh.z_ho, mesh.npoin)
    end
    
    open("./COORDS_HO_faces.dat", "w") do f
        
        ip  = tot_linear_poin + tot_edges_internal_nodes + 1
        for iface_g = 1:mesh.nfaces

            #GGNS numbering
            @show ip1 = mesh.conn_unique_faces[iface_g][1]
            @show ip2 = mesh.conn_unique_faces[iface_g][3]
            @show ip3 = mesh.conn_unique_faces[iface_g][4]
            @show ip4 = mesh.conn_unique_faces[iface_g][2]
            
            x1, y1, z1 = mesh.x[ip1], mesh.y[ip1], mesh.z[ip1]
            x2, y2, z2 = mesh.x[ip2], mesh.y[ip2], mesh.z[ip2]
            x3, y3, z3 = mesh.x[ip3], mesh.y[ip3], mesh.z[ip3]
            x4, y4, z4 = mesh.x[ip4], mesh.y[ip4], mesh.z[ip4]
            
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
                    
                    mesh.z_ho[ip] =  (z1*(1 - ξ)*(1 - ζ)*0.25
		                      + z2*(1 + ξ)*(1 - ζ)*0.25
		                      + z3*(1 + ξ)*(1 + ζ)*0.25
		                      + z4*(1 - ξ)*(1 + ζ)*0.25)

                    @printf(f, " %.6f %.6f %.6f %d\n", mesh.x_ho[ip],  mesh.y_ho[ip], mesh.z_ho[ip], ip)
                    
	            ip = ip + 1
                end
            end
        end 
    end #file
    
    println(" # POPULATE GRID with SPECTRAL NODES ............................ FACES DONE")
    
end


function  add_high_order_nodes_volumes!(mesh::St_mesh, lgl::St_lgl)

    if (mesh.nop < 2) return end
    
    
    x1, y1, z1 = Float64(0.0), Float64(0.0), Float64(0.0)
    x2, y2, z2 = Float64(0.0), Float64(0.0), Float64(0.0)
    x3, y3, z3 = Float64(0.0), Float64(0.0), Float64(0.0)
    x4, y4, z4 = Float64(0.0), Float64(0.0), Float64(0.0)
    x5, y5, z5 = Float64(0.0), Float64(0.0), Float64(0.0)
    x6, y6, z6 = Float64(0.0), Float64(0.0), Float64(0.0)
    x7, y7, z7 = Float64(0.0), Float64(0.0), Float64(0.0)
    x8, y8, z8 = Float64(0.0), Float64(0.0), Float64(0.0)
    
    ξ::typeof(lgl.ξ[1]) = 0.0
    η::typeof(lgl.ξ[1]) = 0.0
    ζ::typeof(lgl.ξ[1]) = 0.0

    ngl                      = mesh.nop + 1
    tot_linear_poin          = mesh.npoin_linear
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    tot_faces_internal_nodes = mesh.nfaces*(ngl-2)*(ngl-2)
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)*(ngl-2)*(ngl-2)

    #Increase number of grid points from linear count to titak high-order points
    mesh.npoin = tot_linear_poin + tot_edges_internal_nodes + tot_faces_internal_nodes + tot_vol_internal_nodes
    
    if length(mesh.x_ho) < mesh.npoin
        resize!(mesh.x_ho, mesh.npoin)
    end
    if length(mesh.y_ho) < mesh.npoin
        resize!(mesh.y_ho, mesh.npoin)
    end
    if length(mesh.z_ho) < mesh.npoin
        resize!(mesh.z_ho, mesh.npoin)
    end
    
    open("./COORDS_HO_vol.dat", "w") do f

        ip  = tot_linear_poin + tot_edges_internal_nodes + tot_faces_internal_nodes + 1
        for iel_g = 1:mesh.nelem

            #CGNS numbering
            @show ip1 = mesh.cell_node_ids[iel_g][2]
            @show ip2 = mesh.cell_node_ids[iel_g][6]
            @show ip3 = mesh.cell_node_ids[iel_g][8]
            @show ip4 = mesh.cell_node_ids[iel_g][4]
            @show ip5 = mesh.cell_node_ids[iel_g][1]
            @show ip6 = mesh.cell_node_ids[iel_g][5]
            @show ip7 = mesh.cell_node_ids[iel_g][7]
            @show ip8 = mesh.cell_node_ids[iel_g][3]
            
            
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
		        
                        @printf(f, " %.6f %.6f %.6f %d\n", mesh.x_ho[ip],  mesh.y_ho[ip], mesh.z_ho[ip], ip)
                        
	                ip = ip + 1
                    end
                end
            end 
        end
    end #file
    
    println(" # POPULATE GRID with SPECTRAL NODES ............................ VOLUMES DONE")
    
end

#=
function mod_mesh_cgns_ordering!(cell_node_ids::Table{Int64,Vector{Int64},Vector{Int64}})

    nelem     = Int64(size(cell_node_ids, 1))
    nnodes_el = Int64(size(cell_node_ids[1], 1))
    
    temp1 = Int64(1)
    temp2 = Int64(1)
    temp3 = Int64(1)
    temp4 = Int64(1)
    temp5 = Int64(1)
    temp6 = Int64(1)
    temp7 = Int64(1)
    temp8 = Int64(1)
    
    @info " before "
    for iel = 1:1
        
        @info cell_node_ids[iel][1], " ", cell_node_ids[iel][2], " ",
        cell_node_ids[iel][3], " ", cell_node_ids[iel][4], " ",
        cell_node_ids[iel][5], " ", cell_node_ids[iel][6], " ",
        cell_node_ids[iel][7], " ", cell_node_ids[iel][8]
    end
    @info " CGNS"
    for iel = 1:1
        
        temp1 = copy(cell_node_ids[iel][8]);
	temp2 = copy(cell_node_ids[iel][6]);
	temp3 = copy(cell_node_ids[iel][5]);
	temp4 = copy(cell_node_ids[iel][7]);
        temp5 = copy(cell_node_ids[iel][4]);
        temp6 = copy(cell_node_ids[iel][2]);
        temp7 = copy(cell_node_ids[iel][1]);
        temp8 = copy(cell_node_ids[iel][3]);

        @info temp1, " t1 ? cell[8] ", cell_node_ids[iel][8];
        @info temp2, " t1 ? cell[6] ", cell_node_ids[iel][6];
        @info temp3, " t1 ? cell[5] ", cell_node_ids[iel][5];
        @info temp4, " t1 ? cell[7] ", cell_node_ids[iel][7];
        @info temp5, " t1 ? cell[4] ", cell_node_ids[iel][4];
        @info temp6, " t1 ? cell[2] ", cell_node_ids[iel][2];
        @info temp7, " t1 ? cell[1] ", cell_node_ids[iel][1];
        @info temp8, " t1 ? cell[3] ", cell_node_ids[iel][3];
        @info " === "
        
	#Rewrite cell_node_ids
	cell_node_ids[iel][1] = copy(temp1);
	cell_node_ids[iel][2] = copy(temp2);
	cell_node_ids[iel][3] = copy(temp3);
	cell_node_ids[iel][4] = copy(temp4);
	cell_node_ids[iel][5] = copy(temp5);
	cell_node_ids[iel][6] = copy(temp6);
	cell_node_ids[iel][7] = copy(temp7);
	cell_node_ids[iel][8] = copy(temp8);

        @info temp1, " t1 ? cell[1] ", cell_node_ids[iel][1];
        @info temp2, " t1 ? cell[2] ", cell_node_ids[iel][2];
        @info temp3, " t1 ? cell[3] ", cell_node_ids[iel][3];
        @info temp4, " t1 ? cell[4] ", cell_node_ids[iel][4];
        @info temp5, " t1 ? cell[5] ", cell_node_ids[iel][5];
        @info temp6, " t1 ? cell[6] ", cell_node_ids[iel][6];
        @info temp7, " t1 ? cell[7] ", cell_node_ids[iel][7];
        @info temp8, " t1 ? cell[8] ", cell_node_ids[iel][8];
        
    end
     
    @info " After reorder "
    for iel = 1:1
        
        @info cell_node_ids[iel][1], " ", cell_node_ids[iel][2], " ",
        cell_node_ids[iel][3], " ", cell_node_ids[iel][4], " ",
        cell_node_ids[iel][5], " ", cell_node_ids[iel][6], " ",
        cell_node_ids[iel][7], " ", cell_node_ids[iel][8]
    end
    @info " "
    
end
=#

#=function mod_mesh_build_edges_faces!(mesh::St_mesh)
    
    if mesh.nsd == 3
        mesh.NEDGES_EL  = 12
        mesh.NFACES_EL  = 6
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 4
    elseif mesh.nsd == 2
        mesh.NEDGES_EL  = 4
        mesh.NFACES_EL  = 1
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 4
    elseif mesh.nsd == 1
        mesh.NEDGES_EL  = 1
        mesh.NFACES_EL  = 0
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 0
    else
        error( " WRONG NSD: This is not theoretical physics: we only handle 1, 2, or 3 dimensions!")
    end
    
    mesh.conn_edge_el = Array{Int64, 3}(undef,  mesh.nelem, mesh.NEDGES_EL, mesh.EDGE_NODES)
    mesh.conn_face_el = Array{Int64, 3}(undef,  mesh.nelem, mesh.NFACES_EL, mesh.FACE_NODES)
    mesh.face_in_elem = Array{Int64, 3}(undef,  mesh.nelem, mesh.NFACES_EL, mesh.FACE_NODES)
    conn_face_el_sort = Array{Int64, 3}(undef,  mesh.nelem, mesh.NEDGES_EL, mesh.EDGE_NODES)
  
    if (mesh.nsd == 1)
    
        #mesh.conn_ho = Array{Int64, 4}(undef,  mesh.nelem, 1)        
    elseif (mesh.nsd == 2)
        @error " Only 3D is currently coded. 2D is not currently supported."
        #mesh.conn_ho = Array{Int64, 4}(undef,  mesh.nelem, mesh.nop+1, mesh.nop+1)
    elseif (mesh.nsd == 3) 
        mesh.conn_ho = Array{Int64, 4}(undef,  mesh.nelem, mesh.nop+1, mesh.nop+1, mesh.nop+1)
    end
    
    populate_conn_edge_el!(mesh)
    populate_conn_face_el!(mesh)
    
    #local sorting for comparison done later
    #conn_face_el_sort = copy(mesh.conn_face_el)
    #sort!(conn_face_el_sort, dims = 3)

    #count internal and boundary nfaces_int, nfaces_bdy
    #populate_face_in_elem!(mesh.face_in_elem, mesh.nelem, mesh.NFACES_EL, conn_face_el_sort)

    #add high order nodes to all unique edges and faces
end
=#
