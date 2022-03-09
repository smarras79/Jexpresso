using Gridap
using Gridap.Arrays: Table
using GridapGmsh
using LinearAlgebra
using Revise

export St_mesh

export mod_mesh_build_mesh!
export mod_mesh_read_gmsh!

const VERTEX_NODES = Int32(1)
const EDGE_NODES   = Int32(2)
const FACE_NODES   = Int32(4)

Base.@kwdef mutable struct St_mesh{TInt, TFloat}
    
    x::Union{Array{TFloat}, Missing} = zeros(2)
    y::Union{Array{TFloat}, Missing} = zeros(2)
    z::Union{Array{TFloat}, Missing} = zeros(2)

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
    npoin::Union{TInt, Missing} = 1
    
    nsd::Union{TInt, Missing} = 1
    nop::Union{TInt, Missing} = 4

    #low and high order connectivity tables
    cell_node_ids::Table{Int32,Vector{Int32},Vector{Int32}} = Gridap.Arrays.Table(zeros(nelem), zeros(npoin))
    cell_node_ids_ho::Table{Int32,Vector{Int32},Vector{Int32}} = Gridap.Arrays.Table(zeros(nelem), zeros(npoin))
    
    conn_edge_el = Array{Int32, 3}(undef,  nelem, 12, 2)
    conn_face_el = Array{Int32, 3}(undef,  nelem, 6,  4)
    face_in_elem = Array{Int32, 3}(undef,  nelem, 6, 2)
    
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
    
    model      = GmshDiscreteModel(gmsh_filename)
    mesh.npoin = length(model.grid.node_coordinates)
    mesh.nelem = length(model.grid.cell_node_ids)

    mesh.cell_node_ids = model.grid.cell_node_ids
    
    @info model.grid.cell_node_ids
    #@info length(model.grid.cell_node_ids)
    
    resize!(mesh.x, mesh.npoin)
    resize!(mesh.y, mesh.npoin)
    resize!(mesh.z, mesh.npoin)

    for ip = 1:mesh.npoin
        mesh.x[ip] = model.grid.node_coordinates[ip][1]
        mesh.y[ip] = model.grid.node_coordinates[ip][2]
        mesh.z[ip] = model.grid.node_coordinates[ip][3]
    end

    mod_mesh_build_edges_faces!(mesh)
    
    #writevtk(model,"gmsh_grid")
end

function mod_mesh_build_edges_faces!(mesh::St_mesh)

    #### Convert to CGNS numbering
    #### not correct; it doesnt modify mesh.cell_node_
    #### @assert "not working yet"
    ###mod_mesh_cgns_ordering!(mesh.cell_node_ids)
    
    if mesh.nsd == 3
        NEDGES_EL = 12
        NFACES_EL = 6
        EDGE_NODES = 2
        FACE_NODES = 4
    elseif mesh.nsd == 2
        NEDGES_EL = 4
        NFACES_EL = 1
        EDGE_NODES = 2
        FACE_NODES = 4
    elseif mesh.nsd == 1
        NEDGES_EL = 1
        NFACES_EL = 0
        EDGE_NODES = 2
        FACE_NODES = 0
    else
        error( " WRONG NSD: This is not theoretical physics: we only handle 1, 2, or 3 dimensions!")
    end
    
    mesh.conn_edge_el = Array{Int32, 3}(undef,  mesh.nelem, NEDGES_EL, EDGE_NODES)
    mesh.conn_face_el = Array{Int32, 3}(undef,  mesh.nelem, NFACES_EL, FACE_NODES)
    mesh.face_in_elem = Array{Int32, 3}(undef,  mesh.nelem, NFACES_EL, FACE_NODES)    
    mesh.face_in_elem = Array{Int32, 3}(undef,  mesh.nelem, NFACES_EL, 2)
    
    conn_face_el_sort = Array{Int32, 3}(undef,  mesh.nelem, NEDGES_EL, EDGE_NODES)
    
    populate_conn_edge_el!(mesh)
    populate_conn_face_el!(mesh)

    #local sorting for comparison done later
    conn_face_el_sort = copy(mesh.conn_face_el)
    sort!(conn_face_el_sort, dims = 3)
    
    nint_faces = populate_face_in_elem!(mesh.face_in_elem, mesh.nelem, NFACES_EL, conn_face_el_sort)
    println(" # Ninternal faces: ", nint_faces)
    
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
        mesh.conn_face_el[iel,1,1] = mesh.cell_node_ids[iel][1]
        mesh.conn_face_el[iel,1,2] = mesh.cell_node_ids[iel][3]
        mesh.conn_face_el[iel,1,3] = mesh.cell_node_ids[iel][4]
        mesh.conn_face_el[iel,1,4] = mesh.cell_node_ids[iel][2]

        mesh.conn_face_el[iel,3,1] = mesh.cell_node_ids[iel][5]
        mesh.conn_face_el[iel,3,2] = mesh.cell_node_ids[iel][6]
        mesh.conn_face_el[iel,3,3] = mesh.cell_node_ids[iel][8]
        mesh.conn_face_el[iel,3,4] = mesh.cell_node_ids[iel][7]
        
        mesh.conn_face_el[iel,2,1] = mesh.cell_node_ids[iel][1]
        mesh.conn_face_el[iel,2,2] = mesh.cell_node_ids[iel][2]
        mesh.conn_face_el[iel,2,3] = mesh.cell_node_ids[iel][6]
        mesh.conn_face_el[iel,2,4] = mesh.cell_node_ids[iel][5]
        
        mesh.conn_face_el[iel,4,1] = mesh.cell_node_ids[iel][3]
        mesh.conn_face_el[iel,4,2] = mesh.cell_node_ids[iel][7]
        mesh.conn_face_el[iel,4,3] = mesh.cell_node_ids[iel][8]
        mesh.conn_face_el[iel,4,4] = mesh.cell_node_ids[iel][4]
        
        mesh.conn_face_el[iel,5,1] = mesh.cell_node_ids[iel][2]
        mesh.conn_face_el[iel,5,2] = mesh.cell_node_ids[iel][4]
        mesh.conn_face_el[iel,5,3] = mesh.cell_node_ids[iel][8]
        mesh.conn_face_el[iel,5,4] = mesh.cell_node_ids[iel][6]
        
        mesh.conn_face_el[iel,6,1] = mesh.cell_node_ids[iel][1]
        mesh.conn_face_el[iel,6,2] = mesh.cell_node_ids[iel][5]
        mesh.conn_face_el[iel,6,3] = mesh.cell_node_ids[iel][7]
        mesh.conn_face_el[iel,6,4] = mesh.cell_node_ids[iel][3]
        
    end
    
end #populate_face_el

function populate_face_in_elem!(face_in_elem::Array{Int32, 3}, nelem, NFACES_EL, conn_face_el_sort::Array{Int32, 3})

    iface = Int32(0)
    for iel = 1:nelem
	for ifac = 1:NFACES_EL
	    for jel = iel:nelem
		for jfac = 1:NFACES_EL
		    
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
    
    nint_faces = Int32(iface)
    return nint_faces
end


function mod_mesh_cgns_ordering!(cell_node_ids::Table{Int32,Vector{Int32},Vector{Int32}})

    nelem     = Int32(size(cell_node_ids, 1))
    nnodes_el = Int32(size(cell_node_ids[1], 1))
    
    temp1 = Int32(1)
    temp2 = Int32(1)
    temp3 = Int32(1)
    temp4 = Int32(1)
    temp5 = Int32(1)
    temp6 = Int32(1)
    temp7 = Int32(1)
    temp8 = Int32(1)
    
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
