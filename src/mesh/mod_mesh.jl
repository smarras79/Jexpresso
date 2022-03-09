using Gridap
using Gridap.Arrays: Table
using GridapGmsh
using LinearAlgebra
using Revise

export St_mesh

export mod_mesh_build_mesh!
export mod_mesh_read_gmsh!

const NPOIN_EL3D  = Int32(8)
const NEDGES_EL3D = Int32(12)
const NFACES_EL3D = Int32(6)

const NPOIN_EL2D  = Int32(4)
const NEDGES_EL2D = Int32(4)
const NFACES_EL2D = Int32(0)

const NPOIN_EL1D  = Int32(2)
const NEDGES_EL1D = Int32(1)
const NFACES_EL1D = Int32(0)

const VERTEX = Int32(1)
const EDGE   = Int32(2)

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
    
    
    conn_edge_el = Array{Int32, 3}(undef,  mesh.nelem, NEDGES_EL3D, EDGE)
    @info size(conn_edge_el)
    @info typeof(conn_edge_el)
    for iel = 1:2 #mesh.nelem
        
	# Edges bottom face:
	iedg_el = 1
        conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][1]
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][2]       
	iedg_el = 2;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][2];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][3];
	iedg_el = 3;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][3];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][4];
	iedg_el = 4;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][4];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][1];

        #Edges top face
	iedg_el = 5;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][5];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][6];
	iedg_el = 6;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][6];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][7];
	iedg_el = 7;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][7];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][8];
	iedg_el = 8;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][8];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][5];
	
	#Vertical edges
	iedg_el = 9;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][1];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][5];
	iedg_el = 10;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][2];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][6];
	iedg_el = 11;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][3];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][7];
	iedg_el = 12;
	conn_edge_el[iel,iedg_el,1] = mesh.cell_node_ids[iel][4];
	conn_edge_el[iel,iedg_el,2] = mesh.cell_node_ids[iel][8];
        
        #@info " mesh.cell_node_ids[iel][1] = ",  iel, " ", conn_edge_el[iel,iedg_el,1]
        
     end
    
end

function mod_mesh_!(mesh::St_mesh)

end

#end #module
