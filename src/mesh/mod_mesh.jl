using Gridap
using Gridap.Arrays: Table
using GridapGmsh
using Revise

export St_mesh

export mod_mesh_build_mesh!
export mod_mesh_read_gmsh!


mutable struct St_mesh{TInt, TFloat}
    
    x::Array{TFloat}
    y::Array{TFloat}
    z::Array{TFloat}

    xmin::TFloat; xmax::TFloat
    ymin::TFloat; ymax::TFloat
    zmin::TFloat; zmax::TFloat
    
    npx::TInt
    npy::TInt
    npz::TInt
    
    npoin::TInt
    
    nsd::TInt
    nop::TInt

end #St_mesh

mutable struct St_conn{}
    
    cell_node_ids_ho::Table{Int32,Vector{Int32},Vector{Int32}}
    
end #St_conn

function mod_mesh_build_mesh!(mesh::St_mesh)

    T = TFloat
    Δx::T=0
    Δy::T=0
    Δz::T=0
    nsd = mesh.nsd
    
    if (nsd == 1)
        mesh.npy = mesh.npz = 1
        Δx = abs(mesh.xmax - mesh.xmin)/(mesh.npx - 1)
        Δy = 0.0
        Δz = 0.0
    elseif (nsd == 2)
        mesh.npz = 1
        Δx = abs(mesh.xmax - mesh.xmin)/(mesh.npx - 1)
        Δy = abs(mesh.ymax - mesh.ymin)/(mesh.npy - 1)
        Δz = 0
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
    
    model = GmshDiscreteModel(gmsh_filename)
    mesh.npoin = length(model.grid.node_coordinates)
    @info model.grid.cell_node_ids #CONNN QUI
    npoin = mesh.npoin
    @info typeof(model.grid.cell_node_ids)
    conn = typeof(model.grid.cell_node_ids)
    @info model.grid.cell_node_ids[1][1:4]
    
    resize!(mesh.x, mesh.npoin)
    resize!(mesh.y, mesh.npoin)
    resize!(mesh.z, mesh.npoin)

    for ip = 1:npoin
        mesh.x[ip] = model.grid.node_coordinates[ip][1]
        mesh.y[ip] = model.grid.node_coordinates[ip][2]
        mesh.z[ip] = model.grid.node_coordinates[ip][3]
    end
        
    writevtk(model,"gmsh_grid")
end

#end #module
