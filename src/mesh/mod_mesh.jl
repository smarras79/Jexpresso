using Revise

export St_mesh
export mod_mesh_build_mesh2d!

mutable struct St_mesh{TInt,TFloat}
    
    x::Array{TFloat}
    y::Array{TFloat}
    z::Array{TFloat}

    xmin::TFloat; xmax::TFloat
    ymin::TFloat; ymax::TFloat
    zmin::TFloat; zmax::TFloat
    
    npx::TInt
    npy::TInt
    npz::TInt

end #St_mesh

    
function mod_mesh_build_mesh2d!(mesh::St_mesh)
    
    Δx = abs(mesh.xmax - mesh.xmin)/(mesh.npx - 1)
    Δy = abs(mesh.ymax - mesh.ymin)/(mesh.npy - 1)
    Δz = abs(mesh.zmax - mesh.zmin)/(mesh.npz - 1)
    
    for i = 1:mesh.npx
        mesh.x[i] = (i - 1)*Δx
    end
    for j = 1:mesh.npy
        mesh.y[j] = (j - 1)*Δy
    end
    for k = 1:mesh.npz
        mesh.z[k] = (k - 1)*Δz
    end
    
end #build_mesh2d!

#end #module

