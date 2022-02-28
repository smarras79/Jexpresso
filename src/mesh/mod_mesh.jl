using Revise

export St_mesh
export mod_mesh_build_mesh!

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

#end #module

