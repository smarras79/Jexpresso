module Mod_grid

export St_grid

struct St_grid{T}
    x::T
    y::T
    z::T

    xmin::T, xmax::T
    ymin::T, ymax::T
    zmin::T, zmax::T

    npx::Int8
    npy::Int8
    npy::Int8
    
end

function build_grid2d(mesh::St_grid)

    mesh.x = range(mesh.xmin, mesh.xmax, length = mesh.npx)
    mesh.y = range(mesh.ymin, mesh.ymax, length = mesh.npy)
    mesh.z = range(mesh.zmin, mesh.zmax, length = mesh.npz)

    
end
