#---------------------------------------------------------------------------------
# write_vtk_sphere.jl
#
# VTK output of the HIGH-ORDER spherical-shell grid built by
# src/kernel/mesh/sphere_mesh.jl.
#
# Each spectral element (ngl × ngl LGL nodes) is written as (ngl-1)² linear
# VTK_QUAD sub-cells whose corners ARE the LGL nodes. This is the same
# convention as write_vtk_grid_only() for the flat grids, and it is the one
# that lets you SEE the actual node distribution in ParaView (a
# VTK_LAGRANGE_QUADRILATERAL would be re-tessellated by ParaView and would
# hide exactly the thing we want to inspect).
#
# What is written, and why it is the useful thing to look at when debugging
# the numbering of a CLOSED shell:
#
#   point data
#     ip         global node index — colour by it and the seams between gmsh
#                panels must show a CONTINUOUS field: a jump means the two
#                sides of the seam carry different node numbers, i.e. the
#                shell was silently torn open.
#     node_type  0 = linear vertex, 1 = edge node, 2 = element interior.
#     lon, lat   degrees.
#     radius     |x|; must be flat at R everywhere.
#   cell data
#     iel        owning element
#     panel      gmsh physical/entity tag (the cubed-sphere panel)
#
# Also provided: write_vtk_sphere_edges(), which writes ONLY the unique-edge
# node lines. Opened together with the surface it makes a torn seam obvious
# at a glance.
#---------------------------------------------------------------------------------

export write_vtk_sphere_grid
export write_vtk_sphere_edges


function write_vtk_sphere_grid(mesh::St_mesh_sphere,
                               file_name::String,
                               OUTPUT_DIR::String;
                               verbose = true)

    if !isdir(OUTPUT_DIR)
        mkpath(OUTPUT_DIR)
    end

    ngl    = mesh.ngl
    nelem  = mesh.nelem
    nsub   = nelem*(ngl-1)*(ngl-1)

    cells    = [MeshCell(VTKCellTypes.VTK_QUAD, Int64[1, 2, 3, 4]) for _ = 1:nsub]
    cell_iel = Vector{Float64}(undef, nsub)
    cell_pan = Vector{Float64}(undef, nsub)

    isel = 1
    for iel = 1:nelem
        for j = 1:ngl-1
            for i = 1:ngl-1
                #
                # Counter-clockwise as seen from OUTSIDE the sphere, because
                # the elements were re-oriented that way in
                # orient_elements_outward!. Keeping that order here means the
                # VTK normals also point outward and ParaView's "backface
                # culling"/lighting shows a lit sphere, not a black one.
                #
                ip1 = mesh.connijk[iel, i,   j  ]
                ip2 = mesh.connijk[iel, i+1, j  ]
                ip3 = mesh.connijk[iel, i+1, j+1]
                ip4 = mesh.connijk[iel, i,   j+1]

                cells[isel]    = MeshCell(VTKCellTypes.VTK_QUAD, Int64[ip1, ip2, ip3, ip4])
                cell_iel[isel] = Float64(iel)
                cell_pan[isel] = Float64(mesh.elem_tag[iel])
                isel += 1
            end
        end
    end

    fout_name = string(OUTPUT_DIR, "/", file_name)

    vtkf = vtk_grid(fout_name,
                    mesh.x[1:mesh.npoin],
                    mesh.y[1:mesh.npoin],
                    mesh.z[1:mesh.npoin],
                    cells, compress = false)

    vtkf["ip",        VTKPointData()] = Float64.(collect(1:mesh.npoin))
    vtkf["node_type", VTKPointData()] = Float64.(mesh.node_type)
    vtkf["lon",       VTKPointData()] = mesh.lon .* (180.0/π)
    vtkf["lat",       VTKPointData()] = mesh.lat .* (180.0/π)
    vtkf["radius",    VTKPointData()] = sqrt.(mesh.x.^2 .+ mesh.y.^2 .+ mesh.z.^2)

    vtkf["iel",   VTKCellData()] = cell_iel
    vtkf["panel", VTKCellData()] = cell_pan

    out = vtk_save(vtkf)

    verbose && println(" # Wrote high-order spherical shell grid: ", fout_name, ".vtu  (",
                       mesh.npoin, " nodes, ", nsub, " sub-cells)")

    return out
end


#---------------------------------------------------------------------------------
# The unique-edge skeleton: one polyline of ngl nodes per unique edge.
#
# Colour it by `iedge` and every seam of the closed shell must be covered by a
# SINGLE edge — two overlapping polylines with different ids there mean the
# panel-seam vertices were never merged.
#---------------------------------------------------------------------------------
function write_vtk_sphere_edges(mesh::St_mesh_sphere,
                                file_name::String,
                                OUTPUT_DIR::String;
                                verbose = true)

    if !isdir(OUTPUT_DIR)
        mkpath(OUTPUT_DIR)
    end

    ngl    = mesh.ngl
    nseg   = mesh.nedges*(ngl-1)
    cells  = [MeshCell(VTKCellTypes.VTK_LINE, Int64[1, 2]) for _ = 1:nseg]
    cid    = Vector{Float64}(undef, nseg)

    is = 1
    for ie = 1:mesh.nedges
        for l = 1:ngl-1
            cells[is] = MeshCell(VTKCellTypes.VTK_LINE,
                                 Int64[mesh.poin_in_edge[ie, l], mesh.poin_in_edge[ie, l+1]])
            cid[is]   = Float64(ie)
            is += 1
        end
    end

    fout_name = string(OUTPUT_DIR, "/", file_name)

    vtkf = vtk_grid(fout_name,
                    mesh.x[1:mesh.npoin],
                    mesh.y[1:mesh.npoin],
                    mesh.z[1:mesh.npoin],
                    cells, compress = false)

    vtkf["iedge", VTKCellData()]  = cid
    vtkf["ip",    VTKPointData()] = Float64.(collect(1:mesh.npoin))

    out = vtk_save(vtkf)

    verbose && println(" # Wrote spherical shell edge skeleton:  ", fout_name, ".vtu  (",
                       mesh.nedges, " unique edges)")

    return out
end
