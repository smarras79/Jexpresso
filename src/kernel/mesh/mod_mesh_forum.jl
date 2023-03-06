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
using ElasticArrays
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


include("../bases/basis_structs.jl")

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
    npoin_el::Union{TInt, Missing} = 1 #Total number of points in the reference element
    
    NNODES_EL::Union{TInt, Missing}  =  2^nsd
    NEDGES_EL::Union{TInt, Missing}  = 12
    NFACES_EL::Union{TInt, Missing}  =  6
    EDGE_NODES::Union{TInt, Missing} =  2
    FACE_NODES::Union{TInt, Missing} =  4

    
    #low and high order connectivity tables
    cell_node_ids::Table{Int64,Vector{Int64},Vector{Int64}} = Gridap.Arrays.Table(zeros(nelem), zeros(8))
    cell_node_ids_ho::Table{Int64,Vector{Int64},Vector{Int64}} = Gridap.Arrays.Table(zeros(nelem), zeros(8))
    
    conn_ptr       = ElasticArray{Int64}(undef, nelem)    
    conn           = ElasticArray{Int64}(undef, ngl*nelem)
    conn_unique_edges = ElasticArray{Int64}(undef,  1, 2)
    conn_unique_faces = ElasticArray{Int64}(undef,  1, 4)

    conn_edge_L2G     = ElasticArray{Int64}(undef, 1, NEDGES_EL, nelem)
    conn_face_L2G     = ElasticArray{Int64}(undef, 1, NFACES_EL, nelem)
    
    conn_edge_el      = ElasticArray{Int64}(undef, 2, NEDGES_EL, nelem)
    conn_face_el      = ElasticArray{Int64}(undef, 4, NFACES_EL, nelem)
    face_in_elem      = ElasticArray{Int64}(undef, 2, NFACES_EL, nelem)
    
end

#=
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

end=#

function mod_mesh_read_gmsh!(mesh::St_mesh, gmsh_filename::String)

    #
    # Read GMSH grid from file
    #
    model    = GmshDiscreteModel(gmsh_filename, renumber=true)
    topology = get_grid_topology(model)
    mesh.nsd = num_cell_dims(model)
    
    #d_to_num_dfaces = [num_vertices(model), num_edges(model), num_cells(model)]
    #@show labels = FaceLabeling(d_to_num_dfaces)
    #@show get_face_entity(labels,0) .= get_isboundary_face(model,0) .+ 1
    #@show get_face_entity(labels,1) .= get_isboundary_face(model,1) .+ 1
    #@show get_face_entity(labels,2) .= get_isboundary_face(model,2) .+ 1
    
    if mesh.nsd == 3
        mesh.NNODES_EL  = 8
        mesh.NEDGES_EL  = 12
        mesh.NFACES_EL  = 6
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 4
    elseif mesh.nsd == 2
        mesh.NNODES_EL  = 4
        mesh.NEDGES_EL  = 4
        mesh.NFACES_EL  = 1
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 4
    elseif mesh.nsd == 1
        mesh.NNODES_EL  = 2
        mesh.EL_NODES   = 2
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
    #=
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
    
    ngl                     = mesh.nop + 1
    tot_linear_poin         = mesh.npoin_linear
    
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    tot_faces_internal_nodes = mesh.nfaces*(ngl-2)^(mesh.nsd-1)
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)^(mesh.nsd)
    
    el_edges_internal_nodes = mesh.NEDGES_EL*(ngl-2)
    el_faces_internal_nodes = mesh.NFACES_EL*(ngl-2)^(mesh.nsd-1)
    el_vol_internal_nodes   = (ngl-2)^(mesh.nsd)
    
    #Update number of grid points from linear count to total high-order points
    mesh.npoin = tot_linear_poin + tot_edges_internal_nodes + tot_faces_internal_nodes + tot_vol_internal_nodes
    
    if (mesh.nop > 1)
        println(" # GMSH HIGH-ORDER GRID PROPERTIES")
        println(" # N. edges internal points   : ", tot_edges_internal_nodes)
        println(" # N. faces internal points   : ", tot_faces_internal_nodes)
        println(" # N. volumes internal points : ", tot_vol_internal_nodes)
        println(" # N. total high order points : ", mesh.npoin)
        println(" # GMSH HIGH-ORDER GRID PROPERTIES ...................... END")
    end
    
    
    #
    # Resize (using resize! from ElasticArrays) as needed
    #
    resize!(mesh.x, (mesh.npoin_linear))
    resize!(mesh.y, (mesh.npoin_linear))
    resize!(mesh.z, (mesh.npoin_linear))

    resize!(mesh.conn_edge_L2G, (1, mesh.NEDGES_EL, mesh.nelem))
    resize!(mesh.conn_face_L2G, (1, mesh.NFACES_EL, mesh.nelem))
    
    resize!(mesh.conn_edge_el,  (2, mesh.NEDGES_EL, mesh.nelem))
    resize!(mesh.conn_face_el,  (4, mesh.NFACES_EL, mesh.nelem))

    mesh.npoin_el = mesh.NNODES_EL + el_edges_internal_nodes + el_faces_internal_nodes + el_vol_internal_nodes
    resize!(mesh.conn, (mesh.npoin_el*mesh.nelem))
    resize!(mesh.conn_ptr, (mesh.nelem))
    =#
    #
    # Connectivity matrices
    #
    mesh.cell_node_ids     = model.grid.cell_node_ids
    mesh.conn_unique_faces = get_face_nodes(model, FACE) #faces --> 4 nodes
    mesh.conn_unique_edges = get_face_nodes(model, EDGE) #edges --> 2 nodes
    
   #= mesh.conn = reshape(mesh.conn, mesh.npoin_el, mesh.nelem)
    if (mesh.nsd == 1)
        #mesh.conn = reshape(mesh.conn, mesh.ngl, mesh.nelem)
    elseif (mesh.nsd == 2)
        #mesh.conn = reshape(mesh.conn, mesh.ngl, mesh.ngl - 4,  mesh.nelem)
    elseif (mesh.nsd == 3)
        #mesh.conn = reshape(mesh.conn, mesh.ngl, mesh.ngl, mesh.ngl, mesh.nelem)
        for iel = 1:mesh.nelem
            
            mesh.conn[1, iel] = mesh.cell_node_ids[iel][2]
            mesh.conn[2, iel] = mesh.cell_node_ids[iel][6]
            mesh.conn[3, iel] = mesh.cell_node_ids[iel][8]
            mesh.conn[4, iel] = mesh.cell_node_ids[iel][4]
            mesh.conn[5, iel] = mesh.cell_node_ids[iel][1]
            mesh.conn[6, iel] = mesh.cell_node_ids[iel][5]
            mesh.conn[7, iel] = mesh.cell_node_ids[iel][7]
            mesh.conn[8, iel] = mesh.cell_node_ids[iel][3]
            
        end
    end=#


#
# Add high-order points to edges, faces, and elements (volumes)
#
# initialize LGL struct and buyild Gauss-Lobatto-xxx points
Legendre = St_Legendre{Float64}(0.0, 0.0, 0.0, 0.0)
lgl      = St_lgl{Float64}(zeros(mesh.nop+1),
                           zeros(mesh.nop+1))
build_lgl!(Legendre, lgl, mesh.nop)

#Edges
#populate_conn_edge_el!(mesh)
@time add_high_order_nodes_edges!(mesh, lgl)
end

#=
function populate_conn_edge_el!(mesh::St_mesh)
    @info " AA"
    for iel = 1:mesh.nelem

        #
        # CGNS numbering
        #
        ip1 = mesh.cell_node_ids[iel][2]
        ip2 = mesh.cell_node_ids[iel][6]
        ip3 = mesh.cell_node_ids[iel][8]
        ip4 = mesh.cell_node_ids[iel][4]
        ip5 = mesh.cell_node_ids[iel][1]
        ip6 = mesh.cell_node_ids[iel][5]
        ip7 = mesh.cell_node_ids[iel][7]
        ip8 = mesh.cell_node_ids[iel][3]

        # Edges bottom face:
        iedg_el = 1
        mesh.conn_edge_el[1, iedg_el, iel] = ip1
        mesh.conn_edge_el[2, iedg_el, iel] = ip2
        iedg_el = 2
        mesh.conn_edge_el[1, iedg_el, iel] = ip2
        mesh.conn_edge_el[2, iedg_el, iel] = ip3
        iedg_el = 3
        mesh.conn_edge_el[1, iedg_el, iel] = ip3
        mesh.conn_edge_el[2, iedg_el, iel] = ip4
        iedg_el = 4
        mesh.conn_edge_el[1, iedg_el, iel] = ip4
        mesh.conn_edge_el[2, iedg_el, iel] = ip1

        #Vertical edges
        iedg_el = 5
        mesh.conn_edge_el[1, iedg_el, iel] = ip1
        mesh.conn_edge_el[2, iedg_el, iel] = ip5
        iedg_el = 6
        mesh.conn_edge_el[1, iedg_el, iel] = ip2
        mesh.conn_edge_el[2, iedg_el, iel] = ip6
        iedg_el = 7
        mesh.conn_edge_el[1, iedg_el, iel] = ip3
        mesh.conn_edge_el[2, iedg_el, iel] = ip7
        iedg_el = 8
        mesh.conn_edge_el[1, iedg_el, iel] = ip4
        mesh.conn_edge_el[2, iedg_el, iel] = ip8

        #Edges top face
        iedg_el = 9
        mesh.conn_edge_el[1, iedg_el, iel] = ip5
        mesh.conn_edge_el[2, iedg_el, iel] = ip6
        iedg_el = 10
        mesh.conn_edge_el[1, iedg_el, iel] = ip6
        mesh.conn_edge_el[2, iedg_el, iel] = ip7
        iedg_el = 11
        mesh.conn_edge_el[1, iedg_el, iel] = ip7
        mesh.conn_edge_el[2, iedg_el, iel] = ip8
        iedg_el = 12
        mesh.conn_edge_el[1, iedg_el, iel] = ip8
        mesh.conn_edge_el[2, iedg_el, iel] = ip5

    end
    @info " BB"
    #@info "CONN EDGE " size(mesh.conn_edge_el)

end #populate_edge_el!
=#
#=
function populate_conn_face_el!(mesh::St_mesh)

for iel = 1:mesh.nelem

#
# CGNS numbering
#
ip1 = mesh.cell_node_ids[iel][2]
ip2 = mesh.cell_node_ids[iel][6]
ip3 = mesh.cell_node_ids[iel][8]
ip4 = mesh.cell_node_ids[iel][4]
ip5 = mesh.cell_node_ids[iel][1]
ip6 = mesh.cell_node_ids[iel][5]
ip7 = mesh.cell_node_ids[iel][7]
ip8 = mesh.cell_node_ids[iel][3]

#
# Local faces node connectivity:
# i.e. what nodes belong to a given local face in iel:
#
face_el = 1
mesh.conn_face_el[1, face_el, iel] = ip1
mesh.conn_face_el[2, face_el, iel] = ip4
mesh.conn_face_el[3, face_el, iel] = ip3
mesh.conn_face_el[4, face_el, iel] = ip2

face_el = 2
mesh.conn_face_el[1, face_el, iel] = ip1
mesh.conn_face_el[2, face_el, iel] = ip2
mesh.conn_face_el[3, face_el, iel] = ip6
mesh.conn_face_el[4, face_el, iel] = ip5

face_el = 3
mesh.conn_face_el[1, face_el, iel] = ip2
mesh.conn_face_el[2, face_el, iel] = ip3
mesh.conn_face_el[3, face_el, iel] = ip7
mesh.conn_face_el[4, face_el, iel] = ip6

face_el = 4
mesh.conn_face_el[1, face_el, iel] = ip3
mesh.conn_face_el[2, face_el, iel] = ip4
mesh.conn_face_el[3, face_el, iel] = ip8
mesh.conn_face_el[4, face_el, iel] = ip7

face_el = 5
mesh.conn_face_el[1, face_el, iel] = ip1
mesh.conn_face_el[2, face_el, iel] = ip5
mesh.conn_face_el[3, face_el, iel] = ip8
mesh.conn_face_el[4, face_el, iel] = ip4

face_el = 6
mesh.conn_face_el[1, face_el, iel] = ip5
mesh.conn_face_el[2, face_el, iel] = ip6
mesh.conn_face_el[3, face_el, iel] = ip7
mesh.conn_face_el[4, face_el, iel] = ip8

end

end #populate_face_el
=#
#=
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

face_in_elem[1, ifac, iel] = iel
face_in_elem[2, ifac, iel] = jel

face_in_elem[1, jfac, jel] = jel
face_in_elem[2, jfac, jel] = iel

iface = iface + 1;
end
end
end
end
end     
nfaces_int = Int64(iface)

end=#

function  add_high_order_nodes!(mesh::St_mesh) end

function  add_high_order_nodes_edges!(mesh::St_mesh, lgl::St_lgl)
    @info " CC"
    if (mesh.nop < 2) return end

    
    conn_edge_poin::Array{Int64, 2}  = zeros(mesh.nedges, mesh.ngl)
    @info typeof(mesh.conn_unique_edges) #This has type Table (As defined in the Gridap --- Verdugo paper)
    
    #
    # Populate mesh.conn_edge_L2G
    #
    cache_unique_edges = array_cache(mesh.conn_unique_edges) # allocation here
    for iedge_g = 1:mesh.nedges

        ai = getindex!(cache_unique_edges, mesh.conn_unique_edges, iedge_g)
        ai[2]
        
        for iel = 1:mesh.nelem
            for iedge_el = 1:mesh.NEDGES_EL
                
                #ai = getindex!(cache_unique_edges, mesh.conn_unique_edges, iedge_g)
                #ai[2]

                #ai = getindex!(cache_unique_edges, mesh.conn_unique_edges, iedge_g)
                #ai[2]

                ip1, ip2 = mesh.conn_unique_edges[iedge_g][1], mesh.conn_unique_edges[iedge_g][2]

                if(issetequal([ai[1], ai[2]], [ip1, ip2]))
                    #if(ai[2] == ai[1])
                    mesh.conn_edge_L2G[1, iedge_el, iel] = iedge_g
                end
                
            end
        end
    end
    #end
    @info " DD"
    
    return 
    
end
