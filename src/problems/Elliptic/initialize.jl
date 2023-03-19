using WriteVTK

include("../AbstractProblems.jl")
include("../../kernel/globalStructs.jl")
include("../../kernel/mesh/mesh.jl")
include("../../io/plotting/jeplots.jl")

function initialize(SD::NSD_1D, ET::Elliptic, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    nothing
end


function initialize(SD::NSD_2D, ET::Elliptic, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    println(" # Initialize fields for ∇²(q) = f........................")
        
    ngl  = mesh.nop + 1
    nsd  = mesh.nsd
    neqs = 1    
    q    = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, neqs)
    
    test_case = "giraldo.12.14"
    if (test_case == "giraldo.12.14")

        c = 2.0
        xc, yc = 0.0, 0.0 #(maximum(mesh.x) + minimum(mesh.x))/2, (maximum(mesh.y) + minimum(mesh.y))/2
        
        for iel_g = 1:mesh.nelem
            for i=1:ngl
                for j=1:ngl

                    ip = mesh.connijk[i,j,iel_g];
                    #@info "iel: " iel_g, i,j, ip
                    x, y = mesh.x[ip], mesh.y[ip];
                    
                    q.qn[ip,1] = sinpi(c*(x - xc))*sinpi(c*(y - yc))
                    
                end
            end
        end
    end
    println(" # Initialize fields for ∇²(q) = f........................ DONE")
    
    println(" # Write q(init) to file ........................ ")
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:mesh.npoin]
    vtk_grid("./points", mesh.x, mesh.y, q.qn[:,1], cells) do vtk
        vtk["qinit", VTKPointData()] = q.qn[:,1]
    end
    println(" # Write q(init) to file ........................ DONE")
    
    return q
end
