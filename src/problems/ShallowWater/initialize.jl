include("../AbstractProblems.jl")

include("../../kernel/globalStructs.jl")
include("../../kernel/mesh/mesh.jl")
include("../../io/plotting/jeplots.jl")

function initialize(SD::NSD_2D, PT::ShallowWater, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
    
    """
    @info " Initialize fields for system of Linear Conservation Laws ........................ "
    
    ngl   = mesh.ngl
    nsd   = mesh.nsd
    nelem = mesh.nelem
    npoin = mesh.npoin
    neqs  = 3

    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, neqs,TFloat)
    
    #Cone properties:
       
    @info "Constant height and no flow shallow water" 
    for iel_g = 1:mesh.nelem
        for i=1:ngl
            for j=1:ngl

                ip = mesh.connijk[i,j,iel_g]
                x  = mesh.x[ip]
                y  = mesh.y[ip]
                q.qn[ip,1] = 0.75 + exp(-4*(abs(x)^2 + abs(y)^2))/4                                    #H
                q.qn[ip,2] = 0.001 * q.qn[ip,1]                                    #Hu
                q.qn[ip,3] = 0.001 * q.qn[ip,1]                                   #Hv 
                
            end
        end
    end
    
    @info "Initialize fields for system of Shallow Water equations ........................ DONE"
    
    return q
end
