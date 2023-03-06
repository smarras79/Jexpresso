include("../AbstractProblems.jl")

include("../../kernel/globalStructs.jl")
include("../../kernel/mesh/mesh.jl")
include("../../io/plotting/jeplots.jl")

function initialize(ET::Burgers, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    @info " Initialize fields for Burgers ........................ "
        
    ngl = mesh.nop + 1
    nsd = mesh.nsd
    
    q = St_SolutionVars{TFloat}(zeros(mesh.npoin, 3),  # qn+1
                                zeros(mesh.npoin, 3),  # qn
                                zeros(mesh.npoin, 3),  # qn-1
                                zeros(1, 1),           # qn-2
                                zeros(1, 1),           # qn-3
                                zeros(mesh.npoin, 3),  # qe
                                zeros(1, 1, 1, 1),     # qnelⁿ
                                zeros(1, 1, 1, 1),     # Fⁿ⁺¹
                                zeros(1, 1, 1, 1),     # Gⁿ⁺¹
                                zeros(1, 1, 1, 1))     # Hⁿ⁺¹

    #Cone properties:
    for iel_g = 1:mesh.nelem
        for i=1:ngl
            for j=1:ngl
                
                ip = mesh.connijk[i,j,iel_g]
                x  = mesh.x[ip]
                y  = mesh.y[ip]
                
                q.qn[ip,1] = sin(π*x)
                q.qe[ip,1] = q.qn[ip,1]
                q.qe[ip,2] = 0.0 #constant
                q.qe[ip,3] = 0.0 #constant

                q.qnm1[ip,1] = q.qn[ip,1]
                
            end
        end
    end
        
    #------------------------------------------
    # Plot initial condition:
    # Notice that I scatter the points to
    # avoid sorting the x and q which would be
    # becessary for a smooth curve plot.
    #------------------------------------------
    title = string( "Tracer: initial condition")
    jcontour(mesh.x, mesh.y, q.qn[:,1], title, string(OUTPUT_DIR, "/INIT.png"))
    
    @info " Initialize fields for Burgers ........................ DONE"
    
    return q
end
