include("../AbstractProblems.jl")

include("../../kernel/globalStructs.jl")
include("../../kernel/mesh/mesh.jl")
include("../../io/plotting/jeplots.jl")

function initialize(ET::SW, mesh::St_mesh, inputs::Dict, TFloat)

    @info " Initialize fields for AdvDiff ........................ "
    
    ngl = mesh.nop + 1
    nsd = mesh.nsd
    q   = St_SolutionVars{TFloat,Int8}(zeros(mesh.npoin, nsd+1),
                                       zeros(mesh.npoin, nsd+1),
                                       zeros(1, nsd+1),
                                       zeros(1, nsd+1),
                                       zeros(1, nsd+1),
                                       zeros(1, nsd+1),
                                       zeros(ngl, ngl, mesh.nelem, nsd+1))
    
    #Cone properties:
    σ = 32.0
    (xc, yc) = (-0.5, 0.0)
    #(xc, yc) = (-0.0, 0.0)
    
    for iel_g = 1:mesh.nelem
        for i=1:ngl
            for j=1:ngl

                ip = mesh.connijk[i,j,iel_g]
                x  = mesh.x[ip]
                y  = mesh.y[ip]
                
                q.qn[ip,1] = exp(-σ*((x - xc)*(x - xc) + (y - yc)*(y - yc)))
                q.qn[ip,2] = +y
                q.qn[ip,3] = -x

                q.qnel[i,j,iel_g,1] = q.qn[ip,1]
                q.qnel[i,j,iel_g,2] = q.qn[ip,2]
                q.qnel[i,j,iel_g,3] = q.qn[ip,3]
            end
        end
    end
    
    #------------------------------------------
    # Plot initial condition:
    # Notice that I scatter the points to
    # avoid sorting the x and q which would be
    # becessary for a smooth curve plot.
    #------------------------------------------
    jcontour(mesh.x, mesh.y, q.qn[:,1], "Initial conditions: tracer")

    @info " Initialize fields for AdvDiff ........................ DONE"
    
    return q
end
