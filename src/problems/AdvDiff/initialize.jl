include("../AbstractProblems.jl")
include("../../kernel/globalStructs.jl")
include("../../kernel/mesh/mesh.jl")
include("../../io/plotting/jeplots.jl")

function initialize(ET::AdvDiff, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    @info " Initialize fields for AdvDiff ........................ "
        
    ngl  = mesh.nop + 1
    nsd  = mesh.nsd
    neqs = 1
    
    q = allocate_q(mesh.nelem, mesh.npoin, ngl, neqs)

    test_case = "kopriva.5.3.5"
    #test_case = "giraldo.15.8"
    if (test_case == "kopriva.5.3.5")
        #Cone properties:
        ν = inputs[:νx] 
        if ν == 0.0
            ν = 0.01
        end
        σ = 1.0/ν
        (xc, yc) = (0.5, 0.5)
        
        for iel_g = 1:mesh.nelem
            for i=1:ngl
                for j=1:ngl

                    ip = mesh.connijk[i,j,iel_g]
                    x  = mesh.x[ip]
                    y  = mesh.y[ip]

                    q.qn[ip,1] = exp(-σ*((x - xc)*(x - xc) + (y - yc)*(y - yc)))
                    q.qe[ip,1] = q.qn[ip,1]
                    q.qn[ip,2] = 0.8 #constant
                    q.qn[ip,3] = 0.8 #constant

                    q.qnm1[ip,1] = q.qn[ip,1]                    
                    #q.qnel[i,j,iel_g,1] = q.qn[ip,1]
                end
            end
        end
    elseif (test_case == "giraldo.15.8")
        σ = 32.0
        (xc, yc) = (-0.5, 0.0)
                
        for iel_g = 1:mesh.nelem
            for i=1:ngl
                for j=1:ngl

                    ip = mesh.connijk[i,j,iel_g]
                    x  = mesh.x[ip]
                    y  = mesh.y[ip]

                    q.qn[ip,1] = exp(-σ*((x - xc)*(x - xc) + (y - yc)*(y - yc)))
                    q.qe[ip,1] = q.qn[ip,1]             
                    q.qu[ip,2] = +y #constant
                    q.qu[ip,3] = -x #constant

                    #q.qnel[i,j,iel_g,1] = q.qn[ip,1]
                end
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
    
    @info " Initialize fields for AdvDiff ........................ DONE"
    
    return q
end
