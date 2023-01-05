include("../AbstractProblems.jl")
include("../../kernel/mesh/mesh.jl")
include("../../io/plotting/jeplots.jl")

mutable struct St_SolutionVectors{TFloat}

    qnp1::Array{TFloat} #qⁿ⁺¹
    qn::Array{TFloat}   #qⁿ
    qnm1::Array{TFloat} #qⁿ⁻¹
    qnm2::Array{TFloat} #qⁿ⁻²
    qnm3::Array{TFloat} #qⁿ⁻³
    qe::Array{TFloat}   #qexact    
    qnel::Array{TFloat} #qⁿ[ngl,ngl,ngl,nelem]
    
    #Finv ::Array{TFloat} #Inviscid flux
    #Fvisc::Array{TFloat} #Viscous flux
    
end

function initialize(ET::Wave1D, mesh::St_mesh, inputs::Dict, TFloat)

    q = St_SolutionVectors{TFloat}(zeros(mesh.npoin),
                                   zeros(mesh.npoin),
                                   zeros(mesh.npoin),
                                   zeros(mesh.npoin),
                                   zeros(mesh.npoin),
                                   zeros(mesh.npoin),
                                   zeros(1, 1))

    ngl = mesh.nop + 1
    for iel_g = 1:mesh.nelem
        for l=1:ngl
            ip       = mesh.conn[l, iel_g]
            
            x        = mesh.x[ip]
            q.qn[ip] = exp(-64.0*x*x)
        end
    end
    #------------------------------------------
    # Plot initial condition:
    # Notice that I scatter the points to
    # avoid sorting the x and q which would be
    # becessary for a smooth curve plot.
    #------------------------------------------
    #plt = scatter() #Clear plot
    #display(scatter!(mesh.x, q.qn))
    
    return q
end


function initialize(ET::Adv2D, mesh::St_mesh, inputs::Dict, TFloat)

    @info " Initialize fields for Adv2D ........................ "
    
    ngl = mesh.nop + 1
    nsd = mesh.nsd
    q = St_SolutionVectors{TFloat}(zeros(mesh.npoin, nsd+1),
                                   zeros(mesh.npoin, nsd+1),
                                   zeros(1, nsd+1),
                                   zeros(1, nsd+1),
                                   zeros(1, nsd+1),
                                   zeros(1, nsd+1),
                                   zeros(ngl, ngl, mesh.nelem, nsd+1) )

    #Cone properties:
    σ = 32.0
    (xc, yc) = (-0.5, 0.0)
    
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

    @info " Initialize fields for Adv2D ........................ DONE"
    
    return q
end
