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
