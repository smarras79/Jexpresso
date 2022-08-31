include("../mesh/mod_mesh.jl")
include("./plotting/jeplots.jl")

mutable struct St_SolutionVectors{TFloat}

    qn::Array{TFloat}   #qⁿ
    qnm1::Array{TFloat} #qⁿ⁻¹
    qnm2::Array{TFloat} #qⁿ⁻²
    qnm3::Array{TFloat} #qⁿ⁻³
    
    qe::Array{TFloat}  #qexact
    
end

function mod_initialize_initialize(mesh::St_mesh, inputs::Dict, TFloat)

    q = St_SolutionVectors{TFloat}(zeros(mesh.npoin),
                                   zeros(mesh.npoin),
                                   zeros(mesh.npoin),
                                   zeros(mesh.npoin),
                                   zeros(mesh.npoin))
    
    if (inputs[:problem] === "wave1d")

        ngl = mesh.nop + 1
        for iel_g = 1:mesh.nelem
            for l=1:ngl
                ip = mesh.conn[l, iel_g]
                x = mesh.x[ip]
                
                q.qn[ip]  = exp(-64.0*x*x)
            end
        end

        #------------------------------------------
        # Plot initial condition:
        # Notice that I scatter the points to
        # avoid sorting the x and q which would be
        # becessary for a smooth curve plot.
        #------------------------------------------
        scatter_curve(mesh.x, q.qn, "")
        
    end
    
    return q
end
