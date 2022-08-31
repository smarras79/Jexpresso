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

    if (inputs[:problem] === "test")
        
        @inbounds for i = 1:2:length(q.qn)-1
            q.qn[i]  = 1.0     # odd i
            q.qn[i+1] = 2.0   # even i
        end
        if isodd(length(q.qn))
            q.qn[end] = 0
        end

    elseif (inputs[:problem] === "wave1d")
        
        @inbounds for i = 1:length(q.qn)
            q.qn[i]  = exp(-64.0*mesh.x[i]*mesh.x[i])
        end
        
    end

    
    return q
end
