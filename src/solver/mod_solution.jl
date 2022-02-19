include("../mesh/mod_mesh.jl")

using LinearAlgebra

export St_solution
export mod_solution_initial_conditions!

mutable struct St_solution{TInt,TFloat}
        
    q    ::Array{TFloat} #q(∀ vars, x̅(npoin)) at n
    q1   ::Array{TFloat} #q(∀ vars, x̅(npoin)) at n-1
    qq2  ::Array{TFloat} #q(∀ vars, x̅(npoin)) at n-2
    
    Finv ::Array{TFloat} #Inviscid flux
    Fvisc::Array{TFloat} #Viscous flux
    
    J ::Array{TFloat} #Jacobian matrix

end #St_mesh
  
function mod_solution_initial_conditions!(mesh::St_mesh, qinit::St_solution, npx::Int64, npy::Int64, npz::Int64, problem="burgers1d")

    println(problem)

    npoin = npx*npy*npz
    if (problem == "burgers1d")
        
        @show typeof(problem)
        for i = 1:npoin
            x = mesh.x[i]            
            qinit.q[1, i] = sin(2π*x)
        end
    end
end

