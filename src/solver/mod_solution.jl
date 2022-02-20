include("../mesh/mod_mesh.jl")

#using LinearAlgebra
#using PlotlyJS, DataFrames

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
    
    npoin = npx*npy*npz
    if (problem == "burgers1d")
        
        for i = 1:npoin
            x = mesh.x[i]
            qinit.q[1, i] = sin(2π*x) + 0.5*sin(π*x)
        end
    end

    #plot(scatter(mesh.x, qinit.q[1,:], mode="markers"))
    (display(plot(cos, 0, 2π, mode="lines", Layout(title="cos(t)"))))
end

