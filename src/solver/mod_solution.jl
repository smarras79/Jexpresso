include("../mesh/mod_mesh.jl")

using Revise
using LinearAlgebra

#Plots
using Plots, DataFrames
plotlyjs()


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
  
function mod_solution_initial_conditions!(mesh::St_mesh,
                                          qsol::St_solution,
                                          nsd,
                                          npx, npy, npz, 
                                          problem="burgers1d")
    
    if (problem == "burgers1d" || problem == "burgers")

        if (nsd == 1)
            Δx = (mesh.xmax - mesh.xmin)/(npx-1)
            for i = 1:npx
                x = mesh.x[i]
                qsol.q[1, i] = sin(2π*x) + 0.5*sin(π*x)
            end
        end
    end

    
    plt = plot(mesh.x, qinit.q[1,:], w = 3)
    plot(scatter!(mesh.x, zeros(length(mesh.x)), x=:sepal_width, y=:sepal_length, mode="markers"))
    savefig("~/Work/Codes/jexpresso/figs/initial_conditions.png")
end


function mod_solution_solveODE!(mesh::St_mesh,
                                qsol::St_solution,
                                nsd,
                                npx, npy, npz, 
                                problem="burgers1d")

    
    
    
end

