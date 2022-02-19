#module Solution

using LinearAlgebra
export St_solution
export initial_conditions

mutable struct St_solution{TInt,TFloat}
    
    nvars::TInt
    npoin::TInt
    
    q    ::Matrix{TFloat} #q(∀ vars, x̅(npoin)) at n
    q1   ::Matrix{TFloat} #q(∀ vars, x̅(npoin)) at n-1
    qq2  ::Matrix{TFloat} #q(∀ vars, x̅(npoin)) at n-2
    
    Finv ::Matrix{TFloat} #Inviscid flux
    Fvisc::Matrix{TFloat} #Viscous flux
    
    J ::Matrix{TFloat} #Jacobian matrix

end #St_mesh
  
function initial_conditions!(mesh::St_mesh, qinit::St_solution, npx::Int64, npy::Int64, npz::Int64; case="burgers1d")

    if ( cmp(case, "burgers1d") )
        for i = 1:npx
            x = mesh.x[i]
            qinit.q[i, 1, i] = sin(2π*x)
        end
    end
end
#end
