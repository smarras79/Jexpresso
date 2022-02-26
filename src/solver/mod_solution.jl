include("../mesh/mod_mesh.jl")

using DiffEqOperators, DifferentialEquations
using LinearAlgebra
using Revise
using SparseMatricesCSR

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
                                          problem="burgers1d")
    
    nsd = mesh.nsd
    npx,  npy, npz = mesh.npx,  mesh.npy,  mesh.npz
    xmin, xmax     = mesh.xmin, mesh.xmax
    ymin, ymax     = mesh.ymin, mesh.ymax
    zmin, zmax     = mesh.zmin, mesh.zmax
    
    if (problem == "burgers1d" || problem == "burgers")

        if (nsd == 1)
            qsol.q[1, 1] = 0
            qsol.q[1, npx] = 0
            Δx = (mesh.xmax - mesh.xmin)/(npx-1)
            for i::TInt = 2:npx-1
                x = mesh.x[i]
                qsol.q[1, i] = sin(2*π*x) + 0.5*sin(π*x)
            end
        end
    end

end


function mod_solution_solveODE!(mesh::St_mesh,
                                qsol::St_solution,
                                problem="burgers1d")


    q0::Typeof(qsol.q)
    
    npx = mesh.npx, npy = mesh.npy, npz = mesh.npz
    xmin, xmax = mesh.xmin, mesh.xmax
    ymin, ymax = mesh.ymin, mesh.ymax
    zmin, zmax = mesh.zmin, mesh.zmax
    
    # Define a problem
    p = (1.0,2.0,1.5,1.25) # a,b,c,d
    
    ϵ  = 1.0e-2
    Δx = (xmax - xmin)/(npx-1)
    a  = ϵ/(Δx*Δx)
    
    x = mesh.x
    qsol.q[1, 1]   = 0
    qsol.q[1, npx] = 0
    Δx = (xmax - xmin)/(npx-1)
    q0 = sin(2*π*x) + 0.5*sin(π*x)

   #= for i = 1:npx
        for j = 1:npx
        A[j, i] = 
    end
    
        u = qsol.q
    #    for i = 2:npx
    #        rhs[i] = (ϵ/Δx)*(u[i+1] - 2.0*u[i] + u[i-1]) #- (0.25/Δx)*(u^2[i+1] - u^2[i-1])
    #    end

    =#
    
    # Plot the solution using the plot recipe
    #plot(sol,title="All Plots.jl Attributes are Available")
    
  #=  f = function (du,u,p,t) # Define f as an in-place update into du

        a,b,c,d = p
        ϵ = 1.0e-2
        Δx = (mesh.xmax - mesh.xmin)/(npx-1)
        
        for i=2:npx-1
            RHS[i] = ϵ*()
            du = a*u[1] - b*u[1]*u[2]
        end
    end
    =#
end
