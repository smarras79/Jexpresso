using PlotlyJS

include("./initialize.jl")
include("./rhs.jl")

include("../../kernel/abstractTypes.jl")
include("../../kernel/timeIntegration/TimeIntegrators.jl")
include("../../io/plotting/jeplots.jl")


function time_loop(TD::RK5,
                   SD::NSD_2D,
                   QT::Inexact,
                   PT::Adv2D,
                   mesh::St_mesh,
                   metrics::St_metrics,
                   basis, ω,
                   qp,
                   M,
                   Nt, Δt,
                   inputs::Dict,
                   BCT, 
                   T)
     
    for it = 1:Nt
        time = it*Δt 
        rk!(qp; TD, SD, QT, PT,
            mesh, metrics, basis, ω, M, Δt, inputs, BCT, time, T)
        
    end
    #title = string("solution at final step ", Nt)
    #jcontour(mesh.x, mesh.y, qp.qn[:,1], title)
    
end
