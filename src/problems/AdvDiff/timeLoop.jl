using PlotlyJS

include("../../kernel/abstractTypes.jl")
include("../../kernel/timeIntegration/TimeIntegrators.jl")
include("../../io/plotting/jeplots.jl")


function time_loop!(TD::RK5,
                    SD::NSD_2D,
                    QT,
                    PT::Adv2D,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp,
                    M,
                    Nt, Δt,
                    inputs::Dict, 
                    T)
    it = 0
    t  = inputs[:tinit]
    t0 = t

    #Diffusion coefficient
    for it = 1:Nt

        @printf "Solution at t=%.6f sec\n" t
        t = t0 + Δt
        t0 = t
        
        rk!(qp; TD, SD, QT, PT,
            mesh, metrics, basis, ω, M, Δt, inputs, T)
        
    end
    
end
