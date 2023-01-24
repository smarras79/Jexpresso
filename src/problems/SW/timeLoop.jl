using PlotlyJS

include("../../kernel/abstractTypes.jl")
include("../../kernel/timeIntegration/TimeIntegrators.jl")
include("../../io/plotting/jeplots.jl")

function time_loop!(TD::RK5,
                    SD::NSD_2D,
                    QT,
                    PT::SW,
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

    it_interval = inputs[:diagnostics_interval]
    for it = 1:Nt

        if mod(it, it_interval) == 0 || it == Nt
            @printf "   Solution at t = %.6f sec\n" t
            @printf "      min(q) = %.6f\n" minimum(qp.qn[:,1])
            @printf "      max(q) = %.6f\n" maximum(qp.qn[:,1])
        end
        t = t0 + Δt
        t0 = t
        
        rk!(qp; TD, SD, QT, PT,
            mesh, metrics, basis, ω, M, Δt, inputs, T)
        
    end
    
end
