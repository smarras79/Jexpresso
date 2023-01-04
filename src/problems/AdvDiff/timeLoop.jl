using PlotlyJS

include("../../io/plotting/jeplots.jl")
include("../../kernel/abstractTypes.jl")
include("../../kernel/timeIntegration/TimeIntegrators.jl")
include("./rhs.jl")

function rk(TD::RK, rkcoeff::RK_Integrator, RHS::Array)
    nothing   
end

function time_loop(SD::NSD_2D,
                   QT::Inexact,
                   PT::Adv2D,
                   mesh::St_mesh,
                   metrics::St_metrics,
                   basis, ω,
                   qp,
                   M,
                   Nt, Δt,
                   T)
    
    RK = RK_Integrator{T}(zeros(T,5),zeros(T,5),zeros(T,5))
    buildRK5Integrator!(RK)

    
    #trace = PlotlyJS.contour(;z=qp.qn[:,1], x=mesh.x, y=mesh.y)
    #data = [trace]
    #pl = plot(data)
    p = jcontour(mesh.x, mesh.y, qp.qn[:,1], "Initial conditions: tracer")
    for it = 1:Nt
        
        dq = zeros(mesh.npoin)
        qe = zeros(mesh.ngl)
        for s = 1:length(RK.a)
            
            #
            # rhs[ngl,ngl,nelem]
            #
            rhs_el = build_rhs(SD, QT, PT, qp, basis.ψ, basis.dψ, ω, mesh, metrics)

            #
            # RHS[npoin] = DSS(rhs)
            #
            RHS = DSSijk_rhs(SD, QT, rhs_el, mesh.connijk, mesh.nelem, mesh.npoin, mesh.nop, T)
            RHS .= RHS./M
            for I=1:mesh.npoin
                dq[I] = RK.a[s]*dq[I] + Δt*RHS[I]
                qp.qn[I,1] = qp.qn[I,1] + RK.b[s]*dq[I]
            end
            
            #
            # B.C.: solid wall
            #
            #qp[1] = 0.0
            #qp[mesh.npoin_linear] = 0.0

        end #stages

        #plot
        #title = @printf("solution at tstep %d", it)
        #layout = Layout(;title=title)
        #trace = PlotlyJS.contour(;z=qp.qn[:,1], x=mesh.x, y=mesh.y)
        #data = [trace]
        #sleep(0.05)
        #display(react!(pl, data, layout))

    end
    title = string("solution at final step ", Nt)
    jcontour(mesh.x, mesh.y, qp.qn[:,1], title)
    
end
