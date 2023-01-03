include("../AbstractTypes.jl")

mutable struct RK_Integrator{TFloat}
  a::Array{TFloat}
  b::Array{TFloat}
  c::Array{TFloat}
end

function buildRK3Integrator!(RK::RK_Integrator)
  RK.a=[0.0 -5/9 -153/128]
  RK.b=[0.0 1/3 3/4]
  RK.c=[1/3 15/16 8/15]
end


function buildRK5Integrator!(RK::RK_Integrator)
    
    RK.a = [TFloat(0),
           TFloat(-567301805773)  / TFloat(1357537059087),
           TFloat(-2404267990393) / TFloat(2016746695238),
           TFloat(-3550918686646) / TFloat(2091501179385),
           TFloat(-1275806237668) / TFloat(842570457699 )]

    RK.b = [TFloat(1432997174477) / TFloat(9575080441755 ),
           TFloat(5161836677717) / TFloat(13612068292357),
           TFloat(1720146321549) / TFloat(2090206949498 ),
           TFloat(3134564353537) / TFloat(4481467310338 ),
           TFloat(2277821191437) / TFloat(14882151754819)]

    RK.c = [TFloat(0),
           TFloat(1432997174477) / TFloat(9575080441755),
           TFloat(2526269341429) / TFloat(6820363962896),
           TFloat(2006345519317) / TFloat(3224310063776),
           TFloat(2802321613138) / TFloat(2924317926251)]

end

function rk(RK::RK_Integrator, RHS::Array{T}) where > {T}

end

function time_advance(SD::NSD_2D, QT::Inexact(), mesh::St_mesh, metrics::St_metrics, T)

    #
    # ALGO 5.6 FROM GIRALDO: GLOBAL VERSION WITH SOLID-WALL B.C. AS A FIRST TEST
    #
    RK = RK_Integrator{Tc}(zeros(T,5),zeros(T,5),zeros(T,5))
    buildRK5Integrator!(RK)
    for it = 1:Nt
        
        dq = zeros(mesh.npoin);
        qe = zeros(mesh.ngl);
        for s = 1:length(RK.a)
            
            #
            # rhs[ngl,ngl,nelem]
            #
            rhs_el = build_rhs(SD, QT, Adv2D(), qp, basis.ψ, basis.dψ, ω, mesh, metrics)
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
        #@info size(RHS) size(mesh.x) size(mesh.y)
                
        clf()
        #title = string(" RHS for N=", Nξ, " & ", QT_String, " integration")        
        #frhs = PyPlot.tricontourf(mesh.x, mesh.y, qp.qn[:,1], levels=30)
        #PyPlot.colorbar(frhs)        
        #plt[:show]()
        display(PyPlot.tricontourf(mesh.x, mesh.y, qp.qn[:,1], levels=30))
        #plt[:show]()
    end
    
end
