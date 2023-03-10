using DifferentialEquations
using LinearAlgebra
#using DiffEqBase
#using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler
#import SciMLBase

include("../abstractTypes.jl")
include("../infrastructure/element_matrices.jl")
include("../../io/plotting/jeplots.jl")

mutable struct RK_Integrator{TFloat}
  a::Array{TFloat}
  b::Array{TFloat}
  c::Array{TFloat}
end

function buildRKIntegrator(RKtype::RK3, TFloat)

    RKcoef = RK_Integrator{TFloat}(zeros(TFloat,3),zeros(TFloat,3),zeros(TFloat,3))
    
    RKcoef.a=[0.0 -5/9 -153/128]
    RKcoef.b=[0.0 1/3 3/4]
    RKcoef.c=[1/3 15/16 8/15]
    
    return RKcoef
end


function buildRKIntegrator!(RKtype::RK5, TFloat)

    RKcoef = RK_Integrator{TFloat}(zeros(TFloat,5),zeros(TFloat,5),zeros(TFloat,5))
    
    RKcoef.a = [TFloat(0),
                TFloat(-567301805773)  / TFloat(1357537059087),
                TFloat(-2404267990393) / TFloat(2016746695238),
                TFloat(-3550918686646) / TFloat(2091501179385),
                TFloat(-1275806237668) / TFloat(842570457699 )]

    RKcoef.b = [TFloat(1432997174477) / TFloat(9575080441755 ),
                TFloat(5161836677717) / TFloat(13612068292357),
                TFloat(1720146321549) / TFloat(2090206949498 ),
                TFloat(3134564353537) / TFloat(4481467310338 ),
                TFloat(2277821191437) / TFloat(14882151754819)]

    RKcoef.c = [TFloat(0),
                TFloat(1432997174477) / TFloat(9575080441755),
                TFloat(2526269341429) / TFloat(6820363962896),
                TFloat(2006345519317) / TFloat(3224310063776),
                TFloat(2802321613138) / TFloat(2924317926251)]

    return RKcoef

end

function rhs!(RHSn, q, params, time)

    T   = Float64
    TD  = params.TD
    SD  = params.SD
    QT  = params.QT
    PT  = params.PT
    BCT = params.BCT
    neqns   = params.neqns
    basis   = params.basis
    mesh    = params.mesh
    metrics = params.metrics
    inputs  = params.inputs
    ω       = params.ω
    M       = params.M
    
    
    
    #
    # rhs[ngl,ngl,nelem]
    #
    rhs_el      = build_rhs(SD, QT, PT, neqns, q.qn, basis.ψ, basis.dψ, ω, mesh, metrics, T)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, neqns, q.qn, basis.ψ, basis.dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
    
    apply_boundary_conditions!(rhs_el, q.qn, mesh, inputs, SD,QT,metrics,basis.ψ,basis.dψ, ω,time,BCT,neqns)
    
    #
    # RHS[npoin] = DSS(rhs)
    #
    for ieqn=1:neqns
        RHSn = DSSijk_rhs(SD,
                          rhs_el[:,:,:,ieqn], # + inputs[:δvisc]*rhs_diff_el[:,:,:,ieqn],
                          mesh.connijk,
                          mesh.nelem, mesh.npoin, mesh.nop,
                          T)
        #RHSn = RHSn + inputs[:νx]*L*q.qn[:,ieqn]
        divive_by_mass_matrix!(RHSn, M, QT)
    end

    return RHSn
end

function time_loopnew!(TD,
                    SD,
                    QT,
                    PT,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp,
                    M, L, 
                    Nt, Δt,
                    neqns, 
                    inputs::Dict,
                    BCT,
                    OUTPUT_DIR::String,
                    T)
    it = 0
    t  = inputs[:tinit]
    t0 = t
    
    it_interval = inputs[:diagnostics_interval]
    it_diagnostics = 1
    
    params = (;TD, SD, QT, PT, BCT, neqns, mesh, metrics, basis, ω, inputs, M)
    tspan  = (FT(0), Nt*FT(Δt))

    RHSn     = zeros(T, mesh.npoin)
    RHSnm1   = copy(RHSn)
    RHSnm2   = copy(RHSnm1)
    qp.qnm1 .= qp.qn
    for it = 1:Nt
        if (mod(it, it_interval) == 0 || it == Nt)
            @printf "   Solution at t = %.6f sec\n" t
            for ieq = 1:neqns
                @printf "      min/max(q[%d]) = %.6f %.6f\n" ieq minimum(qp.qn[:,ieq]) maximum(qp.qn[:,ieq])
            end
        end
        
        rk!(qp, RHSn, RHSnm1, RHSnm2; TD, SD, QT, PT,
            mesh, metrics, basis, ω, M, L, Δt, neqns, inputs, BCT, time=t, T)

        
        title = @sprintf("Tracer: final solution at t=%6.4f", t)
        jcontour(SD, mesh.x, mesh.y, qp.qn[:,1], title, string(OUTPUT_DIR, "/it.", it_diagnostics, ".png"))
        it_diagnostics = it_diagnostics + 1
        t = t0 + Δt
        t0 = t
        
    end

    return
end

function time_loop!(TD,
                    SD,
                    QT,
                    PT,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp,
                    M, L, 
                    Nt, Δt,
                    neqns, 
                    inputs::Dict,
                    BCT,
                    OUTPUT_DIR::String,
                    T)
    it = 0
    t  = inputs[:tinit]
    t0 = t
    
    plot_at_times = [0.25, 0.5, 1.0, 1.5]
   
    it_interval = inputs[:diagnostics_interval]
    it_diagnostics = 1

    #
    # RK for first 2 steps
    #
    RHSn     = zeros(T, mesh.npoin)
    #RHSnm1   = copy(RHSn)
    #RHSnm2   = copy(RHSnm1)
    #qp.qnm1 .= qp.qn
    for it = 1:Nt
        if (mod(it, it_interval) == 0 || it == Nt)
            @printf "   Solution at t = %.6f sec\n" t
            for ieq = 1:neqns
                @printf "      min/max(q[%d]) = %.6f %.6f\n" ieq minimum(qp.qn[:,ieq]) maximum(qp.qn[:,ieq])
            end
            
            title = @sprintf "Tracer: final solution at t=%6.4f" t
            jcontour(SD, mesh.x, mesh.y, qp.qn[:,1], title, string(OUTPUT_DIR, "/it.", it_diagnostics, ".png"))
            it_diagnostics = it_diagnostics + 1
            
            t = t0 + Δt
            t0 = t
        end
            
        rk!(qp, RHSn, 0.0, 0.0; TD, SD, QT, PT, mesh, metrics, basis, ω, M, L, Δt, neqns, inputs, BCT, time=t, T)
  
    end
    
    #Plot final solution
    title = @sprintf "Tracer: final solution at t=%6.4f" inputs[:tend]
    jcontour(SD, mesh.x, mesh.y, qp.qn[:,1], title, string(OUTPUT_DIR, "/it", "end", ".png"))
    
end
 
       

function rknew!(q::St_SolutionVars,
             RHSn,
             RHSnm1,
             RHSnm2;
             TD,
             SD,
             QT,
             PT,
             mesh::St_mesh,
             metrics::St_metrics,
             basis, ω,
             M, L, Δt,
             neqns, 
             inputs::Dict,
             BCT,
             time,
             T)
    
   nothing
end


function rk!(q::St_SolutionVars,
             RHSn,
             RHSnm1,
             RHSnm2;
             TD,
             SD,
             QT,
             PT,
             mesh::St_mesh,
             metrics::St_metrics,
             basis, ω,
             M, L, Δt,
             neqns, 
             inputs::Dict,
             BCT,
             time,
             T)
    
    dq     = zeros(mesh.npoin, neqns)
    RKcoef = buildRKIntegrator!(TD, T)
    
    for s = 1:length(RKcoef.a)
      
        #
        # RHS[npoin] = DSS(rhs)
        #
        for ieqn=1:neqns
              
            #
            # rhs[ngl,ngl,nelem]
            #
            rhs_el      = build_rhs(SD, QT, PT, neqns, q.qn, basis.ψ, basis.dψ, ω, mesh, metrics, T)
            rhs_diff_el = build_rhs_diff(SD, QT, PT, neqns, q.qn, basis.ψ, basis.dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
            
            apply_boundary_conditions!(rhs_el, q.qn, mesh, inputs, SD,QT,metrics,basis.ψ,basis.dψ, ω,time,BCT,neqns)
            
            RHSn = DSSijk_rhs(SD,
                              rhs_el + inputs[:δvisc]*inputs[:νx]*rhs_diff_el,
                              mesh.connijk,
                              mesh.nelem, mesh.npoin, mesh.nop,
                              T)
            #RHSn = RHSn + inputs[:νx]*L*q.qn[:,ieqn]
            divive_by_mass_matrix!(RHSn, M, QT)
            
            for I=1:mesh.npoin
                dq[I, ieqn] = RKcoef.a[s]*dq[I, ieqn] + Δt*RHSn[I]
                q.qn[I, ieqn] = q.qn[I, ieqn] + RKcoef.b[s]*dq[I, ieqn]
            end

            apply_periodicity!(rhs_el, q.qn, mesh, inputs, SD,QT,metrics,basis.ψ,basis.dψ, ω,time,BCT,neqns)     
        end
    end #stages
    
end


function bdf2!(q::St_SolutionVars,
               RHSn,
               RHSnm1,
               RHSnm2;
               TD,
               SD,
               QT,
               PT,
               mesh::St_mesh,
               metrics::St_metrics,
               basis, ω,
               M, Δt,
               neqns, 
               inputs::Dict,
               BCT,
               time,
               T)
    
    dq     = zeros(mesh.npoin, neqns)
    RKcoef = buildRKIntegrator!(TD, T)
    
    #
    # rhs[ngl,ngl,nelem]
    #
    rhs_el      = build_rhs(SD, QT, PT, neqns, q.qn, basis.ψ, basis.dψ, ω, mesh, metrics, T)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, neqns, q.qn, basis.ψ, basis.dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
    
    apply_boundary_conditions!(rhs_el, q.qn, mesh, inputs, SD,QT,metrics,basis.ψ,basis.dψ, ω,time,BCT,neqns)

    #
    # RHS[npoin] = DSS(rhs) ∀ equation
    #
    for ieqn=1:neqns
        RHSn = DSSijk_rhs(SD,
                         rhs_el[:,:,:,ieqn] ,# + inputs[:δvisc]*rhs_diff_el[:,:,:,ieqn],
                         mesh.connijk,
                         mesh.nelem, mesh.npoin, mesh.nop,
                         T)
        divive_by_mass_matrix!(RHSn, M, QT)

        for ip = 1:mesh.npoin
            q.qn[ip,ieqn] = (18.0*q.qn[ip,ieqn] - 9.0*q.qnm1[ip,ieqn] + 2.0*q.qnm2[ip,ieqn])/11.0 + (6.0*Δt/11.0)*(3.0*RHSn[ip,ieqn] - 3.0*RHSnm1[ip,ieqn] + RHSnm2[ip,ieqn])
        end
    end
    apply_periodicity!(rhs_el,q.qn, mesh, inputs, SD,QT,metrics,basis.ψ,basis.dψ, ω,time,BCT,neqns)
    
end
