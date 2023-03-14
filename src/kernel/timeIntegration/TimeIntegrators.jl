import Plots
using LinearAlgebra
using DiffEqBase
using OrdinaryDiffEq
using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler
import SciMLBase


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


function rhs_old!(du, u, params, t)

    T       = Float64
    TD      = params.TD
    SD      = params.SD
    QT      = params.QT
    PT      = params.PT
    BCT     = params.BCT
    neqns   = params.neqns
    basis   = params.basis
    mesh    = params.mesh
    metrics = params.metrics
    inputs  = params.inputs
    ω       = params.ω
    M       = params.M
    
    u[mesh.npoin_linear] = 0.0
    #
    # rhs[ngl,ngl,nelem]
    #
    rhs_el      = build_rhs(SD, QT, PT, neqns, u, basis.ψ, basis.dψ, ω, mesh, metrics, T)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, neqns, u, basis.ψ, basis.dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
    
    apply_boundary_conditions!(SD, rhs_el, u, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω,t, BCT, neqns)
    
    du = DSSijk_rhs(SD, rhs_el + inputs[:δvisc]*rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, mesh.nop, T)

    divive_by_mass_matrix!(du, M, QT)
    
    apply_periodicity!(SD, rhs_el, u, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, t, BCT, neqns)
    
    @info "" maximum(u) maximum(du)
    #error("sasa")
    return du
end


function rhs!(du, u, params, time)

    #SD::NSD_1D, QT::Inexact, PT::Wave1D, mesh::St_mesh, metrics::St_metrics, M, el_mat, u)
    T       = Float64
    TD      = params.TD
    SD      = params.SD
    QT      = params.QT
    PT      = params.PT
    BCT     = params.BCT
    neqns   = params.neqns
    basis   = params.basis
    mesh    = params.mesh
    metrics = params.metrics
    inputs  = params.inputs
    ω       = params.ω
    M       = params.M
    el_mat  = params.el_mat

    RHS = build_rhs(SD, QT, PT, neqns, u, basis.ψ, basis.dψ, ω, mesh, metrics, M, el_mat, time, T)
    
    du .= RHS

    
    return du #This is already DSSed
end



function time_loop!(TD,
                    SD,
                    QT,
                    PT,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp::St_SolutionVars,
                    M, el_mat,
                    Nt, Δt,
                    neqns, 
                    inputs::Dict,
                    BCT,
                    OUTPUT_DIR::String,
                    T)
    
    #
    # ODE
    #
    u = zeros(T, mesh.npoin);
    u .= qp.qn[:,1];
    params = (; T, TD, SD, QT, PT, BCT, neqns, basis, ω, mesh, metrics, inputs, M, el_mat)
    tspan = (inputs[:tinit], inputs[:tend])
    
    prob = ODEProblem(rhs!,
                      u,
                      tspan,
                      params);

    alg = RK4()     #RK4() WORK like a jewl
    #alg = SSPRK53()
    #alg = (SSPRK104(), #1
    #       SSPRK53(),  #2
    #       Tsit5(),    #3
    #       RK4())      #4
    println(" # Solving ODE with ................................" , string(alg), "\n")
    @info " " Δt inputs[:tinit] inputs[:tend] alg
    
    sol = solve(prob,
                alg,
                dt = Δt,
                saveat = range(T(0.), Nt*T(Δt), length=inputs[:ndiagnostics_outputs]),
                progress = true,
                progress_message = (dt, u, p, t) -> t)
    println(" # Solving ODE with    ................................", string(alg), " DONE\n")

    #Plot
    for iout = 1: inputs[:ndiagnostics_outputs]
        title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
        jcontour(SD, mesh.x, mesh.y, sol.u[iout], title, string(OUTPUT_DIR, "/it.", iout, ".png"))
    end
    
  #=  p1 = Plots.scatter()
    p1 = Plots.scatter( p1, mesh.x, sol.u[1], label = "1", markershape = :diamond)
    p1 = Plots.scatter!(p1, mesh.x, sol.u[2], label = "2", markershape = :circle)
    p1 = Plots.scatter!(p1, mesh.x, sol.u[3], label = "3", markershape = :square)
    p1 = Plots.scatter!(p1, mesh.x, sol.u[4], label = "4", markershape = :star)
    p1 = Plots.scatter!(p1, mesh.x, sol.u[end], label = "5", markershape = :diamond)
    p1 = Plots.scatter!(p1, title = " q₁")
    Plots.plot(p1)=#
    
    #Plot final solution
    #    title = @sprintf "Tracer: final solution at t=%6.4f" inputs[:tend]
    #    jcontour(SD, mesh.x, mesh.y, qp.qn[:,1], title, string(OUTPUT_DIR, "/it", "end", ".png"))
    
end

function time_loop_old!(TD,
                      SD,
                      QT,
                      PT,
                      mesh::St_mesh,
                      metrics::St_metrics,
                      basis, ω,
                      qp::St_SolutionVars,
                      M, el_mat, 
                      Nt, Δt,
                      neqns, 
                      inputs::Dict,
                      BCT,
                      OUTPUT_DIR::String,
                      T)

    it_interval    = inputs[:diagnostics_interval]
    it_diagnostics = 1

    #
    # RK
    #
    t = inputs[:tinit]
    for it = 1:Nt
        if (mod(it, it_interval) == 0 || it == Nt)
            @printf "   Solution at t = %.6f sec\n" t
            for ieq = 1:length(qp.qn[1,:])
                @printf "      min/max(q[%d]) = %.6f %.6f\n" ieq minimum(qp.qn[:,ieq]) maximum(qp.qn[:,ieq])
            end
            
            title = @sprintf "Tracer: final solution at t=%6.4f" t
            jcontour(SD, mesh.x, mesh.y, qp.qn[:,1], title, string(OUTPUT_DIR, "/it.", it_diagnostics, ".png"))
            it_diagnostics = it_diagnostics + 1
            
        end
        
        rk!(qp, TD, SD, QT, PT, mesh, metrics, basis, ω, M, el_mat, Δt, neqns, inputs, BCT, t, T)
  
        t += Δt
    end
    
    #Plot final solution
    title = @sprintf "Tracer: final solution at t=%6.4f" inputs[:tend]
    jcontour(SD, mesh.x, mesh.y, qp.qn[:,1], title, string(OUTPUT_DIR, "/it", "end", ".png"))
    
end

function rk!(q::St_SolutionVars,
             TD,
             SD,
             QT,
             PT,
             mesh::St_mesh,
             metrics::St_metrics,
             basis, ω,
             M, el_mat, Δt,
             neqns, 
             inputs::Dict,
             BCT,
             time,
             T)

    
    params = (;T, TD, SD, QT, PT, BCT, neqns, basis, mesh, metrics, inputs, ω, M)
    
    RHS    = zeros(mesh.npoin)
    dq     = zeros(mesh.npoin, neqns)
    RKcoef = buildRKIntegrator!(TD, T)

    #
    # RHS[npoin] = DSS(rhs)
    #
    for ieqn=1:neqns
        
        for s = 1:length(RKcoef.a)
            
            # RHS = rhs!(RHS, q.qn, params, time)
            RHS = rhs!(SD, QT, Wave1D(), mesh, metrics, M, el_mat, q.qn)
            
            for I=1:mesh.npoin
                dq[I, ieqn] = RKcoef.a[s]*dq[I, ieqn] + Δt*RHS[I]
                q.qn[I, ieqn] = q.qn[I, ieqn] + RKcoef.b[s]*dq[I, ieqn]
            end

        end #stages
    end
    
end
