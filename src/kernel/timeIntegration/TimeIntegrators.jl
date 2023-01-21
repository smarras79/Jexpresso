using PrettyTables

include("../infrastructure/element_matrices.jl")

#
# Time discretization
#
abstract type AbstractTime end
struct RK <: AbstractTime end
struct RK3 <: AbstractTime end
struct RK5 <: AbstractTime end

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

function rk!(q::St_SolutionVectors;
             TD::RK5,
             SD::NSD_2D,
             QT,
             PT::AdvDiff,
             mesh::St_mesh,
             metrics::St_metrics,
             basis, ω,
             M, Δt,
             inputs::Dict,
             T)
    
    dq     = zeros(mesh.npoin)    
    RKcoef = buildRKIntegrator!(TD, T)
    
    for s = 1:length(RKcoef.a)
        
        #
        # rhs[ngl,ngl,nelem]
        #
        rhs_el      =      build_rhs(SD, QT, PT, q, basis.ψ, basis.dψ, ω,         mesh, metrics, T)
        rhs_diff_el = build_rhs_diff(SD, QT, PT, q, basis.ψ, basis.dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
        
        #
        # RHS[npoin] = DSS(rhs)
        #
        RHS = DSSijk_rhs(SD, rhs_el + rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, mesh.nop, T)
        divive_by_mass_matrix!(RHS, M, QT)
        
        for I=1:mesh.npoin
            dq[I] = RKcoef.a[s]*dq[I] + Δt*RHS[I]
            q.qn[I,1] = q.qn[I,1] + RKcoef.b[s]*dq[I]
        end
        
        #
        # B.C.
        #
        apply_boundary_conditions!(q, mesh, inputs, SD)
        
    end #stages

    #return qp
    
end
