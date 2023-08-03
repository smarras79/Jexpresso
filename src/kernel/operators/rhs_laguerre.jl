include("rhs_laguerre_diff.jl")
function _build_rhs(SD::NSD_2D, QT::Inexact, PT, qp::Array, neqs, basis1, basis2, ω1, ω2
                    mesh::St_mesh, metrics1::St_metrics, metrics2::St_metrics, M, De, Le, time, inputs, Δt, deps, T; qnm1=zeros(Float64,1,1), qnm2=zeros(Float64,1,1), μ=zeros(Float64,1,1))
    
    F      = zeros(mesh.ngl,mesh.ngl, neqs)
    G      = zeros(mesh.ngl,mesh.ngl, neqs)
    S      = zeros(mesh.ngl,mesh.ngl, neqs)
    rhs_el = zeros(mesh.ngl,mesh.ngl, mesh.nelem, neqs)
    qq     = zeros(mesh.npoin,neqs)   
    
    #ωJe = zeros(mesh.ngl,mesh.ngl)
    
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] .= 0.0 .+ view(qp, idx+1:i*mesh.npoin)
    end
    ωJe = zeros(mesh.ngl,mesh.ngl)
    
    
    F      = zeros(T, mesh.ngl, mesh.ngl, mesh.nelem_semi_inf, neqs)
    G      = zeros(T, mesh.ngl, mesh.ngl, mesh.nelem_semi_inf, neqs)
    rhs_el = zeros(T, mesh.ngl, mesh.ngl, mesh.nelem_semi_inf, neqs)

   
    for iel=1:mesh.nelem_semi_inf

        for j=1:mesh.ngr, i=1:mesh.ngl
            ip = mesh.connijk_lag[i,j,iel]

            user_flux!(@view(F[i,j,1:neqs]), @view(G[i,j,1:neqs]), SD, @view(qq[ip,1:neqs]), mesh; neqs=neqs)
            if (inputs[:lsource] == true)
                user_source!(@view(S[i,j,1:neqs]), @view(qq[ip,1:neqs]), mesh.npoin; neqs=neqs)
            end
        end
        ωJe[:,:] .= @view(metrics.ωJe[:,:,iel])
        
        for ieq = 1:neqs
           
            for j=1:mesh.ngr, i=1:mesh.ngl
                
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0
                for k = 1:mesh.ngl
                    dFdξ += basis2.dψ[k,i]*F[k,j,ieq]
                    dFdη += basis2.dψ[k,j]*F[i,k,ieq]
                    
                    dGdξ += basis1.dψ[k,i]*G[k,j,ieq]
                    dGdη += basis1.dψ[k,j]*G[i,k,ieq]
                end

                dFdx = dFdξ*metrics2.dξdx[i,j,iel] + dFdη*metrics2.dηdx[i,j,iel]
                dGdy = dGdξ*metrics2.dξdy[i,j,iel] + dGdη*metrics2.dηdy[i,j,iel]
                rhs_el[i,j,iel,ieq] -= ωJe[i,j]*(dFdx + dGdy)  - ωJe[i,j]*S[i,j,ieq] #gravity
                
            end
        end
    end
    
    
    apply_boundary_conditions!(SD, rhs_el, qq, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, Δt*(floor(time/Δt)), neqs)
    RHS = DSS_rhs_laguerre(SD, rhs_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)

    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qp[idx+1:i*mesh.npoin] .= qq[:,i]
    end
    
    if (inputs[:lvisc] == true)
        
        if (lowercase(inputs[:visc_model]) === "dsgs")
            
            if (rem(time, Δt) == 0 && time > 0.0)
                qnm1 .= qnm2
                qnm2 .= qq
            end
            
            compute_viscosity!(μ, SD, PT, qq, qnm1, qnm2, RHS, Δt, mesh, metrics, T)
        else
            μ[:] .= inputs[:νx]
        end
        rhs_diff_el = build_rhs_diff(SD, QT, PT, qp, neqs, basis1, basis2, ω1, ω2, inputs, mesh, metrics1, metrics2, μ, T;)
        RHS .= RHS .+ DSS_rhs_laguerre(SD, rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqs, mesh.nop, T)
    end
    
    divive_by_mass_matrix!(RHS, M, QT,neqs)
    
    return RHS
end


#--------------------------------------------------------------------------------------------------------------------------------------------------
# CompEuler:
#--------------------------------------------------------------------------------------------------------------------------------------------------
#
# Optimized (more coud possibly be done)
#
function build_rhs(SD::NSD_2D, QT::Inexact, PT::CompEuler, qp::Array, neqs, basis1, basis2, ω1, ω2,
                   mesh::St_mesh, metrics1::St_metrics, metrics2::St_metrics, M, De, Le, time, inputs, Δt, deps, T; qnm1=zeros(Float64,1,1), qnm2=zeros(Float64,1,1), μ=zeros(Float64,1,1))
    
    RHS = _build_rhs(SD, QT, PT, qp, neqs, basis1, basis2, ω1, ω2, mesh, metrics1, metrics2, M, De, Le, time, inputs, Δt, deps, T; qnm1=qnm1, qnm2=qnm2, μ=μ)
    
    return RHS
    
end

#--------------------------------------------------------------------------------------------------------------------------------------------------
# AdvDiff
#--------------------------------------------------------------------------------------------------------------------------------------------------
function build_rhs(SD::NSD_2D, QT::Inexact, PT::AdvDiff, qp::Array, neqs, basis1, basis2, ω1, ω2, mesh::St_mesh, metrics1::St_metrics, metrics2::St_metrics, M, De, Le, time, inputs, Δt, deps, T;
                   qnm1=zeros(Float64,1,1), qnm2=zeros(Float64,1,1), μ=zeros(Float64,1,1))

    RHS = _build_rhs(SD, QT, PT, qp, neqs, basis1, basis2, ω1, ω2, mesh, metrics1, metrics2, M, De, Le, time, inputs, Δt, deps, T; qnm1=qnm1, qnm2=qnm2, μ=μ)
    
    return RHS
    
end

function build_rhs(SD::NSD_2D, QT::Exact, PT::AdvDiff, qp::Array, neqs, basis1, basis2, ω1, ω2, mesh::St_mesh, metrics1::St_metrics, metrics2::St_metrics, M, De, Le, time, inputs, Δt, deps, T;
                   qnm1=zeros(Float64,1,1), qnm2=zeros(Float64,1,1), μ=zeros(Float64,1,1))
    nothing
end

#--------------------------------------------------------------------------------------------------------------------------------------------------
# LinearCLaw
#--------------------------------------------------------------------------------------------------------------------------------------------------
function build_rhs(SD::NSD_2D, QT::Inexact, PT::LinearCLaw, qp::Array, neqs, basis, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, deps, T;
                   qnm1=zeros(Float64,1,1), qnm2=zeros(Float64,1,1), μ=zeros(Float64,1,1))    
    
    
    RHS = _build_rhs(SD, QT, PT, qp, neqs, basis, ω, mesh, metrics, M, De, Le, time, inputs, Δt, deps, T; qnm1=qnm1, qnm2=qnm2, μ=μ)
    
    return RHS
end

#--------------------------------------------------------------------------------------------------------------------------------------------------
# Source terms:
#--------------------------------------------------------------------------------------------------------------------------------------------------
function build_rhs_source(SD::NSD_2D,
                          QT::Inexact,
                          q::Array,
                          mesh::St_mesh,
                          M::AbstractArray, #M is sparse for exact integration
                          T)

    S = user_source(q, mesh, T)
    
    return M.*S    
end

function build_rhs_source(SD::NSD_2D,
                          QT::Exact,
                          q::Array,
                          mesh::St_mesh,
                          M::Matrix, #M is sparse for exact integration
                          T)

    S = user_source(q, mesh, T)
    
    return M*S   
end
