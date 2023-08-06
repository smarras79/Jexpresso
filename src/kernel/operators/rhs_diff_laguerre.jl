#---------------------------------------------------------------------------
# Fetch equations name to access the user_rhs functions
#---------------------------------------------------------------------------
if (length(ARGS) === 1) #equations
    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/user_source.jl"))
        user_source_dir = string(@__DIR__, "/../../equations/", ARGS[1], "/user_source.jl")
    else
        user_source_dir = "../../fallbacks/source.jl"
    end
elseif (length(ARGS) === 2)  #equations/equations_case_name
    user_flux_dir   = string("../../equations/", ARGS[1], "/", ARGS[2], "/user_flux.jl")
    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/", ARGS[2], "/user_source.jl"))
        user_source_dir = string(@__DIR__, "/../../equations/", ARGS[1], "/", ARGS[2], "/user_source.jl")
    else
        @info " user_source.jl not defined. The fallback ../../fallbacks/source.jl will be used."
        user_source_dir =  "../../fallbacks/source.jl"
    end
end
include(user_flux_dir)
include(user_source_dir)
include("../ArtificialViscosity/DynSGS.jl")
#---------------------------------------------------------------------------
#
# AdvDiff
#

function build_rhs_diff(SD::NSD_2D, QT::Inexact, PT::AdvDiff, qp::Array, nvars, basis1, basis2, ω1, ω2, inputs,  mesh::St_mesh, metrics1::St_metrics, metrics2::St_metrics, μ, T;)

    N = mesh.ngl - 1
    
    qnel = zeros(mesh.ngl,mesh.ngr,mesh.nelem_semi_inf)
    
    rhsdiffξ_el = zeros(mesh.ngl,mesh.ngr,mesh.nelem_semi_inf)
    rhsdiffη_el = zeros(mesh.ngl,mesh.ngr,mesh.nelem_semi_inf)
    
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem_semi_inf

        for j=1:mesh.ngr, i=1:mesh.ngl
            m = mesh.connijk[i,j,iel]            
            qnel[i,j,iel] = qp[m,1]
        end
        
        for k = 1:mesh.ngl, l = 1:mesh.ngr
            ωJkl = ω[k]*ω[l]*metrics2.Je[k, l, iel]
            
            dqdξ = 0.0
            dqdη = 0.0
            for i = 1:mesh.ngl
                dqdξ = dqdξ + basis1.dψ[i,k]*qnel[i,l,iel]
                dqdη = dqdη + basis2.dψ[i,l]*qnel[k,i,iel]
            end
            
            dqdx = dqdξ*metrics2.dξdx[k,l,iel] + dqdη*metrics2.dηdx[k,l,iel]
            dqdy = dqdξ*metrics2.dξdy[k,l,iel] + dqdη*metrics2.dηdy[k,l,iel]
            
            ∇ξ∇q_kl = metrics2.dξdx[k,l,iel]*dqdx + metrics2.dξdy[k,l,iel]*dqdy
            ∇η∇q_kl = metrics2.dηdx[k,l,iel]*dqdx + metrics2.dηdy[k,l,iel]*dqdy
            
            for i = 1:mesh.ngl
                hll,     hkk     =  basis2.ψ[l,l],  basis1.ψ[k,k]
                dhdξ_ik = basis1.dψ[i,k]

                rhsdiffξ_el[i,l,iel] -= ωJkl*dhdξ_ik*hll*∇ξ∇q_kl
            end
            for i = 1:mesh.ngr
                hll,     hkk     =  basis2.ψ[l,l],  basis1.ψ[k,k]
                dhdη_il = basis2.dψ[i,l]

                rhsdiffη_el[k,i,iel] -= ωJkl*hkk*dhdη_il*∇η∇q_kl
            end   
        end
    end

    return (rhsdiffξ_el*inputs[:νx] + rhsdiffη_el*inputs[:νy])
    
end

#
# LinearCLaw
#
function build_rhs_diff(SD::NSD_2D, QT, PT::LinearCLaw, qp, neqs, basis1, basis2, ω1, ω2, inputs,  mesh::St_mesh, metrics1::St_metrics, metrics2::St_metrics, μ, T;)
    
    N = mesh.ngl - 1

    qnel = zeros(mesh.ngl,mesh.ngr,mesh.nelem_semi_inf, neqs)

    rhsdiffξ_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem_semi_inf, neqs)
    rhsdiffη_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem_semi_inf, neqs)
    qq = zeros(mesh.npoin,neqs)

    #
    # qp[1:npoin]         <-- qq[1:npoin, "p"]
    # qp[npoin+1:2npoin]  <-- qq[1:npoin, "u"]
    # qp[2npoin+1:3npoin] <-- qq[1:npoin, "v"]
    #
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] = qp[idx+1:i*mesh.npoin]
    end
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem_semi_inf

        for j=1:mesh.ngr, i=1:mesh.ngl
            m = mesh.connijk_lag[i,j,iel]
            qnel[i,j,iel,1:neqs] = qq[m,1:neqs]
        end

        for k = 1:mesh.ngl, l = 1:mesh.ngr
            ωJkl = ω[k]*ω[l]*metrics2.Je[k, l, iel]

            for ieq = 1:neqs
                dqdξ = 0.0
                dqdη = 0.0
                for i = 1:mesh.ngl
                    dqdξ = dqdξ + basis1.dψ[i,k]*qnel[i,l,iel,ieq]
                    dqdη = dqdη + basis2.dψ[i,l]*qnel[k,i,iel,ieq]
                end
                dqdx = dqdξ*metrics2.dξdx[k,l,iel] + dqdη*metrics2.dηdx[k,l,iel]
                dqdy = dqdξ*metrics2.dξdy[k,l,iel] + dqdη*metrics2.dηdy[k,l,iel]

                ∇ξ∇q_kl = metrics2.dξdx[k,l,iel]*dqdx + metrics2.dξdy[k,l,iel]*dqdy
                ∇η∇q_kl = metrics2.dηdx[k,l,iel]*dqdx + metrics2.dηdy[k,l,iel]*dqdy
                
                for i = 1:mesh.ngl
                  hll,     hkk     =  basis2.ψ[l,l],  basis1.ψ[k,k]
                  dhdξ_ik = basis1.dψ[i,k]

                  rhsdiffξ_el[i,l,iel,ieq] -= ωJkl*dhdξ_ik*hll*∇ξ∇q_kl
                end
                for i = 1:mesh.ngr
                  hll,     hkk     =  basis2.ψ[l,l],  basis1.ψ[k,k]
                  dhdη_il = basis2.dψ[i,l]

                  rhsdiffη_el[k,i,iel,ieq] -= ωJkl*hkk*dhdη_il*∇η∇q_kl
                end  
            end
        end
     end

    return (rhsdiffξ_el*inputs[:νx] + rhsdiffη_el*inputs[:νy])

end

#
# CompEuler

function build_rhs_diff(SD::NSD_2D, QT, PT::CompEuler, qp, neqs, basis1, basis2, ω1, ω2,inputs, mesh::St_mesh, metrics1::St_metrics, metrics2::St_metrics, μ, T; qoutauxi=zeros(1,1))
    
    ρel = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf)
    uel = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf)
    vel = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf)
    Tel = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf)

    rhsdiffξ_el = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf, neqs)
    rhsdiffη_el = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf, neqs)
    
    qq = zeros(mesh.npoin, neqs)
    
    #
    # qp[1:npoin]         <-- qq[1:npoin, "ρ"]
    # qp[npoin+1:2npoin]  <-- qq[1:npoin, "ρu"]
    # qp[2npoin+1:3npoin] <-- qq[1:npoin, "ρE"]
    #
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        qq[:,i] = qp[idx+1:i*mesh.npoin]
    end
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    if (inputs[:case] === "rtb")
        δenergy = 0.0
    else
        δenergy = 1.0
    end

    
    γ = 1.4
    Pr = 0.1
    for iel=1:mesh.nelem_semi_inf
        
        for j=1:mesh.ngr, i=1:mesh.ngl
            m = mesh.connijk_lag[i,j,iel]
            
            ρel[i,j,iel] = qq[m,1]
            uel[i,j,iel] = qq[m,2]/ρel[i,j,iel]
            vel[i,j,iel] = qq[m,3]/ρel[i,j,iel]
            
            Tel[i,j,iel] = qq[m,4]/ρel[i,j,iel] - δenergy*0.5*(uel[i,j,iel]^2 + vel[i,j,iel]^2)
        end    
        #ν = Pr*μ[iel]/maximum(ρel[:,:,iel])
        #κ = Pr*μ[iel]/(γ - 1.0)
        #ν = μ[iel]#10.0
        #κ = μ[iel]#10.0
        
        ν = 0.0
        μ[iel] = inputs[:νx]
        κ = μ[iel]
        for l = 1:mesh.ngr, k = 1:mesh.ngl
            ωJkl = ω[k]*ω[l]*metrics2.Je[k, l, iel]
            
            #for ieq = 1:neqs
            #dqdξ = 0.0
            dρdξ = 0.0
            dudξ = 0.0
            dvdξ = 0.0
            dTdξ = 0.0

            dρdη = 0.0
            dudη = 0.0
            dvdη = 0.0
            dTdη = 0.0
            for i = 1:mesh.ngl
                dρdξ += basis1.dψ[i,k]*ρel[i,l,iel]
                dudξ += basis1.dψ[i,k]*uel[i,l,iel]
                dvdξ += basis1.dψ[i,k]*vel[i,l,iel]
                dTdξ += basis1.dψ[i,k]*Tel[i,l,iel]

                dρdη += basis2.dψ[i,l]*ρel[k,i,iel]
                dudη += basis2.dψ[i,l]*uel[k,i,iel]
                dvdη += basis2.dψ[i,l]*vel[k,i,iel]
                dTdη += basis2.dψ[i,l]*Tel[k,i,iel]
            end
            
            dρdx =       ν*(dρdξ*metrics2.dξdx[k,l,iel] + dρdη*metrics2.dηdx[k,l,iel])
            dudx =  μ[iel]*(dudξ*metrics2.dξdx[k,l,iel] + dudη*metrics2.dηdx[k,l,iel])
            dvdx =  μ[iel]*(dvdξ*metrics2.dξdx[k,l,iel] + dvdη*metrics2.dηdy[k,l,iel])
            dTdx =       κ*(dTdξ*metrics2.dξdx[k,l,iel] + dTdη*metrics2.dηdx[k,l,iel]) #+μ∇u⋅u
          
            dρdy =       ν*(dρdξ*metrics2.dξdy[k,l,iel] + dρdη*metrics2.dηdy[k,l,iel])
            dudy =  μ[iel]*(dudξ*metrics2.dξdy[k,l,iel] + dudη*metrics2.dηdy[k,l,iel])
            dvdy =  μ[iel]*(dvdξ*metrics2.dξdy[k,l,iel] + dvdη*metrics2.dηdy[k,l,iel])
            dTdy =       κ*(dTdξ*metrics2.dξdy[k,l,iel] + dTdη*metrics2.dηdy[k,l,iel]) #+μ∇u⋅u
            
            ∇ξ∇ρ_kl = metrics2.dξdx[k,l,iel]*dρdx + metrics2.dξdy[k,l,iel]*dρdy
            ∇η∇ρ_kl = metrics2.dηdx[k,l,iel]*dρdx + metrics2.dηdy[k,l,iel]*dρdy
            
            ∇ξ∇u_kl = metrics2.dξdx[k,l,iel]*dudx + metrics2.dξdy[k,l,iel]*dudy
            ∇η∇u_kl = metrics2.dηdx[k,l,iel]*dudx + metrics2.dηdy[k,l,iel]*dudy            

            ∇ξ∇v_kl = metrics2.dξdx[k,l,iel]*dvdx + metrics2.dξdy[k,l,iel]*dvdy
            ∇η∇v_kl = metrics2.dηdx[k,l,iel]*dvdx + metrics2.dηdy[k,l,iel]*dvdy

            ∇ξ∇T_kl = metrics2.dξdx[k,l,iel]*dTdx + metrics2.dξdy[k,l,iel]*dTdy
            ∇η∇T_kl = metrics2.dηdx[k,l,iel]*dTdx + metrics2.dηdy[k,l,iel]*dTdy
            
            for i = 1:mesh.ngl
                
                dhdξ_ik, dhdη_il = basis1.dψ[i,k], basis2.dψ[i,l]
                
                rhsdiffξ_el[i,l,iel,1] -= ωJkl*dhdξ_ik*∇ξ∇ρ_kl
                
                rhsdiffξ_el[i,l,iel,2] -= ωJkl*dhdξ_ik*∇ξ∇u_kl
                
                rhsdiffξ_el[i,l,iel,3] -= ωJkl*dhdξ_ik*∇ξ∇v_kl
                
                rhsdiffξ_el[i,l,iel,4] -= ωJkl*dhdξ_ik*∇ξ∇T_kl
                
            end

            for i = 1:mesh.ngr

                dhdξ_ik, dhdη_il = basis1.dψ[i,k], basis2.dψ[i,l]

                rhsdiffη_el[k,i,iel,1] -= ωJkl*dhdη_il*∇η∇ρ_kl

                rhsdiffη_el[k,i,iel,2] -= ωJkl*dhdη_il*∇η∇u_kl

                rhsdiffη_el[k,i,iel,3] -= ωJkl*dhdη_il*∇η∇v_kl

                rhsdiffη_el[k,i,iel,4] -= ωJkl*dhdη_il*∇η∇T_kl

            end
            # end
        end
    end
    
    return (rhsdiffξ_el + rhsdiffη_el)

end
