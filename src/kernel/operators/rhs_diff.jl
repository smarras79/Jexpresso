include("../ArtificialViscosity/DynSGS.jl")
include("./rhs.jl")

#
# CompEuler
#
function build_rhs_diff_work_array!(ρel, uel, vel, Tel)
    fill!(ρel,  zero(Float64))
    fill!(uel,  zero(Float64))
    fill!(vel,  zero(Float64))
    fill!(Tel,  zero(Float64))
end

function uToPrimitives!(iel, ρel, uel, vel, Tel, u, nelem, ngl, npoin, conn, δenergy)

#    for iel=1:nelem
        for j=1:ngl, i=1:ngl
            
            m1 = conn[iel,i,j]
            m2 =   npoin + m1
            m3 = 2*npoin + m1
            m4 = 3*npoin + m1
            
            ρel[i,j] = u[m1]
            #ρel[i,j] = uaux[m,1]

            uel[i,j] = u[m2]/ρel[i,j]
            #uel[i,j] = uaux[m,2]/ρel[i,j]
            
            vel[i,j] = u[m3]/ρel[i,j]
            #vel[i,j] = uaux[m,3]/ρel[i,j]
            
            Tel[i,j] = u[m4]/ρel[i,j] - δenergy*0.5*(uel[i,j]^2 + vel[i,j]^2)
            #Tel[i,j] = uaux[m,4]/ρel[i,j] - δenergy*0.5*(uel[i,j]^2 + vel[i,j]^2)
            
        end
    end
#end

#function build_rhs_diff!(rhs_diff::SubArray{Float64}, SD::NSD_2D, QT, PT::CompEuler, qp, neqs, basis, ω, inputs, mesh::St_mesh, metrics::St_metrics, μ, T; qoutauxi=zeros(1,1))   
function viscous_rhs_el!(rhs_diff_el, uaux, u, ρel, uel, vel, Tel,
                         mesh, metrics, basis, ω, neqs, inputs, SD::NSD_2D)#; νx=0.0, νy=0.0)
    
    build_rhs_diff_work_array!(ρel, uel, vel, Tel)
        
    #
    # qp[1:npoin]         <-- qq[1:npoin, "ρ"]
    # qp[npoin+1:2npoin]  <-- qq[1:npoin, "ρu"]
    # qp[2npoin+1:3npoin] <-- qq[1:npoin, "ρE"]
    #
    #u2uaux!(uaux, u, neqs, mesh.npoin)
    
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    if (inputs[:case] === "rtb")
        δenergy = 0.0
    else
        δenergy = 1.0
    end

  #=  
    μ = inputs[:νx]
    ν = 0.0
    κ = μ
    γ = 1.4
    Pr = 0.1
    =#
    for iel=1:mesh.nelem
        uToPrimitives!(iel, ρel, uel, vel, Tel, u, mesh.nelem, mesh.ngl, mesh.npoin, mesh.connijk, δenergy)
    end
    #=
        #ν = Pr*μ[iel]/maximum(ρel[:,:])
        #κ = Pr*μ[iel]/(γ - 1.0)
        #ν = μ[iel]#10.0
        #κ = μ[iel]#10.0
    =#
    
    #=    for l = 1:params.mesh.ngl
            for k = 1:params.mesh.ngl
                
                ωJkl = params.ω[k]*params.ω[l]*params.metrics.Je[iel,k,l]
                
                dρdξ = 0.0
                dudξ = 0.0
                dvdξ = 0.0
                dTdξ = 0.0

                dρdη = 0.0
                dudη = 0.0
                dvdη = 0.0
                dTdη = 0.0
                for i = 1:params.mesh.ngl
                    dρdξ += params.basis.dψ[i,k]*params.uaux_el[iel,1,i,l]
                    dudξ += params.basis.dψ[i,k]*params.uaux_el[iel,2,i,l]
                    dvdξ += params.basis.dψ[i,k]*params.uaux_el[iel,3,i,l]
                    dTdξ += params.basis.dψ[i,k]*params.uaux_el[iel,4,i,l]

                    dρdη += params.basis.dψ[i,l]*params.uaux_el[iel,1,k,i]
                    dudη += params.basis.dψ[i,l]*params.uaux_el[iel,2,k,i]
                    dvdη += params.basis.dψ[i,l]*params.uaux_el[iel,3,k,i]
                    dTdη += params.basis.dψ[i,l]*params.uaux_el[iel,4,k,i]
                end
                                
                #
                dξdx_kl = params.metrics.dξdx[iel,k,l]
                dξdy_kl = params.metrics.dξdy[iel,k,l]
                dηdx_kl = params.metrics.dηdx[iel,k,l]
                dηdy_kl = params.metrics.dηdy[iel,k,l]
                
                #
                dρdx =  ν*(dρdξ*dξdx_kl + dρdη*dηdx_kl)
                dudx =  μ*(dudξ*dξdx_kl + dudη*dηdx_kl)
                dvdx =  μ*(dvdξ*dξdx_kl + dvdη*dηdx_kl)
                dTdx =  κ*(dTdξ*dξdx_kl + dTdη*dηdx_kl) #+μ∇u⋅u
                
                dρdy =  ν*(dρdξ*dξdy_kl + dρdη*dηdy_kl)
                dudy =  μ*(dudξ*dξdy_kl + dudη*dηdy_kl)
                dvdy =  μ*(dvdξ*dξdy_kl + dvdη*dηdy_kl)
                dTdy =  κ*(dTdξ*dξdy_kl + dTdη*dηdy_kl) #+μ∇u⋅u
                
                ∇ξ∇ρ_kl = dξdx_kl*dρdx + dξdy_kl*dρdy
                ∇η∇ρ_kl = dηdx_kl*dρdx + dηdy_kl*dρdy
                
                ∇ξ∇u_kl = dξdx_kl*dudx + dξdy_kl*dudy
                ∇η∇u_kl = dηdx_kl*dudx + dηdy_kl*dudy            
                ∇ξ∇v_kl = dξdx_kl*dvdx + dξdy_kl*dvdy
                ∇η∇v_kl = dηdx_kl*dvdx + dηdy_kl*dvdy

                ∇ξ∇T_kl = dξdx_kl*dTdx + dξdy_kl*dTdy
                ∇η∇T_kl = dηdx_kl*dTdx + dηdy_kl*dTdy
                
                for i = 1:params.mesh.ngl
                    
                    dhdξ_ik, dhdη_il = params.basis.dψ[i,k], params.basis.dψ[i,l]
                    
                    params.rhs_diffξ_el[iel,i,l,1] -= ωJkl*dhdξ_ik*∇ξ∇ρ_kl
                    params.rhs_diffη_el[iel,k,i,1] -= ωJkl*dhdη_il*∇η∇ρ_kl
                    
                    params.rhs_diffξ_el[iel,i,l,2] -= ωJkl*dhdξ_ik*∇ξ∇u_kl
                    params.rhs_diffη_el[iel,k,i,2] -= ωJkl*dhdη_il*∇η∇u_kl
                    
                    params.rhs_diffξ_el[iel,i,l,3] -= ωJkl*dhdξ_ik*∇ξ∇v_kl
                    params.rhs_diffη_el[iel,k,i,3] -= ωJkl*dhdη_il*∇η∇v_kl
                    
                    params.rhs_diffξ_el[iel,i,l,4] -= ωJkl*dhdξ_ik*∇ξ∇T_kl
                    params.rhs_diffη_el[iel,k,i,4] -= ωJkl*dhdη_il*∇η∇T_kl
                    
                end
            end
        end
    end

    params.rhs_diff_el .= @views (params.rhs_diffξ_el[:,:,:,:] + params.rhs_diffη_el[:,:,:,:])
    =#
end
