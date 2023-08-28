#---------------------------------------------------------------------------
# Fetch equations name to access the user_rhs functions
#---------------------------------------------------------------------------
if (length(ARGS) === 1) #equations
    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/user_source.jl"))
        user_source_dir = string(@__DIR__, "../../equations/", ARGS[1], "/user_source.jl")
    else
        user_source_dir = "../../fallbacks/source.jl"
    end

    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/user_flux.jl"))
        user_flux_dir = string(@__DIR__, "../../equations/", ARGS[1], "/user_flux.jl")
    else
        user_flux_dir = string("../../fallbacks/flux.jl")
    end
    
    include(user_flux_dir)
    include(user_source_dir)
    
elseif (length(ARGS) === 2)  #equations/equations_case_name   
    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/", ARGS[2], "/user_source.jl"))
        user_source_dir = string(@__DIR__, "/../../equations/", ARGS[1], "/", ARGS[2], "/user_source.jl")
    else
        @info " user_source.jl not defined. The fallback ../../fallbacks/source.jl will be used."
        user_source_dir = "../../fallbacks/source.jl"
    end

    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/", ARGS[2], "/user_flux.jl"))
        user_flux_dir = string(@__DIR__, "/../../equations/", ARGS[1], "/", ARGS[2], "/user_flux.jl")
    else
        @info " user_flux.jl not defined. The fallback ../../fallbacks/flux.jl will be used."
        user_flux_dir = string("../../fallbacks/flux.jl")
    end

    include(user_flux_dir)
    include(user_source_dir)
end

#---------------------------------------------------------------------------
function RHStoDU!(du, RHS, neqs, npoin)
    for i=1:neqs
        idx = (i-1)*npoin
        du[idx+1:i*npoin] = @view RHS[:,i]
    end  
end

function u2uaux!(uaux, u, neqs, npoin)

     for i=1:neqs
        idx = (i-1)*npoin
        uaux[:,i] = view(u, idx+1:i*npoin)
    end
    
end

function rhs!(du, u, params, time)
    
    build_rhs!(@view(params.RHS[:,:]), u, params, time)
    
    RHStoDU!(du, @view(params.RHS[:,:]), params.neqs, params.mesh.npoin)
end

fun_ωJac(x, y, z) = x*y*z

function inviscid_rhs_el!(rhs_el, uaux, u, F, G, S, mesh, metrics, basis, ω, SD::NSD_2D; neqs=1, lsource=false)

    u2uaux!(uaux, u, neqs, mesh.npoin)
    
    #lsource = inputs[:lsource]
    for iel=1:mesh.nelem

        for j=1:mesh.ngl, i=1:mesh.ngl
            ip = mesh.connijk[iel,i,j]
            user_flux!(@view(F[i,j,:]), @view(G[i,j,:]), SD, @view(uaux[ip,:]), mesh; neqs=neqs)
            if lsource
                user_source!(@view(S[i,j,:]), @view(uaux[ip,:]), mesh.npoin; neqs=neqs)
            end
        end
        
        for j=1:mesh.ngl
            for i=1:mesh.ngl
                ωJac = fun_ωJac(ω[i], ω[j], metrics.Je[iel,i,j])
                dξdx_ij = metrics.dξdx[iel,i,j]
                dξdy_ij = metrics.dξdy[iel,i,j]
                dηdx_ij = metrics.dηdx[iel,i,j]
                dηdy_ij = metrics.dηdy[iel,i,j]

                for ieq = 1:neqs
                    
                    dFdξ = 0.0
                    dFdη = 0.0
                    dGdξ = 0.0
                    dGdη = 0.0
                    for k = 1:mesh.ngl
                        dFdξ += basis.dψ[k,i]*F[k,j,ieq]
                        dFdη += basis.dψ[k,j]*F[i,k,ieq]
                        
                        dGdξ += basis.dψ[k,i]*G[k,j,ieq]
                        dGdη += basis.dψ[k,j]*G[i,k,ieq]
                    end

                    dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
                    dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij
                    
                    rhs_el[iel,i,j,ieq] -= ωJac*((dFdx + dGdy)  - S[i,j,ieq])
                end
            end
        end
    end
end

function resetRHSToZero_inviscid!(params)
    fill!(params.rhs_el, zero(params.T))
    fill!(params.RHS,    zero(params.T))
end

function resetRHSToZero_viscous!(params)
    fill!(params.rhs_diff_el,  zero(params.T))
    fill!(params.rhs_diffξ_el, zero(params.T))
    fill!(params.rhs_diffη_el, zero(params.T))
    fill!(params.RHS_visc,     zero(params.T))
end


function uToPrimitives!(iel, ρel, uel, vel, Tel, u, nelem, ngl, npoin, conn)

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
            
            Tel[i,j] = u[m4]/ρel[i,j] #- δenergy*0.5*(uel[i,j]^2 + vel[i,j]^2)
            #Tel[i,j] = uaux[m,4]/ρel[i,j] - δenergy*0.5*(uel[i,j]^2 + vel[i,j]^2)
            
        end
    end
#end


function _build_rhs!(RHS, u, params, time)

    T       = Float64
    SD      = params.SD
    neqs    = params.neqs
    ngl     = params.mesh.ngl
    nelem   = params.mesh.nelem
    npoin   = params.mesh.npoin

    # rhs_el, RHS -> 0.0
    resetRHSToZero_inviscid!(params) 
    
    #
    # Inviscid part:
    #
    inviscid_rhs_el!(params.rhs_el, params.uaux, u,
                     params.F, params.G, params.S,
                     params.mesh, params.metrics,
                     params.basis, params.ω, params.SD;
                     neqs=neqs, lsource=params.inputs[:lsource])
    
    apply_boundary_conditions!(u, params, time)
    
    DSS_rhs!(params.SD, @view(params.RHS[:,:]), @view(params.rhs_el[:,:,:,:]),
             params.mesh.connijk, params.mesh.nelem, params.mesh.npoin, neqs, params.mesh.nop, Float64)

    #
    # Viscous part:
    #
    if (params.inputs[:lvisc] == true)

        resetRHSToZero_viscous!(params)

        μ = params.inputs[:νx]
        ν = 0.0
        κ = μ
        γ = 1.4
        Pr = 0.1
        for iel=1:nelem
            uToPrimitives!(iel, params.ρel, params.uel, params.vel, params.Tel, u, 
                           nelem, ngl, npoin, params.mesh.connijk, 0.0)
            
            for l = 1:params.mesh.ngl
                for k = 1:params.mesh.ngl
                    ωJac = fun_ωJac(params.ω[k], params.ω[l], params.metrics.Je[iel,k,l])
                    
                    dρdξ = 0.0
                    dudξ = 0.0
                    dvdξ = 0.0
                    dTdξ = 0.0

                    dρdη = 0.0
                    dudη = 0.0
                    dvdη = 0.0
                    dTdη = 0.0
                    for i = 1:params.mesh.ngl
                        dρdξ += params.basis.dψ[i,k]*params.ρel[i,l]
                        dudξ += params.basis.dψ[i,k]*params.uel[i,l]
                        dvdξ += params.basis.dψ[i,k]*params.vel[i,l]
                        dTdξ += params.basis.dψ[i,k]*params.Tel[i,l]

                        dρdη += params.basis.dψ[i,l]*params.ρel[k,i]
                        dudη += params.basis.dψ[i,l]*params.uel[k,i]
                        dvdη += params.basis.dψ[i,l]*params.vel[k,i]
                        dTdη += params.basis.dψ[i,l]*params.Tel[k,i]
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
                        
                        params.rhs_diffξ_el[iel,i,l,1] -= ωJac*dhdξ_ik*∇ξ∇ρ_kl
                        params.rhs_diffη_el[iel,k,i,1] -= ωJac*dhdη_il*∇η∇ρ_kl
                        
                        params.rhs_diffξ_el[iel,i,l,2] -= ωJac*dhdξ_ik*∇ξ∇u_kl
                        params.rhs_diffη_el[iel,k,i,2] -= ωJac*dhdη_il*∇η∇u_kl
                        
                        params.rhs_diffξ_el[iel,i,l,3] -= ωJac*dhdξ_ik*∇ξ∇v_kl
                        params.rhs_diffη_el[iel,k,i,3] -= ωJac*dhdη_il*∇η∇v_kl
                        
                        params.rhs_diffξ_el[iel,i,l,4] -= ωJac*dhdξ_ik*∇ξ∇T_kl
                        params.rhs_diffη_el[iel,k,i,4] -= ωJac*dhdη_il*∇η∇T_kl
                        
                    end
                end
            end
        end
        
        # params.rhs_diff_el .= @views (params.rhs_diffξ_el[:,:,:,:] + params.rhs_diffη_el[:,:,:,:])
        
        #= @btime viscous_rhs_el!($params.rhs_diff_el, $params.uaux, $u,
                               $params.ρel, $params.uel, $params.vel, $params.Tel,
                               $params.mesh, $params.metrics,
                               $params.basis, $params.ω, $neqs,
                               $params.inputs, $params.SD)
       =# 
        #$νx=params.inputs[:νx], νy=params.inputs[:νy])
        
        #=DSS_rhs!(params.SD, @view(params.RHS_visc[:,:]), @views(params.rhs_diff_el[:,:,:,:]), 
                 params.mesh.connijk, params.mesh.nelem, params.mesh.npoin,
                 neqs, params.mesh.nop, T)
        
        params.RHS .= params.RHS .+ params.RHS_visc
        =#
    end
    
    divive_by_mass_matrix!(params.RHS, params.M, params.QT, params.neqs)
    
    return params.RHS
end

#--------------------------------------------------------------------------------------------------------------------------------------------------
# CompEuler:
#--------------------------------------------------------------------------------------------------------------------------------------------------
#
# Optimized (more coud possibly be done)
#

function build_rhs!(RHS, u, params, time)   

    _build_rhs!(RHS, u, params, time)
    
end
