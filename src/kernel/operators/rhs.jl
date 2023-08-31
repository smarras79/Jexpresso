#---------------------------------------------------------------------------
# Fetch equations name to access the user_rhs functions
#---------------------------------------------------------------------------
#=
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
=#

#---------------------------------------------------------------------------
# Optimized (more coud possibly be done)
#---------------------------------------------------------------------------
function build_rhs!(RHS, u, params, time)
    #
    # build_rhs()! is called by TimeIntegrators.jl -> time_loop!() via ODEProblem(rhs!, u, tspan, params)
    #
    _build_rhs!(RHS, u, params, time)
    
end

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


function uaux2u!(u, uaux, neqs, npoin)

    for i=1:neqs
        idx = (i-1)*npoin
        for j=1:npoin
            u[idx+j] = uaux[j,i]
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


function uToPrimitives!(ρel, uel, vel, Tel, u, mesh, δtotal_energy, iel)
    
    for j=1:mesh.ngl, i=1:mesh.ngl
        
        m1 = mesh.connijk[iel,i,j]
        m2 = mesh.npoin + m1
        m3 = 2*mesh.npoin + m1
        m4 = 3*mesh.npoin + m1
        
        ρel[i,j] = u[m1]
        uel[i,j] = u[m2]/ρel[i,j]
        vel[i,j] = u[m3]/ρel[i,j]
        Tel[i,j] = u[m4]/ρel[i,j] - δtotal_energy*0.5*(uel[i,j]^2 + vel[i,j]^2)
    end
    
end

function rhs!(du, u, params, time)
    
    build_rhs!(@view(params.RHS[:,:]), u, params, time)
    
    RHStoDU!(du, @view(params.RHS[:,:]), params.neqs, params.mesh.npoin)
end


function _build_rhs!(RHS, u, params, time)

    T       = Float64
    SD      = params.SD
    neqs    = params.neqs
    ngl     = params.mesh.ngl
    nelem   = params.mesh.nelem
    npoin   = params.mesh.npoin
    
    # rhs_el, RHS -> 0.0
    resetRHSToZero_inviscid!(params) 
    
    #-----------------------------------------------------------------------------------
    # Inviscid rhs:
    #-----------------------------------------------------------------------------------
    inviscid_rhs_el!(params.rhs_el, params.uaux, u,
                     params.F, params.G, params.S,
                     params.mesh, params.metrics,
                     params.basis, params.ω, params.SD,
                     neqs, true)
    
    apply_boundary_conditions!(u, params.uaux, time,
                               params.mesh, params.metrics, params.basis,
                               params.rhs_el, params.ubdy,
                               params.ω, SD, neqs, params.inputs)
    
    DSS_rhs!(SD, @view(params.RHS[:,:]), @view(params.rhs_el[:,:,:,:]), params.mesh, nelem, ngl, neqs)
    
    #-----------------------------------------------------------------------------------
    # Viscous rhs:
    #-----------------------------------------------------------------------------------
    if (params.inputs[:lvisc] == true)
        
        resetRHSToZero_viscous!(params) 
        
        viscous_rhs_el!(params.ρel, params.uel, params.vel, params.Tel,
                        params.rhs_diff_el, params.rhs_diffξ_el, params.rhs_diffη_el, 
                        u,
                        params.mesh, params.metrics, params.basis,
                        params.visc_coeff, params.inputs,
                        params.ω, neqs, SD)
        
        DSS_rhs!(SD, @view(params.RHS_visc[:,:]), @view(params.rhs_diff_el[:,:,:,:]), params.mesh, nelem, ngl, neqs)

        params.RHS .= params.RHS .+ params.RHS_visc
        
    end
    
    divive_by_mass_matrix!(@view(params.RHS[:,:]), @view(params.M[:]), params.QT, neqs)
    
end


fun_ωJac(x, y, z) = x*y*z

function inviscid_rhs_el!(rhs_el, uaux, u, F, G, S, mesh, metrics, basis, ω, SD::NSD_2D, neqs, lsource)

    u2uaux!(uaux, u, neqs, mesh.npoin)
    
    for iel=1:mesh.nelem

        for j=1:mesh.ngl, i=1:mesh.ngl
            ip = mesh.connijk[iel,i,j]
            user_flux!(@view(F[i,j,:]), @view(G[i,j,:]), SD, @view(uaux[ip,:]), mesh; neqs=neqs)
            if lsource
                user_source!(@view(S[i,j,:]), @view(uaux[ip,:]), mesh.npoin; neqs=neqs)
            end
        end
        
      _expansion_inviscid!(@view(rhs_el[iel,:,:,:]), metrics, basis,
                           @view(F[:,:,:]), @view(G[:,:,:]), @view(S[:,:,:]),
                           ω, mesh.ngl, mesh.npoin, neqs, iel)
        
    end

    
end

function viscous_rhs_el!(ρel, uel, vel, Tel,
                         rhs_diff_el, rhs_diffξ_el, rhs_diffη_el,
                         u,
                         mesh, metrics, basis, visc_coeff, inputs, 
                         ω, neqs, SD::NSD_2D)
    
    for iel=1:mesh.nelem
        _expansion_visc!(rhs_diff_el, rhs_diffξ_el, rhs_diffη_el, ρel, uel, vel, Tel,
                         u, ω, mesh, basis, metrics, visc_coeff, inputs, iel)
    end
    
    rhs_diff_el .= @views (rhs_diffξ_el[:,:,:,:] .+ rhs_diffη_el[:,:,:,:])
end

function _expansion_inviscid!(rhs_el, metrics, basis, F, G, S,
                                ω, ngl, npoin, neqs, iel)


    for ieq = 1:neqs
        for j=1:ngl
            for i=1:ngl
                ωJac = ω[i]*ω[j]*metrics.Je[iel,i,j]
                
                dFdξ = 0.0
                dFdη = 0.0
                dGdξ = 0.0
                dGdη = 0.0

                for k = 1:ngl
                    dFdξ += basis.dψ[k,i]*F[k,j,ieq]
                    dFdη += basis.dψ[k,j]*F[i,k,ieq]
                    
                    dGdξ += basis.dψ[k,i]*G[k,j,ieq]
                    dGdη += basis.dψ[k,j]*G[i,k,ieq]
                end
                dξdx_ij = metrics.dξdx[iel,i,j]
                dξdy_ij = metrics.dξdy[iel,i,j]
                dηdx_ij = metrics.dηdx[iel,i,j]
                dηdy_ij = metrics.dηdy[iel,i,j]
                
                dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
                dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij

                auxi = ωJac*((dFdx + dGdy)  - S[i,j,ieq])
                rhs_el[i,j,ieq] -= auxi
            end
        end
    end
   
end

function _expansion_visc!(rhs_diff_el, rhs_diffξ_el, rhs_diffη_el, ρel, uel, vel, Tel, u, ω, mesh, basis, metrics, visc_coeff, inputs, iel)
    
    uToPrimitives!(ρel, uel, vel, Tel, u, mesh, inputs[:δtotal_energy], iel)
    
    for l = 1:mesh.ngl
        for k = 1:mesh.ngl
            ωJac = fun_ωJac(ω[k], ω[l], metrics.Je[iel,k,l])
            
            dρdξ = 0.0
            dudξ = 0.0
            dvdξ = 0.0
            dTdξ = 0.0
            
            dρdη = 0.0
            dudη = 0.0
            dvdη = 0.0
            dTdη = 0.0
            for i = 1:mesh.ngl
                
                dρdξ += basis.dψ[i,k]*ρel[i,l]
                dudξ += basis.dψ[i,k]*uel[i,l]
                dvdξ += basis.dψ[i,k]*vel[i,l]
                dTdξ += basis.dψ[i,k]*Tel[i,l]

                dρdη += basis.dψ[i,l]*ρel[k,i]
                dudη += basis.dψ[i,l]*uel[k,i]
                dvdη += basis.dψ[i,l]*vel[k,i]
                dTdη += basis.dψ[i,l]*Tel[k,i]
            end
            dξdx_kl = metrics.dξdx[iel,k,l]
            dξdy_kl = metrics.dξdy[iel,k,l]
            dηdx_kl = metrics.dηdx[iel,k,l]
            dηdy_kl = metrics.dηdy[iel,k,l]
            
            #dρdx =  visc_coeff.νρ*(dρdξ*dξdx_kl + dρdη*dηdx_kl)
            dudx =  visc_coeff.μx*(dudξ*dξdx_kl + dudη*dηdx_kl)
            dvdx =  visc_coeff.μx*(dvdξ*dξdx_kl + dvdη*dηdx_kl)
            dTdx =  visc_coeff.κ *(dTdξ*dξdx_kl + dTdη*dηdx_kl) #+μ∇u⋅u
            
            #dρdy =  visc_coeff.νρ*(dρdξ*dξdy_kl + dρdη*dηdy_kl)
            dudy =  visc_coeff.μy*(dudξ*dξdy_kl + dudη*dηdy_kl)
            dvdy =  visc_coeff.μy*(dvdξ*dξdy_kl + dvdη*dηdy_kl)
            dTdy =  visc_coeff.κ *(dTdξ*dξdy_kl + dTdη*dηdy_kl) #+μ∇u⋅u
            
            #∇ξ∇ρ_kl = dξdx_kl*dρdx + dξdy_kl*dρdy
            #∇η∇ρ_kl = dηdx_kl*dρdx + dηdy_kl*dρdy
            
            ∇ξ∇u_kl = dξdx_kl*dudx + dξdy_kl*dudy
            ∇η∇u_kl = dηdx_kl*dudx + dηdy_kl*dudy            
            ∇ξ∇v_kl = dξdx_kl*dvdx + dξdy_kl*dvdy
            ∇η∇v_kl = dηdx_kl*dvdx + dηdy_kl*dvdy

            ∇ξ∇T_kl = dξdx_kl*dTdx + dξdy_kl*dTdy
            ∇η∇T_kl = dηdx_kl*dTdx + dηdy_kl*dTdy

            for i = 1:mesh.ngl
                dhdξ_ik, dhdη_il = basis.dψ[i,k], basis.dψ[i,l]
                
                #rhs_diffξ_el[iel, i,l,1] -= ωJac*dhdξ_ik*∇ξ∇ρ_kl
                #rhs_diffη_el[iel, k,i,1] -= ωJac*dhdη_il*∇η∇ρ_kl
                
                rhs_diffξ_el[iel, i,l,2] -= ωJac*dhdξ_ik*∇ξ∇u_kl
                rhs_diffη_el[iel, k,i,2] -= ωJac*dhdη_il*∇η∇u_kl
                
                rhs_diffξ_el[iel, i,l,3] -= ωJac*dhdξ_ik*∇ξ∇v_kl
                rhs_diffη_el[iel, k,i,3] -= ωJac*dhdη_il*∇η∇v_kl
                
                rhs_diffξ_el[iel, i,l,4] -= ωJac*dhdξ_ik*∇ξ∇T_kl
                rhs_diffη_el[iel, k,i,4] -= ωJac*dhdη_il*∇η∇T_kl
            end
           
        end
    end
    
end

