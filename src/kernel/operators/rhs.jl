#---------------------------------------------------------------------------
# Fetch equations name to access the user_rhs functions
#---------------------------------------------------------------------------
user_source_dir = "../../fallbacks/source.jl"
user_flux_dir   = "../../fallbacks/flux.jl"
if (length(ARGS) === 1) #equations
    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/user_source.jl"))
        user_source_dir = string(@__DIR__, "../../equations/", ARGS[1], "/user_source.jl")
    else
        user_source_dir = "../../fallbacks/source.jl"
    end

    if isfile(string(@__DIR__, "/../../equations/", ARGS[1], "/user_flux.jl"))
        user_flux_dir = string(@__DIR__, "../../equations/", ARGS[1], "/user_flux.jl")
    else
        user_flux_dir = "../../fallbacks/flux.jl"
    end
    
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
        user_flux_dir = "../../fallbacks/flux.jl"
    end
end
include(user_flux_dir)
include(user_source_dir)
#---------------------------------------------------------------------------
function rhs!(du, u, params, time)
    
    RHS = build_rhs(u, params, time)
    
    for i=1:params.neqs
        idx = (i-1)*params.mesh.npoin
        du[idx+1:i*params.mesh.npoin] = @view RHS[:,i]
    end  
    return du #This is already DSSed
end


function inviscid_rhs_el!(F, G, S, rhs_el, uaux, u, SD::NSD_2D, mesh, metrics, basis, ω; neqs, lsource=false)
    
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        uaux[:,i] = view(u, idx+1:i*mesh.npoin)
    end
    
    lsource = inputs[:lsource]
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
                ωJac = ω[i]*ω[j]*metrics.Je[iel,i,j]
                    
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
                    rhs_el[iel,i,j,ieq] -= ωJac*((dFdx + dGdy)  - S[i,j,ieq]) #gravity
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


function _build_rhs(u, params, time)

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
    inviscid_rhs_el!(params.F, params.G, params.S,
                     params.rhs_el, params.uaux, u,
                     params.SD, params.mesh, params.metrics,
                     params.basis, params.ω;
                     neqs, lsource=params.inputs[:lsource])
    
    apply_boundary_conditions!(u, params, time)
    
    DSS_rhs!(params.SD, @view(params.RHS[:,:]), @view(params.rhs_el[:,:,:,:]),
             params.mesh.connijk, params.mesh.nelem, params.mesh.npoin, neqs, params.mesh.nop, Float64)

    #
    # Viscous part:
    #
    #=if (inputs[:lvisc] == true)

        resetRHSToZero_viscous!(params)
        
        #=if (lowercase(params.inputs[:visc_model]) === "dsgs")
            
            if (rem(time, Δt) == 0 && time > 0.0)
                qnm1 .= qnm2
                qnm2 .= uaux
            end
            
            compute_viscosity!(params.qp.μ, params.SD, params.PT, params.uaux,
    params.u.qnm1, params.u.qnm2, params.RHS, params.Δt, mesh, metrics, T)
        else=#
        #    μ[:] .= inputs[:νx]
        ##end
       # for i=1:neqs
       #     idx = (i-1)*params.mesh.npoin
       #     for j=1:params.mesh.npoin
       #         u[idx+j] = params.uaux[j,i]
       #     end
       # end
=#
        @btime viscous_rhs_el!($params.rhs_diff_el, $params.uaux, $u,
                               $params.SD, $params.mesh, $params.metrics,
                               $params.basis, $params.ω, $neqs, $params.inputs)
    #$νx=params.inputs[:νx], νy=params.inputs[:νy])
  #=      
        DSS_rhs!(params.SD, @view(params.RHS_visc[:,:]), @views(params.rhs_diff_el[:,:,:,:]), 
                 params.mesh.connijk, params.mesh.nelem, params.mesh.npoin,
                 neqs, params.mesh.nop, T)
        
        params.RHS .= params.RHS .+ params.RHS_visc
        
    end=#
    
    divive_by_mass_matrix!(params.RHS, params.M, params.QT, params.neqs)
    
    return params.RHS
end

#--------------------------------------------------------------------------------------------------------------------------------------------------
# CompEuler:
#--------------------------------------------------------------------------------------------------------------------------------------------------
#
# Optimized (more coud possibly be done)
#

function build_rhs(u, params, time)   

    RHS = _build_rhs(u, params, time)

    return RHS
    
end
