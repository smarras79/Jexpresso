#---------------------------------------------------------------------------
# Optimized (more coud possibly be done)
#---------------------------------------------------------------------------
@kernel function _build_rhs_lag_gpu_2D_v0!(RHS, u, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, dψ_lag, ω, ω_lag, Minv, flux, source, ngl, ngr, neq, PhysConst)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i_ngl = il[1]
    @inbounds i_ngr = il[2]
    @inbounds ip = connijk[ie,i_ngl,i_ngr]

    DIM_1 = @uniform @groupsize()[1]
    DIM_2 = @uniform @groupsize()[2]
    ### define and populate flux array as shared memory then make sure blocks are synchronized
    F = @localmem eltype(RHS) (DIM_1+1, DIM_2+1)
    G = @localmem eltype(RHS) (DIM_1+1, DIM_2+1)
    S = @localmem eltype(RHS) (DIM_1+1, DIM_2+1)

    uip = @view(u[ip,1:neq])
    @inbounds flux[ie, i_ngl, i_ngr, :] .= user_flux(uip,PhysConst)
    
    @inbounds source[ie, i_ngl, i_ngr, :] .= user_source(uip,x[ip],y[ip],PhysConst)

    @synchronize()
    ### do numerical integration
    for ieq =1:neq
        @inbounds F[i_ngl,i_ngr] = flux[ie, i_ngl, i_ngr, ieq]
        @inbounds G[i_ngl,i_ngr] = flux[ie, i_ngl, i_ngr, neq+ieq]
        @inbounds S[i_ngl,i_ngr] = source[ie, i_ngl, i_ngr, ieq]
        @synchronize()
        dFdξ = zero(Float32)
        dFdη = zero(Float32)
        dGdξ = zero(Float32)
        dGdη = zero(Float32)

        for k=1:ngl
            @inbounds dFdξ += dψ[k,i_ngl]*F[k, i_ngr]
            @inbounds dGdξ += dψ[k,i_ngl]*G[k, i_ngr]
        end
        for k=1:ngr
            @inbounds dFdη += dψ_lag[k,i_ngr]*F[i_ngl, k]
            @inbounds dGdη += dψ_lag[k,i_ngr]*G[i_ngl, k]
        end

        @inbounds dξdx_ij = dξdx[ie,i_ngl,i_ngr]
        @inbounds dξdy_ij = dξdy[ie,i_ngl,i_ngr]
        @inbounds dηdx_ij = dηdx[ie,i_ngl,i_ngr]
        @inbounds dηdy_ij = dηdy[ie,i_ngl,i_ngr]

        dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
        dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij

    ### Adding to rhs, DSS and division by the mass matrix can all be done in one combined step
        @inbounds KernelAbstractions.@atomic RHS[ip,ieq] -= ω[i_ngl]*ω_lag[i_ngr]*Je[ie,i_ngl,i_ngr]*((dFdx + dGdy)- S[i_ngl, i_ngr])* Minv[ip]
    end
end


@kernel function _build_rhs_visc_lag_gpu_2D_v0!(RHS_diff, rhs_diffξ_el, rhs_diffη_el, u, uprimitive, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, dψ_lag, ω, ω_lag, Minv, visc_coeff, ngl, ngr, neq, PhysConst)

    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i_ngl = il[1]
    @inbounds i_ngr = il[2]
    @inbounds ip = connijk[ie,i_ngl,i_ngr]

    DIM_1 = @uniform @groupsize()[1]
    DIM_2 = @uniform @groupsize()[2]
    U = @localmem eltype(RHS_diff) (DIM_1+1, DIM_2+1)

    @inbounds uprimitive[ie, i_ngl, i_ngr, 1:neq] .= uToPrimitives_lag!(@view(u[ip,1:neq]))

    @inbounds ωJac = ω[i_ngl]*ω_lag[i_ngr]*Je[ie,i_ngl,i_ngr]

    for ieq=1:neq

        @inbounds U[i_ngl,i_ngr] = uprimitive[ie,i_ngl,i_ngr,ieq]

        @synchronize()

        dqdξ = zero(Float32)
        dqdη = zero(Float32)

        for ii = 1:ngl
            @inbounds dqdξ += dψ[ii,i_ngl]*U[ii,i_ngr]
        end
        for ii = 1:ngr
            @inbounds dqdη += dψ_lag[ii,i_ngr]*U[i_ngl,ii]
        end

        @inbounds dξdx_kl = dξdx[ie,i_ngl,i_ngr]
        @inbounds dξdy_kl = dξdy[ie,i_ngl,i_ngr]
        @inbounds dηdx_kl = dηdx[ie,i_ngl,i_ngr]
        @inbounds dηdy_kl = dηdy[ie,i_ngl,i_ngr]

        auxi = dqdξ*dξdx_kl + dqdη*dηdx_kl
        @inbounds dqdx = visc_coeff[ieq]*auxi

        auxi = dqdξ*dξdy_kl + dqdη*dηdy_kl
        @inbounds dqdy = visc_coeff[ieq]*auxi

        ∇ξ∇u_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
        ∇η∇u_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac

        for i = 1:ngl
            @inbounds dhdξ_ik = dψ[i,i_ngl]
            @inbounds KernelAbstractions.@atomic rhs_diffξ_el[ie,i,i_ngr,ieq] -= dhdξ_ik * ∇ξ∇u_kl
        end
        for i = 1:ngr
            @inbounds dhdη_il = dψ_lag[i,i_ngr]
            @inbounds KernelAbstractions.@atomic rhs_diffη_el[ie,i_ngl,i,ieq] -= dhdη_il * ∇η∇u_kl
        end
        @synchronize()
        @inbounds KernelAbstractions.@atomic RHS_diff[ip,ieq] += (rhs_diffξ_el[ie,i_ngl,i_ngr,ieq] + rhs_diffη_el[ie,i_ngl,i_ngr,ieq])*Minv[ip]
    end

end

function uToPrimitives_lag!(u)

    return Float32(u[1]), Float32(u[2]/u[1]), Float32(u[3]/u[1]), Float32(u[4]/u[1])
end