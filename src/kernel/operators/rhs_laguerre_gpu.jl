#---------------------------------------------------------------------------
# Optimized (more coud possibly be done)
#---------------------------------------------------------------------------
@kernel function _build_rhs_lag_gpu_2D_v0!(RHS, u, qe, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, dψ_lag, ω, ω_lag, Minv, flux, source, 
        ngl, ngr, neq, PhysConst, xmax, xmin, ymax, ymin,lpert)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i_ngl = il[1]
    @inbounds i_ngr = il[2]
    @inbounds ip = connijk[ie,i_ngl,i_ngr]

    T = eltype(RHS)

    DIM_1 = @uniform @groupsize()[1]
    DIM_2 = @uniform @groupsize()[2]
    ### define and populate flux array as shared memory then make sure blocks are synchronized
    F = @localmem eltype(RHS) (DIM_1+1, DIM_2+1)
    G = @localmem eltype(RHS) (DIM_1+1, DIM_2+1)
    S = @localmem eltype(RHS) (DIM_1+1, DIM_2+1)

    uip = @view(u[ip,1:neq])
    qeip = @view(qe[ip,1:neq+1])
    @inbounds flux[ie, i_ngl, i_ngr, :] .= user_flux_gpu(uip,qeip,PhysConst,lpert)
    
    @inbounds source[ie, i_ngl, i_ngr, :] .= user_source_gpu(uip,qeip,x[ip],y[ip],PhysConst, xmax, xmin, ymax, ymin,lpert)

    @synchronize()
    ### do numerical integration
    for ieq =1:neq
        @inbounds F[i_ngl,i_ngr] = flux[ie, i_ngl, i_ngr, ieq]
        @inbounds G[i_ngl,i_ngr] = flux[ie, i_ngl, i_ngr, neq+ieq]
        @inbounds S[i_ngl,i_ngr] = source[ie, i_ngl, i_ngr, ieq]
        @synchronize()
        dFdξ = zero(T)
        dFdη = zero(T)
        dGdξ = zero(T)
        dGdη = zero(T)

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

        @inbounds dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
        @inbounds dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij

    ### Adding to rhs, DSS and division by the mass matrix can all be done in one combined step
        @inbounds KernelAbstractions.@atomic RHS[ip,ieq] -= ω[i_ngl]*ω_lag[i_ngr]*Je[ie,i_ngl,i_ngr]*((dFdx + dGdy)- S[i_ngl, i_ngr])* Minv[ip]
    end
end


@kernel function _build_rhs_visc_lag_gpu_2D_v0!(RHS_diff, rhs_diffξ_el, rhs_diffη_el, u, qe, uprimitive, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je,
        dψ, dψ_lag, ω, ω_lag, Minv, visc_coeff, ngl, ngr, neq, PhysConst, lpert)

    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i_ngl = il[1]
    @inbounds i_ngr = il[2]
    @inbounds ip = connijk[ie,i_ngl,i_ngr]

    T = eltype(RHS_diff)

    DIM_1 = @uniform @groupsize()[1]
    DIM_2 = @uniform @groupsize()[2]
    U = @localmem eltype(RHS_diff) (DIM_1+1, DIM_2+1)

    @inbounds uprimitive[ie, i_ngl, i_ngr, 1:neq] .= user_primitives_gpu(@view(u[ip,1:neq]),@view(qe[ip,1:neq]),lpert)

    @inbounds ωJac = ω[i_ngl]*ω_lag[i_ngr]*Je[ie,i_ngl,i_ngr]

    for ieq=1:neq

        @inbounds U[i_ngl,i_ngr] = uprimitive[ie,i_ngl,i_ngr,ieq]

        @synchronize()

        dqdξ = zero(T)
        dqdη = zero(T)

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

@kernel function apply_boundary_conditions_lag_gpu!(uaux,u,qe,x,y,t,connijk,qbdy,ngl,ngr,neq,npoin,nelem_semi_inf,lperiodic,lpert)

    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_ngl = il[1]
    i_ngr = il[2]
    T = eltype(uaux)
    if (i_ngr == ngr)
        ip = connijk[ie,i_ngl,ngr]
        @inbounds qbdy[ie, i_ngl, i_ngr, 1:neq] .= T(1234567)
        @inbounds qbdy[ie, i_ngl, i_ngr, 1:neq] .= user_bc_dirichlet_gpu(@view(uaux[ip,:]),@view(qe[ip,:]),x[ip],y[ip],t,T(0.0),T(1.0),@view(qbdy[ie, i_ngl, i_ngr, :]),lpert)
        for ieq =1:neq
            if !(qbdy[ie, i_ngl, i_ngr, ieq] == T(1234567)) && !(qbdy[ie, i_ngl, i_ngr, ieq] == uaux[ip,ieq])
            # if use the commented line in CUDA, somehow get errors
            # @inbounds KernelAbstractions.@atomic u[(ieq-1)*npoin+ip] = qbdy[iedge, ik, ieq] 
                @inbounds KernelAbstractions.@atomic u[(ieq-1)*npoin+ip] = qbdy[ie, i_ngl, i_ngr, ieq]
            end
        end
    end
    
    if !(lperiodic)
        if (ie == Int32(1))
            ip = connijk[ie,1,i_ngr]
            @inbounds qbdy[ie, i_ngl, i_ngr, 1:neq] .= T(1234567)
            @inbounds qbdy[ie, i_ngl, i_ngr, 1:neq] .= user_bc_dirichlet_gpu(@view(uaux[ip,:]),@view(qe[ip,:]),x[ip],y[ip],t,T(-1.0),T(0.0),@view(qbdy[ie, i_ngl, i_ngr, :]),lpert)
            for ieq =1:neq
                if !(qbdy[ie, i_ngl, i_ngr, ieq] == T(1234567)) && !(qbdy[ie, i_ngl, i_ngr, ieq] == uaux[ip,ieq])
                # if use the commented line in CUDA, somehow get errors
                # @inbounds KernelAbstractions.@atomic u[(ieq-1)*npoin+ip] = qbdy[iedge, ik, ieq] 
                    @inbounds KernelAbstractions.@atomic u[(ieq-1)*npoin+ip] = qbdy[ie, i_ngl, i_ngr, ieq]
                end
            end
        end

        if (ie == nelem_semi_inf)
            ip = connijk[ie,ngl,i_ngr]
            @inbounds qbdy[ie, i_ngl, i_ngr, 1:neq] .= T(1234567)
            @inbounds qbdy[ie, i_ngl, i_ngr, 1:neq] .= user_bc_dirichlet_gpu(@view(uaux[ip,:]),@view(qe[ip,:]),x[ip],y[ip],t,T(1.0),T(0.0),@view(qbdy[ie, i_ngl, i_ngr, :]),lpert)
            for ieq =1:neq
                if !(qbdy[ie, i_ngl, i_ngr, ieq] == T(1234567)) && !(qbdy[ie, i_ngl, i_ngr, ieq] == uaux[ip,ieq])
                # if use the commented line in CUDA, somehow get errors
                # @inbounds KernelAbstractions.@atomic u[(ieq-1)*npoin+ip] = qbdy[iedge, ik, ieq]
                    @inbounds KernelAbstractions.@atomic u[(ieq-1)*npoin+ip] = qbdy[ie, i_ngl, i_ngr, ieq]
                end
            end
        end
    end


end

#=function uToPrimitives_lag!(u,qe)
    T = eltype(u)
    return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T((u[4]+qe[4])/(u[1]+qe[1]) - qe[4]/qe[1])
end=#
