@kernel function _build_rhs_gpu_v0!(RHS, u, connijk, dψ, ω, M, ngl)
    ig = @index(Group, Linear)
    il = @index(Local, Linear)
    ip = connijk[ig,il,1]
    
    KernelAbstractions.@atomic RHS[ip] = Float32(0.0)
    DIM = @uniform @groupsize()[1]
    ### define and populate flux array as shared memory then make sure blocks are synchronized
    F = @localmem eltype(RHS) (DIM+1,1)
    F[il] = Float32(1.0)*u[ip] #user_flux(u[ip])
    @synchronize()
    ### do numerical integration
    dFdxi = zero(Float32)
    for k=1:ngl
        dFdxi += dψ[k,il]*F[k]
    end
    ### Adding to rhs, DSS and division by the mass matrix can all be done in one combined step
    KernelAbstractions.@atomic RHS[ip] -= ω[il]*dFdxi/ M[ip]
end

@kernel function _build_rhs_gpu_2D_v0!(RHS, u, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, ω, Minv, flux, source, ngl, neq, PhysConst)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    ip = connijk[ie,i_x,i_y]

    DIM = @uniform @groupsize()[1]
    ### define and populate flux array as shared memory then make sure blocks are synchronized
    F = @localmem eltype(RHS) (DIM+1,DIM+1)
    G = @localmem eltype(RHS) (DIM+1,DIM+1)
    S = @localmem eltype(RHS) (DIM+1,DIM+1)

    flux[ie, i_x, i_y, :] .= user_flux(@view(u[ip,1:neq]),PhysConst)
    
    source[ie, i_x, i_y, :] .= user_source(@view(u[ip,1:neq]),x[ip],y[ip],PhysConst)
    @synchronize()
    ### do numerical integration
    for ieq =1:neq
        F[i_x,i_y] = flux[ie, i_x, i_y, ieq]
        G[i_x,i_y] = flux[ie, i_x, i_y, neq+ieq]
        S[i_x,i_y] = source[ie, i_x, i_y, ieq]
        @synchronize()
        dFdξ = zero(Float32)
        dFdη = zero(Float32)
        dGdξ = zero(Float32)
        dGdη = zero(Float32)

        for k=1:ngl
            dFdξ += dψ[k,i_x]*F[k,i_y]
            dFdη += dψ[k,i_y]*F[i_x,k]
            dGdξ += dψ[k,i_x]*G[k,i_y]
            dGdη += dψ[k,i_y]*G[i_x,k]
        end

        dξdx_ij = dξdx[ie,i_x,i_y]
        dξdy_ij = dξdy[ie,i_x,i_y]
        dηdx_ij = dηdx[ie,i_x,i_y]
        dηdy_ij = dηdy[ie,i_x,i_y]

        dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
        dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij

    ### Adding to rhs, DSS and division by the mass matrix can all be done in one combined step
        KernelAbstractions.@atomic RHS[ip,ieq] -= ω[i_x]*ω[i_y]*Je[ie,i_x,i_y]*((dFdx + dGdy)- S[i_x,i_y])* Minv[ip]
    end
end

@kernel function _build_rhs_diff_gpu_2D_v0!(RHS_diff, rhs_diffξ_el, rhs_diffη_el, u, uprimitive, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, ω, Minv, visc_coeff, ngl, neq, PhysConst)

    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    ip = connijk[ie,i_x,i_y]

    DIM = @uniform @groupsize()[1]
    U = @localmem eltype(RHS_diff) (DIM,DIM)
    
    uprimitive[ie, i_x, i_y, :] .= uToPrimitives_gpu(@view(u[ip,1:neq]))
    
    ωJac = ω[i_x]*ω[i_y]*Je[ie,i_x,i_y]

    for ieq=1:neq

        U[i_x,i_y] = uprimitive[ie,i_x,i_y,ieq]
        
        @synchronize()
        
        dqdξ = zero(Float32)
        dqdη = zero(Float32)

        for ii = 1:ngl
            dqdξ += dψ[ii,i_x]*U[ii,i_y]
            dqdη += dψ[ii,i_y]*U[i_x,ii]
        end
       
        dξdx_kl = dξdx[ie,i_x,i_y]
        dξdy_kl = dξdy[ie,i_x,i_y]
        dηdx_kl = dηdx[ie,i_x,i_y]
        dηdy_kl = dηdy[ie,i_x,i_y]
        
        auxi = dqdξ*dξdx_kl + dqdη*dηdx_kl
        dqdx = visc_coeff[ieq]*auxi

        auxi = dqdξ*dξdy_kl + dqdη*dηdy_kl
        dqdy = visc_coeff[ieq]*auxi
        
        ∇ξ∇u_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
        ∇η∇u_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac

        for i = 1:ngl
            dhdξ_ik = dψ[i,i_x]
            dhdη_il = dψ[i,i_y]

            KernelAbstractions.@atomic rhs_diffξ_el[ie,i,i_y,ieq] -= dhdξ_ik * ∇ξ∇u_kl
            KernelAbstractions.@atomic rhs_diffη_el[ie,i_x,i,ieq] -= dhdη_il * ∇η∇u_kl
        end
        @synchronize()
        KernelAbstractions.@atomic RHS_diff[ip,ieq] += (rhs_diffξ_el[ie,i_x,i_y,ieq] + rhs_diffη_el[ie,i_x,i_y,ieq])*Minv[ip]
    end

end


@kernel function apply_boundary_conditions_gpu!(uaux,u,x,y,t,nx,ny,poin_in_bdy_edge,qbdy,ngl,neq,npoin)

    iedge = @index(Group, Linear)
    ik = @index(Local, Linear)
    ip = poin_in_bdy_edge[iedge,ik]
    
    qbdy[iedge,ik,:] .= 1234567
    qbdy[iedge,ik,:] .= user_bc_dirichlet(@view(uaux[ip,:]),x[ip],y[ip],t,nx[iedge,ik],ny[iedge,ik],@view(qbdy[iedge,ik,:]))
    for ieq =1:neq
        if !(qbdy[iedge,ik,ieq] == 1234567) && !(qbdy[iedge,ik,ieq] == uaux[ip,ieq])
            KernelAbstractions.@atomic u[(ieq-1)*npoin+ip] = qbdy[iedge, ik, ieq] 
        end
    end
end

@kernel function utouaux_gpu!(u,uaux,npoin,neq)
    id = @index(Global, NTuple)
    ip = id[1]
    ieq = id[2]
    idx = (ieq-1)*npoin + ip
    uaux[ip,ieq] = u[idx]
end

function uToPrimitives_gpu(u)

    return Float32(u[1]), Float32(u[2]/u[1]), Float32(u[3]/u[1]), Float32(u[4]/u[1])
end

