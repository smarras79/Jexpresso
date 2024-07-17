@kernel function _build_rhs_gpu_v0!(RHS, u, uaux, qe, x, t, connijk, dψ, ω, Minv, flux, source, PhysConst, xmax, xmin, n_x, neq, lpert, lperiodic_1d, npoin_linear, npoin)
    ie = @index(Group, Linear)
    i = @index(Local, Linear)
    ip = connijk[ie,i,1]
   
    T = eltype(RHS)
    KernelAbstractions.@atomic RHS[ip] = T(0.0)
    DIM = @uniform @groupsize()[1]
    ### define and populate flux array as shared memory then make sure blocks are synchronized
    F = @localmem eltype(RHS) (DIM+1,1)
    S = @localmem eltype(RHS) (DIM+1,1)
    
    uip = @view(uaux[ip,1:neq])
    qeip = @view(qe[ip,1:neq+1])
    
    if (ip == 1 || ip == npoin_linear) && !(lperiodic_1d)
        @inbounds uaux[ip,1:neq] .= user_bc_dirichlet_gpu(uip,qeip,x[ip],t,lpert)
        for ieq =1:neq

            idx = ip + (ieq-1)*npoin
            @inbounds u[idx] = uaux[ip,ieq]
        end
    end
    @inbounds flux[ie, i, :] .= user_flux_gpu(uip, qeip, PhysConst, lpert)


    @inbounds source[ie, i, :] .= user_source_gpu(uip, qeip, x[ip], PhysConst, xmax, xmin, lpert)

    @synchronize()
    ### do numerical integration
    for ieq = 1:neq
        F[i] = flux[ie, i, ieq]
        S[i] = source[ie, i, ieq]
        dFdξ = zero(T)
        for k=1:n_x
            dFdξ += dψ[k,i]*F[k]
        end
    ### Adding to rhs, DSS and division by the mass matrix can all be done in one combined step
        KernelAbstractions.@atomic RHS[ip,ieq] -= ω[i]*(dFdξ - S[i])* Minv[ip]
    end
end

@kernel function _build_rhs_gpu_2D_v0!(RHS, u, qe, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, ω, Minv, flux, source, ngl, neq, PhysConst, xmax, xmin, ymax, ymin, lpert)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i_x = il[1]
    @inbounds i_y = il[2]
    @inbounds ip = connijk[ie,i_x,i_y]

    T = eltype(RHS)
    DIM = @uniform @groupsize()[1]
    ### define and populate flux array as shared memory then make sure blocks are synchronized
    F = @localmem eltype(RHS) (DIM+1,DIM+1)
    G = @localmem eltype(RHS) (DIM+1,DIM+1)
    S = @localmem eltype(RHS) (DIM+1,DIM+1)

    uip = @view(u[ip,1:neq])
    qeip = @view(qe[ip,1:neq+1])
    @inbounds flux[ie, i_x, i_y, :] .= user_flux_gpu(uip,qeip,PhysConst,lpert)
    
    
    @inbounds source[ie, i_x, i_y, :] .= user_source_gpu(uip,qeip,x[ip],y[ip],PhysConst, xmax, xmin, ymax, ymin,lpert)

    @synchronize()
    ### do numerical integration
    for ieq =1:neq
        @inbounds F[i_x,i_y] = flux[ie, i_x, i_y, ieq]
        @inbounds G[i_x,i_y] = flux[ie, i_x, i_y, neq+ieq]
        @inbounds S[i_x,i_y] = source[ie, i_x, i_y, ieq]
        @synchronize()
        dFdξ = zero(T)
        dFdη = zero(T)
        dGdξ = zero(T)
        dGdη = zero(T)

        for k=1:ngl
            @inbounds dFdξ += dψ[k,i_x]*F[k,i_y]
            @inbounds dFdη += dψ[k,i_y]*F[i_x,k]
            @inbounds dGdξ += dψ[k,i_x]*G[k,i_y]
            @inbounds dGdη += dψ[k,i_y]*G[i_x,k]
        end

        @inbounds dξdx_ij = dξdx[ie,i_x,i_y]
        @inbounds dξdy_ij = dξdy[ie,i_x,i_y]
        @inbounds dηdx_ij = dηdx[ie,i_x,i_y]
        @inbounds dηdy_ij = dηdy[ie,i_x,i_y]

        dFdx = dFdξ*dξdx_ij + dFdη*dηdx_ij
        dGdy = dGdξ*dξdy_ij + dGdη*dηdy_ij

    ### Adding to rhs, DSS and division by the mass matrix can all be done in one combined step
        @inbounds KernelAbstractions.@atomic RHS[ip,ieq] -= ω[i_x]*ω[i_y]*Je[ie,i_x,i_y]*((dFdx + dGdy)- S[i_x,i_y])* Minv[ip]
    end
end

@kernel function _build_rhs_gpu_3D_v0!(RHS, u, qe, x, y, z, connijk, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, Je, dψ, ω, Minv, flux, source, ngl, neq, PhysConst, 
        xmax, xmin, ymax, ymin, zmax, zmin, lpert)
    
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i_x = il[1]
    @inbounds i_y = il[2]
    @inbounds i_z = il[3]
    @inbounds ip = connijk[ie,i_x,i_y,i_z]

    T = eltype(RHS)
    DIM = @uniform @groupsize()[1]
    ### define and populate flux array as shared memory then make sure blocks are synchronized
    F = @localmem eltype(RHS) (DIM+1,DIM+1,DIM+1)
    G = @localmem eltype(RHS) (DIM+1,DIM+1,DIM+1)
    H = @localmem eltype(RHS) (DIM+1,DIM+1,DIM+1)
    S = @localmem eltype(RHS) (DIM+1,DIM+1,DIM+1)

    uip = @view(u[ip,1:neq])
    qeip = @view(qe[ip,1:neq+1])
    @inbounds flux[ie, i_x, i_y, i_z,:] .= user_flux_gpu(uip,qeip,PhysConst,lpert)

    @inbounds source[ie, i_x, i_y, i_z,:] .= user_source_gpu(uip,qeip,x[ip],y[ip],z[ip],PhysConst, xmax, xmin, ymax, ymin, zmax, zmin,lpert)

    @synchronize()
    ### do numerical integration
    for ieq =1:neq
        @inbounds F[i_x,i_y,i_z] = flux[ie, i_x, i_y, i_z, ieq]
        @inbounds G[i_x,i_y,i_z] = flux[ie, i_x, i_y, i_z, neq+ieq]
        @inbounds H[i_x,i_y,i_z] = flux[ie, i_x, i_y, i_z, 2*neq+ieq]
        @inbounds S[i_x,i_y,i_z] = source[ie, i_x, i_y, i_z, ieq]
        @synchronize()
        dFdξ = zero(T)
        dFdη = zero(T)
        dFdζ = zero(T)
        dGdξ = zero(T)
        dGdη = zero(T)
        dGdζ = zero(T)
        dHdξ = zero(T)
        dHdη = zero(T)
        dHdζ = zero(T)

        for k=1:ngl
            @inbounds dFdξ += dψ[k,i_x]*F[k,i_y,i_z]
            @inbounds dFdη += dψ[k,i_y]*F[i_x,k,i_z]
            @inbounds dFdζ += dψ[k,i_z]*F[i_x,i_y,k]
            
            @inbounds dGdξ += dψ[k,i_x]*G[k,i_y,i_z]
            @inbounds dGdη += dψ[k,i_y]*G[i_x,k,i_z]
            @inbounds dGdζ += dψ[k,i_z]*G[i_x,i_y,k]

            @inbounds dHdξ += dψ[k,i_x]*H[k,i_y,i_z]
            @inbounds dHdη += dψ[k,i_y]*H[i_x,k,i_z]
            @inbounds dHdζ += dψ[k,i_z]*H[i_x,i_y,k]
            
        end
        @synchronize()

        @inbounds dξdx_ijk = dξdx[ie,i_x,i_y,i_z]
        @inbounds dξdy_ijk = dξdy[ie,i_x,i_y,i_z]
        @inbounds dξdz_ijk = dξdz[ie,i_x,i_y,i_z]
        
        @inbounds dηdx_ijk = dηdx[ie,i_x,i_y,i_z]
        @inbounds dηdy_ijk = dηdy[ie,i_x,i_y,i_z]
        @inbounds dηdz_ijk = dηdz[ie,i_x,i_y,i_z]

        @inbounds dζdx_ijk = dζdx[ie,i_x,i_y,i_z]
        @inbounds dζdy_ijk = dζdy[ie,i_x,i_y,i_z]
        @inbounds dζdz_ijk = dζdz[ie,i_x,i_y,i_z]

        dFdx = dFdξ*dξdx_ijk + dFdη*dηdx_ijk + dFdζ*dζdx_ijk
        dGdy = dGdξ*dξdy_ijk + dGdη*dηdy_ijk + dGdζ*dζdy_ijk
        dHdz = dHdξ*dξdz_ijk + dHdη*dηdz_ijk + dHdζ*dζdz_ijk


    ### Adding to rhs, DSS and division by the mass matrix can all be done in one combined step
        @inbounds KernelAbstractions.@atomic RHS[ip,ieq] -= ω[i_x]*ω[i_y]*ω[i_z]*Je[ie,i_x,i_y,i_z]*((dFdx + dGdy + dHdz)- S[i_x,i_y,i_z])* Minv[ip]
        @synchronize()
    end
end


@kernel function _build_rhs_diff_gpu_3D_v0!(RHS_diff, rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el, u, qe, uprimitive, x, y, z, connijk, 
        dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, Je, dψ, ω, Minv, visc_coeff, ngl, neq, PhysConst, lpert)

    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i_x = il[1]
    @inbounds i_y = il[2]
    @inbounds i_z = il[3]
    @inbounds ip = connijk[ie,i_x,i_y,i_z]

    T = eltype(RHS_diff)
    DIM = @uniform @groupsize()[1]
    U = @localmem eltype(RHS_diff) (DIM+1,DIM+1,DIM+1)
    
    @inbounds uprimitive[ie, i_x, i_y, i_z, 1:neq] .= user_primitives_gpu(@view(u[ip,1:neq]),@view(qe[ip,1:neq]),lpert)
    
    @inbounds ωJac = ω[i_x]*ω[i_y]*ω[i_z]*Je[ie,i_x,i_y,i_z]

    for ieq=1:neq

        @inbounds U[i_x,i_y,i_z] = uprimitive[ie,i_x,i_y,i_z,ieq]
        
        @synchronize()
        
        dqdξ = zero(T)
        dqdη = zero(T)
        dqdζ = zero(T)

        for ii = 1:ngl
            @inbounds dqdξ += dψ[ii,i_x]*U[ii,i_y,i_z]
            @inbounds dqdη += dψ[ii,i_y]*U[i_x,ii,i_z]
            @inbounds dqdζ += dψ[ii,i_z]*U[i_x,i_y,ii]
        end
       
        @inbounds dξdx_klm = dξdx[ie,i_x,i_y,i_z]
        @inbounds dξdy_klm = dξdy[ie,i_x,i_y,i_z]
        @inbounds dξdz_klm = dξdz[ie,i_x,i_y,i_z]
        @inbounds dηdx_klm = dηdx[ie,i_x,i_y,i_z]
        @inbounds dηdy_klm = dηdy[ie,i_x,i_y,i_z]
        @inbounds dηdz_klm = dηdz[ie,i_x,i_y,i_z]
        @inbounds dζdx_klm = dζdx[ie,i_x,i_y,i_z]
        @inbounds dζdy_klm = dζdy[ie,i_x,i_y,i_z]
        @inbounds dζdz_klm = dζdz[ie,i_x,i_y,i_z]

        
        auxi = dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm
        @inbounds dqdx = visc_coeff[ieq]*auxi

        auxi = dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm
        @inbounds dqdy = visc_coeff[ieq]*auxi 
       
        auxi = dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm
        @inbounds dqdz = visc_coeff[ieq]*auxi

        ∇ξ∇u_klm = (dξdx_klm*dqdx + dξdy_klm*dqdy + dξdz_klm*dqdz)*ωJac
        ∇η∇u_klm = (dηdx_klm*dqdx + dηdy_klm*dqdy + dηdz_klm*dqdz)*ωJac
        ∇ζ∇u_klm = (dζdx_klm*dqdx + dζdy_klm*dqdy + dζdz_klm*dqdz)*ωJac

        for i = 1:ngl
            @inbounds dhdξ_ik = dψ[i,i_x]
            @inbounds dhdη_il = dψ[i,i_y]
            @inbounds dhdζ_im = dψ[i,i_z]

            @inbounds KernelAbstractions.@atomic rhs_diffξ_el[ie,i,i_y,i_z,ieq] -= dhdξ_ik * ∇ξ∇u_klm
            @inbounds KernelAbstractions.@atomic rhs_diffη_el[ie,i_x,i,i_z,ieq] -= dhdη_il * ∇η∇u_klm
            @inbounds KernelAbstractions.@atomic rhs_diffζ_el[ie,i_x,i_y,i,ieq] -= dhdζ_im * ∇ζ∇u_klm
        end
        @synchronize()
        @inbounds KernelAbstractions.@atomic RHS_diff[ip,ieq] += (rhs_diffξ_el[ie,i_x,i_y,i_z,ieq] + rhs_diffη_el[ie,i_x,i_y,i_z,ieq] + rhs_diffζ_el[ie,i_x,i_y,i_z,ieq])*Minv[ip]
        #@inbounds KernelAbstractions.@atomic RHS_diff[ip,ieq] += (rhs_diffξ_el[ie,i_x,i_y,i_z,ieq] + rhs_diffζ_el[ie,i_x,i_y,i_z,ieq])*Minv[ip]
    end

end

@kernel function _build_viscosity_nu_gpu_3D!(u, qe, uprimitive, nu, x, y, z, connijk, 
    dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, Je, dψ, ω, Minv, visc_coeff, ngl, neq, PhysConst, lpert)

    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i_x = il[1]
    @inbounds i_y = il[2]
    @inbounds i_z = il[3]
    @inbounds ip = connijk[ie,i_x,i_y,i_z]

    T = eltype(RHS_diff)
    DIM = @uniform @groupsize()[1]
    U = @localmem eltype(RHS_diff) (DIM+1,DIM+1,DIM+1)

    @inbounds uprimitive[ie, i_x, i_y, i_z, 1:neq] .= user_primitives_gpu(@view(u[ip,1:neq]),@view(qe[ip,1:neq]),lpert)

    @inbounds ωJac = ω[i_x]*ω[i_y]*ω[i_z]*Je[ie,i_x,i_y,i_z]


    @inbounds U[i_x,i_y,i_z] = uprimitive[ie,i_x,i_y,i_z,2]
    @inbounds V[i_x,i_y,i_z] = uprimitive[ie,i_x,i_y,i_z,3]
    @inbounds W[i_x,i_y,i_z] = uprimitive[ie,i_x,i_y,i_z,4]
    
    @synchronize()
    
    dudξ = zero(T)
    dudη = zero(T)
    dudζ = zero(T)

    dvdξ = zero(T)
    dvdη = zero(T)
    dvdζ = zero(T)

    dwdξ = zero(T)
    dwdη = zero(T)
    dwdζ = zero(T)

    for ii = 1:ngl
        @inbounds dudξ += dψ[ii,i_x]*U[ii,i_y,i_z]
        @inbounds dudη += dψ[ii,i_y]*U[i_x,ii,i_z]
        @inbounds dudζ += dψ[ii,i_z]*U[i_x,i_y,ii]


        @inbounds dvdξ += dψ[ii,i_x]*V[ii,i_y,i_z]
        @inbounds dvdη += dψ[ii,i_y]*V[i_x,ii,i_z]
        @inbounds dvdζ += dψ[ii,i_z]*V[i_x,i_y,ii]


        @inbounds dwdξ += dψ[ii,i_x]*W[ii,i_y,i_z]
        @inbounds dwdη += dψ[ii,i_y]*W[i_x,ii,i_z]
        @inbounds dwdζ += dψ[ii,i_z]*W[i_x,i_y,ii]
    end

    @inbounds dξdx_klm = dξdx[ie,i_x,i_y,i_z]
    @inbounds dξdy_klm = dξdy[ie,i_x,i_y,i_z]
    @inbounds dξdz_klm = dξdz[ie,i_x,i_y,i_z]
    @inbounds dηdx_klm = dηdx[ie,i_x,i_y,i_z]
    @inbounds dηdy_klm = dηdy[ie,i_x,i_y,i_z]
    @inbounds dηdz_klm = dηdz[ie,i_x,i_y,i_z]
    @inbounds dζdx_klm = dζdx[ie,i_x,i_y,i_z]
    @inbounds dζdy_klm = dζdy[ie,i_x,i_y,i_z]
    @inbounds dζdz_klm = dζdz[ie,i_x,i_y,i_z]

    
    dudx = dudξ*dξdx_klm + dudη*dηdx_klm + dudζ*dζdx_klm
    dvdx = dvdξ*dξdx_klm + dvdη*dηdx_klm + dvdζ*dζdx_klm
    dwdx = dwdξ*dξdx_klm + dwdη*dηdx_klm + dwdζ*dζdx_klm

    dudy = dudξ*dξdy_klm + dudη*dηdy_klm + dudζ*dζdy_klm
    dvdy = dvdξ*dξdy_klm + dvdη*dηdy_klm + dvdζ*dζdy_klm
    dwdy = dwdξ*dξdy_klm + dwdη*dηdy_klm + dwdζ*dζdy_klm

    dudz = dudξ*dξdz_klm + dudη*dηdz_klm + dudζ*dζdz_klm
    dvdz = dvdξ*dξdz_klm + dvdη*dηdz_klm + dvdζ*dζdz_klm
    dwdz = dwdξ*dξdz_klm + dwdη*dηdz_klm + dwdζ*dζdz_klm

    S11 = dudx
    S12 = (dudy + dvdx) * 0.5
    S13 = (dudz + dwdx) * 0.5
    S21 = S12
    S22 = dvdy
    S23 = (dvdz + dwdy) * 0.5
    S31 = S13
    S32 = S23
    S33 = dwdz
    # |Sij|
    Sij = sqrt(2.0 * (S11*S11 + S12*S12 + S13*S13 + S21*S21 + S22*S22 + S23*S23 + S31*S31 + S32*S32 + S33*S33))
    delta2::T = 10000.0

    for ieq=1:neq
        @inbounds nu[ie,i_x,i_y,i_z,ieq] = visc_coeff[ieq] * Sij * delta2
    end
end


@kernel function _build_rhs_diff_gpu_3D_v1!(RHS_diff, rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el, u, qe, uprimitive, x, y, z, connijk, 
    dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, Je, dψ, ω, Minv, visc_coeff, ngl, neq, PhysConst, lpert)

    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i_x = il[1]
    @inbounds i_y = il[2]
    @inbounds i_z = il[3]
    @inbounds ip = connijk[ie,i_x,i_y,i_z]

    T = eltype(RHS_diff)
    DIM = @uniform @groupsize()[1]
    U = @localmem eltype(RHS_diff) (DIM+1,DIM+1,DIM+1)
    V = @localmem eltype(RHS_diff) (DIM+1,DIM+1,DIM+1)
    W = @localmem eltype(RHS_diff) (DIM+1,DIM+1,DIM+1)

    @inbounds uprimitive[ie, i_x, i_y, i_z, 1:neq] .= user_primitives_gpu(@view(u[ip,1:neq]),@view(qe[ip,1:neq]),lpert)

    @inbounds ωJac = ω[i_x]*ω[i_y]*ω[i_z]*Je[ie,i_x,i_y,i_z]


    @inbounds U[i_x,i_y,i_z] = uprimitive[ie,i_x,i_y,i_z,2]
    @inbounds V[i_x,i_y,i_z] = uprimitive[ie,i_x,i_y,i_z,3]
    @inbounds W[i_x,i_y,i_z] = uprimitive[ie,i_x,i_y,i_z,4]
    
    @synchronize()
    
    dudξ = zero(T)
    dudη = zero(T)
    dudζ = zero(T)

    dvdξ = zero(T)
    dvdη = zero(T)
    dvdζ = zero(T)

    dwdξ = zero(T)
    dwdη = zero(T)
    dwdζ = zero(T)

    for ii = 1:ngl
        @inbounds dudξ += dψ[ii,i_x]*U[ii,i_y,i_z]
        @inbounds dudη += dψ[ii,i_y]*U[i_x,ii,i_z]
        @inbounds dudζ += dψ[ii,i_z]*U[i_x,i_y,ii]


        @inbounds dvdξ += dψ[ii,i_x]*V[ii,i_y,i_z]
        @inbounds dvdη += dψ[ii,i_y]*V[i_x,ii,i_z]
        @inbounds dvdζ += dψ[ii,i_z]*V[i_x,i_y,ii]


        @inbounds dwdξ += dψ[ii,i_x]*W[ii,i_y,i_z]
        @inbounds dwdη += dψ[ii,i_y]*W[i_x,ii,i_z]
        @inbounds dwdζ += dψ[ii,i_z]*W[i_x,i_y,ii]
    end
    @synchronize()

    @inbounds dξdx_klm = dξdx[ie,i_x,i_y,i_z]
    @inbounds dξdy_klm = dξdy[ie,i_x,i_y,i_z]
    @inbounds dξdz_klm = dξdz[ie,i_x,i_y,i_z]
    @inbounds dηdx_klm = dηdx[ie,i_x,i_y,i_z]
    @inbounds dηdy_klm = dηdy[ie,i_x,i_y,i_z]
    @inbounds dηdz_klm = dηdz[ie,i_x,i_y,i_z]
    @inbounds dζdx_klm = dζdx[ie,i_x,i_y,i_z]
    @inbounds dζdy_klm = dζdy[ie,i_x,i_y,i_z]
    @inbounds dζdz_klm = dζdz[ie,i_x,i_y,i_z]

    
    dudx = dudξ*dξdx_klm + dudη*dηdx_klm + dudζ*dζdx_klm
    dvdx = dvdξ*dξdx_klm + dvdη*dηdx_klm + dvdζ*dζdx_klm
    dwdx = dwdξ*dξdx_klm + dwdη*dηdx_klm + dwdζ*dζdx_klm

    dudy = dudξ*dξdy_klm + dudη*dηdy_klm + dudζ*dζdy_klm
    dvdy = dvdξ*dξdy_klm + dvdη*dηdy_klm + dvdζ*dζdy_klm
    dwdy = dwdξ*dξdy_klm + dwdη*dηdy_klm + dwdζ*dζdy_klm

    dudz = dudξ*dξdz_klm + dudη*dηdz_klm + dudζ*dζdz_klm
    dvdz = dvdξ*dξdz_klm + dvdη*dηdz_klm + dvdζ*dζdz_klm
    dwdz = dwdξ*dξdz_klm + dwdη*dηdz_klm + dwdζ*dζdz_klm

    S11 = dudx
    S12 = (dudy + dvdx) * 0.5
    S13 = (dudz + dwdx) * 0.5
    S21 = S12
    S22 = dvdy
    S23 = (dvdz + dwdy) * 0.5
    S31 = S13
    S32 = S23
    S33 = dwdz
    # |Sij|
    Sij = sqrt(2.0 * (S11*S11 + S12*S12 + S13*S13 + S21*S21 + S22*S22 + S23*S23 + S31*S31 + S32*S32 + S33*S33))
    delta2::T = 10000.0

    for ieq=1:neq

        @inbounds U[i_x,i_y,i_z] = uprimitive[ie,i_x,i_y,i_z,ieq]
        
        @synchronize()
        
        dqdξ = zero(T)
        dqdη = zero(T)
        dqdζ = zero(T)

        for ii = 1:ngl
            @inbounds dqdξ += dψ[ii,i_x]*U[ii,i_y,i_z]
            @inbounds dqdη += dψ[ii,i_y]*U[i_x,ii,i_z]
            @inbounds dqdζ += dψ[ii,i_z]*U[i_x,i_y,ii]
        end
        @synchronize()
        @inbounds dξdx_klm = dξdx[ie,i_x,i_y,i_z]
        @inbounds dξdy_klm = dξdy[ie,i_x,i_y,i_z]
        @inbounds dξdz_klm = dξdz[ie,i_x,i_y,i_z]
        @inbounds dηdx_klm = dηdx[ie,i_x,i_y,i_z]
        @inbounds dηdy_klm = dηdy[ie,i_x,i_y,i_z]
        @inbounds dηdz_klm = dηdz[ie,i_x,i_y,i_z]
        @inbounds dζdx_klm = dζdx[ie,i_x,i_y,i_z]
        @inbounds dζdy_klm = dζdy[ie,i_x,i_y,i_z]
        @inbounds dζdz_klm = dζdz[ie,i_x,i_y,i_z]

        
        auxi = dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm
        @inbounds dqdx = visc_coeff[ieq] * Sij * delta2 * auxi

        auxi = dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm
        @inbounds dqdy = visc_coeff[ieq] * Sij * delta2 * auxi 
    
        auxi = dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm
        @inbounds dqdz = visc_coeff[ieq] * Sij * delta2 * auxi

        ∇ξ∇u_klm = (dξdx_klm*dqdx + dξdy_klm*dqdy + dξdz_klm*dqdz)*ωJac
        ∇η∇u_klm = (dηdx_klm*dqdx + dηdy_klm*dqdy + dηdz_klm*dqdz)*ωJac
        ∇ζ∇u_klm = (dζdx_klm*dqdx + dζdy_klm*dqdy + dζdz_klm*dqdz)*ωJac

        for i = 1:ngl
            @inbounds dhdξ_ik = dψ[i,i_x]
            @inbounds dhdη_il = dψ[i,i_y]
            @inbounds dhdζ_im = dψ[i,i_z]

            @inbounds KernelAbstractions.@atomic rhs_diffξ_el[ie,i,i_y,i_z,ieq] -= dhdξ_ik * ∇ξ∇u_klm
            @inbounds KernelAbstractions.@atomic rhs_diffη_el[ie,i_x,i,i_z,ieq] -= dhdη_il * ∇η∇u_klm
            @inbounds KernelAbstractions.@atomic rhs_diffζ_el[ie,i_x,i_y,i,ieq] -= dhdζ_im * ∇ζ∇u_klm
        end
        @synchronize()
        @inbounds KernelAbstractions.@atomic RHS_diff[ip,ieq] += (rhs_diffξ_el[ie,i_x,i_y,i_z,ieq] + rhs_diffη_el[ie,i_x,i_y,i_z,ieq] + rhs_diffζ_el[ie,i_x,i_y,i_z,ieq])*Minv[ip]
        #@inbounds KernelAbstractions.@atomic RHS_diff[ip,ieq] += (rhs_diffξ_el[ie,i_x,i_y,i_z,ieq] + rhs_diffζ_el[ie,i_x,i_y,i_z,ieq])*Minv[ip]
    end

end

@kernel function _build_rhs_diff_gpu_2D_v0!(RHS_diff, rhs_diffξ_el, rhs_diffη_el, u, qe, uprimitive, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, ω, Minv, visc_coeff, ngl, neq, PhysConst, lpert)

    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    @inbounds i_x = il[1]
    @inbounds i_y = il[2]
    @inbounds ip = connijk[ie,i_x,i_y]

    T = eltype(RHS_diff)
    DIM = @uniform @groupsize()[1]
    U = @localmem eltype(RHS_diff) (DIM,DIM)

    @inbounds uprimitive[ie, i_x, i_y, 1:neq] .= user_primitives_gpu(@view(u[ip,1:neq]),@view(qe[ip,1:neq]),lpert)

    @inbounds ωJac = ω[i_x]*ω[i_y]*Je[ie,i_x,i_y]

    for ieq=1:neq

        @inbounds U[i_x,i_y] = uprimitive[ie,i_x,i_y,ieq]

        @synchronize()

        dqdξ = zero(T)
        dqdη = zero(T)

        for ii = 1:ngl
            @inbounds dqdξ += dψ[ii,i_x]*U[ii,i_y]
            @inbounds dqdη += dψ[ii,i_y]*U[i_x,ii]
        end

        @inbounds dξdx_kl = dξdx[ie,i_x,i_y]
        @inbounds dξdy_kl = dξdy[ie,i_x,i_y]
        @inbounds dηdx_kl = dηdx[ie,i_x,i_y]
        @inbounds dηdy_kl = dηdy[ie,i_x,i_y]

        auxi = dqdξ*dξdx_kl + dqdη*dηdx_kl
        @inbounds dqdx = visc_coeff[ieq]*auxi

        auxi = dqdξ*dξdy_kl + dqdη*dηdy_kl
        @inbounds dqdy = visc_coeff[ieq]*auxi

        ∇ξ∇u_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy)*ωJac
        ∇η∇u_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy)*ωJac

        for i = 1:ngl
            @inbounds dhdξ_ik = dψ[i,i_x]
            @inbounds dhdη_il = dψ[i,i_y]

            @inbounds KernelAbstractions.@atomic rhs_diffξ_el[ie,i,i_y,ieq] -= dhdξ_ik * ∇ξ∇u_kl
            @inbounds KernelAbstractions.@atomic rhs_diffη_el[ie,i_x,i,ieq] -= dhdη_il * ∇η∇u_kl
        end
        @synchronize()
        @inbounds KernelAbstractions.@atomic RHS_diff[ip,ieq] += (rhs_diffξ_el[ie,i_x,i_y,ieq] + rhs_diffη_el[ie,i_x,i_y,ieq])*Minv[ip]
    end

end


@kernel function apply_boundary_conditions_gpu!(uaux,u,qe,x,y,t,nx,ny,poin_in_bdy_edge,qbdy,ngl,neq,npoin,lpert)

    iedge = @index(Group, Linear)
    ik = @index(Local, Linear)
    @inbounds ip = poin_in_bdy_edge[iedge,ik]
   
    T = eltype(u)
    @inbounds qbdy[iedge,ik,1:neq] .= T(1234567)
    @inbounds qbdy[iedge,ik,1:neq] .= user_bc_dirichlet_gpu(@view(uaux[ip,:]),@view(qe[ip,:]),x[ip],y[ip],t,nx[iedge,ik],ny[iedge,ik],@view(qbdy[iedge,ik,:]),lpert)
    for ieq =1:neq
        if !(qbdy[iedge,ik,ieq] == T(1234567)) && !(qbdy[iedge,ik,ieq] == uaux[ip,ieq])
            # if use the commented line in CUDA, somehow get errors
            # @inbounds KernelAbstractions.@atomic u[(ieq-1)*npoin+ip] = qbdy[iedge, ik, ieq] 
            @inbounds u[(ieq-1)*npoin+ip] = qbdy[iedge, ik, ieq] 
        end
    end
end

@kernel function apply_boundary_conditions_gpu_3D!(uaux,u,qe,x,y,z,t,nx,ny,nz,poin_in_bdy_face,qbdy,ngl,neq,npoin,lpert)

    iface = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    @inbounds ip = poin_in_bdy_face[iface,i_x,i_y]

    T = eltype(u)
    @inbounds qbdy[iface,i_x,i_y,1:neq] .= T(123456.0)
    @inbounds qbdy[iface,i_x,i_y,1:neq] .= user_bc_dirichlet_gpu(@view(uaux[ip,:]),@view(qe[ip,:]),x[ip],y[ip],z[ip],t,nx[iface,i_x,i_y],ny[iface,i_x,i_y],nz[iface,i_x,i_y],
                                                                 @view(qbdy[iface,i_x,i_y,:]),lpert)
    for ieq =1:neq
        if !(qbdy[iface,i_x,i_y,ieq] == T(123456.0)) && !(qbdy[iface,i_x,i_y,ieq] == uaux[ip,ieq])
            @inbounds u[(ieq-1)*npoin+ip] = qbdy[iface, i_x,i_y, ieq]
        end
    end
end

@kernel function utouaux_gpu!(u,uaux,npoin,neq)
    id = @index(Global, NTuple)
    @inbounds ip = id[1]
    @inbounds ieq = id[2]
    idx = (ieq-1)*npoin + ip
    @inbounds uaux[ip,ieq] = u[idx]
end

@kernel function uauxtou_gpu!(u,uaux,npoin,neq)
    id = @index(Global, NTuple)
    @inbounds ip = id[1]
    @inbounds ieq = id[2]
    idx = (ieq-1)*npoin + ip
    @inbounds u[idx] = uaux[ip,ieq]
end

@kernel function RHStodu_gpu!(RHS,du,npoin,neq)
    id = @index(Global, NTuple)
    @inbounds ip = id[1]
    @inbounds ieq = id[2]
    idx = (ieq-1)*npoin + ip
    @inbounds du[idx] = RHS[ip,ieq]
end

@kernel function lin_solve_rhs_gpu_2d!(RHS, u, qe, x, y, neqs)
    ip = @index(Global, Linear)

    uip = @view(u[ip,:])
    qeip = @view(qe[ip,:])

    RHS[ip,:] .= user_source_gpu(uip, qeip, x[ip], y[ip])

end

@kernel function apply_boundary_conditions_gpu_lin_solve!(RHS, A, poin_in_bdy_edge, npoin)

    iedge = @index(Group, Linear)
    ik = @index(Local, Linear)
    @inbounds ip = poin_in_bdy_edge[iedge,ik]

    T = eltype(RHS)
    for jp = 1:npoin
        A[ip,jp] = T(0.0)
    end
    A[ip,ip] = T(1.0)
    RHS[ip,:] .= T(0.0)
end
