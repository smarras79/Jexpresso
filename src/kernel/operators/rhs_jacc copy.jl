using JACC
using Atomix

function build_rhs_jacc!(idx, RHS, u, uaux, qe, x, t, connijk, dψ, ω, Minv, flux, source, PhysConst, xmax, xmin, n_x, neq, lpert, lperiodic_1d, npoin_linear, npoin)
    ie = (idx - 1) ÷ n_x + 1
    i  = (idx - 1) % n_x + 1
    ip = connijk[ie, i, 1]
    T  = eltype(RHS)

    uip  = @view(uaux[ip, 1:neq])
    qeip = @view(qe[ip, 1:neq+1])

    if (ip == 1 || ip == npoin_linear) && !(lperiodic_1d)
        uaux[ip, 1:neq] .= user_bc_dirichlet_gpu(uip, qeip, x[ip], t, lpert)
        for ieq = 1:neq
            uidx = ip + (ieq-1)*npoin
            @inbounds u[uidx] = uaux[ip, ieq]
        end
    end

    @inbounds flux[ie, i, :]   .= [user_flux_gpu(uip, qeip, PhysConst, lpert)...]
    @inbounds source[ie, i, :] .= [user_source_gpu(uip, qeip, x[ip], PhysConst, xmax, xmin, lpert)...]

    for ieq = 1:neq
        dFdξ = zero(T)
        Si   = source[ie, i, ieq]
        for k = 1:n_x
            dFdξ += dψ[k, i] * flux[ie, k, ieq]
        end
        rhs_update = ω[i] * (dFdξ - Si) * Minv[ip]
        Atomix.@atomic RHS[ip, ieq] -= rhs_update
    end
end

function _build_rhs_jacc_2D_v0!(idx, RHS, u, qe, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, ω, Minv, flux, source, ngl, neq, PhysConst, xmax, xmin, ymax, ymin, lpert)
    ngl2 = ngl * ngl
    ie   = (idx - 1) ÷ ngl2 + 1
    rem  = (idx - 1) % ngl2
    i_x  = rem ÷ ngl + 1
    i_y  = rem % ngl + 1
    T    = eltype(RHS)

    @inbounds ip   = connijk[ie, i_x, i_y]
    uip  = @view(u[ip, 1:neq])
    qeip = @view(qe[ip, 1:neq+1])

    @inbounds flux[ie, i_x, i_y, :]   .= user_flux_gpu(uip, qeip, PhysConst, lpert)
    @inbounds source[ie, i_x, i_y, :] .= user_source_gpu(uip, qeip, x[ip], y[ip], PhysConst, xmax, xmin, ymax, ymin, lpert)

    for ieq = 1:neq
        dFdξ = zero(T); dFdη = zero(T)
        dGdξ = zero(T); dGdη = zero(T)
        @inbounds S_ij = source[ie, i_x, i_y, ieq]

        for k = 1:ngl
            @inbounds dFdξ += dψ[k, i_x] * flux[ie, k,   i_y, ieq]
            @inbounds dFdη += dψ[k, i_y] * flux[ie, i_x, k,   ieq]
            @inbounds dGdξ += dψ[k, i_x] * flux[ie, k,   i_y, neq + ieq]
            @inbounds dGdη += dψ[k, i_y] * flux[ie, i_x, k,   neq + ieq]
        end

        @inbounds dξdx_ij = dξdx[ie, i_x, i_y]
        @inbounds dξdy_ij = dξdy[ie, i_x, i_y]
        @inbounds dηdx_ij = dηdx[ie, i_x, i_y]
        @inbounds dηdy_ij = dηdy[ie, i_x, i_y]

        dFdx = dFdξ * dξdx_ij + dFdη * dηdx_ij
        dGdy = dGdξ * dξdy_ij + dGdη * dηdy_ij

        rhs_update = ω[i_x] * ω[i_y] * Je[ie, i_x, i_y] * ((dFdx + dGdy) - S_ij) * Minv[ip]
        Atomix.@atomic RHS[ip, ieq] -= rhs_update
    end
    return nothing
end

function _build_rhs_jacc_3D_v0!(idx, RHS, u, qe, x, y, z, connijk, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, Je, dψ, ω, Minv, flux, source, ngl, neq, PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)
    ngl3 = ngl * ngl * ngl
    ie   = (idx - 1) ÷ ngl3 + 1
    rem  = (idx - 1) % ngl3
    i_x  = rem ÷ (ngl*ngl) + 1
    rem2 = rem % (ngl*ngl)
    i_y  = rem2 ÷ ngl + 1
    i_z  = rem2 % ngl + 1
    T    = eltype(RHS)

    @inbounds ip   = connijk[ie, i_x, i_y, i_z]
    uip  = @view(u[ip, 1:neq])
    qeip = @view(qe[ip, 1:neq+1])

    @inbounds flux[ie, i_x, i_y, i_z, :]   .= user_flux_gpu(uip, qeip, PhysConst, lpert)
    @inbounds source[ie, i_x, i_y, i_z, :] .= user_source_gpu(uip, qeip, x[ip], y[ip], z[ip], PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)

    for ieq = 1:neq
        dFdξ = zero(T); dFdη = zero(T); dFdζ = zero(T)
        dGdξ = zero(T); dGdη = zero(T); dGdζ = zero(T)
        dHdξ = zero(T); dHdη = zero(T); dHdζ = zero(T)
        @inbounds S_ijk = source[ie, i_x, i_y, i_z, ieq]

        for k = 1:ngl
            @inbounds dFdξ += dψ[k, i_x] * flux[ie, k,   i_y, i_z, ieq]
            @inbounds dFdη += dψ[k, i_y] * flux[ie, i_x, k,   i_z, ieq]
            @inbounds dFdζ += dψ[k, i_z] * flux[ie, i_x, i_y, k,   ieq]
            @inbounds dGdξ += dψ[k, i_x] * flux[ie, k,   i_y, i_z, neq+ieq]
            @inbounds dGdη += dψ[k, i_y] * flux[ie, i_x, k,   i_z, neq+ieq]
            @inbounds dGdζ += dψ[k, i_z] * flux[ie, i_x, i_y, k,   neq+ieq]
            @inbounds dHdξ += dψ[k, i_x] * flux[ie, k,   i_y, i_z, 2*neq+ieq]
            @inbounds dHdη += dψ[k, i_y] * flux[ie, i_x, k,   i_z, 2*neq+ieq]
            @inbounds dHdζ += dψ[k, i_z] * flux[ie, i_x, i_y, k,   2*neq+ieq]
        end

        @inbounds dξdx_ijk = dξdx[ie, i_x, i_y, i_z]
        @inbounds dξdy_ijk = dξdy[ie, i_x, i_y, i_z]
        @inbounds dξdz_ijk = dξdz[ie, i_x, i_y, i_z]
        @inbounds dηdx_ijk = dηdx[ie, i_x, i_y, i_z]
        @inbounds dηdy_ijk = dηdy[ie, i_x, i_y, i_z]
        @inbounds dηdz_ijk = dηdz[ie, i_x, i_y, i_z]
        @inbounds dζdx_ijk = dζdx[ie, i_x, i_y, i_z]
        @inbounds dζdy_ijk = dζdy[ie, i_x, i_y, i_z]
        @inbounds dζdz_ijk = dζdz[ie, i_x, i_y, i_z]

        dFdx = dFdξ*dξdx_ijk + dFdη*dηdx_ijk + dFdζ*dζdx_ijk
        dGdy = dGdξ*dξdy_ijk + dGdη*dηdy_ijk + dGdζ*dζdy_ijk
        dHdz = dHdξ*dξdz_ijk + dHdη*dηdz_ijk + dHdζ*dζdz_ijk

        Atomix.@atomic RHS[ip, ieq] -= ω[i_x]*ω[i_y]*ω[i_z]*Je[ie,i_x,i_y,i_z]*((dFdx + dGdy + dHdz) - S_ijk)*Minv[ip]
    end
    return nothing
end

function _build_rhs_jacc_3D_v1!(idx, RHS, u, qe, x, y, z, connijk, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, Je, dψ, ω, Minv, flux, source, ngl, neq, PhysConst, param_set, xmax, xmin, ymax, ymin, zmax, zmin, lpert)
    ngl3 = ngl * ngl * ngl
    ie   = (idx - 1) ÷ ngl3 + 1
    rem  = (idx - 1) % ngl3
    i_x  = rem ÷ (ngl*ngl) + 1
    rem2 = rem % (ngl*ngl)
    i_y  = rem2 ÷ ngl + 1
    i_z  = rem2 % ngl + 1
    T    = eltype(RHS)

    @inbounds ip   = connijk[ie, i_x, i_y, i_z]
    uip  = @view(u[ip, 1:neq])
    qeip = @view(qe[ip, 1:neq+1])

    @inbounds flux[ie, i_x, i_y, i_z, :]   .= user_flux_gpu(uip, qeip, z[ip], PhysConst, param_set, lpert)
    @inbounds source[ie, i_x, i_y, i_z, :] .= user_source_gpu(uip, qeip, x[ip], y[ip], z[ip], PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)

    for ieq = 1:neq
        dFdξ = zero(T); dFdη = zero(T); dFdζ = zero(T)
        dGdξ = zero(T); dGdη = zero(T); dGdζ = zero(T)
        dHdξ = zero(T); dHdη = zero(T); dHdζ = zero(T)
        @inbounds S_ijk = source[ie, i_x, i_y, i_z, ieq]

        for k = 1:ngl
            @inbounds dFdξ += dψ[k, i_x] * flux[ie, k,   i_y, i_z, ieq]
            @inbounds dFdη += dψ[k, i_y] * flux[ie, i_x, k,   i_z, ieq]
            @inbounds dFdζ += dψ[k, i_z] * flux[ie, i_x, i_y, k,   ieq]
            @inbounds dGdξ += dψ[k, i_x] * flux[ie, k,   i_y, i_z, neq+ieq]
            @inbounds dGdη += dψ[k, i_y] * flux[ie, i_x, k,   i_z, neq+ieq]
            @inbounds dGdζ += dψ[k, i_z] * flux[ie, i_x, i_y, k,   neq+ieq]
            @inbounds dHdξ += dψ[k, i_x] * flux[ie, k,   i_y, i_z, 2*neq+ieq]
            @inbounds dHdη += dψ[k, i_y] * flux[ie, i_x, k,   i_z, 2*neq+ieq]
            @inbounds dHdζ += dψ[k, i_z] * flux[ie, i_x, i_y, k,   2*neq+ieq]
        end

        @inbounds dξdx_ijk = dξdx[ie, i_x, i_y, i_z]
        @inbounds dξdy_ijk = dξdy[ie, i_x, i_y, i_z]
        @inbounds dξdz_ijk = dξdz[ie, i_x, i_y, i_z]
        @inbounds dηdx_ijk = dηdx[ie, i_x, i_y, i_z]
        @inbounds dηdy_ijk = dηdy[ie, i_x, i_y, i_z]
        @inbounds dηdz_ijk = dηdz[ie, i_x, i_y, i_z]
        @inbounds dζdx_ijk = dζdx[ie, i_x, i_y, i_z]
        @inbounds dζdy_ijk = dζdy[ie, i_x, i_y, i_z]
        @inbounds dζdz_ijk = dζdz[ie, i_x, i_y, i_z]

        dFdx = dFdξ*dξdx_ijk + dFdη*dηdx_ijk + dFdζ*dζdx_ijk
        dGdy = dGdξ*dξdy_ijk + dGdη*dηdy_ijk + dGdζ*dζdy_ijk
        dHdz = dHdξ*dξdz_ijk + dHdη*dηdz_ijk + dHdζ*dζdz_ijk

        Atomix.@atomic RHS[ip, ieq] -= ω[i_x]*ω[i_y]*ω[i_z]*Je[ie,i_x,i_y,i_z]*((dFdx + dGdy + dHdz) - S_ijk)*Minv[ip]
    end
    return nothing
end

function _build_precipitation_rhs_jacc_3D_v0!(idx, RHS, u, qe, x, y, z, connijk, dξdz, dηdz, dζdz, Je, dψ, ω, Minv, flux_micro, source_micro, ngl, neq, PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert, Pr, Ps, Pg, qi, qn, Tabs, S_micro, MicroConst)
    ngl3 = ngl * ngl * ngl
    ie   = (idx - 1) ÷ ngl3 + 1
    rem  = (idx - 1) % ngl3
    i_x  = rem ÷ (ngl*ngl) + 1
    rem2 = rem % (ngl*ngl)
    i_y  = rem2 ÷ ngl + 1
    i_z  = rem2 % ngl + 1
    T    = eltype(RHS)

    @inbounds ip   = connijk[ie, i_x, i_y, i_z]
    uip  = @view(u[ip, 1:neq])
    qeip = @view(qe[ip, 1:neq+1])

    @inbounds flux_micro[ie, i_x, i_y, i_z, :]   .= precipitation_flux_gpu(uip, qeip, MicroConst, lpert, Pr[ip], Ps[ip], Pg[ip], qi[ip])
    @inbounds source_micro[ie, i_x, i_y, i_z, :] .= precipitation_source_gpu(uip, qeip, lpert, qn[ip], S_micro[ip], PhysConst, MicroConst)

    for ieq = 4:neq
        dHdξ = zero(T); dHdη = zero(T); dHdζ = zero(T)
        @inbounds S_ijk = source_micro[ie, i_x, i_y, i_z, ieq-3]

        for k = 1:ngl
            @inbounds dHdξ += dψ[k, i_x] * flux_micro[ie, k,   i_y, i_z, ieq-3]
            @inbounds dHdη += dψ[k, i_y] * flux_micro[ie, i_x, k,   i_z, ieq-3]
            @inbounds dHdζ += dψ[k, i_z] * flux_micro[ie, i_x, i_y, k,   ieq-3]
        end

        @inbounds dξdz_ijk = dξdz[ie, i_x, i_y, i_z]
        @inbounds dηdz_ijk = dηdz[ie, i_x, i_y, i_z]
        @inbounds dζdz_ijk = dζdz[ie, i_x, i_y, i_z]

        dHdz = dHdξ*dξdz_ijk + dHdη*dηdz_ijk + dHdζ*dζdz_ijk

        Atomix.@atomic RHS[ip, ieq] += ω[i_x]*ω[i_y]*ω[i_z]*Je[ie,i_x,i_y,i_z]*(dHdz + S_ijk)*Minv[ip]
        if (ieq == 6)
            ωn = T(max(T(0), min(T(1), (Tabs[ip]-MicroConst.T00n)/(MicroConst.T0n - MicroConst.T00n))))
            Atomix.@atomic RHS[ip, ieq-1] += ω[i_x]*ω[i_y]*ω[i_z]*Je[ie,i_x,i_y,i_z]*((MicroConst.Lc + ωn*MicroConst.Lf)*dHdz)*Minv[ip]
        end
    end
    return nothing
end

function _build_rhs_diff_jacc_2D_v0!(idx, RHS_diff, rhs_diffξ_el, rhs_diffη_el, u, qe, uprimitive, x, y, connijk, dξdx, dξdy, dηdx, dηdy, Je, dψ, ω, Minv, visc_coeff, ngl, neq, PhysConst, lpert)
    ngl2 = ngl * ngl
    ie   = (idx - 1) ÷ ngl2 + 1
    rem  = (idx - 1) % ngl2
    i_x  = rem ÷ ngl + 1
    i_y  = rem % ngl + 1
    T    = eltype(RHS_diff)

    @inbounds ip   = connijk[ie, i_x, i_y]
    uip  = @view(u[ip, 1:neq])
    qeip = @view(qe[ip, 1:neq])

    @inbounds uprimitive[ie, i_x, i_y, 1:neq] .= user_primitives_gpu(uip, qeip, lpert)
    @inbounds ωJac = ω[i_x] * ω[i_y] * Je[ie, i_x, i_y]

    for ieq = 1:neq
        dqdξ = zero(T); dqdη = zero(T)

        for ii = 1:ngl
            @inbounds dqdξ += dψ[ii, i_x] * uprimitive[ie, ii, i_y, ieq]
            @inbounds dqdη += dψ[ii, i_y] * uprimitive[ie, i_x, ii, ieq]
        end

        @inbounds dξdx_kl = dξdx[ie, i_x, i_y]
        @inbounds dξdy_kl = dξdy[ie, i_x, i_y]
        @inbounds dηdx_kl = dηdx[ie, i_x, i_y]
        @inbounds dηdy_kl = dηdy[ie, i_x, i_y]

        @inbounds dqdx = visc_coeff[ieq] * (dqdξ*dξdx_kl + dqdη*dηdx_kl)
        @inbounds dqdy = visc_coeff[ieq] * (dqdξ*dξdy_kl + dqdη*dηdy_kl)

        ∇ξ∇u_kl = (dξdx_kl*dqdx + dξdy_kl*dqdy) * ωJac
        ∇η∇u_kl = (dηdx_kl*dqdx + dηdy_kl*dqdy) * ωJac

        for i = 1:ngl
            @inbounds dhdξ_ik = dψ[i, i_x]
            @inbounds dhdη_il = dψ[i, i_y]
            Atomix.@atomic rhs_diffξ_el[ie, i, i_y, ieq] -= dhdξ_ik * ∇ξ∇u_kl
            Atomix.@atomic rhs_diffη_el[ie, i_x, i, ieq] -= dhdη_il * ∇η∇u_kl
        end

        @inbounds final_rhs_diff = (rhs_diffξ_el[ie, i_x, i_y, ieq] + rhs_diffη_el[ie, i_x, i_y, ieq]) * Minv[ip]
        Atomix.@atomic RHS_diff[ip, ieq] += final_rhs_diff
    end
    return nothing
end

function _build_rhs_diff_jacc_3D_av!(idx, RHS_diff, rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el, u, qe, uprimitive, x, y, z, connijk, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, Je, dψ, ω, Minv, visc_coeff, ngl, neq, PhysConst, lpert)
    ngl3 = ngl * ngl * ngl
    ie   = (idx - 1) ÷ ngl3 + 1
    rem  = (idx - 1) % ngl3
    i_x  = rem ÷ (ngl*ngl) + 1
    rem2 = rem % (ngl*ngl)
    i_y  = rem2 ÷ ngl + 1
    i_z  = rem2 % ngl + 1
    T    = eltype(RHS_diff)

    @inbounds ip   = connijk[ie, i_x, i_y, i_z]
    @inbounds uprimitive[ie, i_x, i_y, i_z, 1:neq] .= user_primitives_gpu(@view(u[ip, 1:neq]), @view(qe[ip, 1:neq]), lpert)
    @inbounds ωJac = ω[i_x]*ω[i_y]*ω[i_z]*Je[ie, i_x, i_y, i_z]

    for ieq = 1:neq
        dqdξ = zero(T); dqdη = zero(T); dqdζ = zero(T)

        for ii = 1:ngl
            @inbounds dqdξ += dψ[ii, i_x] * uprimitive[ie, ii, i_y, i_z, ieq]
            @inbounds dqdη += dψ[ii, i_y] * uprimitive[ie, i_x, ii, i_z, ieq]
            @inbounds dqdζ += dψ[ii, i_z] * uprimitive[ie, i_x, i_y, ii, ieq]
        end

        @inbounds dξdx_klm = dξdx[ie, i_x, i_y, i_z]
        @inbounds dξdy_klm = dξdy[ie, i_x, i_y, i_z]
        @inbounds dξdz_klm = dξdz[ie, i_x, i_y, i_z]
        @inbounds dηdx_klm = dηdx[ie, i_x, i_y, i_z]
        @inbounds dηdy_klm = dηdy[ie, i_x, i_y, i_z]
        @inbounds dηdz_klm = dηdz[ie, i_x, i_y, i_z]
        @inbounds dζdx_klm = dζdx[ie, i_x, i_y, i_z]
        @inbounds dζdy_klm = dζdy[ie, i_x, i_y, i_z]
        @inbounds dζdz_klm = dζdz[ie, i_x, i_y, i_z]

        @inbounds dqdx = visc_coeff[ieq] * (dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm)
        @inbounds dqdy = visc_coeff[ieq] * (dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm)
        @inbounds dqdz = visc_coeff[ieq] * (dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm)

        ∇ξ∇u_klm = (dξdx_klm*dqdx + dξdy_klm*dqdy + dξdz_klm*dqdz) * ωJac
        ∇η∇u_klm = (dηdx_klm*dqdx + dηdy_klm*dqdy + dηdz_klm*dqdz) * ωJac
        ∇ζ∇u_klm = (dζdx_klm*dqdx + dζdy_klm*dqdy + dζdz_klm*dqdz) * ωJac

        for i = 1:ngl
            @inbounds dhdξ_ik = dψ[i, i_x]
            @inbounds dhdη_il = dψ[i, i_y]
            @inbounds dhdζ_im = dψ[i, i_z]
            Atomix.@atomic rhs_diffξ_el[ie, i, i_y, i_z, ieq] -= dhdξ_ik * ∇ξ∇u_klm
            Atomix.@atomic rhs_diffη_el[ie, i_x, i, i_z, ieq] -= dhdη_il * ∇η∇u_klm
            Atomix.@atomic rhs_diffζ_el[ie, i_x, i_y, i, ieq] -= dhdζ_im * ∇ζ∇u_klm
        end

        Atomix.@atomic RHS_diff[ip, ieq] += (rhs_diffξ_el[ie,i_x,i_y,i_z,ieq] + rhs_diffη_el[ie,i_x,i_y,i_z,ieq] + rhs_diffζ_el[ie,i_x,i_y,i_z,ieq]) * Minv[ip]
    end
    return nothing
end

function _build_rhs_diff_jacc_3D_smag!(idx, RHS_diff, rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el, u, qe, uprimitive, x, y, z, connijk, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, Je, dψ, ω, Minv, visc_coeff, ngl, neq, Δeffective_s, PhysConst, lpert)
    ngl3 = ngl * ngl * ngl
    ie   = (idx - 1) ÷ ngl3 + 1
    rem  = (idx - 1) % ngl3
    i_x  = rem ÷ (ngl*ngl) + 1
    rem2 = rem % (ngl*ngl)
    i_y  = rem2 ÷ ngl + 1
    i_z  = rem2 % ngl + 1
    T    = eltype(RHS_diff)

    @inbounds ip   = connijk[ie, i_x, i_y, i_z]
    @inbounds uprimitive[ie, i_x, i_y, i_z, 1:neq] .= user_primitives_gpu(@view(u[ip, 1:neq]), @view(qe[ip, 1:neq]), lpert)
    @inbounds ωJac = ω[i_x]*ω[i_y]*ω[i_z]*Je[ie, i_x, i_y, i_z]

    dudξ = zero(T); dudη = zero(T); dudζ = zero(T)
    dvdξ = zero(T); dvdη = zero(T); dvdζ = zero(T)
    dwdξ = zero(T); dwdη = zero(T); dwdζ = zero(T)

    for ii = 1:ngl
        @inbounds dudξ += dψ[ii, i_x] * uprimitive[ie, ii, i_y, i_z, 2]
        @inbounds dudη += dψ[ii, i_y] * uprimitive[ie, i_x, ii, i_z, 2]
        @inbounds dudζ += dψ[ii, i_z] * uprimitive[ie, i_x, i_y, ii, 2]
        @inbounds dvdξ += dψ[ii, i_x] * uprimitive[ie, ii, i_y, i_z, 3]
        @inbounds dvdη += dψ[ii, i_y] * uprimitive[ie, i_x, ii, i_z, 3]
        @inbounds dvdζ += dψ[ii, i_z] * uprimitive[ie, i_x, i_y, ii, 3]
        @inbounds dwdξ += dψ[ii, i_x] * uprimitive[ie, ii, i_y, i_z, 4]
        @inbounds dwdη += dψ[ii, i_y] * uprimitive[ie, i_x, ii, i_z, 4]
        @inbounds dwdζ += dψ[ii, i_z] * uprimitive[ie, i_x, i_y, ii, 4]
    end

    @inbounds dξdx_klm = dξdx[ie, i_x, i_y, i_z]
    @inbounds dξdy_klm = dξdy[ie, i_x, i_y, i_z]
    @inbounds dξdz_klm = dξdz[ie, i_x, i_y, i_z]
    @inbounds dηdx_klm = dηdx[ie, i_x, i_y, i_z]
    @inbounds dηdy_klm = dηdy[ie, i_x, i_y, i_z]
    @inbounds dηdz_klm = dηdz[ie, i_x, i_y, i_z]
    @inbounds dζdx_klm = dζdx[ie, i_x, i_y, i_z]
    @inbounds dζdy_klm = dζdy[ie, i_x, i_y, i_z]
    @inbounds dζdz_klm = dζdz[ie, i_x, i_y, i_z]

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
    S12 = (dudy + dvdx) * T(0.5)
    S13 = (dudz + dwdx) * T(0.5)
    S22 = dvdy
    S23 = (dvdz + dwdy) * T(0.5)
    S33 = dwdz
    Sij::T = sqrt(T(2) * (S11*S11 + T(2)*S12*S12 + T(2)*S13*S13 + S22*S22 + T(2)*S23*S23 + S33*S33))
    delta2::T = (T(2) * cbrt(Je[ie, i_x, i_y, i_z]) / (ngl-1))^2

    for ieq = 1:neq
        dqdξ = zero(T); dqdη = zero(T); dqdζ = zero(T)

        for ii = 1:ngl
            @inbounds dqdξ += dψ[ii, i_x] * uprimitive[ie, ii, i_y, i_z, ieq]
            @inbounds dqdη += dψ[ii, i_y] * uprimitive[ie, i_x, ii, i_z, ieq]
            @inbounds dqdζ += dψ[ii, i_z] * uprimitive[ie, i_x, i_y, ii, ieq]
        end

        @inbounds dqdx = visc_coeff[ieq] * Sij * delta2 * (dqdξ*dξdx_klm + dqdη*dηdx_klm + dqdζ*dζdx_klm)
        @inbounds dqdy = visc_coeff[ieq] * Sij * delta2 * (dqdξ*dξdy_klm + dqdη*dηdy_klm + dqdζ*dζdy_klm)
        @inbounds dqdz = visc_coeff[ieq] * Sij * delta2 * (dqdξ*dξdz_klm + dqdη*dηdz_klm + dqdζ*dζdz_klm)

        ∇ξ∇u_klm = (dξdx_klm*dqdx + dξdy_klm*dqdy + dξdz_klm*dqdz) * ωJac
        ∇η∇u_klm = (dηdx_klm*dqdx + dηdy_klm*dqdy + dηdz_klm*dqdz) * ωJac
        ∇ζ∇u_klm = (dζdx_klm*dqdx + dζdy_klm*dqdy + dζdz_klm*dqdz) * ωJac

        for i = 1:ngl
            @inbounds dhdξ_ik = dψ[i, i_x]
            @inbounds dhdη_il = dψ[i, i_y]
            @inbounds dhdζ_im = dψ[i, i_z]
            Atomix.@atomic rhs_diffξ_el[ie, i, i_y, i_z, ieq] -= dhdξ_ik * ∇ξ∇u_klm
            Atomix.@atomic rhs_diffη_el[ie, i_x, i, i_z, ieq] -= dhdη_il * ∇η∇u_klm
            Atomix.@atomic rhs_diffζ_el[ie, i_x, i_y, i, ieq] -= dhdζ_im * ∇ζ∇u_klm
        end

        Atomix.@atomic RHS_diff[ip, ieq] += (rhs_diffξ_el[ie,i_x,i_y,i_z,ieq] + rhs_diffη_el[ie,i_x,i_y,i_z,ieq] + rhs_diffζ_el[ie,i_x,i_y,i_z,ieq]) * Minv[ip]
    end
    return nothing
end

function apply_boundary_conditions_jacc!(idx, uaux, u, qe, x, y, t, nx, ny, poin_in_bdy_edge, qbdy, ngl, neq, npoin, lpert)
    ngl1  = ngl
    iedge = (idx - 1) ÷ ngl1 + 1
    ik    = (idx - 1) % ngl1 + 1
    T     = eltype(u)

    @inbounds ip = poin_in_bdy_edge[iedge, ik]
    @inbounds qbdy[iedge, ik, 1:neq] .= T(1234567)
    @inbounds qbdy[iedge, ik, 1:neq] .= user_bc_dirichlet_gpu(@view(uaux[ip,:]), @view(qe[ip,:]), x[ip], y[ip], t, nx[iedge,ik], ny[iedge,ik], @view(qbdy[iedge,ik,:]), lpert)
    for ieq = 1:neq
        if !(qbdy[iedge,ik,ieq] == T(1234567)) && !(qbdy[iedge,ik,ieq] == uaux[ip,ieq])
            @inbounds u[(ieq-1)*npoin+ip] = qbdy[iedge, ik, ieq]
        end
    end
    return nothing
end

function apply_boundary_conditions_jacc_3D!(idx, uaux, u, qe, x, y, z, t, nx, ny, nz, poin_in_bdy_face, qbdy, ngl, neq, npoin, lpert)
    ngl2  = ngl * ngl
    iface = (idx - 1) ÷ ngl2 + 1
    rem   = (idx - 1) % ngl2
    i_x   = rem ÷ ngl + 1
    i_y   = rem % ngl + 1
    T     = eltype(u)

    @inbounds ip = poin_in_bdy_face[iface, i_x, i_y]
    @inbounds qbdy[iface, i_x, i_y, 1:neq] .= T(123456.0)
    @inbounds qbdy[iface, i_x, i_y, 1:neq] .= user_bc_dirichlet_gpu(@view(uaux[ip,:]), @view(qe[ip,:]), x[ip], y[ip], z[ip], t, nx[iface,i_x,i_y], ny[iface,i_x,i_y], nz[iface,i_x,i_y], @view(qbdy[iface,i_x,i_y,:]), lpert)
    for ieq = 1:neq
        if !(qbdy[iface,i_x,i_y,ieq] == T(123456.0)) && !(qbdy[iface,i_x,i_y,ieq] == uaux[ip,ieq])
            @inbounds u[(ieq-1)*npoin+ip] = qbdy[iface, i_x, i_y, ieq]
        end
    end
    return nothing
end

function saturation_adjustment_jacc_3D!(idx, uaux, qe, z, connijk, neqs, param_set, lpert)
    ngl = size(connijk, 2)
    ngl3 = ngl * ngl * ngl
    ie   = (idx - 1) ÷ ngl3 + 1
    rem  = (idx - 1) % ngl3
    i_x  = rem ÷ (ngl*ngl) + 1
    rem2 = rem % (ngl*ngl)
    i_y  = rem2 ÷ ngl + 1
    i_z  = rem2 % ngl + 1

    @inbounds ip = connijk[ie, i_x, i_y, i_z]
    @inbounds uaux[ip, 1:neqs] .= user_saturation_adjustment(@view(uaux[ip,:]), @view(qe[ip,:]), z[ip], param_set, lpert)
    return nothing
end

function utouaux_jacc!(ip::Int, ieq::Int, u, uaux, npoin, neq)
    idx = (ieq - 1) * npoin + ip
    @inbounds uaux[ip, ieq] = u[idx]
    return nothing
end

function uauxtou_jacc!(ip::Int, ieq::Int, u, uaux, npoin, neq)
    idx = (ieq - 1) * npoin + ip
    @inbounds u[idx] = uaux[ip, ieq]
    return nothing
end

function RHStodu_jacc!(ip::Int, ieq::Int, RHS, du, npoin, neq)
    idx = (ieq - 1) * npoin + ip
    @inbounds du[idx] = RHS[ip, ieq]
    return nothing
end

function apply_boundary_conditions_jacc_lin_solve!(idx::Int, RHS, A, poin_in_bdy_edge, npoin, n_bdy_points)
    T  = eltype(RHS)
    @inbounds ip = poin_in_bdy_edge[idx]
    for jp = 1:npoin
        @inbounds A[ip, jp] = T(0.0)
    end
    @inbounds A[ip, ip] = T(1.0)
    @inbounds RHS[ip, :] .= T(0.0)
    return nothing
end