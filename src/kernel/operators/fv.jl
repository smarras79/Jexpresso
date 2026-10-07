# Finite volumes as order-zero DG: one cell average per element, |K| dq̄/dt = −Σ_f ∫_f F*·n.
# Reuses the DG layout (:nop => 2), face lists, ghost boundaries and numerical_flux!.
# Work arrays, allocated once (concretely typed: no allocation per RHS call)
const _FV_WORK = Ref((zeros(0, 0), [zeros(0) for _ = 1:8]))
function _fv_work(nelem::Int, neqs::Int)
    R, buf = _FV_WORK[]
    if size(R) != (nelem, neqs)
        _FV_WORK[] = (zeros(nelem, neqs), [zeros(neqs) for _ = 1:8])
    end
    return _FV_WORK[]
end

# rhs_el ← each cell's face-flux balance, spread over its nodes ∝ lumped mass,
# so the DG scatter and mass division give dq̄_K/dt at every node of K.
function fv_surface_rhs!(params, uaux, connijk, qe, mesh, time,
                         nelem, ngl, neqs, CL, SVT, nflux, SD::NSD_2D)
    ω    = params.ω
    Minv = params.Minv
    R, buf = _fv_work(nelem, neqs)         # (nelem, neqs) cell residuals; 8 face vectors
    fill!(R, 0.0)
    FL  = buf[1];  GL  = buf[2];  FR  = buf[3];  GR = buf[4]
    FnL = buf[5];  FnR = buf[6];  Fs  = buf[7];  qB = buf[8]

    @inline face_ij(lfid, k) = lfid == 1 ? (1, k)   :
                               lfid == 2 ? (ngl, k) :
                               lfid == 3 ? (k, 1)   : (k, ngl)

    @inbounds for f = 1:length(mesh.dg_face_eL)
        eL  = mesh.dg_face_eL[f];  eR  = mesh.dg_face_eR[f]
        lfL = mesh.dg_face_lfL[f]; lfR = mesh.dg_face_lfR[f]
        rev = mesh.dg_face_revR[f]
        nx  = mesh.dg_face_nx[f];  ny  = mesh.dg_face_ny[f]
        Jf  = mesh.dg_face_Jf[f]
        for k = 1:ngl
            kR = rev ? ngl - k + 1 : k
            iL, jL = face_ij(lfL, k);  iR, jR = face_ij(lfR, kR)
            ipL = connijk[eL, iL, jL];  ipR = connijk[eR, iR, jR]
            qL = @view uaux[ipL, :];    qR = @view uaux[ipR, :]
            _fv_face_flux!(Fs, FL, GL, FR, GR, FnL, FnR, qL, qR, @view(qe[ipL, :]), @view(qe[ipR, :]),
                           mesh, CL, SVT, SD, nx, ny, neqs, ipL, ipR, nflux)
            w = ω[k] * Jf
            for ieq = 1:neqs
                R[eL, ieq] -= w * Fs[ieq]
                R[eR, ieq] += w * Fs[ieq]
            end
        end
    end

    # physical boundaries: exterior trace = ghost state
    @inbounds for f = 1:length(mesh.dg_bfac_e)
        e   = mesh.dg_bfac_e[f];  lf = mesh.dg_bfac_lf[f]
        nx  = mesh.dg_bfac_nx[f]; ny = mesh.dg_bfac_ny[f]
        Jf  = mesh.dg_bfac_Jf[f]; tag = mesh.dg_bfac_tag[f]
        for k = 1:ngl
            i, j = face_ij(lf, k)
            ip   = connijk[e, i, j]
            qL   = @view uaux[ip, :]
            dg_boundary_ghost!(qB, qL, @view(qe[ip, :]), @view(mesh.coords[:, ip]),
                               time, tag, nx, ny, SVT, neqs)
            _fv_face_flux!(Fs, FL, GL, FR, GR, FnL, FnR, qL, qB, @view(qe[ip, :]), @view(qe[ip, :]),
                           mesh, CL, SVT, SD, nx, ny, neqs, ip, ip, nflux)
            w = ω[k] * Jf
            for ieq = 1:neqs
                R[e, ieq] -= w * Fs[ieq]
            end
        end
    end

    @inbounds for e = 1:nelem
        vol = zero(eltype(R))
        for j = 1:ngl, i = 1:ngl
            vol += 1 / Minv[connijk[e, i, j]]
        end
        for j = 1:ngl, i = 1:ngl
            share = (1 / Minv[connijk[e, i, j]]) / vol
            for ieq = 1:neqs
                params.rhs_el[e, i, j, ieq] = share * R[e, ieq]
            end
        end
    end
    return nothing
end

fv_surface_rhs!(params, uaux, connijk, qe, mesh, time, nelem, ngl, neqs, CL, SVT, nflux, SD) =
    error(" # :AD => FV() is implemented in 2D only (got $(typeof(SD))).")

@inline function _fv_face_flux!(Fs, FL, GL, FR, GR, FnL, FnR, qL, qR, qeL, qeR,
                                mesh, CL, SVT, SD, nx, ny, neqs, ipL, ipR, nflux)
    user_flux!(FL, GL, SD, qL, qeL, mesh, CL, SVT; neqs = neqs, ip = ipL)
    user_flux!(FR, GR, SD, qR, qeR, mesh, CL, SVT; neqs = neqs, ip = ipR)
    @inbounds for ieq = 1:neqs
        FnL[ieq] = FL[ieq]*nx + GL[ieq]*ny
        FnR[ieq] = FR[ieq]*nx + GR[ieq]*ny
    end
    λ = max(user_max_wave_speed(qL, qeL, SD, SVT; nx = nx, ny = ny, neqs = neqs),
            user_max_wave_speed(qR, qeR, SD, SVT; nx = nx, ny = ny, neqs = neqs))
    sL, sR = _face_wave_bounds(qL, qR, qeL, qeR, SD, SVT, nx, ny, neqs, nflux)
    numerical_flux!(Fs, FnL, FnR, qL, qR, λ, sL, sR, nx, ny, neqs, nflux)
    return Fs
end

# Lumped-mass L² projection of every column of q onto cell averages (P0)
function fv_cell_averages!(q::AbstractMatrix, connijk, M::AbstractVector, nelem, ngl)
    @inbounds for e = 1:nelem
        vol = zero(eltype(M))
        for j = 1:ngl, i = 1:ngl
            vol += M[connijk[e, i, j]]
        end
        for c = 1:size(q, 2)
            s = zero(eltype(q))
            for j = 1:ngl, i = 1:ngl
                ip = connijk[e, i, j]
                s += M[ip] * q[ip, c]
            end
            s /= vol
            for j = 1:ngl, i = 1:ngl
                q[connijk[e, i, j], c] = s
            end
        end
    end
    return q
end
