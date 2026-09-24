#
# dg_fluxes.jl — Numerical (interface) fluxes for the DG discretization (:AD => DiscGal())
#
# DG couples neighboring elements through a numerical flux evaluated at element
# interfaces from the two face traces (qL, qR). Flux types are subtypes of
# AbstractNumericalFlux (abstractTypes.jl), selected in the case deck via
# :numerical_flux => upwind_flux(), rusanov_flux(), ...
#
# Fluxes are generic over the equation set: they consume the analytic flux
# F(q) from the per-problem user_flux! and a maximum wave speed λ from the
# per-problem wave-speed hook.
#


# Rusanov (local Lax–Friedrichs), generic over the equation set.
# FL,FR = analytic fluxes f(qL),f(qR); qL,qR = interface traces; λ = max wave speed.
function numerical_flux!(Fstar, FL, FR, qL, qR, λ, neqs, ::rusanov_flux)
    @inbounds for ieq = 1:neqs
        Fstar[ieq] = 0.5*(FL[ieq] + FR[ieq]) - 0.5*λ*(qR[ieq] - qL[ieq])
    end
end

# Upwind. For a linear constant-coefficient system with a single characteristic
# speed magnitude (1D advection; the acoustic wave system), upwind == Rusanov
# with λ = that speed. Replace with a characteristic/Roe split when a genuinely
# multi-speed nonlinear system needs directional upwinding (e.g. Euler, SWE).
function numerical_flux!(Fstar, FL, FR, qL, qR, λ, neqs, ::upwind_flux)
    @inbounds for ieq = 1:neqs
        Fstar[ieq] = 0.5*(FL[ieq] + FR[ieq]) - 0.5*λ*(qR[ieq] - qL[ieq])
    end
end

# Strong-form nodal-DG interface term, 1D.
# Volume kernel has already put  -ω[i]*dFdξ + ω[i]*S  into rhs_el. Add the
# boundary correction ±(F_int − F*) at the endpoints; the 1/(ω_i·Je) lift is
# supplied by the later divide_by_mass_matrix! (M_ii = Je·ω_i) — so NO ω/Je here.
function surface_rhs_el!(params, uaux, connijk, qe, mesh,
                         nelem, ngl, neqs, CL, SVT, nflux, SD::NSD_1D)

    lperiodic = params.inputs[:lperiodic_1d]

    # neqs is tiny; these allocate — hoist to params scratch (or reuse two rows
    # of params.F, free after the volume loop) to make this allocation-free.
    FL    = zeros(neqs); FR = zeros(neqs); Fstar = zeros(neqs); Gdum = zeros(neqs)

    nfaces = lperiodic ? nelem : nelem - 1
    for f = 1:nfaces
        eL = f
        eR = (f == nelem) ? 1 : f + 1          # only hits nelem→1 when periodic

        ipL = connijk[eL, ngl, 1]              # eL right-face trace ("−" side)
        ipR = connijk[eR,   1, 1]              # eR left-face  trace ("+" side)

        qL = @view uaux[ipL, :]
        qR = @view uaux[ipR, :]

        user_flux!(FL, Gdum, SD, qL, @view(qe[ipL,:]), mesh, CL, SVT; neqs=neqs, ip=ipL)
        user_flux!(FR, Gdum, SD, qR, @view(qe[ipR,:]), mesh, CL, SVT; neqs=neqs, ip=ipR)

        λ = max(user_max_wave_speed(qL, @view(qe[ipL,:]), SD, SVT; neqs=neqs),
                user_max_wave_speed(qR, @view(qe[ipR,:]), SD, SVT; neqs=neqs))

        numerical_flux!(Fstar, FL, FR, qL, qR, λ, neqs, nflux)

        @inbounds for ieq = 1:neqs
            params.rhs_el[eL, ngl, ieq] += (FL[ieq] - Fstar[ieq])   # right face, n=+1
            params.rhs_el[eR,   1, ieq] -= (FR[ieq] - Fstar[ieq])   # left  face, n=−1
        end
    end
end

# Strong-form nodal-DG interface term, 2D — over the precomputed face list
# (mesh.dg_face_*, built by build_dg_faces_2D!; interior and periodic faces
# are indistinguishable here). Per face node: normal flux Fn = F·nx + G·ny
# on both traces, numerical_flux! in the normal direction, then scatter
#     L:  +ω_k·J_face·(Fn_L − Fn*)      R:  −ω_k·J_face·(Fn_R − Fn*)
# (n points L→R; R's outward normal is −n, and the two minus signs cancel
# into the single leading −, matching the 1D sign pattern). The face measure
# ω_k·J_face is applied HERE; only the volume mass ω_i·ω_j·J_vol is supplied
# by the later divide_by_mass_matrix! — the 1D method adds the bare flux
# difference because a 1D face is a point, and carrying that idiom to 2D
# drops the face quadrature — a silent wrong operator.
function surface_rhs_el!(params, uaux, connijk, qe, mesh,
                         nelem, ngl, neqs, CL, SVT, nflux, SD::NSD_2D)

    ω = params.ω   # 1D LGL weights (length ngl); face quadrature on [-1,1]

    # neqs is tiny; these allocate — same hoist-to-scratch note as the 1D method.
    FL  = zeros(neqs); GL = zeros(neqs); FR = zeros(neqs); GR = zeros(neqs)
    FnL = zeros(neqs); FnR = zeros(neqs); Fstar = zeros(neqs)

    # (lfid, k) → (i, j): slice convention 1=x-min, 2=x-max, 3=y-min, 4=y-max
    @inline face_ij(lfid, k) = lfid == 1 ? (1, k)   :
                               lfid == 2 ? (ngl, k) :
                               lfid == 3 ? (k, 1)   : (k, ngl)

    nfaces = length(mesh.dg_face_eL)
    for f = 1:nfaces
        eL  = mesh.dg_face_eL[f];  eR  = mesh.dg_face_eR[f]
        lfL = mesh.dg_face_lfL[f]; lfR = mesh.dg_face_lfR[f]
        rev = mesh.dg_face_revR[f]
        nx  = mesh.dg_face_nx[f];  ny  = mesh.dg_face_ny[f]
        Jf  = mesh.dg_face_Jf[f]

        for k = 1:ngl
            kR = rev ? ngl - k + 1 : k
            iL, jL = face_ij(lfL, k)
            iR, jR = face_ij(lfR, kR)
            ipL = connijk[eL, iL, jL]
            ipR = connijk[eR, iR, jR]

            qL = @view uaux[ipL, :]
            qR = @view uaux[ipR, :]

            user_flux!(FL, GL, SD, qL, @view(qe[ipL,:]), mesh, CL, SVT; neqs=neqs, ip=ipL)
            user_flux!(FR, GR, SD, qR, @view(qe[ipR,:]), mesh, CL, SVT; neqs=neqs, ip=ipR)

            @inbounds for ieq = 1:neqs
                FnL[ieq] = FL[ieq]*nx + GL[ieq]*ny
                FnR[ieq] = FR[ieq]*nx + GR[ieq]*ny
            end

            λ = max(user_max_wave_speed(qL, @view(qe[ipL,:]), SD, SVT; nx=nx, ny=ny, neqs=neqs),
                    user_max_wave_speed(qR, @view(qe[ipR,:]), SD, SVT; nx=nx, ny=ny, neqs=neqs))

            numerical_flux!(Fstar, FnL, FnR, qL, qR, λ, neqs, nflux)

            @inbounds for ieq = 1:neqs
                params.rhs_el[eL, iL, jL, ieq] += ω[k]  * Jf * (FnL[ieq] - Fstar[ieq])
                params.rhs_el[eR, iR, jR, ieq] -= ω[kR] * Jf * (FnR[ieq] - Fstar[ieq])
            end
        end
    end
    # --- 2:1 mortar faces ------------------------------------------------
    # One parent face against two child halves (mesh.dg_ncfp_*, built by
    # build_dg_faces_2D!). Per half h: evaluate the parent trace at the
    # child's nodes with interp[:,:,h], take the numerical flux there
    # against the child trace, and apply the child correction exactly as on
    # a conforming face. The parent then takes the L2 projection of the two
    # halves' fluxes, project[:,:,1]*F*_1 + project[:,:,2]*F*_2 (project
    # already carries the half-interval factor 1/2), and applies its
    # correction against its OWN trace flux, never the flux of the
    # interpolated states (F(interp*q) != interp*F(q) for a nonlinear flux).
    # The normal points child -> parent, so the parent takes the conforming
    # loop's R-side sign. Conservation: w'*project[:,:,h] = w'/2 and
    # Jfp = 2*Jfc, so the parent absorbs exactly what the children emit.
    npf = length(mesh.dg_ncfp_p)
    if npf > 0
        interp  = params.interp
        project = params.project
        ncol = size(uaux, 2)
        Qp   = zeros(eltype(uaux), ngl, ncol)   # parent trace
        Qs   = zeros(eltype(uaux), ngl, ncol)   # parent trace at one child half's nodes
        Fnp  = zeros(ngl, neqs)                 # parent's own normal flux
        Fsh  = zeros(ngl, neqs, 2)              # numerical flux on each child half
        Fsp  = zeros(ngl, neqs)                 # numerical flux projected to the parent

        for r = 1:npf
            p   = mesh.dg_ncfp_p[r]
            lfp = mesh.dg_ncfp_lfp[r]
            e1  = mesh.dg_ncfp_h1[r]
            nx  = mesh.dg_ncf_nx[e1];  ny = mesh.dg_ncf_ny[e1]
            Jfp = mesh.dg_ncf_Jfp[e1]

            # parent trace, ascending along the slice (the order interp and
            # project assume), and the parent's own normal flux on it
            for k = 1:ngl
                ik, jk = face_ij(lfp, k)
                ipp = connijk[p, ik, jk]
                @inbounds for m = 1:ncol
                    Qp[k, m] = uaux[ipp, m]
                end
                user_flux!(FR, GR, SD, @view(uaux[ipp, :]), @view(qe[ipp,:]), mesh, CL, SVT; neqs=neqs, ip=ipp)
                @inbounds for ieq = 1:neqs
                    Fnp[k, ieq] = FR[ieq]*nx + GR[ieq]*ny
                end
            end

            # the two child halves: parent trace evaluated at the child's
            # nodes, numerical flux against the child trace, child correction
            for h = 1:2
                idx = h == 1 ? mesh.dg_ncfp_h1[r] : mesh.dg_ncfp_h2[r]
                c   = mesh.dg_ncf_c[idx]
                lfc = mesh.dg_ncf_lfc[idx]
                Jfc = mesh.dg_ncf_Jfc[idx]
                @inbounds for m = 1:ncol, b = 1:ngl
                    s = zero(eltype(Qs))
                    for a = 1:ngl
                        s += interp[b, a, h] * Qp[a, m]
                    end
                    Qs[b, m] = s
                end
                for k = 1:ngl
                    ik, jk = face_ij(lfc, k)
                    ipc = connijk[c, ik, jk]
                    qC = @view uaux[ipc, :]
                    qS = @view Qs[k, :]
                    user_flux!(FL, GL, SD, qC, @view(qe[ipc,:]), mesh, CL, SVT; neqs=neqs, ip=ipc)
                    user_flux!(FR, GR, SD, qS, @view(qe[ipc,:]), mesh, CL, SVT; neqs=neqs, ip=ipc)
                    @inbounds for ieq = 1:neqs
                        FnL[ieq] = FL[ieq]*nx + GL[ieq]*ny
                        FnR[ieq] = FR[ieq]*nx + GR[ieq]*ny
                    end
                    λ = max(user_max_wave_speed(qC, @view(qe[ipc,:]), SD, SVT; nx=nx, ny=ny, neqs=neqs),
                            user_max_wave_speed(qS, @view(qe[ipc,:]), SD, SVT; nx=nx, ny=ny, neqs=neqs))
                    numerical_flux!(Fstar, FnL, FnR, qC, qS, λ, neqs, nflux)
                    @inbounds for ieq = 1:neqs
                        params.rhs_el[c, ik, jk, ieq] += ω[k] * Jfc * (FnL[ieq] - Fstar[ieq])
                        Fsh[k, ieq, h] = Fstar[ieq]
                    end
                end
            end

            # gather: L2 projection of the two halves' fluxes onto the parent
            # trace, then the parent correction against its own flux
            @inbounds for ieq = 1:neqs, a = 1:ngl
                s = zero(eltype(Fsp))
                for b = 1:ngl
                    s += project[a, b, 1] * Fsh[b, ieq, 1] + project[a, b, 2] * Fsh[b, ieq, 2]
                end
                Fsp[a, ieq] = s
            end
            for k = 1:ngl
                ik, jk = face_ij(lfp, k)
                @inbounds for ieq = 1:neqs
                    params.rhs_el[p, ik, jk, ieq] -= ω[k] * Jfp * (Fnp[k, ieq] - Fsp[k, ieq])
                end
            end
        end
    end
end
