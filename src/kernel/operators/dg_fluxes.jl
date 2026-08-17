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
# multi-speed nonlinear system needs directional upwinding (Phase 4 / SWE).
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
end