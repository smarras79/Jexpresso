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

# Ghost ("+" side) trace of a physical boundary face, from the case's
# Dirichlet routine.
#
# Under DG a boundary condition is imposed WEAKLY, through the interface
# flux, not by overwriting nodal values (the CG route: see the DiscGal
# no-op in BCs.jl). The exterior state is built by reflecting the interior
# trace about what user_bc_dirichlet! prescribes:
#
#     q⁺ = 2·q_bc − q⁻
#
# so the arithmetic mean of the two traces is exactly q_bc, and the Rusanov
# flux sees the prescribed state at the wall. For a free-slip wall, where
# the case zeroes the normal momentum, this is the textbook mirror state
# (normal momentum negated, everything else kept) — no case has to spell
# the mirror out a second time.
#
# user_bc_dirichlet! writes only the components it wants imposed, leaving the
# rest at the sentinel it was pre-filled with; those are copied from the
# interior trace, which leaves that characteristic free to leave the domain.
const _DG_BC_SENTINEL = 4325789.0

function dg_boundary_ghost!(qR, qL, qe_ip, coords_ip, t, tag, nx, ny, SVT, neqs)
    fill!(qR, _DG_BC_SENTINEL)
    user_bc_dirichlet!(qL, coords_ip, t, tag, qR, nx, ny, qe_ip, SVT)
    @inbounds for ieq = 1:neqs
        qR[ieq] = AlmostEqual(qR[ieq], _DG_BC_SENTINEL) ? qL[ieq] : 2.0*qR[ieq] - qL[ieq]
    end
end

# Strong-form nodal-DG interface term, 1D.
# Volume kernel has already put  -ω[i]*dFdξ + ω[i]*S  into rhs_el. Add the
# boundary correction ±(F_int − F*) at the endpoints; the 1/(ω_i·Je) lift is
# supplied by the later divide_by_mass_matrix! (M_ii = Je·ω_i) — so NO ω/Je here.
#
# NOTE: a non-periodic 1D DG run gets no term at the two domain ends — the
# 1D face list is interior-only. The 2D method below carries the physical
# boundary faces; the 1D equivalent is not implemented.
function surface_rhs_el!(params, uaux, connijk, qe, mesh, time,
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
function surface_rhs_el!(params, uaux, connijk, qe, mesh, time,
                         nelem, ngl, neqs, CL, SVT, nflux, SD::NSD_2D)

    ω = params.ω   # 1D LGL weights (length ngl); face quadrature on [-1,1]

    # neqs is tiny; these allocate — same hoist-to-scratch note as the 1D method.
    FL  = zeros(neqs); GL = zeros(neqs); FR = zeros(neqs); GR = zeros(neqs)
    FnL = zeros(neqs); FnR = zeros(neqs); Fstar = zeros(neqs)
    qB  = zeros(neqs)   # ghost trace of a physical boundary face

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

    #--------------------------------------------------------------------------
    # Physical boundary faces (mesh.dg_bfac_*): one trace from the element,
    # the other from the case's boundary routine via dg_boundary_ghost!. From
    # there the face is an ordinary one — same normal flux, same numerical
    # flux, same +ω_k·J_face·(Fn⁻ − Fn*) scatter as the L side above, and
    # nothing to scatter on the R side because there is no R element.
    # Empty on a fully periodic mesh, so this loop costs nothing there.
    #--------------------------------------------------------------------------
    nbfac = length(mesh.dg_bfac_e)
    for f = 1:nbfac
        e   = mesh.dg_bfac_e[f];  lf = mesh.dg_bfac_lf[f]
        nx  = mesh.dg_bfac_nx[f]; ny = mesh.dg_bfac_ny[f]
        Jf  = mesh.dg_bfac_Jf[f]
        tag = mesh.dg_bfac_tag[f]

        for k = 1:ngl
            i, j = face_ij(lf, k)
            ip   = connijk[e, i, j]

            qL = @view uaux[ip, :]
            dg_boundary_ghost!(qB, qL, @view(qe[ip,:]), @view(mesh.coords[:, ip]),
                               time, tag, nx, ny, SVT, neqs)

            user_flux!(FL, GL, SD, qL, @view(qe[ip,:]), mesh, CL, SVT; neqs=neqs, ip=ip)
            user_flux!(FR, GR, SD, qB, @view(qe[ip,:]), mesh, CL, SVT; neqs=neqs, ip=ip)

            @inbounds for ieq = 1:neqs
                FnL[ieq] = FL[ieq]*nx + GL[ieq]*ny
                FnR[ieq] = FR[ieq]*nx + GR[ieq]*ny
            end

            λ = max(user_max_wave_speed(qL, @view(qe[ip,:]), SD, SVT; nx=nx, ny=ny, neqs=neqs),
                    user_max_wave_speed(qB, @view(qe[ip,:]), SD, SVT; nx=nx, ny=ny, neqs=neqs))

            numerical_flux!(Fstar, FnL, FnR, qL, qB, λ, neqs, nflux)

            @inbounds for ieq = 1:neqs
                params.rhs_el[e, i, j, ieq] += ω[k] * Jf * (FnL[ieq] - Fstar[ieq])
            end
        end
    end
end
