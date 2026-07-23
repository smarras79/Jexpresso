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