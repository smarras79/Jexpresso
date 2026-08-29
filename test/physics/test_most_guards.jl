#=============================================================================
 test_most_guards.jl -- the MOST surface layer: rho in L, the sign of L near
                        neutral, the wind floor, the zeta bounds, and the
                        prescribed-flux closure.

 WHAT WAS WRONG, IN THE ORDER THE TESTS COME
 -------------------------------------------
 1. obukhov_length had no rho. Q_H arrives as a DYNAMIC flux in W/m^2
    (-rho*cp*u_star*theta_star), so the rho*cp that made it dynamic has to be
    in the numerator too. Without it the function returned L/rho: |L| about
    20% small at sea level, every zeta = z/L about 20% large, and the error
    growing with altitude as rho falls.

 2. The near-neutral guard returned +1e6 for EITHER sign of Q_H -- a STABLE
    surface layer reported in the middle of a convective one. It was reached
    whenever u_star got small, which in an LES surface layer happens at some
    node on essentially every step, and it put a sign discontinuity in the
    middle of the map a 20-step undamped Picard iteration was iterating.

 3. Nothing floored the wind speed and nothing bounded zeta. |u| -> 0 at a
    convergence line between thermals sends L ~ u*^2 -> 0 and zeta -> -inf,
    and the drag that comes back is whatever the 20th iterate happened to be.

 4. With :user_heatflux set, BCs.jl discards MOST's own w'theta' (delta_hf = 1)
    and applies rho*user_heatflux -- but theta*, and through it L, psi_m and
    the DRAG, were still built from a free surface node that nothing pins. The
    momentum flux carried a stability correction diagnosed from a heat flux the
    model does not apply.

     julia --project=<env> test/physics/test_most_guards.jl
=============================================================================#

using Test, Printf

const _MOST_SRC = read(joinpath(@__DIR__, "..", "..", "src", "kernel", "physics", "CM_MOST.jl"), String)

# Lift the real definitions rather than restating them -- a test carrying its
# own copy pins the copy, and the whole point here is what the shipped code does.
function _lift_fn!(name)
    for needle in ("@inline function $name(", "function $name(")
        i = findfirst(needle, _MOST_SRC)
        i === nothing && continue
        j = findnext("\nend\n", _MOST_SRC, i[1])
        j === nothing && error("no closing end for $name in CM_MOST.jl")
        include_string(Main, _MOST_SRC[i[1]:j[1]+4], "CM_MOST.jl:$name")
        return
    end
    error("$name not found in CM_MOST.jl -- renamed?")
end

# The Businger-Dyer and guard-rail constants, as one contiguous block.
let i = findfirst("const a_m", _MOST_SRC), j = findfirst("const MOST_L_BIG", _MOST_SRC)
    (i === nothing || j === nothing) && error("MOST constant block not found -- restructured?")
    k = findnext('\n', _MOST_SRC, j[1])
    include_string(Main, _MOST_SRC[i[1]:k], "CM_MOST.jl:constants")
end
for f in ("psi_m", "psi_h", "most_clamp_invL", "most_inv_L",
          "obukhov_length", "_most_iterate", "_most_stress!")
    _lift_fn!(f)
end

Base.@kwdef struct TestPhysConst
    karman::Float64 = 0.4
    cp::Float64     = 1004.0
    g::Float64      = 9.81
end
const PC = TestPhysConst()

say(s) = println(s)

# LESICP2-64x64x60 surface layer.
const Z_REF = 6.91      # first LGL node above the wall
const Z0_M  = 0.1
const Z0_H  = 0.01
const TH    = 300.0     # reference potential temperature [K]
const RHO   = 1.2       # TOTAL density -- under PERT that is uaux[.,1] + qe[.,1]
const WTH   = 0.12      # :user_heatflux, K m/s

say("\n=== 1. obukhov_length carries rho ===")

@testset "obukhov_length: rho" begin
    ustar, tstar = 0.3, -0.05
    Q_H = -RHO * PC.cp * ustar * tstar                       # +18.1 W/m^2, upward

    L = obukhov_length(ustar, TH, Q_H, PC, RHO)
    # The textbook value, written out independently of the implementation.
    L_ref = -RHO * PC.cp * ustar^3 * TH / (PC.karman * PC.g * Q_H)
    @test isapprox(L, L_ref; rtol = 1e-12)

    # It is LINEAR in rho -- that is the whole content of the bug. The old code
    # returned L/rho, i.e. exactly this value at rho = 1.
    @test isapprox(obukhov_length(ustar, TH, Q_H, PC, 2*RHO),
                   2*obukhov_length(ustar, TH, Q_H, PC, RHO); rtol = 1e-12)
    # ... and Q_H itself scales with rho, so the physical L is rho-independent:
    # the two rho's cancel, which is the check that the formula is now
    # self-consistent rather than merely rescaled.
    L2 = obukhov_length(ustar, TH, -2*RHO*PC.cp*ustar*tstar, PC, 2*RHO)
    @test isapprox(L2, L; rtol = 1e-12)
    say(@sprintf("  L = %.2f m   (unstable, u* = %.2f, theta* = %.3f)", L, ustar, tstar))
end

say("\n=== 2. the sign of L, including near neutral ===")

@testset "obukhov_length: sign" begin
    ustar = 0.3
    # Q_H > 0 is an UPWARD, surface-heating, UNSTABLE flux -> L < 0.
    @test obukhov_length(ustar, TH,  50.0, PC, RHO) < 0.0
    @test obukhov_length(ustar, TH, -50.0, PC, RHO) > 0.0

    # THE REGRESSION. The old guard returned +1e6 -- stable -- for both of
    # these. A vanishing Q_H is near-neutral, not stable, and the sign of the
    # residual flux is the one thing that still matters there.
    @test obukhov_length(ustar, TH,  1.0e-9, PC, RHO) < 0.0
    @test obukhov_length(ustar, TH, -1.0e-9, PC, RHO) > 0.0
    @test abs(obukhov_length(ustar, TH, 1.0e-9, PC, RHO)) == MOST_L_BIG

    # u* -> 0 with an unstable flux must not report a stable layer either.
    @test obukhov_length(0.0, TH, 50.0, PC, RHO) < 0.0
    @test obukhov_length(0.0, TH, 50.0, PC, 0.0) < 0.0

    # never NaN, never Inf, never larger than the cap
    for u in (0.0, 1e-30, 0.3, 10.0), q in (-1e3, -1e-9, 0.0, 1e-9, 1e3), r in (0.0, RHO)
        L = obukhov_length(u, TH, q, PC, r)
        @test isfinite(L)
        @test abs(L) <= MOST_L_BIG
    end
    say(@sprintf("  Q_H = +1e-9 W/m^2 -> L = %+.3g m  (old code: +1e6)",
                 obukhov_length(0.3, TH, 1.0e-9, PC, RHO)))
end

@testset "most_inv_L is the reciprocal of obukhov_length" begin
    # The iteration works in 1/L; the two have to agree wherever both are in
    # range, or the reported L and the applied correction are different things.
    for (u, t) in ((0.3, -0.05), (0.5, 0.02), (0.05, -0.3))
        Q_H  = -RHO * PC.cp * u * t
        L    = obukhov_length(u, TH, Q_H, PC, RHO)
        invL = most_inv_L(u, TH, Q_H, PC, RHO)
        abs(L) < MOST_L_BIG || continue
        @test isapprox(invL, 1/L; rtol = 1e-10)
    end
    # and the degenerate inputs answer with the NEUTRAL 0, not with an Inf
    @test most_inv_L(0.0, TH, 50.0, PC, RHO) == 0.0
    @test most_inv_L(0.3, TH, 50.0, PC, 0.0) == 0.0
    @test most_inv_L(0.3, 0.0, 50.0, PC, RHO) == 0.0
end

say("\n=== 3. the wind floor and the zeta bounds ===")

@testset "zeta stays inside its bounds" begin
    # Sweep the whole range an LES surface layer visits, including the calm
    # nodes that used to send zeta to -inf.
    worst = 0.0
    for u in (0.0, 1e-6, 0.01, 0.1, 0.5, 1.0, 5.0, 15.0, 50.0)
        for dθ in (-5.0, -0.5, -0.01, 0.0, 0.01, 0.5, 5.0)
            us, ts, invL, _, _ =
                _most_iterate(u, TH, Z_REF, TH - dθ, Z0_M, Z0_H, PC, RHO, NaN)
            ζ = Z_REF * invL
            @test isfinite(us) && isfinite(ts) && isfinite(invL)
            @test MOST_ZETA_MIN[] - 1e-12 <= ζ <= MOST_ZETA_MAX[] + 1e-12
            @test us >= 0.0
            worst = max(worst, abs(ζ))
        end
    end
    say(@sprintf("  worst |zeta| over the sweep = %.3f  (bounds %.1f .. %.1f)",
                 worst, MOST_ZETA_MIN[], MOST_ZETA_MAX[]))
end

@testset "u* is floored, the stress is not" begin
    # The floor belongs to the DRAG COEFFICIENT. u* must not collapse to zero
    # at a calm node -- free convection over a warmer surface (theta_s above
    # the air, so this is the UNSTABLE branch) still exchanges momentum with
    # the ground ...
    us0, = _most_iterate(0.0,  TH, Z_REF, TH + 1.0, Z0_M, Z0_H, PC, RHO, NaN)
    usm, = _most_iterate(MOST_U_MIN[], TH, Z_REF, TH + 1.0, Z0_M, Z0_H, PC, RHO, NaN)
    @test us0 > 0.0
    @test isapprox(us0, usm; rtol = 1e-12)   # both see the same floored speed

    # ... but the STRESS must still go to zero with the wind, in a direction
    # that is not noise. tau_i = -rho u*^2 u_i / max(|u|, u_min) is linear in
    # u_i below the floor.
    τ = zeros(3)
    _most_stress!(τ, RHO, us0, 0.0, 0.0, 0.0, 0.0)
    @test all(τ .== 0.0)

    prev = 0.0
    for u in (1e-6, 1e-4, 1e-2, 0.05)
        _most_stress!(τ, RHO, us0, u, 0.0, 0.0, u)
        @test all(isfinite, τ)
        @test τ[1] <= 0.0                     # opposes the flow
        @test abs(τ[1]) > prev                # monotone in |u|
        prev = abs(τ[1])
    end
    # linearity below the floor: halving u halves the stress
    _most_stress!(τ, RHO, us0, 0.02, 0.0, 0.0, 0.02); τa = τ[1]
    _most_stress!(τ, RHO, us0, 0.01, 0.0, 0.0, 0.01); τb = τ[1]
    @test isapprox(τa, 2*τb; rtol = 1e-12)

    # above the floor nothing changed: the classic -rho u*^2 * u_i/|u|
    u, v = 3.0, 4.0
    _most_stress!(τ, RHO, 0.3, u, v, 0.0, 5.0)
    @test isapprox(τ[1], -RHO*0.09*(u/5.0); rtol = 1e-12)
    @test isapprox(τ[2], -RHO*0.09*(v/5.0); rtol = 1e-12)
    say(@sprintf("  u* at |u| = 0 is %.4f m/s, and tau there is exactly 0", us0))
end

say("\n=== 4. the prescribed-flux closure ===")

@testset "a prescribed w'theta' comes back unchanged" begin
    # theta* = -w'theta'/u*, so -u* theta* IS the imposed flux, at every wind
    # speed and for either sign. This is what makes the drag's stability
    # correction consistent with the heat flux BCs.jl actually applies.
    for u in (0.0, 0.2, 1.0, 5.0, 15.0), w in (WTH, -WTH, 0.0)
        us, ts, invL, _, _ = _most_iterate(u, TH, Z_REF, TH, Z0_M, Z0_H, PC, RHO, w)
        @test isapprox(-us*ts, w; atol = 1e-14, rtol = 1e-12)
        # and the sign of the stability is the sign of the flux: heating the
        # surface (w > 0) is unstable, i.e. 1/L < 0
        if w > 0
            @test invL < 0.0
        elseif w < 0
            @test invL > 0.0
        else
            @test invL == 0.0
        end
    end
    say("  w'theta' returned to machine precision at every wind speed")
end

@testset "prescribed flux vs. diagnosed: the drag actually differs" begin
    # If the two agreed there would have been nothing to fix. A surface node
    # sitting 0.5 K COOLER than the air (which a free, unpinned theta_sfc does
    # in a convective run) diagnoses a STABLE layer and cuts the drag, while
    # the imposed +0.12 K m/s says the layer is unstable.
    u = 2.0
    us_diag, = _most_iterate(u, TH, Z_REF, TH - 0.5, Z0_M, Z0_H, PC, RHO, NaN)
    us_flux, = _most_iterate(u, TH, Z_REF, TH - 0.5, Z0_M, Z0_H, PC, RHO, WTH)
    @test us_flux > us_diag
    say(@sprintf("  u* = %.4f (diagnosed, reads stable) vs %.4f (imposed flux, unstable)",
                 us_diag, us_flux))
end

say("\n=== 5. limits that must not have moved ===")

@testset "the neutral limit is still the log law" begin
    # w'theta' = 0 -> 1/L = 0 -> psi_m(0) = 0 -> u* = kappa u / ln(z/z0).
    for u in (1.0, 5.0, 15.0)
        us, ts, invL, _, _ = _most_iterate(u, TH, Z_REF, TH, Z0_M, Z0_H, PC, RHO, 0.0)
        @test invL == 0.0
        @test isapprox(us, PC.karman*u/log(Z_REF/Z0_M); rtol = 1e-12)
    end
    # theta_ref == theta_sfc in the diagnosed closure is the same statement
    us, ts, invL, _, _ = _most_iterate(5.0, TH, Z_REF, TH, Z0_M, Z0_H, PC, RHO, NaN)
    @test invL == 0.0
    @test ts == 0.0
    @test isapprox(us, PC.karman*5.0/log(Z_REF/Z0_M); rtol = 1e-12)
end

@testset "unstable raises the drag, stable lowers it" begin
    # The sign of the correction, which is what the +1e6 guard used to get
    # backwards at low wind.
    u = 3.0
    us_n, = _most_iterate(u, TH, Z_REF, TH,       Z0_M, Z0_H, PC, RHO, 0.0)
    us_u, = _most_iterate(u, TH, Z_REF, TH,       Z0_M, Z0_H, PC, RHO,  0.2)
    us_s, = _most_iterate(u, TH, Z_REF, TH,       Z0_M, Z0_H, PC, RHO, -0.2)
    @test us_u > us_n > us_s
    say(@sprintf("  u*: stable %.4f < neutral %.4f < unstable %.4f", us_s, us_n, us_u))
end

@testset "degenerate geometry returns zero fluxes, not NaN" begin
    # z_ref at or below the roughness length, a non-positive theta, a bad rho:
    # a zero flux reaches the RHS harmlessly, a NaN does not.
    for (z, th, r) in ((Z0_M, TH, RHO), (0.5*Z0_M, TH, RHO), (Z_REF, 0.0, RHO),
                       (Z_REF, -1.0, RHO), (Z_REF, TH, 0.0), (Z_REF, TH, -1.0))
        us, ts, invL, _, _ = _most_iterate(5.0, th, z, TH, Z0_M, Z0_H, PC, r, NaN)
        @test us == 0.0 && ts == 0.0 && invL == 0.0
    end
end

@testset "the iteration converges over the operating range" begin
    # Not "it returns something finite" -- it has to actually settle, or the
    # answer is the 20th iterate of an oscillation.
    nconv = 0; ntot = 0
    for u in (0.5, 1.0, 3.0, 8.0, 15.0), w in (-0.05, 0.0, 0.05, 0.12, 0.3)
        _, _, _, _, conv = _most_iterate(u, TH, Z_REF, TH, Z0_M, Z0_H, PC, RHO, w)
        ntot += 1; nconv += conv
    end
    @test nconv == ntot
    say(@sprintf("  converged on %d/%d of the operating range", nconv, ntot))
end

say("\n=== done ===\n")
