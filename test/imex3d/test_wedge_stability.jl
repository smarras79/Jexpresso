#=============================================================================
 test_wedge_stability.jl -- the tableau analysis behind IMEX_ARK(:ARS343),
 recomputed rather than quoted.

     julia --project=<env> test/imex3d/test_wedge_stability.jl

 Needs nothing but LinearAlgebra: no mesh, no MPI, no OrdinaryDiffEq. It is
 pure arithmetic on the Butcher tableaux.

 WHAT IS BEING PROTECTED
 -----------------------
 Two files in this repository make OPPOSITE recommendations about the same
 tableau, and both are right:

     hevi/ark.jl   "ARS343 ... is unusable here"        (for HEVI)
     hevi/imex3d.jl "ARS343, and the reason is ..."     (for the 3D IMEX)

 A reader who finds only one of those will reasonably conclude the other is a
 mistake and "fix" it. The numbers below are the whole argument, and they are
 asserted here so that changing either default without redoing the analysis
 fails a test rather than passing review.

 THE DIFFERENCE IN ONE SENTENCE. HEVI splits the acoustics between the two
 halves, so the explicit half's zE (horizontal sound) and the implicit half's
 zI (vertical sound) are set by INDEPENDENT wavenumbers and the whole
 rectangle is reachable. The 3D IMEX puts all the acoustics implicit, so zE is
 advection and zI is sound AT THE SAME WAVENUMBER -- the reachable set is a
 wedge of opening |v|/c, and ARS343's unstable region (zI ~ 1 with zE ~ 1)
 lies outside it for any subsonic flow.
=============================================================================#

using Test, LinearAlgebra, Printf

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")

# ark.jl is the only file needed, and it names OrdinaryDiffEq types at the
# bottom. Load just the tableau arithmetic by evaluating the file up to the
# algorithm definition -- cheaper and more robust than pulling in the whole
# ODE stack for what is a page of linear algebra.
const ARK_SRC = let s = read(joinpath(SRC, "ark.jl"), String)
    i = findfirst("#-----------------------------------------------------------------------------\n# The OrdinaryDiffEq algorithm", s)
    i === nothing ? s : s[1:first(i)-1]
end
module ARK
    using LinearAlgebra
    # ark.jl opens by resolving OrdinaryDiffEq's internals module. Nothing
    # above the cut uses the result, so an empty stand-in satisfies the lookup
    # without pulling in the ODE stack.
    module OrdinaryDiffEq end
end
Base.include_string(ARK, ARK_SRC, "ark.jl")

using .ARK: ark_tableau, ark_wedge_amplification, ark_wedge_dt_max,
            ark_joint_amplification, ark_imaginary_radius

@testset "IMEX3D wedge stability" begin

#--- 1. the rectangle: why HEVI refuses ARS343 -------------------------------
# On the full rectangle ARS343 amplifies as soon as zE leaves the origin,
# however small the amplification is. This is the measurement that made ARS232
# the HEVI default, and it has to keep coming out this way or the comment in
# ark.jl is wrong.
let tab = ark_tableau(:ARS343)
    r = ark_joint_amplification(tab, 0.5, 20.0; n = 201)
    @test r > 1 + 1e-6
    @printf("  ARS343 on the rectangle zE<=0.5, zI<=20: max|R| = %.6f (> 1, as HEVI found)\n", r)
end
for name in (:ARS232, :ARS443)
    tab = ark_tableau(name)
    r = ark_joint_amplification(tab, 0.5, 20.0; n = 201)
    @test r <= 1 + 1e-6
    @printf("  %-7s on the same rectangle:                max|R| = %.6f\n", String(name), r)
end

#--- 2. the wedge: why the 3D IMEX wants ARS343 ------------------------------
# Neutral out to zE = 2.5 at every Mach number a subsonic atmospheric flow can
# present, and -- the point -- INSENSITIVE to which one, because its unstable
# region is nowhere near the wedge.
for mach in (0.05, 0.1, 0.2, 0.3)
    tab = ark_tableau(:ARS343)
    a, zE, zI = ark_wedge_amplification(tab, 2.5, 2.5/mach*1.5, mach; n = 181)
    @test a <= 1 + 1e-6
    @printf("  ARS343 on the wedge, mach %.2f, zE<=2.5: max|R| = %.6f at (%.2f, %.2f)\n",
            mach, a, zE, zI)
end

#--- 3. and it is the BEST of the family on the wedge ------------------------
# The reverse of the HEVI ranking. Reported as the largest zE each tableau
# stays neutral to, since that is what converts directly into Δt.
function wedge_zE_limit(name; mach = 0.1, hi = 3.0, tol = 1e-3)
    tab = ark_tableau(name)
    ok(z) = ark_wedge_amplification(tab, z, z/mach*1.5, mach; n = 141)[1] <= 1 + 1e-6
    lo = 0.0
    ok(tol) || return 0.0
    while hi - lo > tol
        m = 0.5*(lo + hi)
        ok(m) ? (lo = m) : (hi = m)
    end
    return lo
end
lim = Dict(n => wedge_zE_limit(n) for n in (:ARS232, :ARS343, :ARS443))
for n in (:ARS232, :ARS343, :ARS443)
    tab = ark_tableau(n)
    nrhs = tab.stiffly_accurate ? tab.s - 1 : tab.s
    @printf("  %-7s  wedge zE limit %.2f   over %d RHS/step  =  %.3f per RHS\n",
            String(n), lim[n], nrhs, lim[n]/nrhs)
end
@test lim[:ARS343] > lim[:ARS443] > lim[:ARS232]

#--- 4. ark_wedge_dt_max is consistent with ark_wedge_amplification ----------
# The bisection has a bracketing phase that used to be a fixed cap in the
# joint version and silently returned the cap; check the answer is a real
# limit by evaluating the amplification just below and just above it.
let tab = ark_tableau(:ARS343), re = 1.5, ri = 18.0, mach = 0.06
    dt = ark_wedge_dt_max(tab, re, ri, mach)
    @test 0.0 < dt < 1e8
    below = ark_wedge_amplification(tab, 0.98dt*re, 0.98dt*ri, mach; n = 181)[1]
    above = ark_wedge_amplification(tab, 1.20dt*re, 1.20dt*ri, mach; n = 181)[1]
    @test below <= 1 + 1e-6
    @test above > 1 + 1e-6
    @printf("  ark_wedge_dt_max = %.4f s at %.1f/%.1f 1/s: |R| = %.6f below, %.6f 20%% above\n",
            dt, re, ri, below, above)
end

#--- 5. the explicit imaginary radii, unchanged ------------------------------
# ARS343's radius is the largest in the family. That is the number that reads
# "best" and is the WRONG one for HEVI -- but it is the right one for the
# zI = 0 edge of the wedge, which is where the advected non-acoustic modes sit.
@test isapprox(ark_imaginary_radius(ark_tableau(:ARS343)), 2.828; atol = 1e-2)
@test isapprox(ark_imaginary_radius(ark_tableau(:ARS232)), 1.732; atol = 1e-2)
@test isapprox(ark_imaginary_radius(ark_tableau(:ARS443)), 1.570; atol = 1e-2)

end # testset

println("\nIMEX3D wedge stability: done")
