#=============================================================================
 test_dsgs3d.jl -- the 3D DynSGS coefficient, checked against the properties
 that define it rather than against itself.

 Like test_closures.jl this does NOT load Jexpresso. The kernel of
 compute_dsgs_viscosity!(::SGS_DSGS, ::NSD_3D) (kernel/physics/SGS.jl) is
 transcribed below onto a synthetic one- or two-element mesh, and each test
 asserts something the model MUST do:

   * on a MANUFACTURED residual -- a state whose BDF2 time derivative and
     whose right-hand side are chosen by hand -- the coefficient comes out at
     the value eq. (9) says it should, not merely at a self-consistent one;
   * a solution that satisfies the PDE exactly gets ZERO viscosity, which is
     the entire claim of a residual-based model;
   * the cap min(mu_max, mu_res) binds when, and only when, it should;
   * the denominator is the PERTURBATION spread, so adding a hydrostatic
     reference state to both the solution and qe leaves mu unchanged -- the
     property that makes the model usable on a stratified column at all;
   * mu scales as Delta^2 below the cap and as Delta at it;
   * the floors keep a globally uniform field at mu = 0 instead of 0/eps.

 Run:  julia --project=. test/sgs/test_dsgs3d.jl
=============================================================================#
using Test

# ---- constants matching PhysicalConst -------------------------------------
const GAMMA = 1004.0/718.0
const RAIR  = 287.0
const PREF  = 101200.0
const C0    = (RAIR^GAMMA)/PREF^(GAMMA-1.0)
const PR_T  = 0.7
const SC_T  = 0.7

# ---- transcribed from compute_dsgs_viscosity!(::SGS_DSGS, ::NSD_3D) --------
#
# `q`, `qe`, `q1`, `q2`, `rhs` are (npoin, neqs); `Minv` is (npoin,);
# `elems` is a vector of node-index vectors, one per element (this stands in
# for connijk without needing a mesh); `Δ` is the per-element filter width
# ALREADY divided by nop, i.e. what Δelem_filter[ie]/nop evaluates to.
function dsgs3d(q, q1, q2, qe, rhs, Minv, elems, Δs, Δt;
                C1 = 1.0, C2 = 0.5, lpert = false, neqs = size(q,2))

    npts  = sum(length, elems)
    nelem = length(elems)
    eps   = 1.0e-16

    avg   = zeros(neqs); scale = zeros(neqs); denom = zeros(neqs)
    for el in elems, ip in el, ieq in 1:neqs
        qp = lpert ? q[ip,ieq]            : q[ip,ieq] - qe[ip,ieq]
        qt = lpert ? q[ip,ieq]+qe[ip,ieq] : q[ip,ieq]
        avg[ieq]   += qp
        scale[ieq] += abs(qt)
    end
    avg ./= npts;  scale ./= npts

    for el in elems, ip in el, ieq in 1:neqs
        qp = lpert ? q[ip,ieq] : q[ip,ieq] - qe[ip,ieq]
        denom[ieq] = max(denom[ieq], abs(qp - avg[ieq]))
    end

    ρ̄ = max(scale[1], eps)
    θ̄ = neqs >= 5 ? max(scale[5], eps)/ρ̄ : 300.0
    p̄ = C0*(max(ρ̄*θ̄, 0.0))^GAMMA
    c̄ = sqrt(max(GAMMA*p̄/ρ̄, 0.0))
    FLOOR = 1.0e-3
    denom[1] = max(denom[1] + eps, FLOOR*ρ̄)
    for ieq in 2:min(neqs,4); denom[ieq] = max(denom[ieq] + eps, FLOOR*ρ̄*c̄); end
    neqs >= 5 && (denom[5] = max(denom[5] + eps, FLOOR*ρ̄*θ̄))
    for ieq in 6:neqs; denom[ieq] = max(denom[ieq] + eps, FLOOR*max(scale[ieq],eps)); end

    ν_el = zeros(nelem); μ_el = zeros(nelem)
    for (ie, el) in enumerate(elems)
        Δ = Δs[ie]; ratio = 0.0; uTmx = 0.0; ρ_el = 0.0
        for ip in el
            for ieq in 1:neqs
                # 3q - 4q1 + q2 written as 3(q-q1) - (q1-q2): same number,
                # eps*|dq| of rounding instead of eps*|q|. See SGS.jl.
                dq = 3*(q[ip,ieq] - q1[ip,ieq]) - (q1[ip,ieq] - q2[ip,ieq])
                R  = abs(dq/(2Δt) - Minv[ip]*rhs[ip,ieq])
                ratio = max(ratio, R/denom[ieq])
            end
            ρl = lpert ? q[ip,1] + qe[ip,1] : q[ip,1]
            ρl = ρl > 0 ? ρl : eps
            ul = q[ip,2]/ρl; vl = q[ip,3]/ρl; wl = q[ip,4]/ρl
            θl = lpert ? (q[ip,5]+qe[ip,5])/ρl : q[ip,5]/ρl
            pl = C0*(max(ρl*θl, 0.0))^GAMMA
            cl = sqrt(max(GAMMA*pl/ρl, 0.0))
            uTmx = max(uTmx, sqrt(ul^2+vl^2+wl^2) + cl)
            ρ_el += ρl
        end
        ρ_el /= length(el)
        ν = max(0.0, min(C2*Δ*uTmx, C1*Δ*Δ*ratio))
        ν_el[ie] = ν;  μ_el[ie] = ρ_el*ν
    end
    return (ν = ν_el, μ = μ_el, denom = denom, scale = scale, avg = avg)
end

# SGS_diffusion(::AbstractSGSModel, ::NSD_3D) -- what each slot receives once
# DynSGS has written mu_turb. Transcribed from SGS.jl:1471.
eqvisc(ieq, μ_turb, ρ; mask = 1.0, μ_mol = 1.8e-5, κ_mol = 2.4e-5, ltheta = true) =
    ieq in (2,3,4) ? (μ_mol + μ_turb)*mask :
    ieq == 5       ? (ltheta ? (μ_turb/PR_T)*mask : (ρ*κ_mol + μ_turb/PR_T)*mask) :
                     (ρ*κ_mol + μ_turb/SC_T)*mask

# ---- a small synthetic mesh ------------------------------------------------
# Two elements of eight nodes each, no shared nodes: enough to check that the
# residual norm is per-element while the normalising scales are global.
const NPOIN = 16
const ELEMS = [collect(1:8), collect(9:16)]

"A state at rest on a uniform background: rho = 1.2, theta = 300, u = 0."
function rest_state(; npoin = NPOIN, neqs = 5, ρ = 1.2, θ = 300.0)
    q = zeros(npoin, neqs)
    q[:,1] .= ρ
    q[:,5] .= ρ*θ
    return q
end

@testset "DynSGS 3D" begin

@testset "manufactured residual: mu is exactly C1 Delta^2 R/denom" begin
    # Put ALL the residual in one slot, with a denominator we control, and
    # check the number that comes out is the one eq. (9) predicts -- not just
    # that it is positive.
    Δt = 0.5;  Δ = 40.0
    q  = rest_state();  qe = rest_state()
    # theta perturbation of 1 K on ONE node -> ‖q' - <q'>‖ is set by it.
    q[1,5] += 1.2*1.0
    q1 = copy(q); q2 = copy(q)
    # BDF2 on (q, q1, q2) with q1 == q2 == q is (3-4+1)q/(2dt) = 0, so the
    # whole residual is -Minv*rhs and we can dial it directly.
    rhs = zeros(NPOIN, 5);  Minv = ones(NPOIN)
    # R_5 on node 3 of element 1. Deliberately small: mu_max here is ~7e3
    # m^2/s (see the cap testset), and a residual big enough to saturate it
    # would test the cap rather than eq. (9).
    R5 = 0.5
    rhs[3,5] = -R5

    r = dsgs3d(q, q1, q2, qe, rhs, Minv, ELEMS, [Δ, Δ], Δt)

    # denom[5]: perturbation is 1.2 on node 1 and 0 elsewhere; mean = 1.2/16.
    d5 = max(1.2 - 1.2/16, 1.0e-3*1.2*300.0)
    @test r.denom[5] ≈ d5 rtol = 1e-14
    @test r.ν[1] ≈ 1.0*Δ^2*R5/d5 rtol = 1e-12
    # Element 2 has no residual at all -> EXACTLY zero, not "small". This is
    # what the 3(q-q1)-(q1-q2) form of the BDF2 stencil buys: written as
    # 3q-4q1+q2 the same state gives nu ~ 6e-10 m^2/s from the rounding of
    # 3*1.2-4*1.2+1.2 against a floored denominator.
    @test r.ν[2] == 0.0
    # DYNAMIC out, kinematic in: mu = rho_bar * nu with rho_bar the ELEMENT mean.
    @test r.μ[1] ≈ 1.2*r.ν[1] rtol = 1e-14
end

@testset "an exact solution gets exactly zero viscosity" begin
    # The defining claim of a residual-based model. Give it a state whose BDF2
    # derivative cancels M^-1*rhs node for node.
    Δt = 0.25
    q  = rest_state();  qe = rest_state()
    q[:,2] .= 1.2 .* range(-1, 1, length = NPOIN)     # some structure, so the
    q[:,5] .+= 1.2 .* range(0, 2, length = NPOIN)     # denominators are not floors
    q1 = copy(q); q2 = copy(q)
    Minv = fill(2.0, NPOIN)
    rhs  = zeros(NPOIN, 5)                            # dq/dt = 0 and rhs = 0
    r = dsgs3d(q, q1, q2, qe, rhs, Minv, ELEMS, [40.0, 40.0], Δt)
    @test r.ν == [0.0, 0.0]

    # Now a state that is CHANGING but whose change the RHS accounts for
    # exactly: R = 0 again, even though dq/dt is large.
    q1 = q .- 0.1;  q2 = q .- 0.2                     # BDF2 -> 0.1/(2*0.25)*(3-4*1+1)... 
    dqdt = (3 .* (q .- q1) .- (q1 .- q2)) ./ (2Δt)
    rhs  = dqdt ./ Minv                               # so Minv*rhs == dqdt
    r = dsgs3d(q, q1, q2, qe, rhs, Minv, ELEMS, [40.0, 40.0], Δt)
    # Not exactly zero here and it should not be: `rhs` is BUILT by this test
    # from a division and a multiplication, so Minv*rhs differs from dqdt in
    # the last bit. What matters is that the residual is at eps*|q| and not at
    # anything the model would act on -- 1e-9 m^2/s against a molecular
    # 1.5e-5 and a turbulent O(10).
    @test all(r.ν .< 1e-9)
end

@testset "the cap binds, and only from above" begin
    Δt = 0.5;  Δ = 40.0
    q  = rest_state();  qe = rest_state()
    q[1,5] += 1.2*1.0
    q1 = copy(q); q2 = copy(q)
    Minv = ones(NPOIN)
    d5 = max(1.2 - 1.2/16, 1.0e-3*1.2*300.0)
    # (‖v‖+c)_inf is the ELEMENT maximum, and element 1 owns node 1, whose
    # theta is 301 K rather than 300. Taking c on the background instead is a
    # 0.2% error and the point of the test is that the cap is exact.
    c  = sqrt(GAMMA*C0*(1.2*301.0)^GAMMA/1.2)
    μmax = 0.5*Δ*c                                    # v = 0 here

    # Just under the cap: mu_res governs.
    rhs = zeros(NPOIN,5);  rhs[3,5] = -0.5*μmax*d5/(Δ^2)
    @test dsgs3d(q,q1,q2,qe,rhs,Minv,ELEMS,[Δ,Δ],Δt).ν[1] ≈ 0.5*μmax rtol = 1e-10
    # Far over it: capped at exactly mu_max, never above.
    rhs = zeros(NPOIN,5);  rhs[3,5] = -1.0e4*μmax*d5/(Δ^2)
    @test dsgs3d(q,q1,q2,qe,rhs,Minv,ELEMS,[Δ,Δ],Δt).ν[1] ≈ μmax rtol = 1e-10
    # The cap is enormous in a low-Mach atmosphere, which is the documented
    # reason it is inert here: c ~ 340 m/s, so mu_max ~ 7e3 m^2/s.
    @test 5.0e3 < μmax < 9.0e3
end

@testset "the normalisation sees the perturbation, not the sounding" begin
    # THE 3D-specific change, and the one that decides whether the model does
    # anything at all on a 5 km column. Take a case, then add a strong
    # hydrostatic stratification to BOTH q and qe. The perturbation is
    # untouched, so mu must be untouched -- while a denominator built on the
    # total state would collapse by two orders of magnitude.
    Δt = 0.5;  Δ = 40.0
    q  = rest_state();  qe = rest_state()
    q[1,5] += 1.2*1.0
    q1 = copy(q); q2 = copy(q)
    rhs = zeros(NPOIN,5); rhs[3,5] = -7.0; Minv = ones(NPOIN)
    base = dsgs3d(q, q1, q2, qe, rhs, Minv, ELEMS, [Δ,Δ], Δt)

    # 5 km of sounding: rho 1.2 -> 0.7, theta 300 -> 320.
    ρs = collect(range(1.2, 0.7, length = NPOIN))
    θs = collect(range(300.0, 320.0, length = NPOIN))
    qs = zeros(NPOIN,5); qs[:,1] .= ρs .- 1.2; qs[:,5] .= ρs.*θs .- 1.2*300.0
    q2s  = q  .+ qs;  qes = qe .+ qs
    q1s  = q1 .+ qs;  q2s2 = q2 .+ qs
    strat = dsgs3d(q2s, q1s, q2s2, qes, rhs, Minv, ELEMS, [Δ,Δ], Δt)

    # The theta denominator is unchanged: it is a property of q - qe.
    @test strat.denom[5] ≈ base.denom[5] rtol = 1e-10
    # ... and so, to within the change in the element-mean density, is mu.
    @test strat.ν[1] ≈ base.ν[1] rtol = 1e-10

    # Now the counterfactual: the SAME state normalised on the total field.
    # This is what the 2D path does, and it is why it cannot be used here.
    tot_denom = maximum(abs.(q2s[:,5] .- sum(q2s[:,5])/NPOIN))
    @test tot_denom > 50 * base.denom[5]
end

@testset "scaling in the filter width" begin
    Δt = 0.5
    q  = rest_state();  qe = rest_state()
    q[1,5] += 1.2*1.0
    q1 = copy(q); q2 = copy(q)
    # Small enough that mu_res governs at every width tested -- at the cap the
    # scaling is LINEAR in Delta, not quadratic, and the test would be
    # measuring the wrong branch.
    rhs = zeros(NPOIN,5); rhs[3,5] = -0.1; Minv = ones(NPOIN)
    ν(Δ) = dsgs3d(q,q1,q2,qe,rhs,Minv,ELEMS,[Δ,Δ],Δt).ν[1]
    # Below the cap mu_res = C1 Delta^2 * (...): quadratic.
    @test ν(80.0) ≈ 4*ν(40.0) rtol = 1e-12
    @test ν(20.0) ≈ ν(40.0)/4 rtol = 1e-12
    # C1 is a plain multiplier on that branch.
    r2 = dsgs3d(q,q1,q2,qe,rhs,Minv,ELEMS,[40.0,40.0],Δt; C1 = 2.0)
    @test r2.ν[1] ≈ 2*ν(40.0) rtol = 1e-12
    # At the cap the width scaling is linear instead, since mu_max = C2 Delta
    # (‖v‖+c). Assert the branch change rather than leaving it implicit.
    big = zeros(NPOIN,5); big[3,5] = -1.0e6
    νc(Δ) = dsgs3d(q,q1,q2,qe,big,Minv,ELEMS,[Δ,Δ],Δt).ν[1]
    @test νc(80.0) ≈ 2*νc(40.0) rtol = 1e-12
end

@testset "floors: a uniform field gives 0, not 0/eps" begin
    # Every field of this problem starts horizontally uniform, so every
    # ‖q' - <q'>‖ is exactly zero at t = 0. Without a floor the ratio is
    # garbage; with one it is a clean zero.
    Δt = 0.5
    q  = rest_state();  qe = rest_state()          # q == qe: perturbation is 0
    q1 = copy(q); q2 = copy(q)
    rhs = zeros(NPOIN,5); Minv = ones(NPOIN)
    r = dsgs3d(q,q1,q2,qe,rhs,Minv,ELEMS,[40.0,40.0],Δt)
    @test all(isfinite, r.denom)
    @test r.ν == [0.0, 0.0]
    @test all(r.denom .> 0)
    # The floors are the ones documented: 1e-3 of the natural scale of each.
    ρ̄ = 1.2; θ̄ = 300.0
    c̄ = sqrt(GAMMA*C0*(ρ̄*θ̄)^GAMMA/ρ̄)
    @test r.denom[1] ≈ 1e-3*ρ̄        rtol = 1e-12
    @test r.denom[2] ≈ 1e-3*ρ̄*c̄      rtol = 1e-12
    @test r.denom[5] ≈ 1e-3*ρ̄*θ̄      rtol = 1e-12
    # And a residual against a floored denominator is finite and small, not Inf.
    rhs[3,2] = -1.0e-6
    r = dsgs3d(q,q1,q2,qe,rhs,Minv,ELEMS,[40.0,40.0],Δt)
    @test isfinite(r.ν[1]) && r.ν[1] > 0
end

@testset "the residual max runs over ALL equations" begin
    # A residual in any single slot must be able to drive the coefficient;
    # if one were dropped, the model would be blind to that equation.
    Δt = 0.5
    q  = rest_state();  qe = rest_state()
    q[1,1] += 1.0e-2;  q[1,2] += 1.0;  q[1,5] += 1.2
    q1 = copy(q); q2 = copy(q); Minv = ones(NPOIN)
    for ieq in 1:5
        rhs = zeros(NPOIN,5);  rhs[3,ieq] = -1.0e3
        @test dsgs3d(q,q1,q2,qe,rhs,Minv,ELEMS,[40.0,40.0],Δt).ν[1] > 0
    end
end

@testset "the per-equation split is the standard LES one" begin
    # DynSGS supplies mu_t; the split across equations is SGS_diffusion's, the
    # same one Smagorinsky goes through. That is the whole reason the 3D path
    # writes into the closure struct instead of into visc_coeff_dsgs.
    μt = 12.0;  ρ = 1.2
    @test eqvisc(2, μt, ρ) ≈ 1.8e-5 + μt          # momentum: dynamic, + molecular
    @test eqvisc(3, μt, ρ) == eqvisc(2, μt, ρ)    # isotropic across u, v, w
    @test eqvisc(4, μt, ρ) == eqvisc(2, μt, ρ)
    @test eqvisc(5, μt, ρ) ≈ μt/PR_T              # theta: kappa_t = mu_t/Pr_t
    # The mask switches a slot off completely -- which is how :mu[1] = 0 keeps
    # the scheme strictly conservative in mass.
    @test eqvisc(1, μt, ρ; mask = 0.0) == 0.0
    # :mu[5] = 2 (the Smagorinsky deck's value) really is 2.857x momentum, the
    # factor that makes the theta equation set the parabolic step limit.
    @test eqvisc(5, μt, ρ; mask = 2.0) / eqvisc(2, μt, ρ) ≈ 2/PR_T rtol = 1e-3
end

end
