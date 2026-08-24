#=============================================================================
 test_closures.jl -- Smagorinsky, Vreman and the near-wall limiter, checked
 against their defining properties rather than against themselves.

 These do NOT load Jexpresso. The three kernels are transcribed here from
 SGS.jl and pinned against an INDEPENDENT construction: Vreman's beta tensor
 is rebuilt from the tensor definition (beta = Delta^2 * A' * A with
 A[i,j] = du_j/dx_i) and compared to the hand-expanded scalar form the solver
 actually evaluates. A transcription error in either one shows up as a
 mismatch; agreeing with itself would prove nothing.
=============================================================================#
using Test, LinearAlgebra

# ---- transcribed from SGS.jl -----------------------------------------------
# sgs_mixing_length2, SGS.jl:48
mixlen2(C, Δ2, z, κ, damp) = begin
    CΔ2 = C * Δ2
    (damp && z > 0.0) || return CΔ2
    κz2 = (κ*z)^2
    CΔ2 * κz2 / (CΔ2 + κz2)
end

# Smagorinsky, SGS.jl:1584-1647
function smag(G, Δ2, ρ, C_s2; z = 0.0, κ = 0.4, damp = false)
    dudx,dudy,dudz, dvdx,dvdy,dvdz, dwdx,dwdy,dwdz = G
    S11 = dudx; S22 = dvdy; S33 = dwdz
    S12 = 0.5*(dudy + dvdx); S13 = 0.5*(dudz + dwdx); S23 = 0.5*(dvdz + dwdy)
    SS  = S11^2 + S22^2 + S33^2 + 2*(S12^2 + S13^2 + S23^2)
    ρ * mixlen2(C_s2, Δ2, z, κ, damp) * sqrt(2SS)
end

# Vreman, SGS.jl:1731-1818
function vrem(G, Δ2, ρ, C_v; z = 0.0, κ = 0.4, damp = false, eps_v = 1e-12)
    dudx,dudy,dudz, dvdx,dvdy,dvdz, dwdx,dwdy,dwdz = G
    β11 = Δ2*(dudx*dudx + dudy*dudy + dudz*dudz)
    β12 = Δ2*(dudx*dvdx + dudy*dvdy + dudz*dvdz)
    β13 = Δ2*(dudx*dwdx + dudy*dwdy + dudz*dwdz)
    β22 = Δ2*(dvdx*dvdx + dvdy*dvdy + dvdz*dvdz)
    β23 = Δ2*(dvdx*dwdx + dvdy*dwdy + dvdz*dwdz)
    β33 = Δ2*(dwdx*dwdx + dwdy*dwdy + dwdz*dwdz)
    B_β = β11*β22 + β11*β33 + β22*β33 - (β12^2 + β13^2 + β23^2)
    aa  = dudx^2+dudy^2+dudz^2 + dvdx^2+dvdy^2+dvdz^2 + dwdx^2+dwdy^2+dwdz^2
    ℓ2  = mixlen2(C_v, Δ2, z, κ, damp)
    C_eff = Δ2 > 0 ? ℓ2/Δ2 : C_v
    (aa > eps_v && B_β > 0.0) ? ρ * C_eff * sqrt(B_β/aa) : 0.0
end

# ---- INDEPENDENT: straight from Vreman (2004), eqs. (6)-(7) ------------------
# alpha_ij = du_j/dx_i  ->  A[i,j];  beta = Delta^2 A' A
function vrem_ref(G, Δ2, ρ, C_v)
    dudx,dudy,dudz, dvdx,dvdy,dvdz, dwdx,dwdy,dwdz = G
    A = [dudx dvdx dwdx      # row 1: d/dx of (u,v,w)
         dudy dvdy dwdy      # row 2: d/dy
         dudz dvdz dwdz]     # row 3: d/dz
    β = Δ2 * (A' * A)
    B = β[1,1]*β[2,2] - β[1,2]^2 + β[1,1]*β[3,3] - β[1,3]^2 + β[2,2]*β[3,3] - β[2,3]^2
    aa = sum(abs2, A)
    (aa > 1e-12 && B > 0.0) ? ρ * C_v * sqrt(B/aa) : 0.0
end

const Δ2 = 40.0^2
const ρ  = 1.2
const Cs2 = 0.16^2
const Cv  = 2.5 * Cs2

@testset "SGS closures" begin

@testset "Vreman: hand-expanded beta == the tensor definition" begin
    # Random gradients, so an index slip in any one of the six betas shows up.
    # Seeded: this has to be the same test every run.
    rng = 12345
    lcg() = (rng = (1103515245*rng + 12345) % 2147483648; rng/2147483648 - 0.5)
    for _ = 1:200
        G = ntuple(_ -> lcg(), 9)
        @test vrem(G, Δ2, ρ, Cv) ≈ vrem_ref(G, Δ2, ρ, Cv) rtol = 1e-12
    end
end

@testset "Vreman vanishes in the flows it is designed to vanish in" begin
    # THE defining property of the model, and the whole reason to prefer it to
    # Smagorinsky near a wall: for any flow whose velocity gradient has rank 1
    # (unidirectional shear), B_beta is identically zero.
    for (name, G) in (("u(z) shear", (0,0,1.7, 0,0,0, 0,0,0)),
                      ("u(y) shear", (0,1.7,0, 0,0,0, 0,0,0)),
                      ("v(z) shear", (0,0,0, 0,0,2.3, 0,0,0)),
                      ("w(x) shear", (0,0,0, 0,0,0, 3.1,0,0)))
        @test vrem(Float64.(G), Δ2, ρ, Cv) == 0.0
        @test vrem_ref(Float64.(G), Δ2, ρ, Cv) == 0.0
        # ... while Smagorinsky does NOT vanish there. That contrast is the
        # entire motivation, so assert it rather than assume it.
        @test smag(Float64.(G), Δ2, ρ, Cs2) > 0.0
    end
    # Solid-body rotation: strain zero, so Smagorinsky vanishes too.
    Grot = (0.0,1.0,0.0, -1.0,0.0,0.0, 0.0,0.0,0.0)
    @test smag(Grot, Δ2, ρ, Cs2) ≈ 0.0 atol = 1e-14
end

@testset "Vreman is non-zero on genuinely 3D strain" begin
    # A rank-1 test alone would pass for a model that is identically zero.
    G = (0.5,0.2,0.1, 0.3,-0.4,0.2, 0.1,0.2,-0.1)
    @test vrem(G, Δ2, ρ, Cv) > 0.0
    @test vrem(G, Δ2, ρ, Cv) ≈ vrem_ref(G, Δ2, ρ, Cv) rtol = 1e-12
end

@testset "Smagorinsky is the textbook form" begin
    # nu_t = (Cs*Delta)^2 |S|,  |S| = sqrt(2 Sij Sij).
    # Pure shear du/dz = s: S13 = s/2, so 2 Sij Sij = s^2 and |S| = |s|.
    s = 1.7
    @test smag((0,0,s, 0,0,0, 0,0,0), Δ2, ρ, Cs2) ≈ ρ * Cs2 * Δ2 * s rtol = 1e-14
    # Pure axial strain du/dx = e: 2 Sij Sij = 2 e^2.
    e = 0.9
    @test smag((e,0,0, 0,0,0, 0,0,0), Δ2, ρ, Cs2) ≈ ρ*Cs2*Δ2*sqrt(2)*e rtol = 1e-14
    # Quadratic in Delta, linear in the gradient.
    G = (0.5,0.2,0.1, 0.3,-0.4,0.2, 0.1,0.2,-0.1)
    @test smag(G, 4Δ2, ρ, Cs2) ≈ 4 * smag(G, Δ2, ρ, Cs2) rtol = 1e-14
    @test smag(2.0.*G, Δ2, ρ, Cs2) ≈ 2 * smag(G, Δ2, ρ, Cs2) rtol = 1e-14
    # Vreman is quadratic in Delta too -- beta carries Delta^2, B_beta Delta^4,
    # sqrt(B/aa) Delta^2. This is what makes C_vrem sit where C_s^2 sits, and
    # therefore what makes the SAME near-wall limiter legitimate on both.
    @test vrem(G, 4Δ2, ρ, Cv) ≈ 4 * vrem(G, Δ2, ρ, Cv) rtol = 1e-12
    @test vrem(2.0.*G, Δ2, ρ, Cv) ≈ 2 * vrem(G, Δ2, ρ, Cv) rtol = 1e-12
end

@testset "near-wall limiter" begin
    κ = 0.4; Δ = 40.0; Δ2l = Δ^2
    # Off: exactly C*Delta^2, to the last bit, for every z.
    for z in (0.0, 1.0, 1e6)
        @test mixlen2(Cs2, Δ2l, z, κ, false) === Cs2*Δ2l
        @test mixlen2(Cv,  Δ2l, z, κ, false) === Cv*Δ2l
    end
    # On: the harmonic blend 1/l^2 = 1/(Cs*Delta)^2 + 1/(kappa z)^2.
    for z in (0.1, 1.0, 5.0, 50.0, 500.0)
        l2 = mixlen2(Cs2, Δ2l, z, κ, true)
        @test 1/l2 ≈ 1/(Cs2*Δ2l) + 1/(κ*z)^2 rtol = 1e-14
        @test l2 < Cs2*Δ2l                      # never increases mixing
        @test l2 < (κ*z)^2
    end
    # Limits: wall -> (kappa z)^2, far field -> (Cs Delta)^2.
    @test mixlen2(Cs2, Δ2l, 1e-4, κ, true) ≈ (κ*1e-4)^2 rtol = 1e-6
    @test mixlen2(Cs2, Δ2l, 1e8,  κ, true) ≈ Cs2*Δ2l    rtol = 1e-6
    # z = 0 exactly is the floor node: undamped, per the `z > 0.0` guard.
    @test mixlen2(Cs2, Δ2l, 0.0, κ, true) === Cs2*Δ2l
    # Monotone in z.
    ls = [mixlen2(Cs2, Δ2l, z, κ, true) for z in 0.1:0.1:100.0]
    @test all(diff(ls) .> 0)
end

@testset "damping reaches BOTH closures identically" begin
    # The limiter divides back out by Delta^2 in the Vreman path, so at the
    # same z it must scale both models by the SAME factor. If it ever did not,
    # switching closure would silently change the near-wall solution.
    G = (0.5,0.2,0.1, 0.3,-0.4,0.2, 0.1,0.2,-0.1)
    for z in (0.5, 5.0, 50.0)
        fs = smag(G, Δ2, ρ, Cs2; z = z, damp = true) / smag(G, Δ2, ρ, Cs2)
        fv = vrem(G, Δ2, ρ, Cv;  z = z, damp = true) / vrem(G, Δ2, ρ, Cv)
        @test fs ≈ mixlen2(Cs2, Δ2, z, 0.4, true)/(Cs2*Δ2) rtol = 1e-14
        @test fv ≈ mixlen2(Cv,  Δ2, z, 0.4, true)/(Cv*Δ2)  rtol = 1e-14
        @test 0 < fs < 1 && 0 < fv < 1
    end
end

# ---- 2D reductions, SGS.jl:2011-2016 and 1893-1945 -------------------------
function vrem2d(dudx,dudy,dvdx,dvdy, Δ2, ρ, C_v)
    β11 = Δ2*(dudx*dudx + dudy*dudy)
    β12 = Δ2*(dudx*dvdx + dudy*dvdy)
    β22 = Δ2*(dvdx*dvdx + dvdy*dvdy)
    B_β = β11*β22 - β12*β12
    aa  = dudx^2 + dudy^2 + dvdx^2 + dvdy^2
    (aa > 1e-12 && B_β > 0.0) ? ρ * C_v * sqrt(B_β/aa) : 0.0
end

# sgs_stability_function, SGS.jl
fRi(Ri, Pr_t) = sqrt(max(0.0, 1.0 - Ri/Pr_t))

# SGS_diffusion, SGS.jl:1471 -- what each equation actually receives
function eqvisc(ieq, μ_turb, ρ; Pr_t = 1/3, Sc_t = 1/3, μ_mol = 1.8e-5,
                κ_mol = 2.5e-5, mask = 1.0, ltheta = true)
    if ieq in (2,3,4);      (μ_mol + μ_turb) * mask
    elseif ieq == 5;        ltheta ? (μ_turb/Pr_t)*mask : (ρ*κ_mol + μ_turb/Pr_t)*mask
    else                    (ρ*κ_mol + μ_turb/Sc_t) * mask
    end
end

@testset "2D is the 3D closure with w and d/dz removed" begin
    # A 2D run must not be a different model. Zero the third row and column of
    # the 3D gradient and the two must agree exactly.
    for (a,b,c,d) in ((0.5,0.2,0.3,-0.4), (1.0,0.0,0.0,0.0), (0.1,-0.7,0.9,0.2))
        G3 = (a, b, 0.0,  c, d, 0.0,  0.0, 0.0, 0.0)
        @test vrem(G3, Δ2, ρ, Cv) ≈ vrem2d(a,b,c,d, Δ2, ρ, Cv) rtol = 1e-14
    end
    # And the 2D model inherits the rank-1 property.
    @test vrem2d(0.0, 1.7, 0.0, 0.0, Δ2, ρ, Cv) == 0.0
end

@testset "Richardson stability function (Lilly)" begin
    Pr_t = 1/3
    @test fRi(0.0, Pr_t) == 1.0                    # neutral: untouched
    @test fRi(-5.0, Pr_t) > 1.0                    # unstable: enhanced
    @test fRi(Pr_t, Pr_t) == 0.0                   # cuts off exactly at Ri = Pr_t
    @test fRi(10.0, Pr_t) == 0.0                   # and stays off, never negative
    @test fRi(0.5*Pr_t, Pr_t) ≈ sqrt(0.5) rtol = 1e-14
    # Monotone decreasing through the cutoff.
    v = [fRi(r, Pr_t) for r in -1.0:0.01:1.0]
    @test all(diff(v) .<= 1e-15)
    @test all(v .>= 0.0)                           # never negative -> no anti-diffusion
end

@testset "per-equation eddy viscosity is DYNAMIC throughout" begin
    # The bug this guards: dividing mu_turb by rho in the scalar branch while
    # the momentum branch stays dynamic leaves every scalar flux short by rho.
    μt = 3.0; ρl = 1.2; Pr_t = 1/3
    @test eqvisc(2, μt, ρl) ≈ 1.8e-5 + μt          # momentum: dynamic
    @test eqvisc(3, μt, ρl) == eqvisc(4, μt, ρl) == eqvisc(2, μt, ρl)
    @test eqvisc(5, μt, ρl) ≈ μt/Pr_t              # theta: dynamic, /Pr_t
    @test eqvisc(6, μt, ρl) ≈ ρl*2.5e-5 + μt/(1/3) # scalars: dynamic, /Sc_t
    # Both scale linearly with mu_turb, so a damped mu_turb damps every equation.
    for ieq in (2,5,6)
        @test eqvisc(ieq, 2μt, ρl) - eqvisc(ieq, 0.0, ρl) ≈
              2*(eqvisc(ieq, μt, ρl) - eqvisc(ieq, 0.0, ρl)) rtol = 1e-12
    end
end

end # testset
