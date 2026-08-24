#=============================================================================
 schur.jl -- the scalar Schur-complement (Helmholtz) operator for the 3D
             acoustic stage system, and the preconditioner built from it.

 WHY
 ---
 The 3D IMEX stage solve is preconditioned by a VERTICAL column solve, which
 is exact for the vertical acoustic operator and does nothing at all for the
 horizontal one. Its strength is therefore proportional to the grid's acoustic
 anisotropy h_x/h_z: at 32:1 the horizontal term is a perturbation and GMRES
 eats it in ~8 iterations; at 4:1 -- the production LES mesh -- it is not a
 perturbation and the same solve needs ~54, which costs more than the step
 size it buys. Measured on that mesh: 20.3 Krylov iterations per unit of
 CFL_h = gamma*dt*c/h_x.

 Giraldo's group reports exactly this. In Abdi, Giraldo, Constantinescu, Carr,
 Wilcox & Warburton (Int. J. HPC Appl. 2019; arXiv:1702.04316), section 7:

   "NUMA can benefit from a relative speedup of 5X using IMEX time-integrators
    compared to the explicit RK35 time integrator when using the schur forms
    of IMEX; the standard forms often do not give any speedup over an explicit
    time integrator."

 Our IMEX3D is their "standard form": all five prognostic variables implicit,
 5*Np degrees of freedom, and -- as hevi_column_spectrum measures directly --
 a PURELY IMAGINARY spectrum, which is the worst case for a Krylov method.
 Their Schur form carries one scalar and, in their words, "the condition
 number of the schur form is much better than the standard form (and the
 eigenvalues are all real)".

 WHAT THIS FILE DOES, AND WHAT IT DELIBERATELY DOES NOT
 ------------------------------------------------------
 The stage system (I - lam*A) U = b, with lam = gamma*dt and the blocks as
 verified in test/hevi/test_schur_blocks.jl, is

     rho   + lam*Div[W m]                    = b_rho      (1)
     m     + lam*Grad[P] + lam*g*zhat*rho    = b_m        (2)
     Theta + lam*Div[tb .* W m]              = b_Theta    (3)

 with P = beta.*Theta. Eliminating m leaves a 2x2 system in the SCALARS
 (rho, P) -- 2*Np, not Np. It is Np in the papers because their theta equation
 is ADVECTIVE (u . grad theta0), which makes the coupling pointwise; ours is
 FLUX form, so rho survives the elimination. Dropping the two lam^2*g buoyancy
 couplings leaves the scalar Helmholtz

     H[P] = P./beta - lam^2 * Div[ tb .* W .* Grad[P] ]

 so H is used here as a PRECONDITIONER, not as an exact solve. That is the
 safe direction: GMRES converges to the same answer whatever the
 preconditioner does, so an imperfect H can cost iterations and never
 correctness. It still buys the three things that matter -- coupling in all
 three directions where the column solve has only the vertical, a scalar
 operator, and (see the tests) a real positive spectrum instead of a purely
 imaginary one.

 H IS BUILT OUT OF TWO CALLS TO hevi_apply_A!, ON PURPOSE. A bespoke scalar
 element kernel would be ~5x cheaper and would duplicate the metric handling,
 the free-slip masking, the DSS and the mass division -- every one of which is
 already verified, and every one of which is a place to introduce a bug that
 converges to a plausible wrong answer. This reference form cannot disagree
 with the operator it preconditions, and it is automatically correct in
 parallel because hevi_apply_A! is. An optimised kernel can be added later and
 checked against THIS.
=============================================================================#

"""
    schur_grad!(gx, gy, gz, P, params, op, work)

The assembled weak gradient of the scalar `P`, as the momentum rows of the
operator see it: `(gx, gy, gz) = Grad[P]`.

Read straight off block row 2, which `test_schur_blocks.jl` pins:
`A_rho_u = -Gradx[P]` and likewise for y and z, with `A_rho_w` also carrying
`-g*rho` -- so the momentum slots of `A` applied to a state whose ONLY
non-zero is `Theta = P./beta` are exactly `-Grad[P]`.
"""
function schur_grad!(gx, gy, gz, P, params, op, work)
    V, W = work.V, work.W
    fill!(V, 0.0)
    sθ = op.slot[5]
    @inbounds for ip in eachindex(P)
        V[ip, sθ] = P[ip] / op.beta[ip]        # so that beta*Theta == P
    end
    hevi_apply_A!(W, V, params, op)
    su, sv, sw = op.slot[2], op.slot[3], op.slot[4]
    @inbounds for ip in eachindex(P)
        gx[ip] = -W[ip, su]; gy[ip] = -W[ip, sv]; gz[ip] = -W[ip, sw]
    end
    return nothing
end

"""
    schur_divtb!(out, gx, gy, gz, params, op, work)

`out = Div[ thetabar .* W .* (gx,gy,gz) ]`, from block row 3: `A_Theta =
-Div[tb .* W m]`, so the Theta slot of `A` applied to a state whose only
non-zero is the momentum is exactly minus what we want. The free-slip mask `W`
is applied inside the kernel, which is why it does not appear here.
"""
function schur_divtb!(out, gx, gy, gz, params, op, work)
    V, W = work.V, work.W
    fill!(V, 0.0)
    su, sv, sw = op.slot[2], op.slot[3], op.slot[4]
    @inbounds for ip in eachindex(gx)
        V[ip, su] = gx[ip]; V[ip, sv] = gy[ip]; V[ip, sw] = gz[ip]
    end
    hevi_apply_A!(W, V, params, op)
    sθ = op.slot[5]
    @inbounds for ip in eachindex(out)
        out[ip] = -W[ip, sθ]
    end
    return nothing
end

"""
    SchurWork(npoin, nimp)

Scratch for one `schur_apply!`. Allocated once; `schur_apply!` allocates
nothing.
"""
struct SchurWork{M <: AbstractMatrix{Float64}, V <: AbstractVector{Float64}}
    V::M
    W::M
    gx::V
    gy::V
    gz::V
    d::V
end
SchurWork(npoin::Int, nimp::Int) =
    SchurWork(zeros(Float64, npoin, nimp), zeros(Float64, npoin, nimp),
              zeros(Float64, npoin), zeros(Float64, npoin),
              zeros(Float64, npoin), zeros(Float64, npoin))

"""
    schur_apply!(HP, P, params, op, lam, work)

`HP = H[P] = P./beta - lam^2 * Div[ thetabar .* W .* Grad[P] ]`.

Two applications of the full operator per call. See the header for why that is
the intended cost of the reference form.
"""
function schur_apply!(HP, P, params, op, lam::Real, work::SchurWork)
    schur_grad!(work.gx, work.gy, work.gz, P, params, op, work)
    schur_divtb!(work.d, work.gx, work.gy, work.gz, params, op, work)
    λ2 = Float64(lam)^2
    @inbounds for ip in eachindex(P)
        HP[ip] = P[ip] / op.beta[ip] - λ2 * work.d[ip]
    end
    return HP
end
