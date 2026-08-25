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

# `SchurWork` names `SchurKernel` in a field type, so schur_kernel.jl must be
# loaded before this file is. Jexpresso.jl includes them in that order; the
# standalone tests under test/hevi and test/imex3d each list their own includes
# and several predate the kernel, so they would fail at PARSE time with a bare
# `UndefVarError: SchurKernel` that says nothing about the cause. Pull it in
# here when it is not already there -- inside the package the guard is false and
# nothing is included twice.
isdefined(@__MODULE__, :SchurKernel) ||
    include(joinpath(@__DIR__, "schur_kernel.jl"))

"""
    _schur_kernel(work, op) -> SchurKernel or nothing

The kernel to use for `work` on `op`, or `nothing` to take the reference path.

THE `op.full` CHECK IS NOT DEFENSIVE PADDING. schur_kernel.jl computes the FULL
three-dimensional sweeps -- nine derivative chains for the divergence, three for
the gradient -- and has no vertical-only branch. Handed a `full = false`
operator it would happily return the full-3D answer, and the caller that builds
such an operator is `build_schur_column_precond`: the result would be a
preconditioner band for the wrong operator. That still CONVERGES, just to a
worse iteration count, so it would show up as a performance regression nobody
could explain rather than as a wrong answer. An error is loud; a silent
fallback to the reference form would hide the misuse just as well.
"""
@inline function _schur_kernel(work, op)
    work.kern === nothing && return nothing
    op.full || error("schur: the bespoke kernel is the FULL 3D form and has no ",
                     "vertical-only branch, but it was handed an operator with ",
                     "full = false. Build that SchurState without `kern`.")
    return work.kern
end

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
    kern = _schur_kernel(work, op)
    if kern !== nothing
        return schur_grad_fast!(gx, gy, gz, P, params, op, kern)
    end
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
    # The bespoke scalar sweeps (schur_kernel.jl), or `nothing` for the
    # reference form below. Both are kept, permanently: the reference is what
    # the kernel is CHECKED AGAINST, and a kernel with no independent statement
    # of what it should compute is a kernel nobody can debug. `nothing` here is
    # not a fallback for missing functionality, it selects the slow definition.
    kern::Union{Nothing, SchurKernel}
end
SchurWork(npoin::Int, nimp::Int; kern = nothing) =
    SchurWork(zeros(Float64, npoin, nimp), zeros(Float64, npoin, nimp),
              zeros(Float64, npoin), zeros(Float64, npoin),
              zeros(Float64, npoin), zeros(Float64, npoin), kern)

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

#=============================================================================
 STEP 4 -- the exact reduction, made possible by the advective Theta row.

 With A_Theta = -( tb .* Div[W m] + (W m).grad(tb) ), the stage system is

     rho   + lam*Div[W m]                        = b_rho      (1)
     m     + lam*Grad[P] + lam*g*zhat*rho        = b_m        (2)
     Theta + lam*tb.*Div[W m] + lam*(W m).grad(tb) = b_Theta  (3)

 with P = beta.*Theta. Substituting lam*Div[W m] = b_rho - rho from (1) into
 (3) removes the ONLY derivative of an unknown that appears there:

     Theta - tb.*rho = b_Theta - tb.*b_rho - lam*(W m).grad(tb)

 The left side is (rho*theta)' - tb*rho' = rhobar*theta', Giraldo's theta.
 Call it `q`; it is POINTWISE in m. And since P = beta*(q + tb*rho),

     rho = (P./beta - q) ./ tb                                 POINTWISE

 which is the step the flux form cannot take -- there, eliminating Theta
 introduces Div of the unknown momentum and rho survives into a 2x2 system.

 Putting both into (2) collects every m on the left:

     (I - lam^2*g * zhat (W grad(tb)/tb)^T) m = <known> - lam*Grad[P]
                                                        - lam*g*zhat*P./(beta.*tb)

 The operator on the left is a POINTWISE rank-one update of the identity --
 `zhat` makes it act on the w component alone -- so it inverts by
 Sherman-Morrison at each node, with no solve and no communication. That is
 the matrix Giraldo writes as A = I + lam^2 (g/theta0) (grad theta0)^T.
=============================================================================#

"""
    schur_amat_inv!(mx, my, mz, rx, ry, rz, op, lam)

Apply `A^-1` pointwise, where `A = I + lam^2*g*zhat*(W grad(tb)/tb)^T`.

From the z-momentum row with rho eliminated:

    m_z + lam*Grad_z[P] + lam*g*rho = b_w,   rho = (P./beta - q)./tb
    q   = b_Theta - tb.*b_rho - lam*(W m).grad(tb)

so `lam*g*rho` contributes `+lam^2*g*(W m).grad(tb)./tb` to the LEFT side.

`A = I + u v'` with `u = lam^2*g*zhat` and `v = W grad(tb)/tb`, so
Sherman-Morrison gives `A^-1 = I - u v' / (1 + v'u)` exactly, and `v'u` is the
scalar `lam^2*g*W_z*dtb_dz/tb`. Only the w row is touched, because `u` is
`zhat`.

THE SIGN IS THE WHOLE THING. With a minus the reduction reproduces the full
solve to 5e-07 instead of 1e-15 -- small enough to look like conditioning and
big enough to be wrong, and every block relation still checks out because the
error is downstream of them. It also inverts the physics: `1 + v'u` is
`1 + lam^2*N^2`, which for STABLE stratification is > 1 and can never vanish,
where the wrong sign made it appear singular at lam = 1/N.

For UNSTABLE stratification N^2 < 0 and `1 + lam^2*N^2` genuinely can vanish,
at lam = 1/|N| -- convective overturning reaching the step size. Checked, not
assumed: a silent division by zero would produce a finite, plausible, wrong
answer.
"""
@inline function schur_amat_inv!(mx, my, mz, rx, ry, rz, op::HEVIOperator, lam::Real)
    g  = PhysicalConst{Float64}().g
    λ2g = Float64(lam)^2 * g
    tb = op.thetabar
    wz = op.wall
    @inbounds for ip in eachindex(rx)
        vx = op.wallx[ip] ? 0.0 : op.dtbdx[ip] / tb[ip]
        vy = op.wally[ip] ? 0.0 : op.dtbdy[ip] / tb[ip]
        vz = wz[ip]       ? 0.0 : op.dtbdz[ip] / tb[ip]
        den = 1.0 + λ2g * vz
        # den = 1 + lam^2*N^2. Positive for stable stratification; it can only
        # vanish where N^2 < 0, i.e. convective overturning at the step size.
        abs(den) < 1.0e-10 &&
            error("Schur: the pointwise momentum matrix is singular at node ", ip,
                  " (1 + lam^2*g*dtb/dz/tb = ", den, "). That needs N^2 < 0 -- ",
                  "unstable stratification -- with lam = ", lam, " at the ",
                  "overturning time scale 1/|N|. Lower Delta t.")
        # (I - u v'/den) r, with u = (0,0,lam^2 g)
        s = -λ2g * (vx*rx[ip] + vy*ry[ip] + vz*rz[ip]) / den
        mx[ip] = rx[ip]
        my[ip] = ry[ip]
        mz[ip] = rz[ip] + s
    end
    return nothing
end

"""
    schur_divW!(out, vx, vy, vz, params, op, work)

`out = Div[W v]`, the free-slip-masked divergence, read off block row 1:
`A_rho = -Div[W m]`, so the rho slot of `A` applied to a state whose only
non-zero is the momentum is exactly minus what we want. `W` is applied inside
the kernel, which is why it does not appear here.
"""
function schur_divW!(out, vx, vy, vz, params, op, work)
    kern = _schur_kernel(work, op)
    if kern !== nothing
        return schur_divW_fast!(out, vx, vy, vz, params, op, kern)
    end
    V, W = work.V, work.W
    fill!(V, 0.0)
    su, sv, sw = op.slot[2], op.slot[3], op.slot[4]
    @inbounds for ip in eachindex(vx)
        V[ip, su] = vx[ip]; V[ip, sv] = vy[ip]; V[ip, sw] = vz[ip]
    end
    hevi_apply_A!(W, V, params, op)
    sρ = op.slot[1]
    @inbounds for ip in eachindex(out)
        out[ip] = -W[ip, sρ]
    end
    return nothing
end

"""
    SchurState(npoin, nimp)

Scratch for the reduction. Holds the SchurWork the two assembled operators
need, plus the momentum triple and the pointwise coefficients.
"""
struct SchurState{TF <: AbstractFloat}
    w::SchurWork{Matrix{TF}, Vector{TF}}
    mx::Vector{TF}; my::Vector{TF}; mz::Vector{TF}
    rx::Vector{TF}; ry::Vector{TF}; rz::Vector{TF}
    d::Vector{TF}
    # b_m already combined with the known parts of (3); see schur_setup_rhs!
    cx::Vector{TF}; cy::Vector{TF}; cz::Vector{TF}
    qb::Vector{TF}
    # m_b = A^-1 c, cached by schur_setup_rhs!. NOT the same as c: A^-1 has
    # already been applied, so the P-driven part of the momentum is m - m_b and
    # NOT m - c. Getting that wrong gives a recovery that looks right and is
    # wrong only where grad(tb) is large.
    mbx::Vector{TF}; mby::Vector{TF}; mbz::Vector{TF}
end
function SchurState(npoin::Int, nimp::Int; kern = nothing)
    z() = zeros(Float64, npoin)
    SchurState(SchurWork(npoin, nimp; kern = kern), z(), z(), z(), z(), z(), z(), z(),
               z(), z(), z(), z(), z(), z(), z())
end

"""
    schur_momentum!(st, P, params, op, lam; known)

`m = A^-1 [ c - lam*Grad[P] - lam*g*zhat*P./(beta.*tb) ]`, the momentum written
in terms of the pressure. With `known = false` the `c` term is dropped, which
is what the Helmholtz OPERATOR needs; with `known = true` it is included, which
is what assembling its right-hand side needs. Same code both ways so the two
cannot drift apart.
"""
function schur_momentum!(st::SchurState, P, params, op, lam::Real; known::Bool = false)
    λ = Float64(lam); g = PhysicalConst{Float64}().g
    schur_grad!(st.w.gx, st.w.gy, st.w.gz, P, params, op, st.w)
    tb, β = op.thetabar, op.beta
    @inbounds for ip in eachindex(P)
        bp = P[ip] / (β[ip] * tb[ip])
        st.rx[ip] = -λ * st.w.gx[ip]                + (known ? st.cx[ip] : 0.0)
        st.ry[ip] = -λ * st.w.gy[ip]                + (known ? st.cy[ip] : 0.0)
        st.rz[ip] = -λ * st.w.gz[ip] - λ * g * bp   + (known ? st.cz[ip] : 0.0)
    end
    schur_amat_inv!(st.mx, st.my, st.mz, st.rx, st.ry, st.rz, op, λ)
    return nothing
end

"""
    schur_H!(HP, P, params, op, lam, st)

The scalar Schur operator

    H[P] = P./(beta.*tb) + lam*(W m).grad(tb)./tb + lam*Div[W m],   m = m(P)

i.e. equation (1) with rho written as `(P./beta - q)./tb` and every other
unknown eliminated. ONE scalar field in, one out.
"""
function schur_H!(HP, P, params, op, lam::Real, st::SchurState)
    λ = Float64(lam)
    schur_momentum!(st, P, params, op, λ; known = false)
    schur_divW!(st.d, st.mx, st.my, st.mz, params, op, st.w)
    tb, β = op.thetabar, op.beta
    @inbounds for ip in eachindex(P)
        mx = op.wallx[ip] ? 0.0 : st.mx[ip]
        my = op.wally[ip] ? 0.0 : st.my[ip]
        mz = op.wall[ip]  ? 0.0 : st.mz[ip]
        adv = mx*op.dtbdx[ip] + my*op.dtbdy[ip] + mz*op.dtbdz[ip]
        HP[ip] = P[ip] / (β[ip] * tb[ip]) + λ * adv / tb[ip] + λ * st.d[ip]
    end
    return HP
end

"""
    schur_setup_rhs!(rhs, st, B, params, op, lam)

Turn the 5-field right-hand side `B` of `(I - lam*A) U = B` into the scalar
right-hand side of `H[P] = rhs`, and cache the known parts the recovery needs.

    q_b = b_Theta - tb.*b_rho - lam*(W m_b).grad(tb)
    rhs = b_rho - lam*Div[W m_b] + q_b./tb

with `m_b = A^-1 c`, `c = b_m + lam*g*zhat*(b_Theta - tb.*b_rho)./tb`.
"""
function schur_setup_rhs!(rhs, st::SchurState, B::AbstractMatrix, params, op, lam::Real)
    λ = Float64(lam); g = PhysicalConst{Float64}().g
    sρ, su, sv, sw, sθ = op.slot
    tb = op.thetabar
    # c = b_m + lam*g*zhat*(b_Theta - tb*b_rho)/tb
    @inbounds for ip in eachindex(rhs)
        st.qb[ip] = B[ip, sθ] - tb[ip] * B[ip, sρ]        # b_Theta - tb*b_rho
        st.cx[ip] = B[ip, su]
        st.cy[ip] = B[ip, sv]
        st.cz[ip] = B[ip, sw] + λ * g * st.qb[ip] / tb[ip]
    end
    # m_b = A^-1 c  (P = 0, known = true)
    fill!(st.w.d, 0.0)
    schur_momentum!(st, st.w.d, params, op, λ; known = true)
    schur_divW!(st.d, st.mx, st.my, st.mz, params, op, st.w)
    @inbounds for ip in eachindex(rhs)
        mx = op.wallx[ip] ? 0.0 : st.mx[ip]
        my = op.wally[ip] ? 0.0 : st.my[ip]
        mz = op.wall[ip]  ? 0.0 : st.mz[ip]
        qb = st.qb[ip] - λ * (mx*op.dtbdx[ip] + my*op.dtbdy[ip] + mz*op.dtbdz[ip])
        st.qb[ip] = qb                                     # cached for recovery
        st.mbx[ip] = st.mx[ip]; st.mby[ip] = st.my[ip]; st.mbz[ip] = st.mz[ip]
        rhs[ip] = B[ip, sρ] - λ * st.d[ip] + qb / tb[ip]
    end
    return rhs
end

"""
    schur_recover!(U, P, st, B, params, op, lam)

Rebuild the five prognostic fields from the pressure, undoing the elimination:

    m     = A^-1 [ c - lam*Grad[P] - lam*g*zhat*P./(beta.*tb) ]
    q     = q_b - lam*(W m_P).grad(tb)        (m_P = m - m_b, the P-driven part)
    rho   = (P./beta - q) ./ tb
    Theta = q + tb.*rho = P./beta

`schur_setup_rhs!` must have run first: it leaves `c` and `q_b` in `st`.
"""
function schur_recover!(U::AbstractMatrix, P, st::SchurState, B::AbstractMatrix,
                        params, op, lam::Real)
    λ = Float64(lam)
    sρ, su, sv, sw, sθ = op.slot
    tb, β = op.thetabar, op.beta
    # P-DRIVEN momentum only: known = false gives A^-1 applied to the P terms
    # alone, which is what q needs -- m_b's contribution is already inside q_b.
    schur_momentum!(st, P, params, op, λ; known = false)
    @inbounds for ip in eachindex(P)
        mx = op.wallx[ip] ? 0.0 : st.mx[ip]
        my = op.wally[ip] ? 0.0 : st.my[ip]
        mz = op.wall[ip]  ? 0.0 : st.mz[ip]
        q  = st.qb[ip] - λ * (mx*op.dtbdx[ip] + my*op.dtbdy[ip] + mz*op.dtbdz[ip])
        Θ  = P[ip] / β[ip]
        U[ip, sρ] = (Θ - q) / tb[ip]
        # the FULL momentum is the P-driven part plus the cached m_b
        U[ip, su] = st.mx[ip] + st.mbx[ip]
        U[ip, sv] = st.my[ip] + st.mby[ip]
        U[ip, sw] = st.mz[ip] + st.mbz[ip]
        U[ip, sθ] = Θ
    end
    return U
end
