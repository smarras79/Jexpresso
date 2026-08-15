#---------------------------------------------------------------------------------
# SWsphereDSGS — source terms of the shallow water equations on a spherical shell.
#
# From Marras, Kopera & Giraldo (2015), QJRMS 141: 1727-1739, Eq. (8b), the
# momentum source is
#
#   S_mom = -φ∇φ - f(x × φu) + μx        (+ δν∇²(φu), handled by the viscous path)
#
# with the pressure term -φ∇φ carried in the flux instead (see user_flux.jl), so
# what is left here is the Coriolis force and the Lagrange multiplier:
#
#   S_mom = -f (x × φu) + μ x .
#
#
# CORIOLIS,  f = 2ω z / r²
# -------------------------
# The r² in f is not a typo and does not cancel the usual 2ω sin φ: the cross
# product with the POSITION vector x contributes a factor r, so
#
#   f (x × φu) = (2ω z/r²) · r (n̂ × φu) = 2ω sin(φ_lat) (n̂ × φu),
#
# i.e. the familiar f₀ = 2ω sin(latitude) acting on the tangential momentum,
# written without ever forming a latitude. Being a cross product with x it is
# automatically tangential — x·(x × φu) = 0 — so Coriolis alone never pushes the
# fluid off the shell.
#
#
# THE LAGRANGE MULTIPLIER  μx  — keeping the flow ON the shell
# ------------------------------------------------------------
# The velocity carries three Cartesian components but the fluid lives on a
# two-dimensional surface, so one degree of freedom is redundant and must be
# constrained: u·x = 0 at all times. μx is the constraint force that enforces
# it (paper, section 3.2, after Coté 1988).
#
# μ is USUALLY described as a projection applied to the state after each step —
# paper Eqs (9)-(11):
#
#   (φu)ⁿ⁺¹_c = P (φu)ⁿ⁺¹_u ,   P = I - x xᵀ/r² ,   μ = -(φu)ⁿ⁺¹_u·x / r²
#
# which is what project_momentum_to_sphere! in src/kernel/mesh/sphere_mesh.jl
# does. But μ also has a CLOSED POINTWISE FORM, which is what lets it live in a
# pointwise source callback at all. Requiring the momentum to stay tangential,
#
#   x · ∂(φu)/∂t = 0 ,
#
# and taking x· of the momentum equation:
#
#   * x·(x × φu) = 0                       Coriolis is already tangential
#   * x·∇ₛ(φ²/2) = 0                       a surface gradient is tangential
#   * x·[∇ₛ·(φu⊗u)] = -φ|u|²               the curvature term: differentiating
#                                          the Cartesian components of a
#                                          tangent field along the surface
#                                          produces exactly this normal part
#                                          (∂ˢⱼxᵢ = Pⱼᵢ and u·x = 0)
#
# leaves  φ|u|² + μr² = 0, hence
#
#   μ = -φ |u|² / r²        and        μx = -(φ|u|²/r) n̂ .
#
# The constraint force is CENTRIPETAL, of magnitude φ|u|²/r, pointing at the
# centre of the sphere — the force needed to bend a fluid parcel of "mass" φ
# moving at speed |u| around a circle of radius r. That physical reading is the
# quickest sanity check on the sign: μ must be negative.
#
# BOTH mechanisms are wanted, and they are not redundant:
#
#   * μx here keeps the CONTINUOUS equation consistent, so the exact solution
#     never leaves the shell and the discretization is not fighting a wrong PDE;
#   * project_momentum_to_sphere! removes the normal drift that the DISCRETE
#     operators accumulate anyway (μ above is exact only if the discrete
#     divergence is), and is applied to the state after each stage/step as in
#     the paper.
#
# There is deliberately no user_inputs.jl switch for μx: user_source! receives
# no `inputs`, the same reason ω is a literal below. To run without it — useful
# for seeing how much drift the projection is then left to remove — comment out
# the two `μ*` terms in S[2..4] below.
#
#
# VISCOSITY  δν∇²(φu)
# -------------------
# The paper regularizes with a CONSTANT ν = 1e5 m²/s, switched by δ. That term is
# a second derivative, not a pointwise source, so it belongs to the shell's
# viscous path in src/kernel/operators/sphere_rhs.jl rather than here.
#
# THIS deck does not use a constant ν at all: it sets :visc_model => DSGS_SW(),
# which makes ν residual-based and different in every element at every RK stage.
# The term in the equation above is unchanged — only where its coefficient comes
# from. See user_inputs.jl and the README.
#---------------------------------------------------------------------------------

function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=4, x=0.0, y=0.0, z=0.0,
                      xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0, zmin=0.0, zmax=0.0)

    #
    # Earth's angular velocity, Galewsky et al. (2004) / paper section 3.1.
    # Kept as a literal, like the other ShallowWater cases (TC2/TC3 do the same
    # with f₀), because user_source! receives no `inputs`.
    #
    ω = 7.292e-5

    φ  = q[1]
    φu = q[2]
    φv = q[3]
    φw = q[4]

    r2 = x*x + y*y + z*z

    #
    # Coriolis: -f (x × φu),  f = 2ω z/r²
    #
    f = 2.0*ω*z/r2

    cx = y*φw - z*φv
    cy = z*φu - x*φw
    cz = x*φv - y*φu

    #
    # Lagrange multiplier: μ = -φ|u|²/r², i.e. a centripetal force of magnitude
    # φ|u|²/r directed at the centre of the sphere.
    #
    u2 = (φu*φu + φv*φv + φw*φw)/(φ*φ)
    μ  = -φ*u2/r2

    S[1] = 0.0
    S[2] = -f*cx + μ*x
    S[3] = -f*cy + μ*y
    S[4] = -f*cz + μ*z
end


function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=4, x=0.0, y=0.0, z=0.0,
                      xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0, zmin=0.0, zmax=0.0)

    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, z=z,
                 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
end


function user_source_gpu(q, qe, x, y, z, PhysConst, lpert)
    T = eltype(q)

    ω  = T(7.292e-5)

    φ  = q[1]
    φu = q[2]
    φv = q[3]
    φw = q[4]

    r2 = x*x + y*y + z*z
    f  = T(2.0)*ω*z/r2

    cx = y*φw - z*φv
    cy = z*φu - x*φw
    cz = x*φv - y*φu

    u2 = (φu*φu + φv*φv + φw*φw)/(φ*φ)
    μ  = -φ*u2/r2

    return T(0.0), T(-f*cx + μ*x), T(-f*cy + μ*y), T(-f*cz + μ*z)
end
