#---------------------------------------------------------------------------------
# quasistatic_solid.jl — quasi-static, small-strain solid mechanics on
# Jexpresso's spectral-element discretisation: linear VISCOELASTIC solids and
# solids with QUADRATURE-POINT HISTORY (J2 plasticity; the same storage serves
# any strain-driven model, e.g. a shape-memory alloy).
#
# Selected by `:lquasistatic_solid => true` in user_inputs.jl. The driver then
# skips params_setup and the time integrator: there is no inertia and no
# M dq/dt = RHS. Time — or a pseudo-time load parameter — enters only through
# the material's memory and the loads, and every step solves equilibrium
#
#     R(uⁿ⁺¹) = f_int(uⁿ⁺¹; Hⁿ) − f_ext(tⁿ⁺¹) = 0,        u = ū on Γ_D,
#
# with Hⁿ the history of the material at the last converged step.
#
# Cases: problems/Elasticity/viscoelastic_cylinder2d (Artioli, Beirão da Veiga,
# Lovadina & Sacco, Comput. Mech. 60 (2017), Sec. 4.2) and
# problems/Elasticity/plastic_cylinder2d.
#
# NOTHING HERE IS A NEW DISCRETISATION. Every spatial operator is Jexpresso's
# own, on the grid, basis, metric terms and LGL quadrature of sem_setup:
#
#   stiffness       build_gradient_product_matrices / build_tangent_stiffness_matrices
#                   (element_matrices.jl, the layout of build_laplace_matrix)
#                   assembled by DSS_laplace_sparse
#   strain          _∇f! (operators.jl), Jexpresso's element gradient
#   internal force  the weak divergence ∫∇ψ·σ, contracted with dψ exactly as the
#                   viscous operator _expansion_visc! does, then DSS_rhs!
#   tractions       compute_segment_integral! + DSS_segment_integral!
#                   (surface_integral.jl), with the tangent-based edge Jacobian
#                   and normal of build_metric_terms! — curved-edge aware
#   body force      the lumped LGL mass matrix sem.matrix.M, as standard_linsolve!
#   nodal fields    DSS_rhs! + sem.matrix.M (stress, history, gradients)
#   output          write_output(…, VTK()) through the case's user_uout!
#
# What this file adds is the material models and the time/load stepping.
#
# MATERIAL MODELS — the "black box" of Artioli et al. Sec. 2.2
# ------------------------------------------------------------
# At every quadrature point a model maps the new strain and the converged
# history to the stress, the consistent tangent and the new history:
#
#     (εⁿ⁺¹, Hⁿ, Δt)  ↦  (σⁿ⁺¹, D = ∂σⁿ⁺¹/∂εⁿ⁺¹, Hⁿ⁺¹)      solid_constitutive_update!
#
# Plane strain throughout: ε_zz = 0, the material point is three-dimensional,
# so deviatoric tensors carry four components (xx, yy, zz, xy) and σ_zz is
# computed and output. D is the in-plane tangent in Voigt form with
# ENGINEERING shear: rows/columns (xx, yy, xy), γxy = 2εxy. The case picks the
# model in user_solid_material() (user_flux.jl):
#
#   model = :viscoelastic (default)   (; E, ν, g0, g, τ)
#       Generalized Maxwell in the deviator, elastic in volume (paper Sec. 3.1):
#           σ = K tr(ε) 1 + 2G [g₀ e + Σₘ gₘ qₘ],   q̇ₘ + qₘ/τₘ = ė,   e = dev ε.
#       τₘ is the paper's λ⁽ᵐ⁾, the RELAXATION TIME — not the Lamé constant.
#       History per point: eⁿ and every qₘ, 4 components each.
#
#   model = :j2_plasticity            (; E, ν, σy0, Hiso, Hkin)
#       von Mises plasticity with linear isotropic and kinematic hardening
#       (paper Sec. 3.2; Simo & Hughes 1998, Box 3.2):
#           s = 2G (e − eᵖ),  α = Hkin eᵖ,  X = s − α,
#           f = ‖X‖ − √(2/3) (σy0 + Hiso ēᵖ),   ėᵖ = ζ̇ X/‖X‖,   dēᵖ/dt = √(2/3) ζ̇.
#       Radial return; the consistent tangent of S&H (3.3.x). History per
#       point: eᵖ (4 components) and ēᵖ.
#
# A further model is one more `struct <: AbstractSolidModel` with methods
# nhistory, islinear, solid_history_fields and solid_constitutive_update!.
# A shape-memory alloy (paper Sec. 3.3) keeps the transformation strain eᵗʳ
# where J2 keeps eᵖ: the storage, the Newton loop and the output are the same.
#
# TIME INTEGRATION OF THE VISCOELASTIC MODEL — the exponential recursion
#
#     qₘⁿ⁺¹ = dₘ qₘⁿ + γₘ (eⁿ⁺¹ − eⁿ),        dₘ = exp(−Δt/τₘ)                (*)
#
#   "exponential"     γ = (τ/Δt)(1 − e^{−Δt/τ})   strain LINEAR in time over the
#      (default)        step: Zienkiewicz, Taylor & Fox (2014) and FEAP, which
#                       is what Artioli et al. ran. Second order.
#   "midpoint"        γ = e^{−Δt/2τ}               Simo & Hughes Box 10.1. 2nd order.
#   "backward_euler"  d = γ = 1/(1 + Δt/τ)         first order — the control.
#
# TWO WAYS TO CARRY THE HISTORY  (:solid_history)
#
#   "quadrature"   the general path: H lives at every quadrature point; each
#                  step is a Newton iteration with the consistent tangent,
#                      K_T Δu = f_ext − f_int(u),   K_T = Σ_e ∫ Bᵀ D B,
#                  refactorised every iteration for a nonlinear model.
#                  The only path for plasticity.
#   "nodal"        linear viscoelasticity only, and its default. Because (*) is
#                  linear and the discrete strain is linear in the nodal
#                  displacement, qₘ at every quadrature point is exactly
#                  dev ε(wₘ) of a NODAL field obeying wₘⁿ⁺¹ = dₘwₘⁿ + γₘ(uⁿ⁺¹ − uⁿ).
#                  Each step is then sparse algebra on two assembled operators,
#                      (K A_vol + G̃ A_dev) uⁿ⁺¹ = f − G A_dev Σₘ gₘ(dₘwₘⁿ − γₘuⁿ),
#                      G̃ = G(g₀ + Σₘ gₘγₘ),  v·A_vol u = ∫ div v div u,
#                      v·A_dev z = ∫ ε(v) : 2 dev ε(z),
#                  with one factorisation per Δt. The two paths agree to round-off
#                  (test/test_viscoelastic_cylinder.jl).
#
# THE INITIAL STEP. With :linitial_elastic_solve => true (default) the loads
# are on at t = tinit: a Δt = 0 step (d = γ = 1) gives the instantaneous
# response before the time steps — a Heaviside load. With false the run starts
# unloaded and the loads enter during the FIRST increment, which (*) reads as a
# ramp over [tinit, tinit+Δt] (FEAP; Artioli et al. Fig. 5.2).
#
# BOUNDARY CONDITIONS come from the case through one hook (user_bc.jl):
#
#     user_solid_bc(tag::String, x, y, t, nx, ny) -> NamedTuple
#
# called at every node of every boundary edge with that edge's gmsh tag and
# OUTWARD unit normal. Fields: `ux = value` / `uy = value` prescribe that
# displacement component (a roller is ONE component); `traction = (tx, ty)`
# prescribes σ·n; `(;)` is a traction-free surface. Constrained DOFs are
# eliminated — the solve is on the free block. Their set is fixed at tinit;
# the values may depend on t (displacement control).
#
# THE OTHER HOOKS
#
#     user_solid_material()           (user_flux.jl)    see MATERIAL MODELS
#     user_solid_body_force(x, y, t)  (user_source.jl)  (bx, by) per unit volume
#     user_solid_exact(x, y, t)       (user_analytic.jl, OPTIONAL)
#                                     (; ux, uy[, sxx, syy, sxy, szz])
#     initialize(...)                 qvars: SOLID_FIELDS, then any of the
#                                     extra fields below; displacement zero
#     user_uout!(...)                 the VTK fields, from qvars and the exact
#                                     solution at each node
#
# OUTPUT — what can be plotted
#   * VTK at every :diagnostics_at_times time (:outformat => "vtk"), through
#     the standard writer: iter_<k>.pvtu + simulation.pvd. The fields are the
#     case's qvars: SOLID_FIELDS (ux, uy, sxx, syy, sxy, szz) and any of
#       "svm"    von Mises equivalent stress √(3/2 s:s)
#       "umag"   |u|
#       the model's history fields (quadrature path), projected onto the nodes:
#       J2: ep_bar, ep_xx, ep_yy, ep_zz, ep_xy, alpha_xx, alpha_yy, alpha_zz,
#       alpha_xy;  viscoelastic: q<m>_xx, …
#     :vtk_vectors => [("u", "ux", "uy")] adds the displacement as a VECTOR, so
#     ParaView's Warp By Vector shows the deformed solid in one step.
#   * deformed-<field>-it<k>.png and deformed.gif with :lplot_deformed => true:
#     the deformed SEM grid (magnified by :plot_deformed_scale) coloured by
#     :plot_deformed_field — plot_deformed_mesh in src/io/plotting/jeplots.jl.
#   * solid_history.csv — always: every step, the displacement and the extra
#     fields at the :solid_probes points, the exact solution there, error norms,
#     the reaction force on every Dirichlet tag and the Newton iterations.
#
# LIMITS: 2D (NSD_2D) plane strain, serial, CPU, inexact (LGL-collocated)
# quadrature, no Laguerre/semi-infinite elements, small strains.
#---------------------------------------------------------------------------------

using WriteVTK      # src/io/write_output.jl loads it too

const SOLID_INTEGRATORS = (:exponential, :midpoint, :backward_euler)

# The fields the solver always hands to write_output: the case's initialize()
# declares these first in qvars, then any extra fields it wants written.
const SOLID_FIELDS = ["ux", "uy", "sxx", "syy", "sxy", "szz"]

# Plane-strain Voigt forms (rows/columns xx, yy, xy; engineering shear):
# K·SOLID_DVOL + G·SOLID_DDEV is the isotropic elastic D (λ = K − 2G/3, μ = G).
const SOLID_DVOL = [1.0 1.0 0.0; 1.0 1.0 0.0; 0.0 0.0 0.0]
const SOLID_DDEV = [4/3 -2/3 0.0; -2/3 4/3 0.0; 0.0 0.0 1.0]

"""
    maxwell_step_coefficients(Δt, τ, integrator) -> (d, γ)

Coefficients of one step of `qⁿ⁺¹ = d qⁿ + γ (eⁿ⁺¹ − eⁿ)` for the Maxwell
internal variable `q̇ + q/τ = ė`. `Δt == 0` is the instantaneous (glassy) jump,
`d = γ = 1`, whatever the integrator.
"""
function maxwell_step_coefficients(Δt::Real, τ::Real, integrator::Symbol)
    Δt == 0 && return (1.0, 1.0)
    x = Float64(Δt)/Float64(τ)
    if integrator === :exponential
        return (exp(-x), -expm1(-x)/x)          # (1 − e^{−x})/x, accurate at small x
    elseif integrator === :midpoint
        return (exp(-x), exp(-x/2))
    elseif integrator === :backward_euler
        return (1.0/(1.0 + x), 1.0/(1.0 + x))
    end
    error(" # ERROR quasistatic_solid.jl: unknown :viscoelastic_integrator ", integrator,
          ". Use one of ", join(string.(SOLID_INTEGRATORS), ", "), ".")
end

#=================================================================================
 Material models
=================================================================================#
abstract type AbstractSolidModel end

struct SolidViscoelastic <: AbstractSolidModel
    E::Float64; ν::Float64; K::Float64; G::Float64
    g0::Float64; g::Vector{Float64}; τ::Vector{Float64}     # τ: relaxation times
    integrator::Symbol
end

struct SolidJ2Plasticity <: AbstractSolidModel
    E::Float64; ν::Float64; K::Float64; G::Float64
    σy0::Float64; Hiso::Float64; Hkin::Float64
end

nhistory(m::SolidViscoelastic) = 4*(length(m.g) + 1)       # [eⁿ ; q₁ ; … ; q_M], 4 components each
nhistory(m::SolidJ2Plasticity) = 5                         # [eᵖ (xx, yy, zz, xy) ; ēᵖ]
islinear(::SolidViscoelastic)  = true
islinear(::SolidJ2Plasticity)  = false

"""
    solid_history_fields(model) -> Vector{Pair{String,Function}}

The plottable quantities a model keeps at every quadrature point, each as
`name => h -> value` of the history vector `h` of one point.
"""
function solid_history_fields(m::SolidJ2Plasticity)
    f = Pair{String,Function}["ep_bar" => (h -> h[5])]
    for (c, cn) in enumerate(("xx", "yy", "zz", "xy"))
        push!(f, string("ep_", cn) => (h -> h[c]))
        push!(f, string("alpha_", cn) => (h -> m.Hkin*h[c]))      # back stress α = Hkin eᵖ
    end
    return f
end

function solid_history_fields(m::SolidViscoelastic)
    f = Pair{String,Function}[]
    for mm = 1:length(m.g), (c, cn) in enumerate(("xx", "yy", "zz", "xy"))
        push!(f, string("q", mm, "_", cn) => (h -> h[4*mm + c]))
    end
    return f
end

#
# The case's material, checked and turned into a model.
#
function solid_model_from_case(inputs)

    isdefined(@__MODULE__, :user_solid_material) ||
        error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid => true needs ",
              "user_solid_material() in the case's user_flux.jl (see the header of quasistatic_solid.jl).")

    m    = user_solid_material()
    kind = Symbol(get(m, :model, :viscoelastic))
    E    = Float64(m.E)
    ν    = Float64(m.ν)
    E > 0          || error(" # ERROR quasistatic_solid.jl: E = $E must be positive.")
    -1.0 < ν < 0.5 || error(" # ERROR quasistatic_solid.jl: ν = $ν must lie in (-1, 1/2).")
    K = E/(3.0*(1.0 - 2.0*ν))
    G = E/(2.0*(1.0 + ν))

    if kind === :viscoelastic
        g0 = Float64(get(m, :g0, 1.0))
        g  = Float64[x for x in get(m, :g, Float64[])]
        τ  = Float64[x for x in get(m, :τ, Float64[])]
        length(g) == length(τ) ||
            error(" # ERROR quasistatic_solid.jl: the Prony series has ", length(g),
                  " weights g but ", length(τ), " relaxation times τ.")
        all(>(0), τ) || error(" # ERROR quasistatic_solid.jl: every relaxation time τ must be positive, got $τ.")
        (g0 >= 0 && all(>=(0), g)) ||
            error(" # ERROR quasistatic_solid.jl: Prony weights must be non-negative, got g0 = $g0, g = $g.")
        g0 + sum(g; init=0.0) > 0 || error(" # ERROR quasistatic_solid.jl: g0 + Σg must be positive.")
        integrator = Symbol(lowercase(string(get(inputs, :viscoelastic_integrator, "exponential"))))
        integrator in SOLID_INTEGRATORS ||
            error(" # ERROR quasistatic_solid.jl: :viscoelastic_integrator => \"", integrator,
                  "\" is not one of ", join(string.(SOLID_INTEGRATORS), ", "), ".")
        return SolidViscoelastic(E, ν, K, G, g0, g, τ, integrator)

    elseif kind === :j2_plasticity
        σy0  = Float64(m.σy0)
        Hiso = Float64(get(m, :Hiso, 0.0))
        Hkin = Float64(get(m, :Hkin, 0.0))
        σy0 > 0 || error(" # ERROR quasistatic_solid.jl: the initial yield stress σy0 = $σy0 must be positive.")
        (Hiso >= 0 && Hkin >= 0) ||
            error(" # ERROR quasistatic_solid.jl: hardening moduli must be non-negative (Hiso = $Hiso, Hkin = $Hkin).")
        return SolidJ2Plasticity(E, ν, K, G, σy0, Hiso, Hkin)
    end

    error(" # ERROR quasistatic_solid.jl: user_solid_material() asks for model = :", kind,
          " — use :viscoelastic or :j2_plasticity.")
end

"""
    solid_constitutive_update!(σ, D, hnew, model, εxx, εyy, εxy, hold, δt)

The material point: from the in-plane strain (tensor components, ε_zz = 0) at
the new time, the converged history `hold` and the step `δt`, write the stress
`σ = (σxx, σyy, σzz, σxy)`, the consistent in-plane tangent `D` (Voigt,
engineering shear) and the new history `hnew`. `hold` is not modified, so the
update can be re-evaluated at every Newton iteration.
"""
function solid_constitutive_update!(σ, D, hnew, m::SolidViscoelastic, εxx, εyy, εxy, hold, δt)
    K, G = m.K, m.G
    tr = εxx + εyy                                      # ε_zz = 0
    e1 = εxx - tr/3; e2 = εyy - tr/3; e3 = -tr/3; e4 = εxy
    s1 = m.g0*e1; s2 = m.g0*e2; s3 = m.g0*e3; s4 = m.g0*e4
    Gt = G*m.g0
    @inbounds for mm = 1:length(m.g)
        d, γ = maxwell_step_coefficients(δt, m.τ[mm], m.integrator)
        o  = 4*mm
        q1 = d*hold[o+1] + γ*(e1 - hold[1]); q2 = d*hold[o+2] + γ*(e2 - hold[2])
        q3 = d*hold[o+3] + γ*(e3 - hold[3]); q4 = d*hold[o+4] + γ*(e4 - hold[4])
        hnew[o+1] = q1; hnew[o+2] = q2; hnew[o+3] = q3; hnew[o+4] = q4
        s1 += m.g[mm]*q1; s2 += m.g[mm]*q2; s3 += m.g[mm]*q3; s4 += m.g[mm]*q4
        Gt += G*m.g[mm]*γ
    end
    hnew[1] = e1; hnew[2] = e2; hnew[3] = e3; hnew[4] = e4
    p = K*tr
    σ[1] = p + 2G*s1; σ[2] = p + 2G*s2; σ[3] = p + 2G*s3; σ[4] = 2G*s4
    @inbounds for j = 1:3, i = 1:3
        D[i,j] = K*SOLID_DVOL[i,j] + Gt*SOLID_DDEV[i,j]
    end
    return nothing
end

function solid_constitutive_update!(σ, D, hnew, m::SolidJ2Plasticity, εxx, εyy, εxy, hold, δt)
    K, G = m.K, m.G
    tr = εxx + εyy                                      # ε_zz = 0
    e1 = εxx - tr/3; e2 = εyy - tr/3; e3 = -tr/3; e4 = εxy
    p1, p2, p3, p4, pbar = hold[1], hold[2], hold[3], hold[4], hold[5]

    # Elastic trial: s = 2G(e − eᵖ), relative stress X = s − α, α = Hkin eᵖ.
    s1 = 2G*(e1 - p1); s2 = 2G*(e2 - p2); s3 = 2G*(e3 - p3); s4 = 2G*(e4 - p4)
    x1 = s1 - m.Hkin*p1; x2 = s2 - m.Hkin*p2; x3 = s3 - m.Hkin*p3; x4 = s4 - m.Hkin*p4
    nx = sqrt(x1^2 + x2^2 + x3^2 + 2.0*x4^2)            # ‖X‖, Frobenius (xy counted twice)
    R  = sqrt(2.0/3.0)*(m.σy0 + m.Hiso*pbar)            # radius of the elastic domain
    f  = nx - R

    if f <= 1.0e-12*R                                   # elastic step
        hnew[1] = p1; hnew[2] = p2; hnew[3] = p3; hnew[4] = p4; hnew[5] = pbar
        @inbounds for j = 1:3, i = 1:3
            D[i,j] = K*SOLID_DVOL[i,j] + G*SOLID_DDEV[i,j]
        end
    else                                                # radial return
        n1 = x1/nx; n2 = x2/nx; n3 = x3/nx; n4 = x4/nx
        Hm = 2G + m.Hkin + (2.0/3.0)*m.Hiso
        Δγ = f/Hm
        s1 -= 2G*Δγ*n1; s2 -= 2G*Δγ*n2; s3 -= 2G*Δγ*n3; s4 -= 2G*Δγ*n4
        hnew[1] = p1 + Δγ*n1; hnew[2] = p2 + Δγ*n2
        hnew[3] = p3 + Δγ*n3; hnew[4] = p4 + Δγ*n4
        hnew[5] = pbar + sqrt(2.0/3.0)*Δγ
        # Consistent tangent: K 1⊗1 + 2Gθ I_dev − 2Gθ̄ n⊗n (Simo & Hughes 1998, Box 3.2)
        θ  = 1.0 - 2G*Δγ/nx
        θb = 2G/Hm - (1.0 - θ)
        v  = (n1, n2, n4)
        @inbounds for j = 1:3, i = 1:3
            D[i,j] = K*SOLID_DVOL[i,j] + G*θ*SOLID_DDEV[i,j] - 2G*θb*v[i]*v[j]
        end
    end
    p = K*tr
    σ[1] = p + s1; σ[2] = p + s2; σ[3] = p + s3; σ[4] = s4
    return nothing
end

#=================================================================================
 SEM operators — all built from Jexpresso's element matrices, gradient and DSS
=================================================================================#

#
# DSS_laplace_sparse sizes its result from the largest index it stored, and it
# drops entries below eps: an off-diagonal block such as Sxy can come out
# smaller than npoin x npoin. Square it up before the blocks are stacked.
#
function _solid_npoin_square(S::SparseMatrixCSC, n::Integer)
    size(S) == (n, n) && return S
    I, J, V = findnz(S)
    return sparse(I, J, V, n, n)
end

"""
    solid_sem_operators(sem) -> (A_vol, A_dev)

The two constant-coefficient plane-strain operators, from Jexpresso's element
matrices (build_gradient_product_matrices) and sparse DSS (DSS_laplace_sparse),
in the flat layout u = [ux; uy]: `v·A_vol u = ∫ div v div u`,
`v·A_dev z = ∫ ε(v):2dev ε(z)`.
"""
function solid_sem_operators(sem)
    mesh = sem.mesh
    N    = mesh.nop
    Lxx, Lxy, Lyy = build_gradient_product_matrices(mesh.SD, sem.basis.ψ, sem.basis.dψ, sem.ω,
                                                    mesh.nelem, mesh, sem.metrics, N, N, Float64)
    n   = Int(mesh.npoin)
    Sxx = _solid_npoin_square(DSS_laplace_sparse(mesh, Lxx), n)
    Sxy = _solid_npoin_square(DSS_laplace_sparse(mesh, Lxy), n)
    Syy = _solid_npoin_square(DSS_laplace_sparse(mesh, Lyy), n)
    Syx = sparse(transpose(Sxy))

    A_vol = [Sxx  Sxy;
             Syx  Syy]
    A_dev = [(4/3)*Sxx + Syy    Syx - (2/3)*Sxy;
             Sxy - (2/3)*Syx    Sxx + (4/3)*Syy]
    return A_vol, A_dev
end

"""
    solid_tangent_matrix(Dq, sem) -> K_T

The global tangent stiffness Σ_e ∫ Bᵀ D B for a tangent `Dq[1:3,1:3,k,l,iel]`
that varies from quadrature point to quadrature point:
build_tangent_stiffness_matrices + DSS_laplace_sparse, flat layout [ux; uy].
"""
function solid_tangent_matrix(Dq, sem)
    mesh = sem.mesh
    N    = mesh.nop
    Kxx, Kxy, Kyx, Kyy = build_tangent_stiffness_matrices(mesh.SD, sem.basis.ψ, sem.basis.dψ, sem.ω,
                                                          mesh.nelem, mesh, sem.metrics, N, N, Dq, Float64)
    n = Int(mesh.npoin)
    return [_solid_npoin_square(DSS_laplace_sparse(mesh, Kxx), n)  _solid_npoin_square(DSS_laplace_sparse(mesh, Kxy), n);
            _solid_npoin_square(DSS_laplace_sparse(mesh, Kyx), n)  _solid_npoin_square(DSS_laplace_sparse(mesh, Kyy), n)]
end

"""
    solid_qp_strains!(εq, u, sem, f_el, ∇ux_el, ∇uy_el)

The in-plane strain (εxx, εyy, εxy — tensor shear) at every quadrature point,
from Jexpresso's element gradient _∇f! (which returns ωJ ∇f at the LGL points;
with collocation the quadrature point (k,l) is the node (k,l)).
"""
function solid_qp_strains!(εq, u, sem, f_el, ∇ux_el, ∇uy_el)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    npoin   = Int(mesh.npoin)
    ω       = sem.ω
    for iel = 1:mesh.nelem
        for j = 1:ngl, i = 1:ngl
            f_el[i, j] = u[mesh.connijk[iel, i, j]]
        end
        _∇f!(∇ux_el, f_el, ngl, sem.basis.dψ, ω, metrics.Je,
             metrics.dξdx, metrics.dξdy, metrics.dηdx, metrics.dηdy, iel, CL(), sem.QT, mesh.SD, ContGal())
        for j = 1:ngl, i = 1:ngl
            f_el[i, j] = u[npoin + mesh.connijk[iel, i, j]]
        end
        _∇f!(∇uy_el, f_el, ngl, sem.basis.dψ, ω, metrics.Je,
             metrics.dξdx, metrics.dξdy, metrics.dηdx, metrics.dηdy, iel, CL(), sem.QT, mesh.SD, ContGal())
        for j = 1:ngl, i = 1:ngl
            ωJ = ω[i]*ω[j]*metrics.Je[iel, i, j]
            εq[1, i, j, iel] = ∇ux_el[iel, i, j, 1]/ωJ
            εq[2, i, j, iel] = ∇uy_el[iel, i, j, 2]/ωJ
            εq[3, i, j, iel] = 0.5*(∇ux_el[iel, i, j, 2] + ∇uy_el[iel, i, j, 1])/ωJ
        end
    end
    return nothing
end

"""
    solid_internal_force!(fint, σq, sem, rhs_el, tmp)

f_int = Σ_e ∫ ∇ψ · σ dΩ for a quadrature-point stress σq = (σxx, σyy, σzz, σxy):
the weak divergence of the two rows of σ, contracted with dψ exactly as the
viscous operator _expansion_visc! (rhs.jl) does, then DSS_rhs!.
"""
function solid_internal_force!(fint, σq, sem, rhs_el, tmp)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    npoin   = Int(mesh.npoin)
    ω       = sem.ω
    dψ      = sem.basis.dψ
    fill!(rhs_el, 0.0)
    @inbounds for iel = 1:mesh.nelem
        for l = 1:ngl, k = 1:ngl
            ωJ  = ω[k]*ω[l]*metrics.Je[iel, k, l]
            ξx  = metrics.dξdx[iel, k, l]; ξy = metrics.dξdy[iel, k, l]
            ηx  = metrics.dηdx[iel, k, l]; ηy = metrics.dηdy[iel, k, l]
            sxx = σq[1, k, l, iel]; syy = σq[2, k, l, iel]; sxy = σq[4, k, l, iel]
            # contravariant fluxes of the x row (σxx, σxy) and of the y row (σxy, σyy)
            Fξx = ωJ*(ξx*sxx + ξy*sxy); Fηx = ωJ*(ηx*sxx + ηy*sxy)
            Fξy = ωJ*(ξx*sxy + ξy*syy); Fηy = ωJ*(ηx*sxy + ηy*syy)
            for i = 1:ngl
                rhs_el[iel, i, l, 1] += dψ[i, k]*Fξx
                rhs_el[iel, k, i, 1] += dψ[i, l]*Fηx
                rhs_el[iel, i, l, 2] += dψ[i, k]*Fξy
                rhs_el[iel, k, i, 2] += dψ[i, l]*Fηy
            end
        end
    end
    fill!(tmp, 0.0)
    DSS_rhs!(tmp, rhs_el, mesh.connijk, mesh.nelem, ngl, 2, mesh.SD, ContGal())
    @views fint[1:npoin]         .= tmp[:, 1]
    @views fint[npoin+1:2*npoin] .= tmp[:, 2]
    return fint
end

"""
    solid_qp_to_nodes!(out, value, sem, rhs_el1, tmp1)

Nodal values of a quadrature-point quantity `value(k, l, iel)`: Σ_e ωJ value / M
over the elements sharing a node (DSS_rhs! and the lumped mass matrix).
"""
function solid_qp_to_nodes!(out, value, sem, rhs_el1, tmp1)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    ω       = sem.ω
    @inbounds for iel = 1:mesh.nelem, l = 1:ngl, k = 1:ngl
        rhs_el1[iel, k, l, 1] = ω[k]*ω[l]*metrics.Je[iel, k, l]*value(k, l, iel)
    end
    fill!(tmp1, 0.0)
    DSS_rhs!(tmp1, rhs_el1, mesh.connijk, mesh.nelem, ngl, 1, mesh.SD, ContGal())
    out .= @view(tmp1[:, 1]) ./ sem.matrix.M
    return out
end

"""
    solid_nodal_gradient!(grad, field, sem, f_el, ∇f_el)

The SEM nodal gradient of a nodal field: _∇f! (ωJ-weighted element gradient),
DSS_rhs!, and the lumped mass matrix. `grad` is npoin x 2.
"""
function solid_nodal_gradient!(grad, field, sem, f_el, ∇f_el)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    for iel = 1:mesh.nelem
        for j = 1:ngl, i = 1:ngl
            f_el[i, j] = field[mesh.connijk[iel, i, j]]
        end
        _∇f!(∇f_el, f_el, ngl, sem.basis.dψ, sem.ω, metrics.Je,
             metrics.dξdx, metrics.dξdy, metrics.dηdx, metrics.dηdy,
             iel, CL(), sem.QT, mesh.SD, ContGal())
    end
    fill!(grad, 0.0)
    DSS_rhs!(grad, ∇f_el, mesh.connijk, mesh.nelem, ngl, 2, mesh.SD, ContGal())
    grad ./= sem.matrix.M
    return grad
end

#
# Boundary loads and constraints at time t. Tractions are integrated with
# Jexpresso's segment integral (ω times the curved-edge Jacobian Jef, then the
# DSS onto the nodes); Dirichlet components are collected as dof => value, in
# the flat layout (ux of node ip is dof ip, uy is npoin + ip), and dof => tag.
#
function solid_boundary_conditions!(f, dirichlet::Dict{Int,Float64}, dof_tag::Dict{Int,String},
                                    t, sem, S_edge, F_edge, S_flux)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    npoin   = Int(mesh.npoin)

    empty!(dirichlet)
    fill!(S_edge, 0.0)
    fill!(S_flux, 0.0)
    for iedge = 1:mesh.nedges_bdy
        tag = mesh.bdy_edge_type[iedge]
        tag === nothing && continue
        tag = String(tag)
        fill!(F_edge, 0.0)
        for k = 1:ngl
            ip = mesh.poin_in_bdy_edge[iedge, k]
            x  = mesh.coords[1, ip]; y = mesh.coords[2, ip]
            bc = user_solid_bc(tag, x, y, t, metrics.nx[iedge, k], metrics.ny[iedge, k])
            if haskey(bc, :traction)
                F_edge[k, 1], F_edge[k, 2] = bc.traction
            end
            haskey(bc, :ux) && _solid_set_dirichlet!(dirichlet, dof_tag, ip,         Float64(bc.ux), tag, x, y)
            haskey(bc, :uy) && _solid_set_dirichlet!(dirichlet, dof_tag, npoin + ip, Float64(bc.uy), tag, x, y)
        end
        compute_segment_integral!(S_edge, F_edge, sem.ω, metrics.Jef, iedge, ngl)
    end
    DSS_segment_integral!(S_flux, S_edge, nothing, mesh.nedges_bdy, ngl, mesh.connijk,
                          mesh.poin_in_bdy_edge, mesh.bdy_edge_in_elem, 2)
    @views f[1:npoin]         .+= S_flux[:, 1]
    @views f[npoin+1:2*npoin] .+= S_flux[:, 2]
    return nothing
end

function _solid_set_dirichlet!(dirichlet, dof_tag, dof, value, tag, x, y)
    old = get(dirichlet, dof, nothing)
    if old !== nothing && abs(old - value) > 1.0e-12*max(1.0, abs(value))
        error(" # ERROR quasistatic_solid.jl: two boundary edges prescribe different values ",
              "(", old, " and ", value, ") on the same displacement component at (", x, ", ", y,
              "), tag \"", tag, "\".")
    end
    dirichlet[dof] = value
    haskey(dof_tag, dof) || (dof_tag[dof] = tag)
    return nothing
end

#
# f += ∫ φ b dΩ with the LGL-collocated (lumped) mass, as standard_linsolve!.
#
function solid_add_body_force!(f, t, sem)
    mesh  = sem.mesh
    npoin = Int(mesh.npoin)
    M     = sem.matrix.M
    for ip = 1:npoin
        bx, by = user_solid_body_force(mesh.coords[1, ip], mesh.coords[2, ip], t)
        f[ip]         += M[ip]*bx
        f[npoin + ip] += M[ip]*by
    end
    return nothing
end

function _solid_check_supported(sem, inputs)
    mesh = sem.mesh
    mesh.SD isa NSD_2D ||
        error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid is implemented for 2D (plane strain) only.")
    sem.QT isa Inexact ||
        error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid needs the LGL-collocated (inexact) ",
              "quadrature; remove :lexact_integration from the deck.")
    inputs[:backend] == CPU() ||
        error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid runs on the CPU backend only.")
    MPI.Comm_size(get_mpi_comm()) == 1 ||
        error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid assembles one global sparse ",
              "matrix and is serial; run it on one MPI rank.")
    mesh.nelem_semi_inf == 0 ||
        error(" # ERROR quasistatic_solid.jl: Laguerre/semi-infinite elements are not supported.")
    for f in (:user_solid_bc, :user_solid_body_force)
        isdefined(@__MODULE__, f) ||
            error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid needs `", f,
                  "` from the case (user_bc.jl / user_source.jl). See the header of ",
                  "src/kernel/solvers/quasistatic_solid.jl.")
    end
    return nothing
end

#
# Factorise a free-DOF block: Cholesky for the symmetric positive definite
# tangents of these models, LU where Cholesky refuses (a perfectly plastic
# zone at the onset of collapse can make the block only semi-definite).
#
function _solid_factorize(A)
    try
        return cholesky(Symmetric(A))
    catch e
        e isa PosDefException || rethrow()
        return lu(A)
    end
end

#=================================================================================
 The driver
=================================================================================#
"""
    quasistatic_solid_solve!(sem, qp, inputs, OUTPUT_DIR, TFloat)

March the quasi-static problem from `:tinit` to `:tend` in steps of `:Δt`.
See the header of this file for the models, the hooks and the output.
"""
function quasistatic_solid_solve!(sem, qp, inputs, OUTPUT_DIR::String, TFloat)

    _solid_check_supported(sem, inputs)

    mesh  = sem.mesh
    ngl   = mesh.ngl
    nq    = ngl
    nelem = Int(mesh.nelem)
    npoin = Int(mesh.npoin)
    ndof  = 2*npoin
    M     = sem.matrix.M

    model    = solid_model_from_case(inputs)
    linitial = get(inputs, :linitial_elastic_solve, true) == true
    tinit    = Float64(inputs[:tinit])
    tend     = Float64(inputs[:tend])
    Δt       = Float64(inputs[:Δt])
    (Δt > 0 && tend > tinit) ||
        error(" # ERROR quasistatic_solid.jl: need Δt > 0 and tend > tinit (got Δt = $Δt, [$tinit, $tend]).")

    # Where the history lives.
    hpath = Symbol(lowercase(string(get(inputs, :solid_history,
                                        model isa SolidViscoelastic ? "nodal" : "quadrature"))))
    hpath in (:nodal, :quadrature) ||
        error(" # ERROR quasistatic_solid.jl: :solid_history => \"", hpath, "\" — use \"nodal\" or \"quadrature\".")
    (hpath === :nodal && !(model isa SolidViscoelastic)) &&
        error(" # ERROR quasistatic_solid.jl: :solid_history => \"nodal\" is for linear viscoelasticity only; ",
              typeof(model), " keeps its history at the quadrature points.")
    newton_maxit = Int(get(inputs, :newton_maxit, 25))
    newton_rtol  = Float64(get(inputs, :newton_rtol, 1.0e-10))

    # user_solid_exact returning `nothing` means "no closed form" (the hook is a
    # module method, so another case's may still be defined in this session).
    has_exact = get(inputs, :_has_analytic, false) == true && isdefined(@__MODULE__, :user_solid_exact) &&
                user_solid_exact(mesh.coords[1, 1], mesh.coords[2, 1], tinit) !== nothing
    compare   = has_exact && linitial

    # Fields written: SOLID_FIELDS, then the case's extras (svm, umag, history).
    hfields = solid_history_fields(model)
    qvars   = String.(qp.qvars)
    (length(qvars) >= length(SOLID_FIELDS) && qvars[1:length(SOLID_FIELDS)] == SOLID_FIELDS) ||
        error(" # ERROR quasistatic_solid.jl: the case's initialize() must declare qvars starting with ",
              SOLID_FIELDS, ", got ", qvars, ".")
    extras = qvars[length(SOLID_FIELDS)+1:end]
    for name in extras
        name in ("svm", "umag") && continue
        any(p -> first(p) == name, hfields) ||
            error(" # ERROR quasistatic_solid.jl: no field \"", name, "\" to write. Available: svm, umag, ",
                  join(first.(hfields), ", "), ".")
        hpath === :quadrature ||
            error(" # ERROR quasistatic_solid.jl: the history field \"", name, "\" needs ",
                  ":solid_history => \"quadrature\".")
    end
    nvar = length(qvars)

    # The run's step times: uniform Δt, the last step shortened to land on tend.
    tsteps = Float64[tinit]
    while tsteps[end] < tend - 1.0e-10*max(1.0, abs(tend))
        push!(tsteps, min(tinit + length(tsteps)*Δt, tend))     # tinit + nΔt: no drift
    end
    abs(tsteps[end] - tend) < 1.0e-9*max(1.0, abs(tend)) && (tsteps[end] = tend)
    nsteps = length(tsteps) - 1

    # Output times: the deck's :diagnostics_at_times (a number or a list).
    _dt = get(inputs, :diagnostics_at_times, tend)
    out_times = _dt isa Number ? Float64[_dt] : Float64[x for x in _dt]
    is_output(t) = any(to -> abs(to - t) <= 1.0e-9*max(1.0, abs(t)), out_times)
    lvtk     = inputs[:outformat] isa VTK
    lplotdef = get(inputs, :lplot_deformed, false) == true

    # Probes: the mesh node nearest each requested point.
    probes = Tuple{String,Int}[]
    for p in get(inputs, :solid_probes, ())
        name, px, py = String(p[1]), Float64(p[2]), Float64(p[3])
        ipbest = argmin([(mesh.coords[1, ip] - px)^2 + (mesh.coords[2, ip] - py)^2 for ip = 1:npoin])
        push!(probes, (name, ipbest))
    end

    #-----------------------------------------------------------------------
    # Report the set-up.
    #-----------------------------------------------------------------------
    println(" # Quasi-static solid (plane strain) .................................")
    @printf("   SEM: nop = %d, %d elements, %d nodes, %d displacement DOFs\n", mesh.nop, nelem, npoin, ndof)
    @printf("   material: E = %g, ν = %g  →  K = %.6g, G = %.6g\n", model.E, model.ν, model.K, model.G)
    if model isa SolidViscoelastic
        if isempty(model.g)
            println("   no Prony terms: linear ELASTIC solid")
        else
            @printf("   viscoelastic, Prony series G(t)/G = %g", model.g0)
            for m = 1:length(model.g)
                @printf(" + %g exp(-t/%g)", model.g[m], model.τ[m])
            end
            println("      (τ = relaxation time, the paper's λ; NOT the Lamé constant)")
            abs(model.g0 + sum(model.g) - 1.0) > 1.0e-12 &&
                @printf("   NOTE g0 + Σg = %.6g ≠ 1: the instantaneous shear modulus is G·%.6g, not G\n",
                        model.g0 + sum(model.g), model.g0 + sum(model.g))
        end
        @printf("   integrator \"%s\"\n", model.integrator)
    else
        @printf("   J2 plasticity: σy0 = %g, Hiso = %g, Hkin = %g%s\n", model.σy0, model.Hiso, model.Hkin,
                model.Hiso == 0 && model.Hkin == 0 ? "   (perfectly plastic)" : "")
    end
    @printf("   history at the %s, %d values per point\n",
            hpath === :nodal ? "NODES (linear viscoelastic shortcut)" : "QUADRATURE POINTS", nhistory(model))
    @printf("   time: [%g, %g], Δt = %g, %d steps\n", tinit, tend, Δt, nsteps)
    println(linitial ?
            "   loads ON at t = tinit: instantaneous elastic solve first (Heaviside load)" :
            "   loads enter during the FIRST increment (FEAP convention: a ramp over [tinit, tinit+Δt])")
    for (name, ip) in probes
        @printf("   probe %-4s → node %d at (%.6g, %.6g)\n", name, ip, mesh.coords[1, ip], mesh.coords[2, ip])
    end
    isempty(extras) || println("   extra fields: ", join(extras, ", "))
    has_exact && println("   exact solution: user_solid_exact (user_analytic.jl)",
                         linitial ? "" : " — NOT compared: it describes a Heaviside load at tinit")

    #-----------------------------------------------------------------------
    # Constraints, work arrays.
    #-----------------------------------------------------------------------
    nedges = Int(mesh.nedges_bdy)
    S_edge = zeros(nedges, ngl, 2)
    F_edge = zeros(ngl, 2)
    S_flux = zeros(npoin, 2)

    u = zeros(ndof)
    u[1:npoin]      .= @view qp.qn[1:npoin, 1]
    u[npoin+1:ndof] .= @view qp.qn[1:npoin, 2]
    maximum(abs, u; init=0.0) == 0.0 ||
        error(" # ERROR quasistatic_solid.jl: the initial displacement from initialize() is not zero. ",
              "The quasi-static path starts from the undeformed, unstressed state: the history ",
              "that produced a prior deformation is not known.")

    fext      = zeros(ndof)
    fint      = zeros(ndof)
    dirichlet = Dict{Int,Float64}()
    dof_tag   = Dict{Int,String}()
    solid_boundary_conditions!(fext, dirichlet, dof_tag, tinit, sem, S_edge, F_edge, S_flux)
    cdofs  = sort!(collect(keys(dirichlet)))
    isfree = trues(ndof); isfree[cdofs] .= false
    fdofs  = findall(isfree)
    dtags  = sort!(unique(values(dof_tag)))
    println("   ", length(cdofs), " constrained DOFs, ", length(fdofs), " free")

    function loads!(t)
        fill!(fext, 0.0)
        solid_boundary_conditions!(fext, dirichlet, dof_tag, t, sem, S_edge, F_edge, S_flux)
        sort!(collect(keys(dirichlet))) == cdofs ||
            error(" # ERROR quasistatic_solid.jl: the set of Dirichlet DOFs changed at t = $t; ",
                  "user_solid_bc may change the VALUES it prescribes in time, not WHICH components.")
        solid_add_body_force!(fext, t, sem)
        return nothing
    end

    f_el   = zeros(ngl, ngl)
    ∇a_el  = zeros(nelem, ngl, ngl, 2)
    ∇b_el  = zeros(nelem, ngl, ngl, 2)
    rhs_el = zeros(nelem, ngl, ngl, 2)
    tmp2   = zeros(npoin, 2)
    rhs1   = zeros(nelem, ngl, ngl, 1)
    tmp1   = zeros(npoin, 1)
    newton_iters = Ref(0)

    #-----------------------------------------------------------------------
    # PATH 1 — linear viscoelasticity with nodal history (see the header).
    #-----------------------------------------------------------------------
    if hpath === :nodal
        A_vol, A_dev = solid_sem_operators(sem)
        nM = length(model.g)
        w  = [zeros(ndof) for _ = 1:nM]
        d  = zeros(nM); γ = zeros(nM)
        hist = zeros(ndof)
        z    = zeros(ndof)
        factor_cache = Dict{Float64,Any}()
        A_fc_cache   = Dict{Float64,Any}()
        system_for = function (Gt)
            get!(factor_cache, Gt) do
                A = model.K*A_vol + Gt*A_dev
                A_fc_cache[Gt] = A[fdofs, cdofs]
                try
                    cholesky(Symmetric(A[fdofs, fdofs]))
                catch e
                    e isa PosDefException || rethrow()
                    error(" # ERROR quasistatic_solid.jl: the stiffness on the free DOFs is not positive ",
                          "definite. The Dirichlet conditions in user_solid_bc probably leave a rigid-body ",
                          "motion (translation or rotation) unconstrained.")
                end
            end
        end
        advance! = function (t, δt)
            Gt = model.G*model.g0
            fill!(hist, 0.0)
            for m = 1:nM
                d[m], γ[m] = maxwell_step_coefficients(δt, model.τ[m], model.integrator)
                Gt += model.G*model.g[m]*γ[m]
                hist .+= model.g[m] .* (d[m] .* w[m] .- γ[m] .* u)     # Σ gₘ (dₘ wₘⁿ − γₘ uⁿ)
            end
            loads!(t)
            f = copy(fext)
            nM > 0 && (f .-= model.G .* (A_dev*hist))
            F    = system_for(Gt)
            uc   = Float64[dirichlet[i] for i in cdofs]
            unew = similar(u)
            unew[fdofs] = F \ (f[fdofs] .- A_fc_cache[Gt]*uc)
            unew[cdofs] = uc
            for m = 1:nM
                w[m] .= d[m] .* w[m] .+ γ[m] .* (unew .- u)
            end
            u .= unew
            # f_int = K A_vol u + G A_dev z (z = g₀u + Σ gₘwₘ): the reactions come from it
            z .= model.g0 .* u
            for m = 1:nM
                z .+= model.g[m] .* w[m]
            end
            fint .= model.K .* (A_vol*u) .+ model.G .* (A_dev*z)
            newton_iters[] = 1
            return nothing
        end
        gux = zeros(npoin, 2); guy = zeros(npoin, 2); gzx = zeros(npoin, 2); gzy = zeros(npoin, 2)
        nodal_stress! = function (S)                   # S[ip, 1:4] = σxx, σyy, σzz, σxy
            solid_nodal_gradient!(gux, @view(u[1:npoin]),      sem, f_el, ∇a_el)
            solid_nodal_gradient!(guy, @view(u[npoin+1:ndof]), sem, f_el, ∇a_el)
            solid_nodal_gradient!(gzx, @view(z[1:npoin]),      sem, f_el, ∇a_el)
            solid_nodal_gradient!(gzy, @view(z[npoin+1:ndof]), sem, f_el, ∇a_el)
            for ip = 1:npoin
                p   = model.K*(gux[ip, 1] + guy[ip, 2])
                exx = gzx[ip, 1]; eyy = gzy[ip, 2]; exy = 0.5*(gzx[ip, 2] + gzy[ip, 1])
                trz = exx + eyy
                S[ip, 1] = p + 2.0*model.G*(exx - trz/3.0)
                S[ip, 2] = p + 2.0*model.G*(eyy - trz/3.0)
                S[ip, 3] = p + 2.0*model.G*(-trz/3.0)
                S[ip, 4] = 2.0*model.G*exy
            end
            return S
        end
        history_to_nodes! = (out, getter) -> error("unreachable: history fields need the quadrature path")

    #-----------------------------------------------------------------------
    # PATH 2 — history at the quadrature points, Newton with the consistent
    # tangent (any model).
    #-----------------------------------------------------------------------
    else
        nh     = nhistory(model)
        H      = zeros(nh, nq, nq, nelem)        # converged history Hⁿ
        Htr    = zeros(nh, nq, nq, nelem)        # history of the current iterate
        εq     = zeros(3, nq, nq, nelem)
        σq     = zeros(4, nq, nq, nelem)
        Dq     = zeros(3, 3, nq, nq, nelem)
        lin_cache = Dict{Float64,Any}()          # linear model: tangent depends on δt only

        evaluate! = function (δt)
            solid_qp_strains!(εq, u, sem, f_el, ∇a_el, ∇b_el)
            @inbounds for iel = 1:nelem, l = 1:nq, k = 1:nq
                solid_constitutive_update!(@view(σq[:, k, l, iel]), @view(Dq[:, :, k, l, iel]),
                                           @view(Htr[:, k, l, iel]), model,
                                           εq[1, k, l, iel], εq[2, k, l, iel], εq[3, k, l, iel],
                                           @view(H[:, k, l, iel]), δt)
            end
            solid_internal_force!(fint, σq, sem, rhs_el, tmp2)
            return nothing
        end

        advance! = function (t, δt)
            loads!(t)
            for dof in cdofs                         # the prescribed values at t
                u[dof] = dirichlet[dof]
            end
            converged = false
            iters = 0
            rnorm = Inf; ref = 0.0
            for it = 1:newton_maxit + 1
                evaluate!(δt)
                r     = fext[fdofs] .- fint[fdofs]
                rnorm = norm(r)
                ref   = max(norm(fext), norm(fint), floatmin(Float64))
                if rnorm <= newton_rtol*ref
                    converged = true
                    break
                end
                it > newton_maxit && break
                iters += 1
                F = if islinear(model)
                    get!(() -> _solid_factorize(solid_tangent_matrix(Dq, sem)[fdofs, fdofs]), lin_cache, δt)
                else
                    _solid_factorize(solid_tangent_matrix(Dq, sem)[fdofs, fdofs])
                end
                u[fdofs] .+= F \ r
            end
            converged ||
                error(@sprintf(" # ERROR quasistatic_solid.jl: Newton did not converge at t = %g in %d iterations (|r|/|f| = %.3e). ",
                               t, newton_maxit, rnorm/ref),
                      "Reduce :Δt, or raise :newton_maxit. A perfectly plastic body past its limit load ",
                      "cannot carry a prescribed load: use displacement control there.")
            H .= Htr                                  # accept the step
            newton_iters[] = iters
            return nothing
        end
        nodal_stress! = function (S)
            for c = 1:4
                solid_qp_to_nodes!(@view(S[:, c]), (k, l, iel) -> σq[c, k, l, iel], sem, rhs1, tmp1)
            end
            return S
        end
        history_to_nodes! = (out, getter) ->
            solid_qp_to_nodes!(out, (k, l, iel) -> getter(@view(H[:, k, l, iel])), sem, rhs1, tmp1)
    end

    #-----------------------------------------------------------------------
    # Fields: SOLID_FIELDS and the extras, in the flat layout write_output reads.
    #-----------------------------------------------------------------------
    fields = zeros(npoin*nvar)
    Snod   = zeros(npoin, 4)
    hbuf   = zeros(npoin)
    function compute_fields!()
        nodal_stress!(Snod)
        @views fields[1:npoin]          .= u[1:npoin]
        @views fields[npoin+1:2npoin]   .= u[npoin+1:ndof]
        @views fields[2npoin+1:3npoin]  .= Snod[:, 1]
        @views fields[3npoin+1:4npoin]  .= Snod[:, 2]
        @views fields[4npoin+1:5npoin]  .= Snod[:, 4]
        @views fields[5npoin+1:6npoin]  .= Snod[:, 3]
        for (iv, name) in enumerate(extras)
            dst = @view fields[(length(SOLID_FIELDS) + iv - 1)*npoin+1:(length(SOLID_FIELDS) + iv)*npoin]
            if name == "svm"
                for ip = 1:npoin
                    sxx, syy, szz, sxy = Snod[ip, 1], Snod[ip, 2], Snod[ip, 3], Snod[ip, 4]
                    pm = (sxx + syy + szz)/3
                    dst[ip] = sqrt(1.5*((sxx - pm)^2 + (syy - pm)^2 + (szz - pm)^2 + 2.0*sxy^2))
                end
            elseif name == "umag"
                for ip = 1:npoin
                    dst[ip] = hypot(u[ip], u[npoin + ip])
                end
            else
                getter = last(hfields[findfirst(p -> first(p) == name, hfields)])
                history_to_nodes!(hbuf, getter)
                dst .= hbuf
            end
        end
        return fields
    end
    field_index(name) = findfirst(==(name), qvars)

    #-----------------------------------------------------------------------
    # Diagnostics: error norms (LGL mass-matrix norms), history file, VTK, PNG.
    #-----------------------------------------------------------------------
    qexact = zeros(npoin, nvar + 1)
    function exact_fields!(t)
        has_exact || return qexact
        for ip = 1:npoin
            ex = user_solid_exact(mesh.coords[1, ip], mesh.coords[2, ip], t)
            qexact[ip, 1] = ex.ux; qexact[ip, 2] = ex.uy
            if haskey(ex, :sxx)
                qexact[ip, 3] = ex.sxx; qexact[ip, 4] = ex.syy
                qexact[ip, 5] = ex.sxy; qexact[ip, 6] = ex.szz
            end
        end
        return qexact
    end
    function error_norms()
        compare || return (NaN, NaN, NaN, NaN)
        emax = 0.0; e2 = 0.0; n2 = 0.0; es2 = 0.0; ns2 = 0.0
        for ip = 1:npoin
            dx = fields[ip] - qexact[ip, 1]; dy = fields[npoin + ip] - qexact[ip, 2]
            emax = max(emax, hypot(dx, dy))
            e2  += M[ip]*(dx^2 + dy^2)
            n2  += M[ip]*(qexact[ip, 1]^2 + qexact[ip, 2]^2)
            ds = (fields[2npoin + ip] - qexact[ip, 3], fields[3npoin + ip] - qexact[ip, 4],
                  fields[4npoin + ip] - qexact[ip, 5], fields[5npoin + ip] - qexact[ip, 6])
            es2 += M[ip]*(ds[1]^2 + ds[2]^2 + 2.0*ds[3]^2 + ds[4]^2)
            ns2 += M[ip]*(qexact[ip, 3]^2 + qexact[ip, 4]^2 + 2.0*qexact[ip, 5]^2 + qexact[ip, 6]^2)
        end
        return (emax, sqrt(e2), sqrt(e2/max(n2, eps())), ns2 > 0 ? sqrt(es2/ns2) : NaN)
    end
    # Reaction on each Dirichlet tag: the force the constraint applies, f_int − f_ext.
    function reactions()
        R = Dict(tag => [0.0, 0.0] for tag in dtags)
        for dof in cdofs
            c = dof <= npoin ? 1 : 2
            R[dof_tag[dof]][c] += fint[dof] - fext[dof]
        end
        return R
    end

    hist_path = joinpath(OUTPUT_DIR, "solid_history.csv")
    hio = open(hist_path, "w")
    header = String["t"]
    for (name, _) in probes
        push!(header, string(name, "_ux"), string(name, "_uy"))
        compare && push!(header, string(name, "_ux_exact"), string(name, "_uy_exact"))
        for e in extras
            push!(header, string(name, "_", e))
        end
    end
    for tag in dtags
        push!(header, string("R_", tag, "_x"), string("R_", tag, "_y"))
    end
    push!(header, "newton_iters")
    compare && append!(header, ["err_u_max", "err_u_L2", "err_u_L2_rel", "err_sigma_L2_rel"])
    println(hio, join(header, ","))

    pvd_path = joinpath(OUTPUT_DIR, "simulation.pvd")
    lvtk && init_pvd_file(pvd_path)
    uaux = zeros(npoin, nvar + 1)
    iout = 0
    anim = lplotdef ? Plots.Animation() : nothing
    defname  = string(get(inputs, :plot_deformed_field, "umag"))
    defscale = Float64(get(inputs, :plot_deformed_scale, 1.0))
    lplotdef && defname != "umag" && field_index(defname) === nothing &&
        error(" # ERROR quasistatic_solid.jl: :plot_deformed_field => \"", defname,
              "\" is not among the written fields ", qvars, " (or \"umag\").")

    function record!(t, n)
        compute_fields!()
        exact_fields!(t)
        errs = error_norms()
        R    = reactions()
        row  = Float64[t]
        for (_, ip) in probes
            push!(row, u[ip], u[npoin + ip])
            compare && push!(row, qexact[ip, 1], qexact[ip, 2])
            for e in extras
                push!(row, fields[(field_index(e) - 1)*npoin + ip])
            end
        end
        for tag in dtags
            append!(row, R[tag])
        end
        push!(row, newton_iters[])
        compare && append!(row, collect(errs))
        println(hio, join((@sprintf("%.15e", v) for v in row), ","))
        flush(hio)

        (is_output(t) || n == 0 || n == nsteps) || return nothing

        s = @sprintf("   t = %-9.5g", t)
        for (name, ip) in probes
            s *= @sprintf("  %s: ux = %.8f", name, u[ip])
            compare && (s *= @sprintf(" (exact %.8f)", qexact[ip, 1]))
        end
        hpath === :quadrature && !islinear(model) && (s *= @sprintf("  [%d Newton]", newton_iters[]))
        compare && (s *= @sprintf("  |e|max = %.3e  e_L2,rel = %.3e", errs[1], errs[3]))
        println(s)

        if is_output(t)
            iout += 1
            if lvtk
                write_output(mesh.SD, fields, uaux, t, iout, mesh, nothing,
                             nothing, nothing, 0.0, 0.0, 0.0,
                             OUTPUT_DIR, inputs, qp.qvars, qp.qoutvars, inputs[:outformat];
                             nvar = nvar, qexact = qexact)
                append_pvd_entry(pvd_path, t, "iter_$(iout).pvtu")
            end
            if lplotdef
                vals = defname == "umag" ? hypot.(u[1:npoin], u[npoin+1:ndof]) :
                       fields[(field_index(defname) - 1)*npoin+1:field_index(defname)*npoin]
                plt = plot_deformed_mesh(mesh, u[1:npoin], u[npoin+1:ndof], vals, defname, t,
                                         OUTPUT_DIR, iout; scale = defscale,
                                         clims = get(inputs, :plot_deformed_clims, nothing))
                Plots.frame(anim, plt)
            end
        end
        return nothing
    end

    #-----------------------------------------------------------------------
    # March.
    #-----------------------------------------------------------------------
    try
        linitial && advance!(tinit, 0.0)          # instantaneous response to the load at tinit
        record!(tinit, 0)
        for n = 1:nsteps
            advance!(tsteps[n + 1], tsteps[n + 1] - tsteps[n])
            record!(tsteps[n + 1], n)
        end
    finally
        close(hio)
    end

    if lplotdef && length(anim.frames) > 1
        try
            Plots.gif(anim, joinpath(OUTPUT_DIR, "deformed.gif"); fps = 4, show_msg = false)
        catch e
            @warn "quasistatic_solid.jl: could not assemble deformed.gif (the PNG frames are written)" exception = e
        end
    end

    # Hand the final state back through the usual state array.
    for ivar = 1:nvar, ip = 1:npoin
        qp.qn[ip, ivar] = fields[(ivar - 1)*npoin + ip]
    end

    println("   history: ", hist_path)
    lvtk     && println("   fields:  ", pvd_path, " (", iout, " frames)")
    lplotdef && println("   deformed shape: ", joinpath(OUTPUT_DIR, "deformed-$(defname)-it<k>.png"),
                        " and deformed.gif (magnified x", defscale, ")")
    println(" # Quasi-static solid (plane strain) ................................. DONE")

    return (; u, fields, model, tsteps, probes)
end
