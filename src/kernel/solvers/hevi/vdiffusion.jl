#=============================================================================
 vdiffusion.jl -- the linear implicit VERTICAL diffusion operator D.

 WHY THIS EXISTS
 ---------------
 HEVI and the 3D IMEX both remove the acoustic step-size limit. What is left
 is whatever was second, and on a vertically refined LES mesh that is not
 advection -- it is the SGS diffusion of the vertical derivative:

     lambda_visc ~ nu / dz^2

 with dz the smallest LGL gap in the column. On the LESICP2 mesh (10 m
 elements at order 4, so dz_min = 0.1727 * 10 / 1 = 1.7 m at the surface) that
 rate crosses the acoustic one the moment nu_t reaches a few tens of m^2/s,
 which a convective boundary layer does within the first few hundred seconds.
 Before that it is invisible: the startup CFL report is taken on a laminar
 sounding where nu_t is ~0, so it says diffusion never limits dt, and it is
 right -- at t = 0.

 That is a stiffness the acoustic split cannot touch, and it is anisotropic in
 exactly the way the vertical split already is: nu/dz^2 against nu/dx^2 is
 (dx/dz)^2, a factor of 64 on that mesh. So the same column-local machinery
 that swallowed the vertical sound wave can swallow the vertical diffusion,
 and for the same reason -- it is the direction the mesh is refined in.

 WHAT IS IMPLICIT, AND WHAT IS DELIBERATELY NOT
 -----------------------------------------------
 Only the part of the viscous flux divergence that a single column can see.
 On a mesh whose zeta lines are vertical (which is what the column topology
 already requires) the physical viscous flux in z is

     tau_xz = mu (du/dz + dw/dx)        <- second term is horizontal
     tau_yz = mu (dv/dz + dw/dy)        <- second term is horizontal
     tau_zz = 2 mu dw/dz - (2/3) mu div(u)
            = (4/3) mu dw/dz - (2/3) mu (du/dx + dv/dy)
     q_z    = kappa dtheta/dz           <- entirely column-local

 and the first term of each line is column-local while the rest is not. So D
 carries

     d/dz ( mu   d/dz ) on u and v
     d/dz ( 4/3 mu d/dz ) on w
     d/dz ( kappa d/dz ) on theta

 and the cross terms stay explicit. They are not stiff: they scale as
 mu/(dx dz), which is dz/dx = 1/8 of what D carries. Dropping them costs
 nothing in stability and keeps the operator strictly column-local, which is
 what makes it bandable and what makes it share the existing LU.

 Notice what D is block-diagonal in: the column-local part of the stress
 tensor couples each variable only to itself. mass is untouched by the SGS
 closure, and u, v, w, theta each get their own scalar vertical diffusion.
 There is no cross-variable coupling to get wrong.

 CONSISTENCY IS FREE, ACCURACY OF THE COEFFICIENT IS NOT
 -------------------------------------------------------
 The ARK split is computed as a subtraction, f_exp = rhs! - f_imp, so the
 scheme is consistent for ANY D whatsoever -- a wrong D cannot lose or
 double-count physics, it can only fail to remove stiffness. That is the whole
 reason the approximations above are safe to make.

 What it does NOT excuse is a stale mu. The Jacobian left in the explicit part
 is (D_true - D), so if D is built with mu = 0 (which is what the closure
 returns on a laminar initial sounding) then the explicit part still carries
 the entire viscous rate and nothing has been gained. mu therefore has to be
 refreshed from the running solution, which is what `vdiff_refresh!` does and
 why implicit diffusion under a dynamic closure requires :PS linearisation.
 `vdiff_check_linearisation` refuses the combination that would silently do
 nothing.

 WHY THE WEAK FORM, TO THE LETTER
 --------------------------------
 The element kernel below is `_expansion_visc!` (rhs.jl, NSD_3D, ContGal)
 with the xi and eta sweeps deleted and the flux restricted to its
 column-local part. It has to be, term for term: the point of the split is
 that f_exp = rhs! - f_imp CANCELS the stiff term, and a cancellation between
 a weak-form discretisation and a strong-form one is not a cancellation. It
 would leave an O(1) fraction of the viscous rate in the explicit part -- the
 one thing this file exists to remove -- and the run would simply keep
 blowing up, with a stability report that says it should not.

 The weak form carries no surface integral here, exactly as in rhs.jl, so the
 natural boundary condition is zero viscous flux at the floor and the lid. The
 real surface flux is applied by `apply_boundary_conditions_neumann!` further
 down the explicit RHS and therefore stays entirely explicit, which is right:
 it is a forcing, not a stiff eigenvalue.

 WHAT THE MOMENTUM EQUATIONS LEAVE ON THE LATERAL FACES, AND WHY
 ---------------------------------------------------------------
 Measured, not assumed. Comparing D against the explicit viscous RHS on a
 HORIZONTALLY UNIFORM field -- where the entire flux divergence ought to be
 column-local -- the theta equation agrees to 6e-13, and so do u, v and w
 everywhere except on the domain's lateral faces, where they do not agree at
 all. That is not a defect in D. It is the stress tensor:

     u equation:  flux_x = 2 mu du/dx - (2/3) mu div(u) = -(2/3) mu dw/dz
     w equation:  flux_x = mu (du/dz + dw/dx)           =  mu du/dz

 Both are non-zero on a horizontally uniform field, so the explicit form's xi
 sweep does not vanish there. What it contributes contracts to
 sum_k psi'_i(xi_k) w_k = psi_i(1) - psi_i(-1), i.e. a term carried by the two
 xi-endpoint test functions alone -- which DSS cancels exactly between
 neighbouring elements (and around a periodic seam) and cannot cancel where
 the domain ends. Interior agreement is 4e-16 on u and v and 7e-14 on w; the
 whole difference sits on the boundary nodes.

 So on a laterally periodic LES the cancellation in f_exp is complete to
 round-off, and on a walled domain a boundary-face term of size
 mu/(dx dz) stays explicit -- dz/dx of what D carries, which on the target
 mesh is 1/8. It is not what sets the step size, and it could not be made
 implicit here anyway: it is a horizontal flux, and a horizontal flux is not
 something one column can see.

 The 4/3 is the only constant in this file derived on paper rather than
 transcribed, so it is worth saying what it is worth: with it, the interior
 agreement on the w equation is 5e-16; without it, 1e-1.
=============================================================================#

"""
    VerticalDiffusion

The frozen coefficients of the implicit vertical diffusion operator, packed to
match an operator's implicit variable set.

  `mu`  `npoin x nimp`, the diffusion coefficient actually multiplying
        `d(prim)/dz` for that variable -- already carrying the 4/3 on w and
        whatever `:mu` mask and Prandtl/Schmidt number the closure applies.
  `sc`  `npoin x nimp`, the factor converting the CONSERVATIVE deviation this
        operator acts on into the PRIMITIVE the flux is built from: `1/rho`
        for the specific variables, `1` for density itself.
  `dyn` true when the coefficient comes from a running SGS closure and so has
        to be refreshed; false for a deck-constant AV viscosity, which never
        moves and is therefore safe under :RS.

`mu` and `sc` are separate rather than premultiplied because `sc` is applied
per node BEFORE the zeta derivative and `mu` per node AFTER it. Folding them
together would differentiate `rho` along with the state.
"""
struct VerticalDiffusion{TF <: AbstractFloat}
    mu::Matrix{TF}
    sc::Matrix{TF}
    dyn::Bool
end

"""
    vdiff_enabled(inputs) -> Bool

Whether the deck asked for implicit vertical diffusion. One key serves HEVI
and the 3D IMEX because a deck runs one of them, not both, and because the
operator, the band, the LU and the refresh path are literally the same code.
"""
vdiff_enabled(inputs) = get(inputs, :implicit_vdiff, false) == true &&
                        get(inputs, :lvisc, false) == true

"""
    vdiff_vars(inputs, vars) -> Vector{Int}

Widen an implicit variable set so it can carry vertical diffusion.

Diffusion acts on u, v, w and theta, so the horizontal momenta have to join
the implicit set even on a perfectly vertical mesh where the acoustic operator
does not need them. Leaving them out would put two of the four stiff viscous
terms back in the explicit part, which is not a partial win -- the step size
is set by the largest rate that remains, and mu/dz^2 on u is the same number
as mu/dz^2 on theta.

Density joins only if the deck actually diffuses it (`:mu[1] != 0`, which is
the AV path; every SGS closure leaves mass alone).
"""
function vdiff_vars(params, inputs, vars::Vector{Int})
    vdiff_enabled(inputs) || return vars
    out = copy(vars)
    for ieq in (2, 3, 4, 5)
        ieq in out || push!(out, ieq)
    end
    if 1 ∉ out && _vdiff_diffuses_density(params)
        push!(out, 1)
    end
    return sort!(out)
end

function _vdiff_diffuses_density(params)
    vc = params.visc_coeff
    params.sgs isa AbstractSGSModel && return false     # no closure touches mass
    return length(vc) >= 1 && vc[1] != 0.0
end

"""
    build_vertical_diffusion(params, vars) -> VerticalDiffusion

Allocate the coefficient arrays for the implicit variable set `vars` and fill
them from the reference state. At setup the SGS cache is still zero (no `rhs!`
has run), so a dynamic closure starts at its molecular floor -- which is the
correct value for the laminar initial sounding, and moves as soon as the first
`vdiff_refresh!` lands.
"""
function build_vertical_diffusion(params, vars::Vector{Int})
    npoin = Int(params.mesh.npoin)
    nimp  = length(vars)
    vd = VerticalDiffusion{Float64}(zeros(Float64, npoin, nimp),
                                    ones(Float64, npoin, nimp),
                                    params.sgs isa AbstractSGSModel)
    _vdiff_fill!(vd, params, vars, params.qp.qe, nothing)
    return vd
end

"""
    vdiff_refresh!(vd, params, vars, u)

Recompute `mu` and `sc` from the running state. `u` is the flat
`npoin*neqs` state vector; the SGS eddy viscosity is read from the closure's
per-node cache, which the last `rhs!` evaluation filled -- one step stale,
against an eddy-turnover time of hundreds of steps.
"""
vdiff_refresh!(vd::VerticalDiffusion, params, vars::Vector{Int}, u) =
    _vdiff_fill!(vd, params, vars, nothing, u)

# `qe` and `u` are the two ways of naming the state to read rho from: the
# reference-state matrix (npoin x neqs) at setup, the flat solution vector at
# every refresh. Exactly one is non-nothing.
function _vdiff_fill!(vd::VerticalDiffusion, params, vars::Vector{Int}, qe, u)

    npoin = Int(params.mesh.npoin)
    neqs  = Int(params.neqs)
    sgs   = params.sgs
    vc    = params.visc_coeff
    SD    = params.SD
    mu    = vd.mu
    sc    = vd.sc

    ltheta = sgs isa AbstractSGSModel ? sgs.ltheta_eqn : true

    @inbounds for (q, ieq) in enumerate(vars)
        # The extra 4/3 on w is the (2 - 2/3) that survives when the
        # column-local part of tau_zz = 2 mu dw/dz - (2/3) mu div(u) is taken:
        # dw/dz appears in both terms. u and v keep their plain mu because
        # their vertical stress component is mu du/dz with no divergence part.
        fac = (ieq == 4 && sgs isa AbstractSGSModel) ? 4.0 / 3.0 : 1.0
        for ip = 1:npoin
            ρ = qe === nothing ? u[ip] : qe[ip, 1]
            ρ = ρ > 0 ? ρ : 1.0
            mu[ip, q] = fac * _vdiff_coeff(vc, ieq, ρ, ip, sgs, ltheta, SD)
            # rho, if it is diffused at all, is diffused as itself: the AV path
            # builds its flux from uprimitive[1] = rho, not from rho/rho.
            sc[ip, q] = ieq == 1 ? 1.0 : 1.0 / ρ
        end
    end
    return vd
end

# One call site for the coefficient so the implicit operator and the explicit
# RHS can never drift apart: under a closure this IS the function rhs.jl calls.
@inline _vdiff_coeff(vc, ieq, ρ, ip, sgs::AbstractSGSModel, ltheta, SD) =
    SGS_diffusion(vc, ieq, ρ, ip, sgs, ltheta, SD)
# The AV / constant-viscosity path: `:mu` entries are the diffusivities
# themselves, applied per equation to the gradient of that equation's
# primitive, with no stress tensor and no density weighting.
@inline _vdiff_coeff(vc, ieq, ρ, ip, ::Nothing, ltheta, SD) =
    ieq <= length(vc) ? vc[ieq] : 0.0

"""
    vdiff_check_linearisation(inputs, params, scheme)

Refuse `:implicit_vdiff` under a linearisation policy that would never update
the coefficient.

Under a dynamic closure the eddy viscosity at setup is its molecular floor, so
an operator frozen there is, for the purpose it was built for, the zero
operator. The run would pay the full cost of the wider band and get none of
the stability, and the only symptom would be a blow-up at the model time the
boundary layer becomes turbulent -- i.e. exactly the failure implicit
diffusion was switched on to prevent, arriving on schedule and looking like a
bug somewhere else.

A deck-constant AV viscosity does not move, so :RS is fine there and this says
nothing.
"""
function vdiff_check_linearisation(inputs, params, scheme::AbstractString, lin::Symbol)
    vdiff_enabled(inputs) || return nothing
    params.sgs isa AbstractSGSModel || return nothing
    lin === :PS && return nothing
    key = lowercase(scheme) == "imex3d" ? "imex" : lowercase(scheme)
    error(scheme, ": :implicit_vdiff => true needs :", key,
          "_linearization => :PS under a dynamic SGS closure ($(typeof(params.sgs))). ",
          "The implicit diffusion coefficient is read from the closure, which returns ",
          "its molecular floor on the laminar initial state; :RS would freeze it there ",
          "for the whole run and the implicit operator would carry no diffusion at all, ",
          "while still costing the wider band. Set :PS (with :", key,
          "_update_freq, default 5), or set :implicit_vdiff => false.")
end

"""
    _hevi_D_elements!(rhs_el, V, mu, sc, conn, dψ, ω, Je, dζdz, nelem, ngl, nimp)

The element loop of the implicit vertical diffusion operator, behind a typed
function barrier for the same reason as `_hevi_A_elements!` -- `metrics.dζdz`
and `mesh.connijk` are inferred `::Any` off their untyped parent structs, and
reading them inside this loop would make every index a boxed dispatch.

Accumulates INTO `rhs_el`, which `hevi_apply_A!` has already filled with the
acoustic part and has not yet DSSed. Both halves therefore go through one DSS,
one MPI assembly and one mass-matrix division -- the same tail the explicit
code applies to `RHS + RHS_visc`, which is what makes the subtraction in
`f_exp = rhs! - f_imp` cancel to round-off rather than to discretisation
error.

Indices follow `_expansion_visc!`: `m` runs over the quadrature node along
zeta, `i` over the test function that node scatters into, `(k,l)` over the
horizontal node pair that fixes the column. On a vertical-zeta mesh dzeta/dx
and dzeta/dy vanish, so the flux contracts to dzeta/dz twice -- once turning
d/dzeta into d/dz, once projecting the physical flux back onto zeta.
"""
function _hevi_D_elements!(rhs_el, V, mu, sc, conn, dψ, ω, Je, dζdz,
                           nelem::Int, ngl::Int, nimp::Int)

    @inbounds for iel = 1:nelem, l = 1:ngl, k = 1:ngl
        for q = 1:nimp
            for m = 1:ngl
                ip   = conn[iel, k, l, m]
                μ    = mu[ip, q]
                μ == 0.0 && continue
                dzz  = dζdz[iel, k, l, m]
                ωJac = ω[k] * ω[l] * ω[m] * Je[iel, k, l, m]

                # d(primitive)/dzeta at this node, from the conservative
                # deviation scaled to a primitive one node by node.
                dqdζ = 0.0
                for ii = 1:ngl
                    jp = conn[iel, k, l, ii]
                    dqdζ += dψ[ii, m] * V[jp, q] * sc[jp, q]
                end

                # mu * d/dz, projected back onto zeta and weighted.
                ∇ζflux = dzz * (μ * dzz * dqdζ) * ωJac

                for i = 1:ngl
                    rhs_el[iel, k, l, i, q] -= dψ[i, m] * ∇ζflux
                end
            end
        end
    end
    return nothing
end
