#=============================================================================
 operator.jl -- the linear vertical acoustic-gravity operator A, and its
 assembly into one banded matrix per vertical column.

 WHAT IS TREATED IMPLICITLY
 --------------------------
 Only the terms that carry the vertical sound wave, linearised about the
 reference state qe:

     d(dρ)/dt   = -d/dz [ d(ρw) ]
     d(dρw)/dt  = -d/dz [ beta * d(ρθ) ]  -  g * dρ
     d(dρθ)/dt  = -d/dz [ θbar * d(ρw) ]

 with beta = dp/d(ρθ) evaluated on the reference state. Since p = C0 (ρθ)^γ,
 beta = γ p / (ρθ), which is c²/θ. Vertical *advection* stays explicit: on
 this grid |w|/Δz is ~1/800 of c/Δz, so making it implicit would buy nothing
 and would cost a refactorisation every step.

 The operator is affine-free by construction: it acts on the deviation
 (u - qe), so f_imp(qe) = 0 exactly. The reference state is therefore an
 exact steady state of the implicit part, and whatever discrete hydrostatic
 imbalance the code already has stays entirely inside the explicit part,
 untouched. Nothing about the split can create or destroy balance.

 The explicit part is then, by definition,

     f_exp(u) = rhs!(u) - f_imp(u)

 computed as a subtraction rather than as a second, separately-written RHS.
 That costs one extra (cheap) A-application per stage and buys the guarantee
 that no physics can be lost or double-counted in the split -- every source,
 every boundary condition, the SGS closure, the wall model and the sponge
 all remain exactly as the explicit code computes them.

 WHY THE MATRIX IS ASSEMBLED BY PROBING
 --------------------------------------
 The banded matrix is not derived on paper and transcribed. It is extracted
 by applying the operator to a set of structured probe vectors, which is
 possible because A couples a node only to nodes within +-(ngl-1) levels of
 it in the same column. Colouring the levels modulo S = 2(ngl-1)+1 makes each
 probe pick out exactly one band entry per row, so S * nimp applications of
 the (cheap) operator recover the matrix exactly.

 This matters more than it looks. The assembled matrix is then guaranteed to
 be the matrix OF THE OPERATOR THAT IS ACTUALLY APPLIED -- including the DSS
 sum, the inverse mass matrix, the metric terms on a warped mesh and the
 boundary flux treatment. A hand-derived matrix that disagrees with the
 operator in any of those would not be caught by a convergence test; it would
 quietly degrade the implicit solve into a mediocre preconditioner and show
 up only as an unexplained step-size limit.
=============================================================================#

"""
    HEVIOperator

Everything the vertical implicit operator needs, built once at setup.

  `vars`     which solution components are implicit: `(1,4,5)` on a mesh whose
             ζ lines are vertical, `(1,2,3,4,5)` when the mesh is warped and
             dζ/dx, dζ/dy do not vanish
  `beta`     dp/d(ρθ) on the reference state, per node
  `thetabar` reference θ, per node
  `wtop`     per-node flag: node sits on the domain floor or lid, where the
             implicit vertical mass flux is set to zero
  `V`, `R`   work arrays, npoin x nimp: operator input (deviation) and output
"""
struct HEVIOperator{TF <: AbstractFloat, CT}
    vars::Vector{Int}
    nimp::Int
    # slot[ieq] is the position of equation ieq inside the packed nimp-wide
    # layout, or 0 when that equation is not implicit. Resolved once, here,
    # rather than by findfirst inside hevi_apply_A! -- findfirst returns
    # Union{Nothing,Int}, which makes the whole inner loop type-unstable and
    # would box on every node.
    slot::NTuple{5, Int}
    beta::Vector{TF}
    thetabar::Vector{TF}
    wall::Vector{Bool}
    V::Matrix{TF}
    R::Matrix{TF}
    rhs_el::Array{TF,5}
    RHS::Matrix{TF}
    vaux::Vector{TF}
    # per-ζ-line flux scratch, hoisted out of the element loop
    Fv::Matrix{TF}
    Gv::Matrix{TF}
    Hv::Matrix{TF}
    dss_cache::CT
    lwall_flux::Bool
end

"""
    hevi_choose_vars(metrics, nelem, ngl, comm) -> Vector{Int}

Decide whether the horizontal momenta have to join the implicit set.

The ζ-direction flux divergence is `dF/dζ dζ/dx + dG/dζ dζ/dy + dH/dζ dζ/dz`.
On a mesh whose ζ lines are truly vertical the first two vanish and the sound
wave along a column only involves (ρ, ρw, ρθ). On a terrain-following mesh
they do not, the ζ-direction acoustic operator picks up ρu and ρv, and
leaving them out would make the implicit operator an approximation rather
than the exact Jacobian of the split.

Deciding this from the metrics rather than from `:lwarp` means a mesh that is
warped by any other route still gets the right operator.
"""
function hevi_choose_vars(metrics, comm)
    mx = maximum(abs, metrics.dζdx)
    my = maximum(abs, metrics.dζdy)
    mz = maximum(abs, metrics.dζdz)
    if MPI.Comm_size(comm) > 1
        mx = MPI.Allreduce(mx, MPI.MAX, comm)
        my = MPI.Allreduce(my, MPI.MAX, comm)
        mz = MPI.Allreduce(mz, MPI.MAX, comm)
    end
    return (max(mx, my) <= 1.0e-10 * max(mz, eps())) ? [1, 4, 5] : [1, 2, 3, 4, 5]
end

"""
    build_hevi_operator(params, topo, vars; lwall_flux = true) -> HEVIOperator
"""
function build_hevi_operator(params, topo::ColumnTopology, vars::Vector{Int};
                             lwall_flux::Bool = true)

    mesh  = params.mesh
    npoin = Int(mesh.npoin)
    nelem = Int(mesh.nelem)
    ngl   = Int(mesh.ngl)
    nimp  = length(vars)
    TF    = Float64

    # ρ, ρw and ρθ are what carry the vertical sound wave; without all three
    # the operator is not the acoustic subsystem and the buoyancy coupling has
    # nowhere to read δρ from.
    issubset([1, 4, 5], vars) ||
        error("HEVI: the implicit variable set must contain 1 (ρ), 4 (ρw) and 5 (ρθ); ",
              "got $vars. Check :hevi_vars in the deck.")
    all(v -> 1 <= v <= 5, vars) ||
        error("HEVI: :hevi_vars entries must be equation indices in 1:5; got $vars.")
    PhysConst = PhysicalConst{TF}()
    γ = PhysConst.γ

    qe = params.qp.qe
    beta     = zeros(TF, npoin)
    thetabar = zeros(TF, npoin)
    @inbounds for ip = 1:npoin
        ρ  = qe[ip, 1]
        ρθ = qe[ip, 5]
        if ρ > 0 && ρθ > 0
            p = perfectGasLaw_ρθtoP(PhysConst, ρ, ρθ / ρ)
            beta[ip]     = γ * p / ρθ         # dp/d(ρθ) = c^2/θ
            thetabar[ip] = ρθ / ρ
        end
    end

    # Floor and lid nodes, from the column topology: level 1 is the global
    # ground, level nlev the global lid. Everything in between is interior
    # even where it sits on a rank boundary.
    wall = falses(npoin)
    if lwall_flux
        @inbounds for ic = 1:topo.ncol
            for il in (1, topo.nlev)
                ip = topo.node[il, ic]
                ip != 0 && (wall[ip] = true)
            end
        end
    end

    V   = zeros(TF, npoin, nimp)
    R   = zeros(TF, npoin, nimp)
    rhs_el = zeros(TF, nelem, ngl, ngl, ngl, nimp)
    RHS = zeros(TF, npoin, nimp)
    vaux = zeros(TF, npoin)
    cache = setup_assembler(mesh.SD, RHS, mesh.ip2gip, mesh.gip2owner)

    slot = ntuple(ieq -> something(findfirst(==(ieq), vars), 0), 5)

    return HEVIOperator{TF, typeof(cache)}(
        copy(vars), nimp, slot, beta, thetabar, Vector(wall),
        V, R, rhs_el, RHS, vaux,
        zeros(TF, ngl, nimp), zeros(TF, ngl, nimp), zeros(TF, ngl, nimp),
        cache, lwall_flux)
end


"""
    _hevi_A_elements!(rhs_el, V, beta, thetabar, wall, Fv, Gv, Hv,
                      conn, dψ, ω, Je, dζdx, dζdy, dζdz,
                      nelem, ngl, nimp, slot, g)

The element loop of the vertical acoustic operator, behind a typed function
barrier. Every array arrives as a positional argument so this specialises on
its concrete type -- see the comment at the call site for why that matters
enough to be worth a separate function.
"""
function _hevi_A_elements!(rhs_el, V, beta, thetabar, wall, Fv, Gv, Hv,
                           conn, dψ, ω, Je, dζdx, dζdy, dζdz,
                           nelem::Int, ngl::Int, nimp::Int,
                           slot::NTuple{5,Int}, g::Float64)

    sρ, sρu, sρv, sρw, sρθ = slot
    have_h = sρu != 0

    @inbounds for iel = 1:nelem, j = 1:ngl, i = 1:ngl

        # --- linearised fluxes along this ζ line -----------------------------
        for m = 1:ngl
            ip  = conn[iel, i, j, m]
            β   = beta[ip]
            θb  = thetabar[ip]
            dρθ = V[ip, sρθ]
            dρw = V[ip, sρw]
            noflux = wall[ip]

            for q = 1:nimp
                Fv[m, q] = 0.0; Gv[m, q] = 0.0; Hv[m, q] = 0.0
            end
            # vertical: mass flux, pressure gradient, potential-temperature flux
            Hv[m, sρ]  = noflux ? 0.0 : dρw
            Hv[m, sρw] = β * dρθ
            Hv[m, sρθ] = noflux ? 0.0 : θb * dρw
            if have_h
                dρu = V[ip, sρu]; dρv = V[ip, sρv]
                Fv[m, sρ]  = dρu;  Fv[m, sρu] = β * dρθ;  Fv[m, sρθ] = θb * dρu
                Gv[m, sρ]  = dρv;  Gv[m, sρv] = β * dρθ;  Gv[m, sρθ] = θb * dρv
            end
        end

        # --- ζ derivative and the buoyancy source ----------------------------
        for k = 1:ngl
            ip   = conn[iel, i, j, k]
            ωJac = ω[i] * ω[j] * ω[k] * Je[iel, i, j, k]
            dzx  = have_h ? dζdx[iel, i, j, k] : 0.0
            dzy  = have_h ? dζdy[iel, i, j, k] : 0.0
            dzz  = dζdz[iel, i, j, k]

            for q = 1:nimp
                dF = 0.0; dG = 0.0; dH = 0.0
                for m = 1:ngl
                    d   = dψ[m, k]
                    dF += d * Fv[m, q]
                    dG += d * Gv[m, q]
                    dH += d * Hv[m, q]
                end
                div = dF * dzx + dG * dzy + dH * dzz
                S   = (q == sρw) ? -g * V[ip, sρ] : 0.0
                rhs_el[iel, i, j, k, q] -= ωJac * (div - S)
            end
        end
    end
    return nothing
end

"""
    hevi_apply_A!(out, V, params, op)

Apply the linear vertical acoustic operator to the deviation field `V`
(`npoin x nimp`, columns ordered as `op.vars`) and write the result into
`out` (same shape).

The discretisation mirrors `_expansion_inviscid!` term for term -- strong
form, collocation differentiation along ζ, `-ωJ (div - S)` accumulated into
the element RHS, then DSS, the MPI assembly and the inverse mass matrix. It
has to: an implicit operator that is not the exact ζ part of the explicit
discretisation would leave a residual fast mode behind, and the step size
would stay pinned for no visible reason.
"""
function hevi_apply_A!(out::AbstractMatrix, V::AbstractMatrix, params, op::HEVIOperator)

    mesh  = params.mesh
    nelem = Int(mesh.nelem)
    ngl   = Int(mesh.ngl)
    npoin = Int(mesh.npoin)
    nimp  = op.nimp
    met   = params.metrics

    rhs_el = op.rhs_el
    fill!(rhs_el, 0.0)

    #-------------------------------------------------------------------------
    # TYPED FUNCTION BARRIER -- not optional, and not a micro-optimisation.
    #
    # St_metrics and St_mesh are `Base.@kwdef mutable struct`s whose fields
    # carry NO type annotation, so `metrics.dζdz` and `mesh.connijk` are
    # inferred `::Any`. Reading them inside the element loop makes every single
    # index a boxed dynamic dispatch. Hoisting them into local variables does
    # not help: the local is still `::Any`.
    #
    # Measured on the rtb_hevi case: 0.0929 s per application, 3.3x the cost of
    # a full nonlinear rhs! evaluation, and 60x what the same code costs
    # against concretely-typed structs. It made a HEVI step 3.7x an explicit
    # one and turned a 2.5x win into a 1.5x loss.
    #
    # Passing the arrays positionally lets Julia specialise this callee on
    # their concrete runtime types, and the loop compiles as it should. Same
    # device rhs.jl uses for _inviscid_rhs_el_3d! -- see the comment there.
    #-------------------------------------------------------------------------
    _hevi_A_elements!(rhs_el, V, op.beta, op.thetabar, op.wall,
                      op.Fv, op.Gv, op.Hv,
                      mesh.connijk, params.basis.dψ, params.ω,
                      met.Je, met.dζdx, met.dζdy, met.dζdz,
                      nelem, ngl, nimp, op.slot, PhysicalConst{Float64}().g)

    RHS = op.RHS
    fill!(RHS, 0.0)
    DSS_rhs!(RHS, rhs_el, mesh.connijk, nelem, ngl, nimp, params.SD, params.AD)
    op.dss_cache !== nothing && assemble_mpi!(RHS, op.dss_cache)

    for q = 1:nimp
        divide_by_mass_matrix!(@view(RHS[:, q]), op.vaux, params.Minv, nimp, npoin, params.AD)
    end
    out === RHS || copyto!(out, RHS)
    return out
end

"""
    hevi_implicit_rhs!(du, u, params, hevi)

`f_imp` in the ARK sense: the implicit part evaluated at the current state,
written into the flat `du` layout the rest of the code uses.

`du` is zero in every equation that is not implicit, which is what makes
`f_exp = rhs! - f_imp` leave those equations exactly as the explicit code
computed them.
"""
function hevi_implicit_rhs!(du, u, params, hevi)

    op    = hevi.op
    npoin = Int(params.mesh.npoin)
    neqs  = Int(params.neqs)
    qe    = params.qp.qe

    V = op.V
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            V[ip, q] = u[off + ip] - qe[ip, ieq]
        end
    end

    hevi_apply_A!(op.R, V, params, op)

    fill!(du, 0.0)
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            du[off + ip] = op.R[ip, q]
        end
    end
    return du
end
