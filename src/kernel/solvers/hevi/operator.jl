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
  `wallx`    per-node flag: node sits on a free-slip boundary whose normal is
             x, where the implicit x mass flux is set to zero. Only the FULL
             operator reads it -- the vertical one carries no x flux to zero.
  `wally`    the same for y
  `V`, `R`   work arrays, npoin x nimp: operator input (deviation) and output
  `vd`       implicit vertical diffusion coefficients, or `nothing` when
             diffusion is left explicit. See vdiffusion.jl.
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
    wallx::Vector{Bool}
    wally::Vector{Bool}
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
    # `full = false` is the HEVI operator proper: only the ζ-derivative part of
    # the flux divergence, which is what makes it column-local and bandable.
    # `full = true` keeps the ξ and η sweeps as well, giving the COMPLETE
    # linear acoustic-gravity operator in all three directions. That one is not
    # column-local and is never factorised -- it exists to be APPLIED, as the
    # fast operator of an acoustic-substepping split, where the outer step sees
    # no sound at all. Same fluxes, same reference state, same wall treatment,
    # same DSS and mass division; only the derivative sweeps differ.
    full::Bool
    # element-wide flux scratch, allocated only when `full` (ngl^3 x nimp)
    F3::Array{TF,4}
    G3::Array{TF,4}
    H3::Array{TF,4}
    # Implicit vertical diffusion, or `nothing` when the deck leaves diffusion
    # explicit. It rides on the operator rather than beside it because every
    # consumer -- the matrix probing in factorize.jl, the Krylov matvec, the
    # column preconditioner -- has to see ONE operator. A diffusion term that
    # the band did not know about would turn the exact column solve into a
    # mediocre preconditioner and show up only as an iteration count.
    vd::Union{Nothing, VerticalDiffusion{TF}}
    # ADVECTIVE FORM OF THE THETA ROW.
    #
    # Flux form is  A_Theta = -Div[thetabar .* W m].
    # Advective is  A_Theta = -( thetabar .* Div[W m] + (W m) . grad(thetabar) ).
    #
    # The two are IDENTICAL in the continuum -- Div(tb*m) = tb*Div(m) + m.grad(tb)
    # -- and differ only in how the discrete operator is assembled, i.e. by the
    # aliasing of the LGL quadrature on the product. They are NOT the same
    # matrix, and the difference stays explicit in f_exp = rhs! - f_imp, so it
    # has to be measured rather than assumed harmless: see
    # test/hevi/test_theta_advective.jl.
    #
    # Why bother: with the flux form, eliminating Theta introduces Div OF THE
    # UNKNOWN momentum, so the Schur reduction does not close on one scalar and
    # rho survives it. With the advective form the coupling is POINTWISE
    # (grad(thetabar) is a known coefficient), which is exactly what lets
    # Giraldo's reduction collapse to a single pressure unknown.
    #
    # ONLY THE `full` KERNEL HONOURS THIS. `_hevi_A_elements_full!` takes the
    # flag and the gradient arrays; the vertical `_hevi_A_elements!` does not
    # and always applies the FLUX row. So on a `full = false` operator this flag
    # changes nothing about the operator -- its whole effect is to fill
    # `dtbd*` below. That is not a bug to route around: schur_precond.jl WANTS
    # exactly that, because it builds the scalar Helmholtz out of the
    # continuity and momentum rows plus grad(thetabar) and never reads the
    # Theta row at all. It is recorded here because the flag looks like it
    # should change the operator, and silently does not.
    theta_advective::Bool
    # grad(thetabar), per node, computed once at setup. Zero-length when the
    # flux form is in use.
    dtbdx::Vector{TF}
    dtbdy::Vector{TF}
    dtbdz::Vector{TF}
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
    build_hevi_operator(params, topo, vars; lwall_flux = true, full = false,
                        wallx = nothing, wally = nothing) -> HEVIOperator

`wallx` / `wally` are per-node free-slip flags for the lateral boundaries.
They are only meaningful for the FULL operator -- the vertical one carries no
horizontal mass flux to zero -- and default to "no lateral wall", which is
what a laterally periodic case wants. See `imex_lateral_walls`.
"""
function build_hevi_operator(params, topo::ColumnTopology, vars::Vector{Int};
                             lwall_flux::Bool = true, full::Bool = false,
                             wallx::Union{Nothing, AbstractVector{Bool}} = nothing,
                             wally::Union{Nothing, AbstractVector{Bool}} = nothing,
                             vdiff::Bool = false, theta_advective::Bool = false)

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

    # The full operator carries the horizontal acoustic terms, so it needs the
    # horizontal momenta in its variable set; refusing here beats producing an
    # operator that silently drops half the sound wave.
    if full && !issubset([1, 2, 3, 4, 5], vars)
        error("The full acoustic operator needs all five equations in its ",
              "variable set (it carries the horizontal pressure gradient and ",
              "mass flux); got $vars.")
    end

    # Lateral free-slip flags. Length-checked here rather than trusted: a
    # caller handing over a vector built against a different npoin would
    # otherwise silently index out of a neighbouring node's flag.
    wx = wallx === nothing ? falses(npoin) : Vector{Bool}(wallx)
    wy = wally === nothing ? falses(npoin) : Vector{Bool}(wally)
    (length(wx) == npoin && length(wy) == npoin) ||
        error("HEVI/IMEX: the lateral wall flags must be npoin = $npoin long; got ",
              "$(length(wx)) and $(length(wy)).")

    # The implicit variable set has to be able to hold what diffusion acts on
    # before the coefficients are built against it; `vdiff_vars` is what widens
    # it, and a caller that skipped that step would get an operator whose
    # diffusion silently covers only some of the stiff terms.
    vd = nothing
    if vdiff
        missing_v = filter(ieq -> ieq ∉ vars, (2, 3, 4, 5))
        isempty(missing_v) ||
            error("HEVI/IMEX: implicit vertical diffusion needs equations 2, 3, 4 and 5 ",
                  "in the implicit variable set (diffusion acts on u, v, w and theta); ",
                  "got $vars, missing $(collect(missing_v)). Pass the set through ",
                  "vdiff_vars(params, inputs, vars) before building the operator.")
        vd = build_vertical_diffusion(params, copy(vars))
    end

    z4(n) = zeros(TF, n, n, n, nimp)
    op = HEVIOperator{TF, typeof(cache)}(
        copy(vars), nimp, slot, beta, thetabar, Vector(wall), wx, wy,
        V, R, rhs_el, RHS, vaux,
        zeros(TF, ngl, nimp), zeros(TF, ngl, nimp), zeros(TF, ngl, nimp),
        cache, lwall_flux, full,
        full ? z4(ngl) : z4(0), full ? z4(ngl) : z4(0), full ? z4(ngl) : z4(0),
        vd, theta_advective,
        theta_advective ? zeros(TF, npoin) : TF[],
        theta_advective ? zeros(TF, npoin) : TF[],
        theta_advective ? zeros(TF, npoin) : TF[])
    theta_advective && hevi_fill_gradthetabar!(op, params)
    return op
end

"""
    hevi_fill_gradthetabar!(op, params)

`grad(thetabar)` at every node, for the advective Theta row.

Computed by APPLYING THE OPERATOR ITSELF rather than by a second derivative
kernel: block row 2 is `A_rho_u = -Gradx[P]` with `P = beta.*Theta` (pinned in
test_schur_blocks.jl), so feeding a state whose only non-zero is
`Theta = thetabar./beta` makes the momentum slots exactly `-grad(thetabar)`.

That reuses the metric handling, the DSS and the mass division that are
already verified, and guarantees the gradient is the SAME discrete gradient
the rest of the operator uses. A second implementation could disagree with it
in a way that no self-check would catch, because both would be internally
consistent.

The flux term must be OFF while this runs (`theta_advective` is set but the
gradient is still zero, so the advective branch would read zeros); it is called
with the operator's own Theta row inactive because the momentum rows do not
depend on it.
"""
function hevi_fill_gradthetabar!(op::HEVIOperator, params)
    npoin = length(op.thetabar)
    V = zeros(Float64, npoin, op.nimp)
    W = zeros(Float64, npoin, op.nimp)
    sθ = op.slot[5]
    @inbounds for ip = 1:npoin
        V[ip, sθ] = op.thetabar[ip] / op.beta[ip]      # so that beta*Theta == thetabar
    end
    hevi_apply_A!(W, V, params, op)
    su, sv, sw = op.slot[2], op.slot[3], op.slot[4]
    @inbounds for ip = 1:npoin
        op.dtbdx[ip] = -W[ip, su]
        op.dtbdy[ip] = -W[ip, sv]
        op.dtbdz[ip] = -W[ip, sw]
    end
    return op
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
    _hevi_A_elements_full!(rhs_el, V, beta, thetabar, wall, wallx, wally, F, G, H,
                           conn, dψ, ω, Je,
                           dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                           nelem, ngl, nimp, slot, g)

The element loop of the COMPLETE linear acoustic-gravity operator: the same
linearised fluxes as `_hevi_A_elements!`, but contracted with all three
derivative sweeps instead of ζ alone.

It is the term-for-term linearisation of `_expansion_inviscid!` (NSD_3D,
ContGal) restricted to the acoustic subsystem, which is what it has to be --
a fast operator that is not the exact acoustic part of the explicit
discretisation leaves a residual fast mode in the slow part, and the outer
step stays pinned by sound with nothing in the output to say why.

THE LATERAL WALLS ARE PART OF THE OPERATOR, NOT A DETAIL
--------------------------------------------------------
At a free-slip boundary the normal mass flux vanishes, so the normal
component of the linearised momentum must not enter the continuity or the
potential-temperature flux -- exactly the condition the vertical kernel
already imposes at the floor and the lid through `wall`. The pressure term in
the *momentum* flux stays: `p` at a wall is not zero, it is what pushes back.

The vertical kernel never needed the horizontal version because it carries no
horizontal flux. This one does, and leaving it out is not a small error on a
walled domain: the implicit half would let sound run out through the side
walls, the difference would land in `f_exp = rhs! - f_imp` as a stiff
residual at exactly the nodes where `rhs!` projects the normal momentum out,
and the horizontal acoustic step-size limit this whole scheme exists to
remove would come back along the boundary with nothing in the output to say
so.

Behind the same typed function barrier as the vertical kernel, and for the
same reason (see the comment in `hevi_apply_A!`).
"""
function _hevi_A_elements_full!(rhs_el, V, beta, thetabar, wall, wallx, wally, F, G, H,
                                conn, dψ, ω, Je,
                                dξdx, dξdy, dξdz,
                                dηdx, dηdy, dηdz,
                                dζdx, dζdy, dζdz,
                                nelem::Int, ngl::Int, nimp::Int,
                                slot::NTuple{5,Int}, g::Float64,
                                theta_adv::Bool, dtbdx, dtbdy, dtbdz)

    sρ, sρu, sρv, sρw, sρθ = slot

    @inbounds for iel = 1:nelem

        # --- linearised acoustic fluxes over the whole element ---------------
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip  = conn[iel, i, j, k]
            β   = beta[ip]
            θb  = thetabar[ip]
            dρθ = V[ip, sρθ]
            dρu = V[ip, sρu]
            dρv = V[ip, sρv]
            dρw = V[ip, sρw]
            # Free-slip: the normal MASS flux vanishes, the pressure term in
            # the normal MOMENTUM flux does not.
            mx = wallx[ip] ? 0.0 : dρu
            my = wally[ip] ? 0.0 : dρv
            mz = wall[ip]  ? 0.0 : dρw

            for q = 1:nimp
                F[i, j, k, q] = 0.0; G[i, j, k, q] = 0.0; H[i, j, k, q] = 0.0
            end
            # x flux
            F[i, j, k, sρ]  = mx
            F[i, j, k, sρu] = β * dρθ
            F[i, j, k, sρθ] = θb * mx
            # y flux
            G[i, j, k, sρ]  = my
            G[i, j, k, sρv] = β * dρθ
            G[i, j, k, sρθ] = θb * my
            # z flux -- the vertical mass and heat fluxes vanish at floor/lid
            H[i, j, k, sρ]  = mz
            H[i, j, k, sρw] = β * dρθ
            H[i, j, k, sρθ] = θb * mz
        end

        # --- all three derivative sweeps, plus the buoyancy source -----------
        for q = 1:nimp, k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip   = conn[iel, i, j, k]
            ωJac = ω[i] * ω[j] * ω[k] * Je[iel, i, j, k]

            dFξ = 0.0; dFη = 0.0; dFζ = 0.0
            dGξ = 0.0; dGη = 0.0; dGζ = 0.0
            dHξ = 0.0; dHη = 0.0; dHζ = 0.0
            for m = 1:ngl
                dξm = dψ[m, i]; dηm = dψ[m, j]; dζm = dψ[m, k]
                dFξ += dξm * F[m, j, k, q]
                dFη += dηm * F[i, m, k, q]
                dFζ += dζm * F[i, j, m, q]
                dGξ += dξm * G[m, j, k, q]
                dGη += dηm * G[i, m, k, q]
                dGζ += dζm * G[i, j, m, q]
                dHξ += dξm * H[m, j, k, q]
                dHη += dηm * H[i, m, k, q]
                dHζ += dζm * H[i, j, m, q]
            end

            # NOT WORTH SKIPPING THE ZERO FLUXES. F is zero for rho_v/rho_w, G
            # for rho_u/rho_w, H for rho_u/rho_v -- 18 of these 45 chains per
            # node, 40% of the arithmetic, differentiate identically zero data.
            # Branching them out is bitwise-identical and measured 1.03x on a
            # 4x4x20 p=4 mesh: the loop is bound by the strided F[i,m,k,q] and
            # F[i,j,m,q] reads, not by the multiply-adds, so removing flops
            # buys nothing and the branches cost clarity. Left as it is
            # deliberately; the way to speed this up is less memory traffic
            # (fewer implicit fields), not fewer operations.
            dFdx = dFξ * dξdx[iel,i,j,k] + dFη * dηdx[iel,i,j,k] + dFζ * dζdx[iel,i,j,k]
            dGdy = dGξ * dξdy[iel,i,j,k] + dGη * dηdy[iel,i,j,k] + dGζ * dζdy[iel,i,j,k]
            dHdz = dHξ * dξdz[iel,i,j,k] + dHη * dηdz[iel,i,j,k] + dHζ * dζdz[iel,i,j,k]

            S = (q == sρw) ? -g * V[ip, sρ] : 0.0
            rhs_el[iel, i, j, k, q] -= ωJac * ((dFdx + dGdy + dHdz) - S)
        end

        # --- ADVECTIVE THETA ROW, replacing the flux one --------------------
        #
        #   flux:       A_Theta = -Div[ tb .* W m ]
        #   advective:  A_Theta = -( tb .* Div[W m] + (W m) . grad(tb) )
        #
        # The second sweep is over the MASS fluxes (F[..,srho] = mx etc.), so
        # `divm` here is the same divergence the continuity row takes -- which
        # is the point: it makes the Theta row a POINTWISE multiple of the
        # continuity row plus a pointwise advection, and that is what lets the
        # Schur elimination close on one scalar.
        #
        # tb multiplies BEFORE the DSS and the mass division, which is exact:
        # thetabar is a continuous nodal field, so it takes the same value in
        # every element contribution at a node and commutes with the assembly
        # sum. Same identity that makes the buoyancy source exactly -g*rho
        # (verified to 1.8e-16 in test_schur_blocks.jl).
        if theta_adv
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip   = conn[iel, i, j, k]
            ωJac = ω[i] * ω[j] * ω[k] * Je[iel, i, j, k]
            dFξ = 0.0; dFη = 0.0; dFζ = 0.0
            dGξ = 0.0; dGη = 0.0; dGζ = 0.0
            dHξ = 0.0; dHη = 0.0; dHζ = 0.0
            for m = 1:ngl
                dξm = dψ[m, i]; dηm = dψ[m, j]; dζm = dψ[m, k]
                dFξ += dξm * F[m, j, k, sρ]; dFη += dηm * F[i, m, k, sρ]; dFζ += dζm * F[i, j, m, sρ]
                dGξ += dξm * G[m, j, k, sρ]; dGη += dηm * G[i, m, k, sρ]; dGζ += dζm * G[i, j, m, sρ]
                dHξ += dξm * H[m, j, k, sρ]; dHη += dηm * H[i, m, k, sρ]; dHζ += dζm * H[i, j, m, sρ]
            end
            divm = (dFξ * dξdx[iel,i,j,k] + dFη * dηdx[iel,i,j,k] + dFζ * dζdx[iel,i,j,k]) +
                   (dGξ * dξdy[iel,i,j,k] + dGη * dηdy[iel,i,j,k] + dGζ * dζdy[iel,i,j,k]) +
                   (dHξ * dξdz[iel,i,j,k] + dHη * dηdz[iel,i,j,k] + dHζ * dζdz[iel,i,j,k])
            adv = F[i,j,k,sρ]*dtbdx[ip] + G[i,j,k,sρ]*dtbdy[ip] + H[i,j,k,sρ]*dtbdz[ip]
            # overwrite what the flux sweep just wrote into the Theta row
            rhs_el[iel, i, j, k, sρθ] = -ωJac * (thetabar[ip] * divm + adv)
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
    if op.full
        _hevi_A_elements_full!(rhs_el, V, op.beta, op.thetabar,
                               op.wall, op.wallx, op.wally,
                               op.F3, op.G3, op.H3,
                               mesh.connijk, params.basis.dψ, params.ω, met.Je,
                               met.dξdx, met.dξdy, met.dξdz,
                               met.dηdx, met.dηdy, met.dηdz,
                               met.dζdx, met.dζdy, met.dζdz,
                               nelem, ngl, nimp, op.slot, PhysicalConst{Float64}().g,
                               op.theta_advective, op.dtbdx, op.dtbdy, op.dtbdz)
    else
        _hevi_A_elements!(rhs_el, V, op.beta, op.thetabar, op.wall,
                          op.Fv, op.Gv, op.Hv,
                          mesh.connijk, params.basis.dψ, params.ω,
                          met.Je, met.dζdx, met.dζdy, met.dζdz,
                          nelem, ngl, nimp, op.slot, PhysicalConst{Float64}().g)
    end

    # Implicit vertical diffusion accumulates into the SAME rhs_el, before the
    # single DSS / mass division below, so that the two halves of the implicit
    # operator go through exactly the tail the explicit code applies to
    # RHS + RHS_visc. Same barrier discipline as the acoustic kernels above.
    if op.vd !== nothing
        vd = op.vd
        _hevi_D_elements!(rhs_el, V, vd.mu, vd.sc,
                          mesh.connijk, params.basis.dψ, params.ω,
                          met.Je, met.dζdz, nelem, ngl, nimp)
    end

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
