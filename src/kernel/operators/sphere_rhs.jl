#---------------------------------------------------------------------------------
# sphere_rhs.jl
#
# The right-hand side of the shallow water equations on the spherical shell,
# assembled by the SPECTRAL ELEMENT METHOD — no finite differences anywhere.
#
# Continuous Galerkin, strong form, exactly the convention of
# _expansion_inviscid!(…, NSD_2D, ContGal) in src/kernel/operators/rhs.jl, with
# the third flux component and the third metric direction added:
#
#   rhs_el[e,i,j,eq] -= ω_i ω_j J [ (∂F/∂x + ∂G/∂y + ∂H/∂z) - S ]
#
#   ∂F/∂x = dξdx ∂F/∂ξ + dηdx ∂F/∂η          (and likewise for G,y and H,z)
#   ∂F/∂ξ = Σ_k dψ[k,i] F[k,j]               the LGL differentiation matrix
#
# then DIRECT STIFFNESS SUMMATION (the `+=` into the global array, which is what
# makes the solution continuous across elements) and a division by the diagonal
# mass matrix:
#
#   ∂q/∂t = M⁻¹ Σ_e rhs_el .
#
# There is NO boundary term. The shell is a closed manifold: ∫_Γ vanishes
# identically, as the paper notes for the CG case (Marras, Kopera & Giraldo
# 2015, section 4.1, below Eq. 14). The "periodicity" is the topology itself,
# enforced by the node numbering built in sphere_mesh.jl.
#
# The fluxes come from the case's user_flux!, the sources from its user_source!,
# through the ordinary Jexpresso callbacks — this file never hard-codes an
# equation.
#
# ALSO HERE: the ARTIFICIAL DIFFUSION δν∇²q of Eq. (8b), assembled in weak form
# with the SURFACE (Laplace-Beltrami) gradient. See _sphere_visc! below for why
# the flat viscous path of rhs.jl cannot supply it.
#
# ν is either the CONSTANT of the deck's :μ, or — with :visc_model => DSGS_SW()
# — the residual-based Marras-Nazarov DynSGS coefficient, recomputed per element
# and per RK stage from how badly the discrete solution fails to satisfy the
# equations. The model itself is compute_dsgs_viscosity!(::DSGS_SW, ::NSD_2D) in
# src/kernel/physics/SGS.jl; what lives here is the plumbing — the per-element Δ,
# the second element loop the residual forces (see the DynSGS block in
# _sphere_rhs_kernel! for why it cannot be folded into the first), and the
# packaging of ν|e for output. DSGS.md §5 is the write-up; the case is
# problems/ShallowWater/SWsphereDSGS.
#
# ALSO HERE: the modal filter. The paper applies a Boyd-Vandeven filter at every
# step, "chosen to reduce by 5% the highest modes only", because the inviscid
# high-order solution of an advective problem grows spurious grid-scale modes.
# What is implemented is the exponential modal filter, which serves the same
# purpose and has the same knobs; see sphere_filter! below.
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

export St_sphere_params
export build_sphere_params
export sphere_rhs!
export sphere_filter!
export build_sphere_filter_matrix
export sphere_relative_vorticity!
export sphere_dsgs_requested
export build_sphere_element_size
export sphere_dsgs_extra


#---------------------------------------------------------------------------------
# Scratch space, allocated once and reused every stage of every step.
#---------------------------------------------------------------------------------
mutable struct St_sphere_params{TFloat}
    neqs ::Int
    F    ::Array{TFloat, 3}    # ngl × ngl × neqs
    G    ::Array{TFloat, 3}
    H    ::Array{TFloat, 3}
    S    ::Array{TFloat, 3}
    Filt ::Array{TFloat, 2}    # ngl × ngl modal filter (identity when disabled)
    qf   ::Array{TFloat, 3}    # ngl × ngl × neqs filter scratch
    acc  ::Array{TFloat, 2}    # npoin × neqs   filter accumulator
    lfilter ::Bool
    μ    ::Array{TFloat, 1}    # neqs : artificial viscosity, per equation (0 = off)
    lvisc   ::Bool             # any μ[ieq] > 0
    gξ   ::Array{TFloat, 3}    # ngl × ngl × neqs  ν ∇ₛq · a¹  at the LGL nodes
    gη   ::Array{TFloat, 3}    # ngl × ngl × neqs  ν ∇ₛq · a²

    #--- DynSGS (:visc_model => DSGS_SW()). All sized 1 and unused when off.
    ldsgs   ::Bool             # residual-based ν, recomputed every RK stage
    Δelem   ::Array{TFloat, 1} # nelem : shortest element chord, the Δ of the model
    μ_dsgs  ::Array{TFloat, 2} # nelem × neqs : ν|e per equation, this stage
    qnm1    ::Array{TFloat, 2} # npoin × neqs : qⁿ⁻² of the BDF2 residual
    qnm2    ::Array{TFloat, 2} # npoin × neqs : qⁿ⁻¹  (see the roll in sphere_time_loop.jl)
    davg    ::Array{TFloat, 1} # neqs : ⟨q_i⟩ scratch
    ddenom  ::Array{TFloat, 1} # neqs : ‖q_i − ⟨q_i⟩‖_{∞,Ω} scratch
    νel     ::Array{TFloat, 1} # neqs : μ_dsgs[iel,:] unpacked, handed to _sphere_visc_el!
    C1      ::TFloat           # :dsgs_C1, the residual coefficient
    C2      ::TFloat           # :dsgs_C2, the wave-speed cap coefficient
    Δt      ::Base.RefValue{Float64}    # the step, filled in by sphere_time_loop!
    thist   ::Base.RefValue{Float64}    # when the BDF2 history was last rolled
end


function build_sphere_params(mesh::St_mesh, metrics::St_sphere_metrics, inputs;
                             neqs = 4, TF = TFloat)

    ngl   = Int(mesh.ngl)
    nelem = Int(mesh.nelem)
    npoin = Int(mesh.npoin)

    lfilter = get(inputs, :lfilter, false) == true
    Filt    = lfilter ? build_sphere_filter_matrix(metrics.ξ, inputs, TF) :
                        Matrix{TF}(I, ngl, ngl)

    μ     = build_sphere_viscosity(inputs, neqs, TF)
    lvisc = any(>(zero(TF)), μ)
    # `&& lvisc`: a DSGS_SW deck whose multipliers are all zero has asked for a
    # model whose output is discarded. Turning it off here skips the residual
    # pass rather than computing it and multiplying it by nothing.
    ldsgs = sphere_dsgs_requested(inputs) && lvisc

    #
    # DynSGS carries a per-element Δ and a residual history; a constant-ν run
    # needs neither, so everything below is sized 1 when the model is off.
    #
    ndsgs_e = ldsgs ? nelem : 1
    ndsgs_p = ldsgs ? npoin : 1
    Δelem   = ldsgs ? build_sphere_element_size(mesh, TF) : zeros(TF, 1)

    return St_sphere_params{TF}(neqs,
                                zeros(TF, ngl, ngl, neqs),
                                zeros(TF, ngl, ngl, neqs),
                                zeros(TF, ngl, ngl, neqs),
                                zeros(TF, ngl, ngl, neqs),
                                Filt,
                                zeros(TF, ngl, ngl, neqs),
                                zeros(TF, npoin, neqs),
                                lfilter,
                                μ,
                                lvisc,
                                zeros(TF, ngl, ngl, neqs),
                                zeros(TF, ngl, ngl, neqs),
                                ldsgs,
                                Δelem,
                                zeros(TF, ndsgs_e, neqs),
                                zeros(TF, ndsgs_p, neqs),
                                zeros(TF, ndsgs_p, neqs),
                                zeros(TF, neqs),
                                zeros(TF, neqs),
                                zeros(TF, neqs),
                                TF(get(inputs, :dsgs_C1, 1.0)),
                                TF(get(inputs, :dsgs_C2, 0.5)),
                                Ref(0.0),
                                Ref(-1.0e30))
end


#---------------------------------------------------------------------------------
# sphere_dsgs_requested(inputs) -> Bool
#
# `:visc_model => DSGS_SW()` together with `:lvisc => true`. The model is a
# viscosity, so it is off whenever the master viscosity switch is off — the same
# rule build_sphere_viscosity applies to the constant ν.
#---------------------------------------------------------------------------------
function sphere_dsgs_requested(inputs)
    get(inputs, :lvisc, false) == true || return false
    return get(inputs, :visc_model, nothing) isa DSGS_SW
end


#---------------------------------------------------------------------------------
# build_sphere_element_size(mesh, TF) -> Vector, one entry per element
#
# The Δ_elem of the DynSGS model: the SHORTEST of the two element chords,
#
#   d_ξ = |x(ngl,1) − x(1,1)| ,   d_η = |x(1,ngl) − x(1,1)| ,
#
# which is the manifold counterpart of mesh.Δelem in the flat cases (min corner-
# to-corner distance inside the element). Chords, not arcs: the ratio of the two
# is 1 − O((Δ/r)²) ≈ 1 − 1e-4 on a cubed sphere with elements of a few hundred
# km, i.e. far below the modelling uncertainty in C1 itself.
#
# The model then uses Δ = Δ_elem/(N+1), the per-node effective resolution — done
# at the point of use, exactly as the 1D/2D/MHD paths do.
#---------------------------------------------------------------------------------
function build_sphere_element_size(mesh::St_mesh, TF = TFloat)

    nelem = Int(mesh.nelem)
    ngl   = Int(mesh.ngl)
    crd   = mesh.coords
    conn  = mesh.connijk

    Δelem = zeros(TF, nelem)

    @inbounds for iel = 1:nelem
        i00 = conn[iel, 1,   1  ]
        i10 = conn[iel, ngl, 1  ]
        i01 = conn[iel, 1,   ngl]

        dξ = sqrt((crd[1,i10] - crd[1,i00])^2 +
                  (crd[2,i10] - crd[2,i00])^2 +
                  (crd[3,i10] - crd[3,i00])^2)
        dη = sqrt((crd[1,i01] - crd[1,i00])^2 +
                  (crd[2,i01] - crd[2,i00])^2 +
                  (crd[3,i01] - crd[3,i00])^2)

        Δelem[iel] = TF(min(dξ, dη))
    end

    minimum(Δelem) > 0 ||
        error(" # ERROR sphere_rhs.jl: an element of the shell grid has zero size; DynSGS cannot build Δ.")

    return Δelem
end


#---------------------------------------------------------------------------------
# build_sphere_viscosity(inputs, neqs, TF)
#
# The per-equation artificial viscosity ν of Eq. (8b), read from the deck with
# the SAME switches the flat cases use:
#
#   :lvisc            => false   master switch. Off ⇒ ν ≡ 0, whatever :μ says.
#   :μ                          the coefficient, in m²/s. A scalar applies to
#                               every equation listed in :ivisc_equations; a
#                               vector of length neqs is taken per equation and
#                               :ivisc_equations then only selects from it.
#   :ivisc_equations  => 2:neqs which equations are diffused. The paper puts
#                               δν∇²(φu) in the MOMENTUM equation only — the
#                               continuity equation ∂φ/∂t + ∇·(φu) = 0 carries no
#                               diffusion, and diffusing φ leaks mass across the
#                               shell's height gradients for no physical reason.
#                               Hence the default, which is deliberately NOT
#                               "all equations".
#
# UNDER DynSGS (:visc_model => DSGS_SW()) THE SAME TWO SWITCHES MEAN SOMETHING
# SLIGHTLY DIFFERENT, and it is worth being explicit about it because the key is
# spelled the same. The coefficient is no longer a viscosity: the model supplies
# ν|e itself, per element and per stage, and what is read here is the
# DIMENSIONLESS per-equation MULTIPLIER of Marras eq. (10) — the factor that
# scales (or switches off) each equation's share of the one element coefficient.
# Hence the default of 1.0 rather than 0.0: "DynSGS on, nothing rescaled". A
# deck that leaves :μ at the constant-ν value 1e5 by accident would otherwise
# multiply the model's output by 1e5, so the value is range-checked below.
#---------------------------------------------------------------------------------
function build_sphere_viscosity(inputs, neqs::Int, TF = TFloat)

    μ = zeros(TF, neqs)

    get(inputs, :lvisc, false) == true || return μ

    ldsgs    = sphere_dsgs_requested(inputs)
    μin      = get(inputs, :μ, ldsgs ? 1.0 : 0.0)
    ieq_visc = get(inputs, :ivisc_equations, 2:neqs)

    for ieq in ieq_visc
        1 <= ieq <= neqs ||
            error(string(" # ERROR sphere_rhs.jl: :ivisc_equations lists equation ", ieq,
                         ", but this case has ", neqs, " equations."))
        μ[ieq] = TF(μin isa Number ? μin : μin[ieq])
    end

    any(<(zero(TF)), μ) &&
        error(" # ERROR sphere_rhs.jl: :μ must be non-negative; a negative viscosity is anti-diffusive and blows up immediately.")

    #
    # See the note above: under DynSGS :μ is a dimensionless multiplier, so a
    # value carried over from a constant-ν deck (1e5 m²/s) is an error rather
    # than a strong setting. 100 is far above any sensible rescaling of eq. (10)
    # and far below the smallest physical ν anyone writes here.
    #
    if ldsgs && any(>(TF(100.0)), μ)
        error(string(" # ERROR sphere_rhs.jl: :visc_model => DSGS_SW() makes :μ a DIMENSIONLESS per-equation\n",
                     " #   multiplier of the model's own ν|e, and this deck asks for ", maximum(μ), ".\n",
                     " #   That looks like a constant-viscosity value in m²/s left over from a\n",
                     " #   :visc_model-free deck. Use :μ => 1.0 (or a per-equation vector near 1),\n",
                     " #   or drop :visc_model to get the constant-ν path back."))
    end

    return μ
end


#---------------------------------------------------------------------------------
# sphere_rhs!(RHS, q, qe, mesh, metrics, sp, SVT)
#
# RHS  npoin × neqs, overwritten with ∂q/∂t.
# q    npoin × (neqs+1) conservative state (only the first neqs columns are read)
# qe   npoin × (neqs+1) reference state handed to the user callbacks
#---------------------------------------------------------------------------------
#
# TYPED FUNCTION BARRIER — the same pattern as _viscous_rhs_el_2d! in
# src/kernel/operators/rhs.jl, and it is not optional here.
#
# St_mesh is a Base.@kwdef struct whose fields carry NO type annotation, so
# `mesh.coords`, `mesh.connijk`, `mesh.npoin` … are all `::Any`. Reading them inside
# the element loop makes every access a dynamic lookup that boxes its result:
# on the Galewsky run that cost 63 G allocations / 1.1 TiB and roughly 10x the
# runtime. Hoisting them into concretely-typed arguments once per call lets
# Julia specialise the kernel below, and the inner loops become allocation-free.
#
# Anything read per element or per node belongs in the signature. Anything read
# once per call (SD, the user callbacks' `mesh` argument) can stay dynamic.
#
function sphere_rhs!(RHS, q, qe,
                     mesh::St_mesh,
                     metrics::St_sphere_metrics,
                     sp::St_sphere_params,
                     SVT)

    _sphere_rhs_kernel!(RHS, q, qe,
                        mesh.connijk, mesh.coords,
                        Int(mesh.nelem), Int(mesh.ngl), Int(mesh.npoin),
                        mesh, mesh.SD,
                        metrics.Je, metrics.Minv,
                        metrics.dξdx, metrics.dξdy, metrics.dξdz,
                        metrics.dηdx, metrics.dηdy, metrics.dηdz,
                        metrics.dψ, metrics.ω,
                        sp.F, sp.G, sp.H, sp.S, Int(sp.neqs), SVT,
                        sp.lvisc, sp.μ, sp.gξ, sp.gη,
                        sp.ldsgs, sp.μ_dsgs, sp.Δelem,
                        sp.qnm1, sp.qnm2, sp.davg, sp.ddenom, sp.νel,
                        sp.C1, sp.C2, sp.Δt[])
    return RHS
end

function _sphere_rhs_kernel!(RHS::AbstractMatrix{TF}, q, qe,
                             connijk, crd::AbstractMatrix{TF},
                             nelem::Int, ngl::Int, npoin::Int,
                             mesh, SD,
                             Je, Minv::AbstractVector{TF},
                             dξdx, dξdy, dξdz,
                             dηdx, dηdy, dηdz,
                             dψ, ω,
                             F, G, H, S, neqs::Int, SVT,
                             lvisc::Bool, μ::AbstractVector{TF}, gξ, gη,
                             ldsgs::Bool, μ_dsgs::AbstractMatrix{TF},
                             Δelem::AbstractVector{TF},
                             qnm1::AbstractMatrix{TF}, qnm2::AbstractMatrix{TF},
                             davg::AbstractVector{TF}, ddenom::AbstractVector{TF},
                             νel::AbstractVector{TF},
                             C1::TF, C2::TF, Δt::Float64) where {TF}

    fill!(RHS, zero(TF))

    @inbounds for iel = 1:nelem

        #--- fluxes and sources at the element's LGL nodes
        for j = 1:ngl, i = 1:ngl
            ip = connijk[iel,i,j]

            user_flux!(@view(F[i,j,:]), @view(G[i,j,:]), @view(H[i,j,:]),
                       SD, @view(q[ip,:]), @view(qe[ip,:]),
                       mesh, CL(), SVT; neqs = neqs, ip = ip)

            user_source!(@view(S[i,j,:]), @view(q[ip,:]), @view(qe[ip,:]),
                         npoin, CL(), SVT;
                         neqs = neqs,
                         x = crd[1, ip], y = crd[2, ip], z = crd[3, ip])
        end

        #--- surface divergence + direct stiffness summation
        for ieq = 1:neqs
            for j = 1:ngl
                ωj = ω[j]
                for i = 1:ngl

                    Jij  = Je[iel,i,j]
                    ωJac = ω[i]*ωj*Jij

                    dFdξ = dFdη = zero(TF)
                    dGdξ = dGdη = zero(TF)
                    dHdξ = dHdη = zero(TF)
                    for k = 1:ngl
                        dξ_ki = dψ[k,i]
                        dη_kj = dψ[k,j]
                        dFdξ += dξ_ki*F[k,j,ieq];  dFdη += dη_kj*F[i,k,ieq]
                        dGdξ += dξ_ki*G[k,j,ieq];  dGdη += dη_kj*G[i,k,ieq]
                        dHdξ += dξ_ki*H[k,j,ieq];  dHdη += dη_kj*H[i,k,ieq]
                    end

                    dFdx = dFdξ*dξdx[iel,i,j] + dFdη*dηdx[iel,i,j]
                    dGdy = dGdξ*dξdy[iel,i,j] + dGdη*dηdy[iel,i,j]
                    dHdz = dHdξ*dξdz[iel,i,j] + dHdη*dηdz[iel,i,j]

                    RHS[connijk[iel,i,j], ieq] -= ωJac*((dFdx + dGdy + dHdz) - S[i,j,ieq])
                end
            end
        end

        #--- artificial diffusion, δν∇ₛ²q, in weak form
        #
        # CONSTANT ν only. DynSGS needs the assembled inviscid RHS of the WHOLE
        # domain before it can size ν, so it runs in its own pass below.
        #
        (lvisc && !ldsgs) &&
            _sphere_visc_el!(RHS, q, connijk, iel, ngl,
                             Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                             dψ, ω, μ, gξ, gη, neqs)
    end

    #--- DynSGS: size ν from the residual, THEN diffuse with it
    #
    # WHY A SECOND ELEMENT LOOP. The residual R = ∂q/∂t − M⁻¹·RHS_inviscid is
    # what sets ν, and RHS is only the complete inviscid right-hand side once
    # every element has contributed to it (direct stiffness summation). So the
    # order is forced: assemble everything, measure, then diffuse. Folding the
    # viscous term into the loop above would size ν from a partially assembled
    # RHS whose value at a shared node depends on the element numbering.
    #
    # RHS is still in WEAK form here — the M⁻¹ below is the only place it is
    # applied — which is exactly what compute_dsgs_viscosity! expects: it
    # multiplies by Minv itself when forming the residual, and the viscous term
    # it then adds must go in in weak form like everything else.
    #
    if ldsgs
        #
        # Δt is the denominator of the BDF2 residual, so an unset one does not
        # degrade the model, it inverts it: R → ∞ everywhere and μ pins at the
        # C2·Δ·(|u|+c) cap on every element, i.e. uniform first-order-upwind
        # dissipation dressed up as a residual sensor. sphere_time_loop! fills
        # sp.Δt[] before the first stage; anything else calling sphere_rhs!
        # with DynSGS on has to as well.
        #
        Δt > 0 ||
            error(" # ERROR sphere_rhs.jl: DynSGS is on but sp.Δt[] is unset. Set it to the time step before calling sphere_rhs!.")

        compute_dsgs_viscosity!(μ_dsgs, DSGS_SW(), NSD_2D(),
                                q, qnm2, qnm1, RHS, Minv, μ,
                                davg, ddenom, TF(Δt),
                                connijk, Δelem, C1, C2, nelem, ngl)

        @inbounds for iel = 1:nelem
            for ieq = 1:neqs
                νel[ieq] = μ_dsgs[iel,ieq]
            end
            _sphere_visc_el!(RHS, q, connijk, iel, ngl,
                             Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                             dψ, ω, νel, gξ, gη, neqs)
        end
    end

    #--- M⁻¹ : the mass matrix is diagonal under LGL (inexact) integration
    @inbounds for ieq = 1:neqs, ip = 1:npoin
        RHS[ip,ieq] *= Minv[ip]
    end

    return RHS
end


#---------------------------------------------------------------------------------
# _sphere_visc_el!  — the δν∇²(φu) of Marras/Kopera/Giraldo Eq. (8b), one element.
#
# WHY THIS IS NOT THE FLAT VISCOUS PATH
# -------------------------------------
# ∇² on a shell is the LAPLACE-BELTRAMI operator, and it is built from the same
# 3×2 manifold metrics as the divergence in the kernel above:
#
#   ∇ₛq = a¹ ∂q/∂ξ + a² ∂q/∂η ,   a¹ = (dξdx,dξdy,dξdz) , a² = (dηdx,dηdy,dηdz)
#
# The viscous machinery in src/kernel/operators/rhs.jl builds ∇q from the flat
# 2×2 inverse Jacobian (dξdx, dξdy, dηdx, dηdy — no z), so it cannot even be
# CALLED here: the shell's metric has a third component and the flat one has no
# slot for it. Dropping that component is not an approximation, it is a
# different operator — it differentiates along the projection of the shell onto
# the z = 0 plane, which is singular at the equator. Hence a kernel of its own,
# and hence the fact that setting :lvisc => true did nothing at all before this:
# the shell RHS simply never had a viscous branch.
#
# WEAK FORM
# ---------
# The shell is CLOSED, so integration by parts leaves no boundary term at all —
# the same argument the inviscid part makes:
#
#   ∫ ψ ν∇ₛ²q dΩ = -∫ ∇ₛψ · (ν ∇ₛq) dΩ .
#
# This is the form to use, not the strong ν∇ₛ·(∇ₛq): it needs only ONE
# derivative of the (only C⁰-continuous) solution, it is symmetric negative
# semi-definite by construction — so the term can only ever REMOVE energy, never
# add it, which is exactly what "artificial diffusion that stabilises" has to
# mean — and it is exactly what the flat CG path does.
#
# With the tensor-product nodal basis, ∂ψ_{ij}/∂ξ|_{kl} = dψ[i,k] δ_{jl} and
# ∂ψ_{ij}/∂η|_{kl} = δ_{ik} dψ[j,l], the double sum collapses to two lines:
#
#   ∫∇ₛψ_{ij}·(ν∇ₛq) = Σ_k ω_k ω_j J_{kj} dψ[i,k] gξ_{kj}
#                    + Σ_l ω_i ω_l J_{il} dψ[j,l] gη_{il}
#
#   gξ = ν ∇ₛq·a¹ ,  gη = ν ∇ₛq·a²        (evaluated once per node, per equation)
#
# and it is subtracted from RHS with the sign above, i.e. with the SAME sign as
# the divergence, since the kernel writes `RHS -= ω J (∇·F - S)` and this is
# `RHS -= ∫∇ψ·(ν∇q)`.
#
# THE VECTOR LAPLACIAN, AND WHY COMPONENTWISE IS RIGHT HERE
# ---------------------------------------------------------
# q[2:4] are the CARTESIAN components of φu, and this diffuses each of them
# separately with the scalar Laplace-Beltrami operator — literally what "∇²(φu)"
# means when ∇ = (∂x,∂y,∂z)ᵀ, as it is defined in the paper. It is not the
# covariant (Bochner) vector Laplacian: the two differ by curvature terms of
# order ν|φu|/r², i.e. ~1e5 × 8e6 / 4e13 ≈ 2e-8 m/s² against a Coriolis term of
# ~1e-4 × 8e6 ≈ 8e2 in the same units — ten orders of magnitude down, and this
# is an artificial term whose coefficient is itself a modelling choice. The part
# of it that matters, the small NORMAL component the componentwise operator
# leaves behind, is removed exactly by the Lagrange projection that already runs
# after every RK stage.
#---------------------------------------------------------------------------------
function _sphere_visc_el!(RHS::AbstractMatrix{TF}, q, connijk, iel::Int, ngl::Int,
                          Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                          dψ, ω, μ::AbstractVector{TF}, gξ, gη, neqs::Int) where {TF}

    @inbounds for ieq = 1:neqs

        ν = μ[ieq]
        ν > zero(TF) || continue

        #--- ν ∇ₛq at every node of the element, in the contravariant basis
        for j = 1:ngl, i = 1:ngl

            dqdξ = dqdη = zero(TF)
            for k = 1:ngl
                dqdξ += dψ[k,i]*q[connijk[iel,k,j], ieq]
                dqdη += dψ[k,j]*q[connijk[iel,i,k], ieq]
            end

            a1x = dξdx[iel,i,j]; a1y = dξdy[iel,i,j]; a1z = dξdz[iel,i,j]
            a2x = dηdx[iel,i,j]; a2y = dηdy[iel,i,j]; a2z = dηdz[iel,i,j]

            # ∇ₛq, tangential by construction (a¹ and a² are tangent vectors)
            qx = dqdξ*a1x + dqdη*a2x
            qy = dqdξ*a1y + dqdη*a2y
            qz = dqdξ*a1z + dqdη*a2z

            gξ[i,j,ieq] = ν*(qx*a1x + qy*a1y + qz*a1z)
            gη[i,j,ieq] = ν*(qx*a2x + qy*a2y + qz*a2z)
        end

        #--- -∫ ∇ₛψ·(ν∇ₛq), then direct stiffness summation
        for j = 1:ngl, i = 1:ngl

            t = zero(TF)
            for k = 1:ngl
                t += ω[k]*ω[j]*Je[iel,k,j]*dψ[i,k]*gξ[k,j,ieq] +
                     ω[i]*ω[k]*Je[iel,i,k]*dψ[j,k]*gη[i,k,ieq]
            end

            RHS[connijk[iel,i,j], ieq] -= t
        end
    end

    return RHS
end


#---------------------------------------------------------------------------------
# sphere_dsgs_extra(sp, qvars) -> tuple of name => per-ELEMENT vector
#
# The DynSGS coefficients, packaged for the writers. PER ELEMENT, one value
# each, because that is exactly what they are: the model produces one ν|e per
# element and per RK stage. `write_vtk_sphere_grid` takes these as `extra_cell`
# and writes them as VTK CELL data.
#
# They used to be broadcast onto nodes and written as point data, which is how
# a field that tracks the jet closely (measured: r = +0.79 against both |∇ₛφ|
# and |ζ|, and 5.2× more viscosity inside the jet band than outside) came out
# looking like noise: a node on an element boundary can hold only one of its
# elements' values, and the renderer then interpolates between those arbitrary
# picks, speckling every seam.
#
# HOW MANY FIELDS. Marras eq. (10) is one element coefficient scaled per
# equation by the deck's multiplier, so when those multipliers are all equal —
# which is the shipped SWsphereDSGS deck, :μ => 1.0 on all four equations —
# every slot holds the SAME NUMBER and writing four copies of it is noise that
# invites the reader to look for a difference that cannot exist. One field,
# `mu_dsgs`, is written in that case. When the multipliers differ (the paper's
# :ivisc_equations => [2,3,4] leaves slot 1 at zero, say) the per-equation
# fields carry real information and are written individually as
# `mu_dsgs_<var>`, the same names the flat DSGS paths use (DSGS.md §6).
#
# They are the coefficients of the LAST RK stage of the last completed step, not
# a step average. Since the output cadence is hours and the step is ~a minute,
# that is a snapshot of a quantity that varies from stage to stage; compare a
# field against itself over time rather than reading a single value as exact.
#---------------------------------------------------------------------------------
function sphere_dsgs_extra(sp::St_sphere_params, qvars)

    sp.ldsgs || return ()

    neqs = Int(sp.neqs)
    uniform = all(ieq -> sp.μ[ieq] == sp.μ[1], 1:neqs)

    uniform && return ("mu_dsgs" => @view(sp.μ_dsgs[:, 1]),)

    return ntuple(ieq -> string("mu_dsgs_",
                                ieq <= length(qvars) ? String(qvars[ieq]) : string(ieq)) =>
                         @view(sp.μ_dsgs[:, ieq]), neqs)
end


#---------------------------------------------------------------------------------
# build_sphere_filter_matrix(ξ, inputs, TF)
#
# The exponential modal filter, in nodal form: Filt = V Σ V⁻¹ where V is the
# Legendre Vandermonde on the LGL points and Σ damps the high modes,
#
#   σ_k = 1                                        for k ≤ k₀
#   σ_k = 1 - α ((k - k₀)/(N - k₀))^s              for k > k₀
#
# so σ_N = 1 - α: :filter_alpha is literally "how much the TOP mode is reduced
# by". The paper's Boyd-Vandeven filter "chosen to reduce by 5% the highest
# modes only" corresponds to α = 0.05 with a large s and k₀ near N, which is the
# default here.
#
#   :filter_alpha  damping of the highest mode          (default 0.05)
#   :filter_order  exponent s: larger = sharper cutoff  (default 8)
#   :filter_kcut   first damped mode k₀, as a fraction of N (default 2/3)
#---------------------------------------------------------------------------------
function build_sphere_filter_matrix(ξ, inputs, TF = TFloat)

    ngl = length(ξ)
    N   = ngl - 1

    α  = TF(get(inputs, :filter_alpha, 0.05))
    s  = TF(get(inputs, :filter_order, 8))
    k0 = max(0, min(N-1, round(Int, TF(get(inputs, :filter_kcut, 2/3))*N)))

    # Legendre Vandermonde: V[i,k] = P_{k-1}(ξ_i)
    V = zeros(TF, ngl, ngl)
    for i = 1:ngl
        p0 = one(TF); p1 = ξ[i]
        V[i,1] = p0
        ngl >= 2 && (V[i,2] = p1)
        for k = 2:N
            p0, p1 = p1, ((2k-1)*ξ[i]*p1 - (k-1)*p0)/k
            V[i,k+1] = p1
        end
    end

    σ = ones(TF, ngl)
    for k = k0+1:N
        σ[k+1] = one(TF) - α*((TF(k - k0)/TF(N - k0))^s)
    end

    return V * Diagonal(σ) * inv(V)
end


#---------------------------------------------------------------------------------
# sphere_filter!(q, mesh, metrics, sp)
#
# Apply the modal filter element by element (tensor product: rows then columns),
# then restore continuity with a MASS-WEIGHTED direct stiffness average,
#
#   q[ip] = Σ_e ω_i ω_j J q̃ / M[ip] ,
#
# which is the L² projection back onto the continuous space and therefore
# conserves ∫q exactly — the filter must not create or destroy mass.
#---------------------------------------------------------------------------------
function sphere_filter!(q, mesh::St_mesh, metrics::St_sphere_metrics, sp::St_sphere_params)

    sp.lfilter || return q

    # Typed function barrier; see the note above _sphere_rhs_kernel!.
    _sphere_filter_kernel!(q, mesh.connijk, Int(mesh.nelem), Int(mesh.ngl), Int(mesh.npoin),
                           metrics.Je, metrics.Minv, metrics.ω,
                           sp.Filt, sp.acc, sp.qf, sp.S, Int(sp.neqs))
    return q
end

function _sphere_filter_kernel!(q::AbstractMatrix{TF}, connijk,
                                nelem::Int, ngl::Int, npoin::Int,
                                Je, Minv::AbstractVector{TF}, ω,
                                Fm, acc, qf, S, neqs::Int) where {TF}

    fill!(acc, zero(TF))

    @inbounds for iel = 1:nelem

        # filter along ξ, then along η
        for ieq = 1:neqs
            for j = 1:ngl, i = 1:ngl
                t = zero(TF)
                for k = 1:ngl
                    t += Fm[i,k]*q[connijk[iel,k,j], ieq]
                end
                qf[i,j,ieq] = t
            end
        end
        for ieq = 1:neqs
            for i = 1:ngl
                for j = 1:ngl
                    t = zero(TF)
                    for k = 1:ngl
                        t += Fm[j,k]*qf[i,k,ieq]
                    end
                    S[i,j,ieq] = t             # reuse S as scratch
                end
            end
        end

        # mass-weighted accumulation
        for ieq = 1:neqs, j = 1:ngl, i = 1:ngl
            acc[connijk[iel,i,j], ieq] += ω[i]*ω[j]*Je[iel,i,j]*S[i,j,ieq]
        end
    end

    @inbounds for ieq = 1:neqs, ip = 1:npoin
        q[ip,ieq] = acc[ip,ieq]*Minv[ip]
    end

    return q
end


#---------------------------------------------------------------------------------
# sphere_relative_vorticity!(ζ, q, mesh, metrics)
#
# The RELATIVE VORTICITY of the tangential flow,
#
#   ζ = n̂ · (∇ₛ × u) ,   ∇ₛ × u = a¹ × ∂u/∂ξ + a² × ∂u/∂η ,
#
# with u = φu/φ the Cartesian velocity. This is the field the Galewsky test is
# judged on — the paper (and Galewsky et al.) plot relative vorticity at day 6,
# where the barotropic instability has rolled the jet into a train of vortices.
# The height field barely moves by comparison, so h alone makes a developing
# instability almost invisible.
#
# Computed element-wise from the metric terms and then made continuous with the
# same mass-weighted direct stiffness average used by the filter, so ζ is a
# proper nodal field on the shell rather than a discontinuous per-element one.
#---------------------------------------------------------------------------------
function sphere_relative_vorticity!(ζ, q, mesh::St_mesh,
                                    metrics::St_sphere_metrics, sp::St_sphere_params)
    # Typed function barrier; see the note above _sphere_rhs_kernel!. Without it
    # this cost 86 MB per call (mesh.connijk / mesh.nelem / mesh.npoin are ::Any),
    # and it is called at every diagnostics print and every output.
    _sphere_vorticity_kernel!(ζ, q, mesh.connijk,
                              Int(mesh.nelem), Int(mesh.ngl), Int(mesh.npoin),
                              metrics.dξdx, metrics.dξdy, metrics.dξdz,
                              metrics.dηdx, metrics.dηdy, metrics.dηdz,
                              metrics.nx, metrics.ny, metrics.nz,
                              metrics.Je, metrics.Minv, metrics.dψ, metrics.ω,
                              sp.acc)
    return ζ
end

function _sphere_vorticity_kernel!(ζ::AbstractVector{TF}, q::AbstractMatrix{TF},
                                   connijk, nelem::Int, ngl::Int, npoin::Int,
                                   dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                                   nx, ny, nz, Je, Minv::AbstractVector{TF},
                                   dψ, ω, acc::AbstractMatrix{TF}) where {TF}

    fill!(acc, zero(TF))

    @inbounds for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl

            # ∂u/∂ξ and ∂u/∂η of the CARTESIAN velocity
            dξx = dξy = dξz = zero(TF)
            dηx = dηy = dηz = zero(TF)
            for k = 1:ngl
                ipk = connijk[iel,k,j]
                φk  = q[ipk,1]
                dξx += dψ[k,i]*q[ipk,2]/φk; dξy += dψ[k,i]*q[ipk,3]/φk; dξz += dψ[k,i]*q[ipk,4]/φk

                ipl = connijk[iel,i,k]
                φl  = q[ipl,1]
                dηx += dψ[k,j]*q[ipl,2]/φl; dηy += dψ[k,j]*q[ipl,3]/φl; dηz += dψ[k,j]*q[ipl,4]/φl
            end

            a1x, a1y, a1z = dξdx[iel,i,j], dξdy[iel,i,j], dξdz[iel,i,j]
            a2x, a2y, a2z = dηdx[iel,i,j], dηdy[iel,i,j], dηdz[iel,i,j]

            cx = (a1y*dξz - a1z*dξy) + (a2y*dηz - a2z*dηy)
            cy = (a1z*dξx - a1x*dξz) + (a2z*dηx - a2x*dηz)
            cz = (a1x*dξy - a1y*dξx) + (a2x*dηy - a2y*dηx)

            zt = cx*nx[iel,i,j] + cy*ny[iel,i,j] + cz*nz[iel,i,j]

            acc[connijk[iel,i,j], 1] += ω[i]*ω[j]*Je[iel,i,j]*zt
        end
    end

    @inbounds for ip = 1:npoin
        ζ[ip] = acc[ip,1]*Minv[ip]
    end

    return ζ
end
