#---------------------------------------------------------------------------------
# astroJetWuShu2018 — the classical magnetized astrophysical jet.
#
#   K. Wu, C.-W. Shu, "A provably positive discontinuous Galerkin method for
#   multidimensional ideal magnetohydrodynamics", SIAM J. Sci. Comput. 40(5)
#   (2018) B1302-B1329, Example 5.6 ("Astrophysical jets"),
#
# which magnetizes the Mach 800 gas dynamical jet of
#
#   D. S. Balsara, "Self-adjusting, positivity preserving high order schemes for
#   hydrodynamics and magnetohydrodynamics", JCP 231 (2012) 7504-7517
#
# (itself after X. Zhang & C.-W. Shu, JCP 229 (2010) 8918). The same test, with
# the same constants, is the closing example of most later structure-preserving
# MHD papers — e.g. Peng, Sun & Wu, "Structure-preserving oscillation-
# eliminating DG schemes for ideal MHD" (arXiv:2404.16794), §4.2.6.
#
# Run with:
#   julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
#
# All the physical constants live in user_flux.jl; the setup is spelled out in
# README.md and in the header of initialize.jl. READ README.md BEFORE RUNNING:
# in the beam the gas pressure is 5.6 parts per million of the total energy, and
# every published solution of this test uses an explicit positivity-preserving
# limiter that this code does not have.
#
# Overrides, all optional, all read once at include time:
#
#   physics (user_flux.jl / user_source.jl)
#     JEXPRESSO_AJ_BA2      B_a SQUARED, as the paper quotes it (default 200;
#                           the paper's three cases are 200, 2000, 20000, i.e.
#                           plasma β_a = 1e-2, 1e-3, 1e-4)
#     JEXPRESSO_AJ_UJET     injection speed = beam Mach number (default 800)
#     JEXPRESSO_AJ_SMOOTH   nozzle-lip smoothing length (default 0 = top hat)
#     JEXPRESSO_AJ_PFLOOR   flux-only pressure floor (default 1e-6, 0 = off)
#     JEXPRESSO_AJ_GLMCR    Dedner's c_r (default 0.18)
#
#   discretization (this file)
#     JEXPRESSO_AJ_MESH     "40x60" (default) or "100x150"
#     JEXPRESSO_AJ_NOP      polynomial order (default 4)
#     JEXPRESSO_AJ_DT       time step (default: per mesh, see _aj_dt below)
#     JEXPRESSO_AJ_TEND     final time (default 2e-3, the paper's)
#     JEXPRESSO_AJ_DTOUT    output interval (default tend/20 = 1e-4)
#
#   DynSGS (this file)
#     JEXPRESSO_AJ_SENSOR   "residual" (default) or "legacy"
#     JEXPRESSO_AJ_CR       C_R    (default 1.0)
#     JEXPRESSO_AJ_CMAX     C_max  (default 0.5)
#     JEXPRESSO_AJ_CMIN     C_min background floor (default 0.0 — see the note)
#     JEXPRESSO_AJ_CUTOFF   smoothness cutoff on the normalized residual (0)
#     JEXPRESSO_AJ_NORMS    "domain" (default), "rank" or "element"
#     JEXPRESSO_AJ_REL      :dsgs_rel, the physical-scale floor factor (1.0)
#     JEXPRESSO_AJ_NODAL    1 = one ν per NODE instead of per element (only with
#                           JEXPRESSO_AJ_SENSOR=legacy — see the note on the key)
#     JEXPRESSO_AJ_HOLD     :dsgs_hold_steps (default 0 — see the note on the key)
#---------------------------------------------------------------------------------
_aj_mesh()   = String(strip(get(ENV, "JEXPRESSO_AJ_MESH", "40x60")))
_aj_nop()    = something(tryparse(Int,     get(ENV, "JEXPRESSO_AJ_NOP",    "")), 4)
_aj_tend()   = something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_TEND",   "")), 2.0e-3)
# Output interval. The default is tend/20, so that a run on one of the low-Mach
# rungs of README.md §6 (which needs a proportionally longer tend) writes 21
# snapshots and not 800. At the default tend = 2e-3 it is 1e-4, which includes
# the paper's three output times 1e-3, 1.5e-3 and 2e-3.
_aj_dtout()  = something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_DTOUT",  "")), _aj_tend()/20.0)
_aj_sensor() = String(lowercase(strip(get(ENV, "JEXPRESSO_AJ_SENSOR", "residual"))))
_aj_CR()     = something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_CR",     "")), 1.0)
_aj_Cmax()   = something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_CMAX",   "")), 0.5)
_aj_Cmin()   = something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_CMIN",   "")), 0.0)
_aj_cutoff() = something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_CUTOFF", "")), 0.0)
_aj_norms()  = String(lowercase(strip(get(ENV, "JEXPRESSO_AJ_NORMS",  "domain"))))
_aj_rel()    = something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_REL",    "")), 1.0)
_aj_nodal()  = lowercase(strip(get(ENV, "JEXPRESSO_AJ_NODAL", "false"))) in ("1", "true", "yes", "on")
_aj_hold()   = something(tryparse(Int,     get(ENV, "JEXPRESSO_AJ_HOLD",   "")), 0)

# The two meshes that ship with the case, and the time step each one wants.
#
# The element size h, the smallest LGL spacing at :nop => 4 (a gap of
# 1 - sqrt(3/7) = 0.34535 in the reference element, so 0.17267·h), and the
# Courant number of the default Δt against the largest wave speed of the
# problem (c_h; see initialize.jl):
#
#   mesh      h       elements   LGL points   Δx_min    Δt      CFL      steps
#   40x60     0.025    2 400      160x240     4.317e-3  5e-7   0.09-0.11  4 000
#   100x150   0.01    15 000      400x600     1.727e-3  2e-7   0.09-0.11 10 000
#
# The CFL range spans the paper's three magnetizations: c_h goes from 812 at
# B_a = √200 to 920 at B_a = √20000, only 13 %, because it is set by the beam
# speed 800 and not by the field. ONE Δt per mesh therefore covers all three.
#
# 100x150 is the paper's own resolution: it computes the right half
# [0, 0.5] x [0, 1.5] on 200 x 600 cells, i.e. Δ = 2.5e-3, and 100 x 150
# elements at :nop => 4 give 400 x 600 unique points over the FULL width, the
# same spacing. (This case uses the full domain and not the half domain with a
# reflecting axis at x = 0 — see README.md.)
#
# The VISCOUS limit is not the binding one: DynSGS cannot exceed its own cap
# μ_max = C_max Δ (|v| + c_f) = 0.5·(h/5)·c_h, which is 2.03 on 40x60, and
# 0.5·Δx_min²/μ_max = 4.6e-6 there — nine times the advective Δt. On 100x150 it
# is 1.8e-6, nine times again. Raising :dsgs_Cmax or :μ changes that ratio.
if !@isdefined(AJ_MESHES)
    const AJ_MESHES = Dict("40x60"   => ("AJ_40x60.msh",   5.0e-7),
                           "100x150" => ("AJ_100x150.msh", 2.0e-7))
end

function _aj_mesh_entry()
    k = _aj_mesh()
    haskey(AJ_MESHES, k) || error(string(" problems/MHD/astroJetWuShu2018: JEXPRESSO_AJ_MESH=\"", k,
                                        "\" is not one of ", join(sort(collect(keys(AJ_MESHES))), ", "),
                                        ". Set :gmsh_filename in user_inputs.jl by hand for a mesh of your own."))
    return AJ_MESHES[k]
end
_aj_gmsh() = string("./problems/MHD/astroJetWuShu2018/", _aj_mesh_entry()[1])
_aj_dt()   = something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_DT", "")), _aj_mesh_entry()[2])

function user_inputs()

    tend = _aj_tend()

    # The one env combination that is documented to be wrong (see :ldsgs_nodal
    # below): the nodal kernel's assembled residual is blind to
    # under-resolution, so pairing it with the strong-residual sensor removes
    # the sensor rather than refining it.
    if _aj_nodal() && _aj_sensor() == "residual"
        @warn string(" problems/MHD/astroJetWuShu2018: JEXPRESSO_AJ_NODAL=1 with ",
                     "JEXPRESSO_AJ_SENSOR=residual is a BLIND sensor — the nodal kernel reads the ",
                     "assembled residual, which vanishes on an under-resolved solution (DSGS.md §1.2; ",
                     "measured as a 17x regression on CompEuler/shock_circle_M7). Set ",
                     "JEXPRESSO_AJ_SENSOR=legacy alongside it, or leave the nodal form off.")
    end

    inputs = Dict(
        #---------------------------------------------------------------------------
        # 2D magnetized astrophysical jet, ideal GLM-MHD.
        #
        # Domain [-0.5, 0.5] x [0, 1.5], γ = 1.4, static magnetized ambient
        # medium (ρ, p) = (0.1γ, 1) with B = (0, B_a, 0), and a Mach 800 beam of
        # density γ injected through the nozzle {y = 0, |x| <= 0.05}. Final time
        # 2e-3, which is when the paper's figures are drawn and just after the
        # jet head (800 · 2e-3 = 1.6 > 1.5) would have crossed the box had it not
        # been decelerated by the ambient medium it is ploughing into.
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => _aj_dt(),
        :tinit                => 0.0,
        :tend                 => tend,
        :diagnostics_at_times => (0.0:_aj_dtout():tend),   # includes the references' t = 1e-3, 1.5e-3, 2e-3
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,   # GLM ψ-damping only (Dedner mixed cleaning; user_source.jl)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => _aj_nop(),
        #---------------------------------------------------------------------------
        # STABILIZATION.
        #
        # The papers that run this test solve ideal MHD with a DG or finite-volume
        # scheme whose Riemann solver supplies the numerical dissipation AND with
        # an explicit positivity-preserving limiter on top; Wu & Shu state that
        # without the PP limiter "the simulation will break down after several
        # time steps due to nonphysical numerical solutions", and the OEDG paper
        # says the same of its oscillation-eliminating and PP procedures. A
        # collocated continuous-Galerkin SEM has neither a Riemann solver nor a
        # limiter, so here the ONLY thing standing between the Mach 800 beam and
        # a negative pressure is the Marras-Nazarov DynSGS dissipation. Note
        # that the node-wise realizability repair of src/kernel/positivity/
        # (:lpositivity, used by the Mach-7 CompEuler decks) is NOT available:
        # it is scoped to neqs == nsd + 2 and errors on this nine-field state,
        # because the magnetic energy in p = (γ-1)(ρE − KE − ½|B|² − ½ψ²) is
        # outside what it knows. Extending it to the GLM-MHD state is the
        # obvious next step for this case; until then the flux-only pressure
        # floor of user_flux.jl is the whole of the safety net.
        #
        # DynSGS (DSGS.md; Marras, Nazarov & Giraldo, JCP 301 (2015) 77; Dao &
        # Nazarov, J. Sci. Comput. 92:77 (2022), §4.4) sets the eddy viscosity
        # from the LOCAL RESIDUAL of the governing equations rather than from a
        # tuned constant: one kinematic
        #
        #     ν|_e = max(0, min( C_max Δ (|v|+c_f)|_e ,
        #                        C_R Δ² max_i ‖R_i‖_∞,e / ‖q_i-⟨q_i⟩‖_∞,Ω ))
        #
        # per element, Δ = h/(N+1). It is parameter-free — C_R = 1, C_max = 0.5
        # are the paper's own values — and at a discontinuity the normalized
        # residual is O(10²) while smooth flow gives O(10⁻⁶), so the coefficient
        # saturates at the first-order cap exactly on the Mach shock at the jet
        # head and on the beam/cocoon interface, and is ~0 in the quiescent
        # ambient medium. That is the property this case is here to exercise.
        #
        # CONSERVED FORM (:dsgs_conserved => true). The dissipation is
        # ∇·(ν∇q) on all nine conserved variables — see user_primitives.jl for
        # why the physical form (ν on ρ, ρν on u, κ on T, η on B) is the wrong
        # choice for a 10:1 contact whose pressure is a 5.6e-6 residue of the
        # total energy. Same form as the MHD shock tube problems/MHD/brioWu1d.
        #
        # :μ are per-equation multipliers on ν for
        # (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ), all at full strength. In the
        # conserved form :μ[1] is NOT optional: it is the mass diffusion ν∇ρ
        # that keeps the density jump across the shock from ringing (on
        # problems/CompEuler/ffs_step, setting it to zero is the single most
        # destabilizing change measured).
        #
        # :dsgs_gamma MUST equal γ_mhd = 1.4 of user_flux.jl.
        #---------------------------------------------------------------------------
        :lvisc            => true,
        :μ                => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model       => DSGS_MHD(),
        # SENSOR. "residual" (the default, DSGS.md §1.2) is the element-wise
        # strong residual with the stage-consistent time stencil; "legacy" is the
        # pre-September-2026 sensor (the assembled RHS against a fixed BDF2,
        # in effect |∂ₜq|), which the atmospheric decks and
        # orszagTangBormanis2024 were validated with.
        #
        # "residual" is the default HERE for a reason that is specific to this
        # case: it is the only sensor path that carries the Dirichlet-boundary
        # treatment of rhs.jl (_dsgs_bdy_zero!). At a constrained node the
        # element's own RHS is the CONSTRAINT FORCE, not an under-resolution, and
        # feeding it to the sensor drives ν to its cap along the whole boundary —
        # measured on the rising bubble, where it blew the run up. This case has
        # four non-periodic boundaries and a strongly constrained inflow patch, so
        # that protection matters. What it does NOT do is blind the model at the
        # nozzle: only the boundary NODES are zeroed, and an element on the
        # bottom edge has 5 of its 25 nodes there, so the 20 interior ones still
        # set the element's coefficient.
        :dsgs_sensor      => _aj_sensor(),
        :dsgs_CR          => _aj_CR(),
        # C_max = 0.5 is the value both validated DSGS_MHD cases use
        # (brioWu1d, orszagTangBormanis2024) and it is right for THIS kernel:
        # compute_dsgs_viscosity!(::DSGS_MHD) caps at C_max·Δ·(|v|+c_f) with
        # Δ = Δelem/(N+1) ALREADY divided by the polynomial count, so at
        # nop = 4 it is 0.2h against a smallest LGL gap of 0.17267h — a factor
        # 1.16, not the 5.8 that Δ = Δelem would give. (The Mach-7 CompEuler
        # decks run C_max = 0.1 for that reason; the kernels there divide by
        # ngl too, so read their comment with this in mind before porting the
        # 0.1 here.)
        :dsgs_Cmax        => _aj_Cmax(),
        # BACKGROUND FLOOR C_min·Δ·(|v|+c_f), a fraction of the wave-speed cap.
        # 0 here (pure Marras), NOT brioWu1d's 0.06, and the reason is the beam:
        # the floor is proportional to the LOCAL wave speed, which in the beam is
        # 812 rather than the O(1) of a shock tube. C_min = 0.06 would put
        # ν = 0.06·(0.025/5)·812 = 0.24 inside the beam, and √(2νt) over the run
        # is then 0.031 — 60 % of the beam's own half-width 0.05, i.e. the floor
        # alone would smear the beam away. In the quiescent ambient the same
        # C_min is harmless (ν = 0.011, √(2νt) = 0.0068), so if a node-to-node
        # (checkerboard) mode does appear — the one thing the residual sensor is
        # blind to, because the discrete operator returns almost nothing on it —
        # JEXPRESSO_AJ_CMIN=0.005..0.01 is the first lever to try, and the beam
        # profile is what to check afterwards.
        :dsgs_Cmin        => _aj_Cmin(),
        # Smoothness cutoff ratio -> max(0, ratio - cutoff). 0 = off. At a shock
        # the normalized ratio is O(10²), so a cutoff of 1e-3 would be invisible
        # here; it exists for smooth accuracy studies (smoothVortex) and is left
        # off because this case has nothing smooth to protect.
        :dsgs_cutoff      => _aj_cutoff(),
        # STARTUP HOLD OFF (the kernel default is 2), for the same reason
        # problems/CompEuler/ffs_step turns it off: this initial condition is
        # not smooth data that the sensor would misread, it is a Mach 800 beam
        # started IMPULSIVELY against a medium at rest, and the first steps are
        # the most violent of the whole run. Holding ν at zero through them
        # integrates exactly the steps that need dissipation with none at all,
        # and the oscillation planted at the nozzle lip is what the rest of the
        # run then has to carry. (On ffs_step the hold is the difference between
        # reaching t = 8e-3 and dying at 1.46e-3, in its convex corner.)
        #
        # The hold exists because on SMOOTH data the BDF2 seeded from the
        # initial condition makes the residual the whole flux divergence, and
        # the sensor reads a fully resolved field as unresolved everywhere
        # (measured on the smooth vortex: ν at its cap on the first call). That
        # failure mode cannot happen here: this initial condition is UNIFORM, so
        # ∇·F ≡ 0 and the residual is exactly zero everywhere except in the
        # elements the nozzle datum reaches. ν on step one is therefore at the
        # cap at the nozzle and zero in the rest of the domain — which is
        # precisely what is wanted.
        :dsgs_hold_steps  => _aj_hold(),
        :dsgs_gamma       => 1.4,           # = γ_mhd
        :dsgs_Prt         => 1.0,           # unused in the conserved, non-split form
        :dsgs_conserved   => true,          # Laplacian on the conserved variables (user_primitives.jl)
        :dsgs_ref_weight  => false,         # the reference state is uniform, so this would be a no-op
        :dsgs_nazarov_energy => false,      # keep ONE ν on the energy slot (user_primitives.jl)
        :dsgs_nodal_rho   => false,         # (inactive under :dsgs_conserved anyway)
        # :ldsgs_nodal => false (one ν per element) is the default and it MUST
        # stay false with the "residual" sensor above: the nodal kernel reads
        # the ASSEMBLED residual, which the lumped-LGL assembly cancels on an
        # under-resolved solution just as it does on a resolved one (DSGS.md
        # §1.2), so nodal + "residual" is a BLIND sensor. Measured on
        # problems/CompEuler/shock_circle_M7 as a 17x regression (231 steps
        # against 3969). Paired with "legacy" the nodal form is legitimate, and
        # that pairing is what JEXPRESSO_AJ_NODAL=1 is for — the deck warns
        # below if only one of the two is set.
        :ldsgs_nodal      => _aj_nodal(),
        :dsgs_Cl          => 0.4,           # Dao & Nazarov's C_l, used only by the nodal form
        # NORMALIZATION SCOPE. "domain" is the method's own norm over Ω and the
        # only choice that makes the answer independent of the MPI partition
        # (two small Allreduce per RHS call). "rank" normalizes by each
        # subdomain's spread, which on a case like this — where 90 % of the
        # ranks hold nothing but quiescent ambient gas and a few hold the whole
        # jet — is exactly the pathology described in SGS.jl: the quiet ranks
        # normalize by their own floor and apply a different viscosity to the
        # same solution than their neighbours.
        :dsgs_norms       => _aj_norms(),
        # How far below its own physical scale a variable's spread may fall
        # before that scale, and not the spread, normalizes its residual.
        #
        # CAVEAT worth knowing on this case. The B floor in the kernel is
        # √ρ̄ c̄ = √0.14 · √10 = 1.18, which is the field strength of a β ≈ 1
        # plasma — but this problem starts at β_a = 1e-2, so B_a = 14.1 is 12×
        # that (120× at B_a = √20000). At t = 0 the field is uniform, its spread
        # is zero and the floor IS its normalization, so the B equations would be
        # read against a scale well below their own magnitude. What makes it
        # harmless is that the NUMERATOR is zero too at that moment: the injected
        # field is the SAME (0, B_a, 0) as the ambient field, so the nozzle
        # introduces no jump in B at all and the induction residual only grows
        # once the beam starts shearing the field — by which time By has a spread
        # of its own order, which wins the max. Do NOT raise :dsgs_rel to "fix"
        # it: the key scales every floor, so it would desensitize ρ, ρv and E at
        # the same time.
        :dsgs_rel         => _aj_rel(),
        :lrichardson      => false,      # no gravity/stratification in this problem
        # Slot 4 carries the TOTAL ENERGY ρE. ("theta" would be wrong for any
        # shock: ρθ is an entropy variable, conserved across a contact but not
        # across a shock, so the Euler-θ system carries the wrong shock speed no
        # matter how it is stabilized.)
        :energy_equation  => "energy",
        #---------------------------------------------------------------------------
        # No entropy-stable / kinetic-energy-preserving machinery:
        #---------------------------------------------------------------------------
        :lkep              => false,
        :entropy_variables => false,
        #---------------------------------------------------------------------------
        # Mesh. Structured quads on [-0.5,0.5] x [0,1.5]; both meshes put the
        # nozzle lip |x| = 0.05 exactly on an element boundary (0.05/h = 2 and 5),
        # which is what makes the top-hat inflow patch of user_bc.jl land on
        # nodes rather than inside an element. Boundary physical-curve tags:
        # "bottom", "right", "top", "left". See AJ.geo to regenerate or refine.
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => _aj_gmsh(),
        #---------------------------------------------------------------------------
        # Filter parameters.
        #
        # OFF on purpose. The Boyd-Vandeven "erf" filter is Jexpresso's other
        # stabilization mechanism, but it filters the CONSERVATIVE variables
        # INDEPENDENTLY, and on this problem that is not a small error: filtering
        # ρ, ρv and ρE separately perturbs each of the three large terms of
        # p = (γ-1)(ρE - ½ρ|v|² - ½|B|²) by a different amount, and their
        # near-cancellation in the beam means a 1e-5 relative change is enough to
        # take the pressure negative. DynSGS is applied in divergence form to the
        # whole conserved state at once, which does not have that failure mode.
        #---------------------------------------------------------------------------
        :lfilter             => false,
        #---------------------------------------------------------------------------
        # Plotting. ParaView: the 14 output variables of initialize.jl (including
        # log10rho and log10p, which are what the papers plot), the DynSGS
        # coefficient fields mu_dsgs_*, and the two numerical-schlieren fields.
        #
        # :lschlieren gives |∇ρ| and exp(-k|∇ρ|/max|∇ρ|) — the second one, in a
        # REVERSED greyscale, is the "Schlieren image of the density logarithm"
        # the reference figures actually are. Computed at output times only, so it
        # costs nothing in the time loop.
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => false,   # a new output dir per run: this case is meant to be swept
        :lwrite_initial      => true,
        :output_dir          => "./output/",
        :loutput_pert        => false,   # the total state, not the departure from the ambient medium
        :lschlieren          => true,
        :schlieren_k         => 20.0,    # contrast; Hadjadj uses 10-100
        #---------------------------------------------------------------------------
        # AMR off. The mesh is uniform and DynSGS handles what is left
        # under-resolved. NOTE :linitial_refine => true with :init_refine_lvl => n
        # would uniformly h-refine the mesh read above by 2ⁿ per direction through
        # p4est (kernel/mesh/mesh.jl), which is a cheap way to go finer than
        # 100x150 without generating a mesh — but it also changes Δx_min, so Δt
        # has to come down by the same factor.
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 0,
        :ladapt              => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
