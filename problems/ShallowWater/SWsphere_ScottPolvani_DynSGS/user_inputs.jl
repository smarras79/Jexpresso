#---------------------------------------------------------------------------------
# SWsphere_ScottPolvani_DynSGS — forced-dissipative shallow-water turbulence on
# the sphere with giant-planet parameters, stabilised by the DYNAMIC SGS
# (residual-based) viscosity of Marras, Nazarov & Giraldo (2015) instead of the
# constant ν∇ₛ² + modal filter of its sibling SWsphere_ScottPolvani. Everything
# else — grid, equations, planet, forcing, dissipation, outputs — is the same
# deck; see the STABILISATION block and the README. The setup is that of
#
#   Scott, R. K. & Polvani, L. M. (2007), "Forced-dissipative shallow-water
#   turbulence on the sphere and the atmospheric circulation of the giant
#   planets", J. Atmos. Sci. 64, 3158-3176.
#
# SAME GRID, SAME EQUATIONS, SAME DISCRETISATION as problems/ShallowWater/
# SWsphere (the Galewsky jet): the cubed sphere shipped next to this file, the
# conservative Cartesian form q = [φ, φu, φv, φw] of Marras, Kopera & Giraldo
# (2015), continuous-Galerkin spectral elements, SSP-RK3 with the Lagrange
# projection after every stage, the modal filter and the ν∇ₛ² artificial
# diffusion. What changes is the SETUP — a planet instead of the Earth, a fluid
# at rest instead of a balanced jet, and a forcing that keeps injecting energy:
#
#   * the fluid starts AT REST with a uniform depth H, fixed by the planet's
#     Rossby deformation radius L_D = √(gH)/(2Ω)          (initialize.jl);
#   * a random, isotropic, small-scale VORTICITY forcing injects kinetic energy
#     at a prescribed rate ε₀ in a band of spherical-harmonic wavenumbers
#     around n_f                              (src/kernel/operators/sphere_forcing.jl);
#   * a linear large-scale dissipation — Rayleigh friction on the momentum
#     and/or radiative relaxation of the height — lets the energy equilibrate;
#   * Coriolis uses the planet's Ω                              (user_source.jl).
#
# The paper's parameters are NONDIMENSIONAL: lengths in planetary radii a,
# times in rotation periods T = 2π/Ω. Jexpresso integrates in SI units, so the
# block below converts once and the Dict carries only dimensional numbers; the
# nondimensional ones are kept alongside (:sp_*) for the record and the
# printout. The two numbers that characterise a planet in this model are
#
#   Ro  = U/(2aΩ)          the Rossby number of the equilibrated flow
#   L_D = √(gH)/(2Ω a)     the polar deformation radius, in planetary radii
#
# and the paper's section 7 runs — its Figs. 13 and 14 — (values from Cho &
# Polvani 1996b):
#
#   Jupiter         Ro = 0.002   L_D = 0.025     Fig. 14 (a)/(d)
#   Saturn          Ro = 0.013   L_D = 0.025     Fig. 14 (b)/(e)
#   Uranus/Neptune  Ro = 0.06    L_D = 0.1       Fig. 14 (c)/(f)
#
# All three are selectable below (`planet`, or SP_PLANET in the environment).
#
# Ro is an OUTCOME of forcing and dissipation, not an input. With Rayleigh
# friction ν_l the paper's closure is E ≈ ε/(2ν_l) with ε ≈ 0.4 ε₀ the part of
# the input that cascades upscale (section 4a), and E = ½U², so a target Ro is
# turned into an input rate ε₀ = 2ν_l · ½(4πRo)² / 0.4 below. That is a
# scaling estimate, and with ν∇ₛ² instead of the paper's ∇⁸ hyperdiffusion the
# upscale fraction will differ — read the Ro the diagnostics actually report.
#
# WHAT IS NOT THE PAPER. Two of its numerical ingredients are not reproduced:
# the ∇⁸ hyperdiffusion at the small scales (the shell's modal filter + ν∇ₛ²
# stand in for it) and the (-Δ)⁻¹ hypodiffusion it uses for most of section 5
# (a global inverse Laplacian; the Rayleigh friction it shows to give similar
# equilibria, section 4, is used instead). See the README.
#
# RESOLUTION. The shipped grid (10 elements per panel edge, nop = 5) has 200
# nodes around the equator, i.e. it resolves total wavenumbers up to n ≈ 50-60
# (T170 in the paper, with n_f = N/4 ≈ 42). The forcing here is centred at
# n_f = 24 — 8 nodes per wavelength at the equator — which leaves the jets
# (n₀ ≈ 5-17 in the paper's Tables 1-2) room below it. Lower n_f on a coarser
# grid, raise it on a finer one; keep N/n_f ≳ 2-4.
#
# The five sibling user_*.jl in this directory are the REAL implementations (a
# case directory owns its six files; see test/test_case_includes.jl for why a
# one-line `include` of a sibling case is not an option).
#---------------------------------------------------------------------------------

function user_inputs()

    #---------------------------------------------------------------------------
    # THE PROBLEM: which planet. This is the one user choice of the deck;
    # everything dimensional below follows from it.
    #
    # The three simulations of the paper's section 7 — Fig. 13 (zonal-mean
    # wind) and Fig. 14 (potential vorticity, panels a-c; vorticity, panels
    # d-f), parameters from Cho & Polvani (1996b):
    #
    #   :jupiter    Fig. 14 (a)/(d)    Ro = 0.002   L_D = 0.025
    #   :saturn     Fig. 14 (b)/(e)    Ro = 0.013   L_D = 0.025
    #   :neptune    Fig. 14 (c)/(f)    Ro = 0.06    L_D = 0.1
    #
    # The paper labels the third "Uranus/Neptune" and treats the two as one
    # case; :uranus is accepted as an alias of that panel, with Uranus's own
    # a, Ω, g (the dynamics only see Ro and L_D, which are the same).
    #
    # SELECT HERE — a symbol or a string, either works:
    #
    #   planet = :jupiter        planet = :saturn        planet = :neptune
    #
    # or leave the line alone and set SP_PLANET in the environment, which is
    # how one scripts the three panels in a row (it overrides the line):
    #
    #   SP_PLANET=saturn julia --project=. -e 'using Jexpresso; Jexpresso.run_case("ShallowWater", "SWsphere_ScottPolvani")'
    #
    # SP_NROT overrides the run length in rotations, and SP_MESH the grid file,
    # the same way.
    #
    #   a   equatorial radius [m]      Ω   rotation rate [1/s]     g  gravity [m/s²]
    #   LD  deformation radius / a     Ro  target Rossby number    fig  the panel
    #
    # a, Ω, g are the bodies' measured values (equatorial radius, sidereal
    # rotation, equatorial gravity); only L_D and Ro enter the dynamics, since
    # the equations are scale-free in a and Ω and carry φ = gh rather than h.
    #---------------------------------------------------------------------------
    planet = :jupiter

    # the environment override, and :Saturn / "saturn" / "SATURN" all read as :saturn
    haskey(ENV, "SP_PLANET") && (planet = ENV["SP_PLANET"])
    planet = Symbol(lowercase(strip(string(planet))))

    planets = Dict(
        :jupiter => (a = 7.1492e7, Ω = 1.7585e-4, g = 24.79, LD = 0.025, Ro = 0.002, fig = "Fig. 14 (a)/(d), Jupiter"),
        :saturn  => (a = 6.0268e7, Ω = 1.6378e-4, g = 10.44, LD = 0.025, Ro = 0.013, fig = "Fig. 14 (b)/(e), Saturn"),
        :neptune => (a = 2.4764e7, Ω = 1.0834e-4, g = 11.15, LD = 0.1,   Ro = 0.06,  fig = "Fig. 14 (c)/(f), Uranus/Neptune"),
        :uranus  => (a = 2.5559e7, Ω = 1.0124e-4, g =  8.87, LD = 0.1,   Ro = 0.06,  fig = "Fig. 14 (c)/(f), Uranus/Neptune (Uranus's constants)"),
    )
    haskey(planets, planet) ||
        error(" # ERROR user_inputs.jl: unknown planet :", planet,
              "; choose :jupiter, :saturn or :neptune (the three panels of Fig. 14; :uranus is an alias of the third), or set SP_PLANET")
    P = planets[planet]

    a, Ω, g = P.a, P.Ω, P.g
    T       = 2π/Ω                     # one rotation [s]: the paper's unit of time

    #---------------------------------------------------------------------------
    # NONDIMENSIONAL parameters, in the paper's units (a, T).
    #---------------------------------------------------------------------------
    LD_nd    = P.LD                    # deformation radius / a
    Ro_target= P.Ro                    # target Rossby number (sets ε₀ through the closure)
    nul_nd   = 1.0e-4                  # Rayleigh friction, per rotation (paper Fig. 2a: 1e-4)
    nuh_nd   = 0.0                     # radiative relaxation, per rotation (paper Fig. 4a: 1; 0 = off)
    tau_nd   = 10.0                    # forcing decorrelation time, rotations (paper: c_r = 10; 0 = white)
    f_up     = 0.4                     # fraction of ε₀ that cascades upscale (paper, section 4a)
    # (no nu_nd here: the small-scale dissipation is DynSGS, see STABILISATION)

    E_nd     = 0.5*(4π*Ro_target)^2    # ½U² with U = 4π Ro in a/T units
    eps0_nd  = 2.0*nul_nd*E_nd/f_up    # the Rayleigh-friction closure E = 0.4 ε₀/(2ν_l)

    #---------------------------------------------------------------------------
    # DIMENSIONAL conversions. These are what the code reads.
    #---------------------------------------------------------------------------
    H     = (2Ω*LD_nd*a)^2/g           # mean depth from L_D = √(gH)/(2Ω)   [m]
    eps0  = eps0_nd*a^2/T^3            # energy input rate                  [m²/s³]
    nu_r  = nul_nd/T                   # Rayleigh friction                  [1/s]
    nu_h  = nuh_nd/T                   # radiative relaxation               [1/s]
    tau   = tau_nd*T                   # forcing decorrelation time         [s]

    nrot  = parse(Float64, get(ENV, "SP_NROT", "500"))   # length of the run, in rotations

    inputs = Dict(
        :lspherical_shell     => true,
        :lgrid_only           => false,
        :linit_only           => false,
        #---------------------------------------------------------------------------
        # Grid checks (see check_sphere_metrics in src/kernel/mesh/sphere_metrics.jl)
        #---------------------------------------------------------------------------
        :lcheck_grid             => true,
        :lstop_on_bad_grid       => true,
        :lproject_to_sphere      => true,
        # The .msh carries the Earth's radius; the reader rescales the shell to
        # the planet's. Everything metric follows from mesh.radius.
        :sphere_radius        => a,
        :cubed_sphere_map     => :none,        # the grid is already equiangular — see SWsphere
        :sphere_metrics       => :curl_invariant,
        #---------------------------------------------------------------------------
        # THE PLANET (read by initialize.jl, which also hands Ω and g to
        # user_source.jl / user_primitives.jl through their module globals).
        #---------------------------------------------------------------------------
        :sp_planet            => planet,
        :sp_figure            => P.fig,
        :sp_radius            => a,
        :sp_Omega             => Ω,
        :sp_gravity           => g,
        :sp_H                 => H,
        :sp_LD                => LD_nd,
        :sp_Ro_target         => Ro_target,
        :sp_epsilon0_nd       => eps0_nd,
        :sp_nul_nd            => nul_nd,
        :sp_nuh_nd            => nuh_nd,
        :sp_tau_nd            => tau_nd,
        :sp_f_up              => f_up,
        #---------------------------------------------------------------------------
        # FORCING AND LARGE-SCALE DISSIPATION (src/kernel/operators/sphere_forcing.jl)
        #
        #   ε₀   energy input rate; n_f, Δn the forced band |n - n_f| ≤ Δn/2
        #        (paper Eq. 11: Δn = 4, n_f = N/4);
        #   τ    Markovian decorrelation time (paper Eq. 13, c_r = 10 rotations;
        #        0 gives the δ-correlated forcing of Eq. 12);
        #   normalize  renormalise the forcing every step so that it injects
        #        exactly ε₀ (paper Eq. 14). Without it the fixed amplitude of
        #        Eq. (12)-(13) injected only 0.19 ε₀ on Jupiter at τ = 10:
        #        Coriolis and geostrophic adjustment decorrelate the flow from
        #        a slowly varying forcing long before τ;
        #   ν_r  Rayleigh friction on the momentum, ν_h radiative relaxation of
        #        the height towards the rest depth (paper Eq. 15 with δ_l = 0).
        #
        # The seed makes a run reproducible; change it for another realisation.
        #---------------------------------------------------------------------------
        :lsphere_forcing      => true,
        :forcing_epsilon      => eps0,
        :forcing_nf           => 24,
        :forcing_dn           => 4,
        :forcing_tau          => tau,
        :forcing_seed         => 1234,
        :forcing_normalize    => true,
        :rayleigh_friction    => nu_r,
        :radiative_relaxation => nu_h,
        :sphere_Omega         => Ω,            # for the Ro diagnostic and the PV output
        :sphere_gravity       => g,            # for the PV output q = (ζ + f)/h, Fig. 14 (a)-(c)
        #---------------------------------------------------------------------------
        # Keeping the flow ON the shell: the μx source of user_source.jl and the
        # P = I - xxᵀ/r² projection, applied after every RK stage.
        #---------------------------------------------------------------------------
        :llagrange_projection   => true,
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 5,
        #---------------------------------------------------------------------------
        # Mesh: the cubed sphere that ships with this case (600 quads, 602
        # vertices, 10 elements per panel edge) — the SWsphere grid.
        #---------------------------------------------------------------------------
        # SP_MESH in the environment points a run at another .msh (a finer
        # cubed sphere, say) without editing this line, like SP_PLANET.
        :lread_gmsh           => true,
        :gmsh_filename        => get(ENV, "SP_MESH", "./problems/ShallowWater/SWsphere_ScottPolvani_DynSGS/cubed_sphere.msh"),
        #---------------------------------------------------------------------------
        # Time integration. The paper integrates for 10⁴-10⁵ rotations; this deck
        # covers the spin-up (the energy grows linearly to t ≈ 500 rotations,
        # paper Fig. 1) at ~70 steps per rotation. Raise nrot for equilibrium.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK33(),
        :lcfl_dt              => true,           # take Δt from the CFL condition
        :cfl                  => 0.35,
        :tinit                => 0.0,
        :tend                 => nrot*T,
        :ndiagnostics_outputs => 50,             # a VTK dump every nrot/50 rotations
        :ndiagnostics_prints  => 500,            # steps between diagnostic lines
        :max_steps            => 5_000_000,
        :case                 => "swsphere_scottpolvani_dynsgs",
        :SOL_VARS_TYPE        => TOTAL(),
        :lsource              => true,
        #---------------------------------------------------------------------------
        # STABILISATION at the small scales: DynSGS.
        #
        # The sibling case uses a constant ν∇ₛ²(φu) plus the modal filter. Both
        # act everywhere at all times: ν has to be kept tiny so as not to damp
        # the forced band (its damping rate there is ν n_f(n_f+1)/a²), and the
        # filter then carries the grid scale on its own. The paper's own choice
        # is scale-selective (∇⁸ hyperdiffusion).
        #
        # DynSGS (Marras, Nazarov & Giraldo 2015, JCP 301:77-101, as run by the
        # flat cases with :visc_model => DSGS()) is selective in SPACE and TIME
        # instead: per element,
        #
        #   ν_e = max(0, min( C₂ Δ (|u|+√φ) ,  C₁ Δ² max_i ‖R_i‖∞,e/‖q_i-⟨q_i⟩‖∞ ))
        #
        # with R_i the strong-form residual of equation i (BDF2 time derivative
        # over the last two steps minus the assembled RHS) and Δ the LGL node
        # spacing of the element. Where the discrete equations are satisfied —
        # the resolved, cascading scales — ν_e is ~0; where they are not, it
        # rises, up to the wave-speed cap that keeps the explicit step stable.
        # A resting start has zero residual, so the run begins inviscid. The
        # coefficient is applied to the momentum only ([0,1,1,1], Marras
        # Eq. 10) through the same Laplace-Beltrami weak form as the constant
        # ν. Implementation: src/kernel/operators/sphere_dsgs.jl.
        #
        # The filter is OFF: DynSGS is meant to carry the grid scale by
        # itself. Switch it back on if a run shows grid-scale noise in ζ.
        #
        #   :dsgs_C1, :dsgs_C2   the paper's 1.0 and 0.5
        #   :dsgs_floor          floor of the normalisation, as a fraction of
        #                        the natural scales φ̄ and φ̄√φ̄ (see the module
        #                        header for why a resting layer needs one)
        #   :ldsgs_global_norms  domain norms instead of rank-local (costs
        #                        Allreduces every RHS call)
        #
        # The diagnostics print max/mean ν_e and the fraction of elements at
        # the cap every :ndiagnostics_prints steps, and the VTK files carry
        # dsgs_nu (the largest ν_e among the elements at each node).
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => DSGS(),
        :ivisc_equations      => [2, 3, 4],
        :μ                    => [0.0, 1.0, 1.0, 1.0],
        :dsgs_C1              => 1.0,
        :dsgs_C2              => 0.5,
        :dsgs_floor           => 1.0e-3,
        :ldsgs_global_norms   => false,
        :lfilter              => false,
        :filter_alpha         => 0.05,
        :filter_order         => 8,
        :filter_kcut          => 2/3,
        #---------------------------------------------------------------------------
        # Output: VTK with h, (u_zonal, v_merid, w_radial), ζ (Fig. 14 d-f), the
        # potential vorticity q = (ζ + f)/h (Fig. 14 a-c) and the vorticity
        # forcing; plus the paper's other diagnostic — the ZONAL-MEAN zonal
        # velocity against latitude (Fig. 13) — as a text file
        # zonal_mean_NNNN.dat per output.
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        :lzonal_mean          => true,
        :zonal_mean_nbins     => 90,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs

end
