#---------------------------------------------------------------------------------
# fluxEmergenceSon2025DSGS — the flux-emergence problem of Son, Jang & Magara
# (2025, ApJS 277:46) with the SAME physics, mesh, initial and boundary
# conditions, time step and DynSGS constants as problems/MHD/fluxEmergenceSon2025,
# but NO positivity limiter: the Marras-Nazarov DynSGS dissipation alone keeps
# the solution admissible. Nothing is clipped, nothing is added; every term
# is in divergence form, so mass, momentum and energy are conserved.
#
# What is different from the sibling (see README.md and user_primitives.jl):
#   * :ode_solver is the plain CarpenterKennedy2N54 (no stage_limiter!);
#   * :dsgs_ref_weight => true. The conserved-form operator diffuses the
#     RELATIVE departure from the magnetostatic reference state,
#     ∇·(μ ρ_e ∇((q − q_e)/ρ_e)) on (ρ, ρu, ρv, E, ρw), instead of the
#     absolute one ∇·(μ ∇(q − q_e)) of the sibling. Both vanish at rest;
#     only the relative form obeys a maximum principle across the 25×
#     reference jump of the transition region, which is where the sibling's
#     operator lost positivity (the coefficient was at its cap and the run
#     still evacuated the coronal side of the contact — a −10% departure of
#     the dense side, diffused across, is more than the whole coronal
#     density). The magnetic slots keep the absolute form (B_e is not a
#     density-like weight).
# The physics files (initialize.jl, user_flux.jl, user_source.jl, user_bc.jl)
# include the sibling's, so the two cases cannot drift apart.
#---------------------------------------------------------------------------------
function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # 2D flux emergence in a two-temperature solar atmosphere (Parker
        # instability of a horizontal flux sheet), ideal GLM-MHD.
        #
        # D. Son, Y. Jang, T. Magara,
        # "A Comparative Analysis of High-resolution Shock-capturing Schemes
        #  for Two-dimensional Magnetohydrodynamic Simulation of Flux Emergence
        #  in the Solar Atmosphere", ApJS 277:46 (2025), Secs. 2.1-2.2 and 4.
        #
        # Run with (10 MPI ranks):
        #   mpiexec -n 10 julia --project=. src/Jexpresso.jl MHD fluxEmergenceSon2025DSGS
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),   # no stage limiter: DynSGS alone (see top of file)
        # Δt = 2.5e-3 τ₀ is CFL ≈ 0.07 against the initial maximum wave speed
        # (the coronal sound speed, = c_h ≈ 5.05 C_s) and the smallest LGL
        # spacing (0.173 H₀ at :nop => 4 on 1 H₀ elements). The margin is for
        # the late-time flows: the emerged loop reaches V_A ≈ 4-7 C_s and the
        # lateral downflows 4-5 C_s (paper Sec. 4.1), i.e. |v| + c_f ≈ 10 C_s
        # and CFL ≈ 0.15. The paper's own Courant number is 0.23.
        :Δt                   => 7.5e-3,
        :tinit                => 0.0,
        :tend                 => 54.0,  # paper Fig. 5 runs to t = 54 τ₀ (snapshots of Fig. 2 at t = 51 τ₀)
        :diagnostics_at_times => (0.0:1.0:54.0),
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,   # gravity + GLM ψ damping + absorbing layer (user_source.jl)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #
        # :nop => 4 on the 80×35 element grid below (1 H₀ × 1 H₀ elements)
        # gives 320 × 140 unique LGL points, a mean nodal spacing of 0.25 H₀:
        # the coarsest grid that still resolves the 0.5-0.6 H₀ tanh
        # transitions of the flux sheet and of the transition region with ~3
        # points. The paper's coarsest mesh is 300² cells (Δx = 0.27 H₀,
        # Δz = 0.12 H₀). See FE_80x35.geo to refine.
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 4,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #
        # The paper solves IDEAL MHD with fifth-order WENO/TENO finite-volume
        # reconstructions and an HLLD Riemann solver supplying the numerical
        # dissipation. A collocated continuous-Galerkin SEM has none, so the
        # fast/slow/intermediate shocks of the expanding loop (paper Sec. 4.1)
        # must be regularized explicitly: Marras-Nazarov DynSGS, the eddy
        # viscosity set by the LOCAL RESIDUAL of the governing equations, so
        # that it appears at the shocks and stays near zero elsewhere.
        #
        # :μ are per-equation multipliers on the DynSGS coefficient for
        # (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ), all at full strength. The
        # operator runs in its CONSERVED-VARIABLE form (:dsgs_conserved, see
        # user_primitives.jl and kernel/physics/SGS.jl): one kinematic
        # coefficient, a Laplacian on every conserved variable, no τ·u term.
        # The physical form of the Orszag-Tang case (u, v, T primitives) is
        # not usable here: the 25× density drop of the transition region
        # needs mass diffusion to stay positive, and mass diffusion under a
        # T-based energy closure drove p negative within a few τ₀ (measured),
        # while its κ∇T smeared the temperature jump and launched a 0.7 C_s
        # pulse into the corona. C1/C2 are Marras's residual and
        # wave-speed-cap coefficients; dsgs_gamma MUST match γ_mhd = 1.05 of
        # user_flux.jl (the DynSGS wave speed and pressure are built from it).
        #---------------------------------------------------------------------------
        :lvisc            => true,
        :μ                => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model       => DSGS_MHD(),
        :dsgs_C1          => 1.0,
        :dsgs_C2          => 0.5,
        # Background floor of 3% of the wave-speed cap: the residual sensor
        # cannot see a node-to-node mode (the discrete operator returns
        # nearly nothing on it), and in the corona above the rising crest
        # such a mode grew from 0.02 to 0.5 C_s between t = 8 and 11 τ₀ and
        # ended the run (measured). μ_floor = 0.03·Δ·c damps it at ≈ 7/τ₀
        # and spreads a resolved structure by ≈ 1 H₀ over the whole run.
        :dsgs_C0          => 0.03,
        :dsgs_gamma       => 1.05,
        :dsgs_Prt         => 0.7,
        # Stratification variants of the model (kernel/physics/SGS.jl): the
        # residual is normalized per element, not by the domain spread that
        # the 10⁸-times denser photosphere sets (the sensor was blind to the
        # corona and a grid-scale sawtooth grew across the transition
        # region), and the dynamic coefficient uses the nodal density, not
        # the element mean (which over-diffuses the light side of a
        # stratified element by ρ̄/ρ, up to 25 at the transition region,
        # past the explicit viscous limit).
        :dsgs_local_norms => true,
        :dsgs_local_rel   => 1.0,   # floor of the element spread = the local ρ, ρc, ρc², √ρc themselves (see SGS.jl)
        :dsgs_conserved   => true,  # Laplacian on the conserved variables (user_primitives.jl); implies no nodal-ρ scaling
        :dsgs_ref_weight  => true,  # ... of the RELATIVE departure (q − q_e)/ρ_e, coefficient μ·ρ_e (see top of file)
        :dsgs_nodal_rho   => true,  # (inactive with :dsgs_conserved; kept for the physical-form variant)
        :lrichardson      => false,      # gravity enters through user_source.jl, not the SGS closure
        # Slot 4 carries the TOTAL ENERGY ρE: "energy" keeps the kernel's τ·u
        # viscous-work augmentation of the energy equation active.
        :energy_equation  => "energy",
        #---------------------------------------------------------------------------
        # No entropy-stable / kinetic-energy-preserving machinery:
        #---------------------------------------------------------------------------
        :lkep              => false,
        :entropy_variables => false,
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #
        # [0, 80 H₀] × [0, 35 H₀], 80×35 quads, periodic in x, "bottom" and
        # "top" walls (user_bc.jl). The mesh ships with the case so that it
        # runs out of the box; regenerate it with
        #
        #   gmsh -2 problems/MHD/fluxEmergenceSon2025/FE_80x35.geo \
        #        -o problems/MHD/fluxEmergenceSon2025/FE_80x35.msh
        #
        # FE_80x70.msh (0.5 H₀ tall elements, the paper's vertical spacing)
        # also ships: point :gmsh_filename at it and halve :Δt.
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./problems/MHD/fluxEmergenceSon2025/FE_80x35.msh",   # shared with the sibling case
        #---------------------------------------------------------------------------
        # Filter parameters (off: DynSGS provides the dissipation; see the
        # Orszag-Tang case for why the Boyd-Vandeven filter is a poor
        # substitute on the conservative MHD variables).
        #---------------------------------------------------------------------------
        :lfilter             => false,
        #---------------------------------------------------------------------------
        # Output: PNG figures styled after the paper, written directly by the
        # solver at every diagnostics time (gathered on rank 0 under MPI).
        #
        #   ρ-it<n>.png      log₁₀(ρ/ρ₀) on the paper's "jet" scale [-8.5, 0]
        #                    with magnetic field lines (black, isocontours of
        #                    the vector potential A_y) and velocity vectors
        #                    (white, reference arrow = 5 C_s) — paper Fig. 2
        #   <var>-it<n>.png  v (= V_z/C_s), vA (= V_A/C_s), Bx, p, T, β
        #   profile-it<n>.png  vertical profiles at x = X_max/2 of V_z/C_s,
        #                    V_A/C_s, log₁₀(B_x/B₀), log₁₀(ρ/ρ₀) on the axes of
        #                    the paper's Fig. 5, with z_cor = 18 H₀ marked
        #
        # Switch to :outformat => "vtk" for ParaView output (plus the
        # mu_dsgs_<var> DynSGS fields).
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :plot_matrix         => false,        # silent per-variable PNGs, no GR window
        :plot_colormap       => :jet,         # the paper's colormap (Fig. 2)
        :plot_vars           => ["ρ", "v", "vA", "Bx", "p", "T", "β"],
        :plot_log10          => ["ρ", "p", "β"],
        :plot_clims          => Dict("ρ" => (-8.5, 0.0)),   # paper Fig. 2 colorbar
        :plot_fieldlines     => ("Bx", "By"),
        :plot_fieldlines_levels => 40,
        :plot_vectors        => ("u", "v"),
        :plot_vectors_ref    => 5.0,          # paper Fig. 2: reference arrow "= 5.0"
        :plot_vectors_n      => (30, 13),
        :plot_overlay_on     => ["ρ"],        # field lines/vectors only on the density panel
        :plot_xlabel         => "X/H₀",
        :plot_ylabel         => "Z/H₀",
        :plot_time_unit      => " τ₀",
        # The DynSGS coefficient actually applied, log10_μ_dsgs_ρ-it<n>.png:
        # log₁₀ of the kinematic μ (H₀ C_s), floored at 1e-6. In the conserved
        # form every slot carries the same coefficient, so the ρ slot stands
        # for all nine. (:outformat => "vtk" writes all nine as mu_dsgs_<var>.)
        :plot_dsgs           => true,
        :plot_dsgs_vars      => ["ρ"],
        :plot_dsgs_log10     => true,
        :plot_dsgs_floor     => 1.0e-6,
        :plot_profile_x      => 40.0,         # paper Fig. 5: x = X_max/2
        :plot_profile_vars   => ["v", "vA", "Bx", "ρ"],
        :plot_profile_log10  => ["Bx", "ρ"],
        :plot_profile_ylims  => Dict("v" => (0.0, 1.7), "vA" => (0.0, 4.5), "Bx" => (-2.0, 1.0), "ρ" => (-8.5, 1.0)),
        :plot_profile_vlines => [18.0],       # z_cor
        :loverwrite_output   => false,
        :lwrite_initial      => true,
        :output_dir          => "./output/",
        :loutput_pert        => false,
        #---------------------------------------------------------------------------
        # AMR (off)
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :ladapt              => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
