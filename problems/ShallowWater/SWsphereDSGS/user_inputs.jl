function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # SWsphereDSGS — the Galewsky jet on the shell, stabilized by DynSGS.
        #
        # Same equations, same grid and the same initial condition as
        # ShallowWater/SWsphere. What changes is the STABILIZATION, and that is
        # the whole point of the case:
        #
        #   SWsphere       modal filter + constant ν = 2e5 m²/s everywhere
        #   SWsphereDSGS   the residual-based Marras-Nazarov DynSGS model,
        #                  ν = ν(x, t), NO filter, no tuned coefficient
        #
        # DynSGS measures how badly the discrete solution fails to satisfy the
        # PDE and puts viscosity exactly there:
        #
        #   μ_res|e = C1 Δ² max_i ‖R_i‖∞,e / ‖q_i − ⟨q_i⟩‖∞,Ω
        #   μ_max|e = C2 Δ  (‖u‖ + √φ)∞,e
        #   ν|e     = max(0, min(μ_max, μ_res))
        #
        # with R the BDF2 residual of each of the four equations. In smooth,
        # resolved regions R is small and so is ν; at a front it spikes and the
        # model damps there and nowhere else. See DSGS.md for the model and the
        # README next to this file for what it buys on this test.
        #---------------------------------------------------------------------------
        :lspherical_shell     => true,
        :lgrid_only           => false,
        :linit_only           => false,
        #---------------------------------------------------------------------------
        # Grid checks (see check_sphere_metrics in src/kernel/mesh/sphere_metrics.jl)
        #---------------------------------------------------------------------------
        :lcheck_grid             => true,
        :lstop_on_bad_grid       => true,
        :lproject_to_sphere      => true,
        :sphere_radius        => 6.37122e6,
        #---------------------------------------------------------------------------
        # Metric terms of the 2D manifold. See SWsphere/user_inputs.jl and the
        # header of sphere_metrics.jl; unchanged here.
        #---------------------------------------------------------------------------
        :sphere_metrics       => :curl_invariant,
        #---------------------------------------------------------------------------
        # Initial condition: Galewsky, Scott & Polvani (2004), Tellus 56A,
        # 429-440 — barotropically unstable mid-latitude jet. Identical to
        # SWsphere so the two runs are comparable field by field.
        #---------------------------------------------------------------------------
        :lgalewsky_perturbation => true,
        :galewsky_umax          => 80.0,
        :galewsky_phi0          => π/7,
        :galewsky_phi1          => π/2 - π/7,
        :galewsky_Omega         => 7.292e-5,
        :galewsky_gravity       => 9.80616,
        :galewsky_hmean         => 10000.0,
        :galewsky_hhat          => 120.0,
        :galewsky_alpha         => 1.0/3.0,
        :galewsky_beta          => 1.0/15.0,
        :galewsky_phi2          => π/4,
        :galewsky_nquad         => 512,
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
        # vertices, 10 elements per panel edge).
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/SWsphereDSGS/cubed_sphere.msh",
        #---------------------------------------------------------------------------
        # Time integration.
        #
        # :cfl is the SAME 0.35 as SWsphere, deliberately: the point of the
        # comparison is the stabilization, not the step. DynSGS costs nothing
        # here — sphere_cfl_dt takes the min of the wave-speed limit and the
        # diffusive limit Δ²/ν, and ν is bounded by the model's own
        # C2·Δ·(|u|+c) cap, which leaves the diffusive branch an order of
        # magnitude slacker than the gravity waves. The banner prints both.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK33(),
        :lcfl_dt              => true,           # take Δt from the CFL condition
        :cfl                  => 0.35,
        :tinit                => 0.0,
        :tend                 => 6*24*3600,      # 6 days: when the roll-up is judged
        :ndiagnostics_outputs => 24,             # a VTK dump every 6 h
        :ndiagnostics_prints  => 200,            # steps between diagnostic lines
        :case                 => "swsphere",
        :SOL_VARS_TYPE        => TOTAL(),
        :lsource              => true,
        #---------------------------------------------------------------------------
        # STABILIZATION — DynSGS, and nothing else.
        #
        # :visc_model => DSGS_SW() is what selects the residual-based model on
        # the shell. It is a separate tag from DSGS() because the shallow water
        # system has no equation of state: the wave speed in the cap is the
        # GRAVITY-wave speed √φ, and the coefficient is applied to the
        # CONSERVATIVE variables (ν∇ₛ²(φu), the paper's Eq. 8b) rather than to
        # primitives, so it carries no ρ factor. See src/kernel/physics/SGS.jl.
        #
        # :μ MEANS SOMETHING DIFFERENT HERE than in SWsphere. With DSGS_SW it is
        # NOT a viscosity in m²/s — the model supplies ν itself, per element and
        # per RK stage. It is the DIMENSIONLESS per-equation multiplier of
        # Marras eq. (10), i.e. how much of the one element coefficient each
        # equation gets. 1.0 = "as the model says". build_sphere_viscosity
        # rejects values above 100 so a 1e5 left over from a constant-ν deck is
        # an error rather than a 10⁵× overdose.
        #
        # :ivisc_equations => [1,2,3,4] — ALL FOUR, including continuity, which
        # is NOT where the paper puts its constant ν (it uses [2,3,4], momentum
        # only). The reason is the filter: SWsphere runs with the modal filter
        # ON, and the filter damps every field including φ. This deck has no
        # filter, and continuous Galerkin has no upwinding, so with [2,3,4] the
        # geopotential would carry no dissipation of any kind — which is
        # precisely the failure the error message at the end of
        # sphere_time_loop.jl warns about. Set it to [2,3,4] to reproduce the
        # paper's placement; then switch :lfilter back on.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => DSGS_SW(),
        :ivisc_equations      => [1, 2, 3, 4],
        :μ                    => 1.0,
        #---------------------------------------------------------------------------
        # The two model coefficients, Marras/Nazarov/Giraldo (JCP 2015) eq. (8):
        #
        #   :dsgs_C1  on the residual viscosity C1·Δ²·‖R‖/‖q−⟨q⟩‖
        #   :dsgs_C2  on the first-order-upwind cap C2·Δ·(|u|+c)
        #
        # These are the paper's own values, and they are NOT tuned per case —
        # that is what "parameter-free" means for this model: the sensor, not a
        # constant, decides where the viscosity goes. Raise C1 if the run needs
        # more dissipation everywhere the residual is already large; raising C2
        # only lifts the ceiling and does nothing while μ_res governs (which,
        # per the diagnostics line, it does here).
        #---------------------------------------------------------------------------
        :dsgs_C1              => 1.0,
        :dsgs_C2              => 0.5,
        #---------------------------------------------------------------------------
        # NO modal filter. DynSGS is meant to replace it, and leaving it on
        # would make the comparison against SWsphere meaningless — the filter is
        # what keeps that deck alive, so any result here would be its.
        #
        # Switch it on (with :filter_alpha => 0.05) only to run the two
        # mechanisms together, which is a legitimate configuration but a
        # different experiment.
        #---------------------------------------------------------------------------
        :lfilter              => false,
        :filter_alpha         => 0.05,
        :filter_order         => 8,
        :filter_kcut          => 2/3,
        #---------------------------------------------------------------------------
        # Plotting parameters. The VTK files carry, on top of the SWsphere
        # fields, one `mu_dsgs_<var>` point field per equation: the element
        # coefficient the model chose, broadcast to nodes. That is the field to
        # look at to see whether the viscosity is tracking the fronts.
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs

end
