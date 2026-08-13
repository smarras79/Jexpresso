function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # SWsphere — shallow water equations on a spherical shell.
        #
        # Galewsky, Scott & Polvani (2004) barotropically unstable jet, solved
        # with the shallow water equations of Marras, Kopera & Giraldo (2015),
        # QJRMS 141: 1727-1739, Eq. (8), on a cubed-sphere shell.
        #
        # The run: grid -> manifold metrics -> initial condition -> :ode_solver.
        #
        #   :lspherical_shell => true   solve on the shell: manifold metrics,
        #                               surface-divergence RHS and the shell
        #                               time loop replace the flat ones. The GRID is
        #                               read by the ordinary gmsh path, which
        #                               detects a 2D manifold embedded in 3D
        #                               and keeps z (`lmanifold` in mesh.jl).
        #
        #   :lgrid_only       => true   stop right after the grid, without
        #                               touching initialize.jl.
        #   :linit_only       => true   stop after the initial condition,
        #                               without integrating. :lgrid_only wins if
        #                               both are set.
        #---------------------------------------------------------------------------
        :lspherical_shell     => true,
        :lgrid_only           => false,
        :linit_only           => false,
        #---------------------------------------------------------------------------
        # Grid checks (see check_sphere_metrics in src/kernel/mesh/sphere_metrics.jl).
        #
        # A spherical shell is CLOSED: no boundary, every edge shared by exactly
        # two elements, and gmsh files built panel-by-panel routinely carry
        # duplicated nodes along the panel seams. These switches control what is
        # done about that.
        #
        #                             before the topology is built (this is what
        #                             turns "six panels that look joined" into
        #                             one watertight shell)
        #   :lcheck_grid              run the T1..T9 consistency tests
        #   :lstop_on_bad_grid        abort if any of them fails
        #   :lproject_to_sphere       snap every node radially onto the shell
        #---------------------------------------------------------------------------
        :lcheck_grid             => true,
        :lstop_on_bad_grid       => true,
        :lproject_to_sphere      => true,
        # Radius of the shell. Setting it explicitly RESCALES the grid onto a
        # sphere of that radius (a unit-sphere .msh is blown up to the Earth);
        # comment it out to take the radius from the gmsh file instead, as the
        # mean |x| over its nodes.
        #
        # Pinned here to the Galewsky Earth radius so the test is the published
        # one whatever the .msh happens to carry.
        :sphere_radius        => 6.37122e6,
        #---------------------------------------------------------------------------
        # Initial condition: Galewsky, Scott & Polvani (2004), Tellus 56A,
        # 429-440 — barotropically unstable mid-latitude jet.
        #
        # All of these are optional: initialize.jl falls back to the published
        # values shown. They are spelled out so a parameter sweep never needs to
        # touch the source.
        #
        #   the jet          u_max, and the band φ ∈ [φ0, φ1]
        #   the balance      Ω, g, and the global mean depth h_mean that fixes
        #                    the constant of integration h₀
        #   the trigger      ĥ cos φ exp[-(λ/α)²] exp[-((φ2-φ)/β)²]
        #   :galewsky_nquad  Simpson intervals for the balance integral. 512
        #                    gives ~5e-8 m on a 10 km depth; the answer is
        #                    unchanged to 10 digits from 128.
        #
        # Set :lgalewsky_perturbation => false for the BALANCED control run: an
        # exact steady solution, which is the natural first test of the scheme.
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
        # Keeping the flow ON the shell — Marras, Kopera & Giraldo (2015),
        # QJRMS 141: 1727-1739, section 3.2 (after Coté 1988).
        #
        # The velocity is a full 3-D Cartesian vector on a 2-D surface, so one
        # degree of freedom is redundant and must be constrained: u·x = 0.
        # Two complementary mechanisms:
        #
        #   the μx SOURCE          in user_source.jl, with the closed form
        #                          μ = -φ|u|²/r² (a centripetal force). Keeps
        #                          the CONTINUOUS equation consistent. It has no
        #                          switch here because user_source! receives no
        #                          `inputs` — comment the two μ lines out in
        #                          user_source.jl to run without it.
        #
        #   :llagrange_projection  the P = I - xxᵀ/r² projection of Eq.(9)-(11),
        #                          applied to the momentum. Removes the normal
        #                          drift the DISCRETE operators accumulate. In a
        #                          real run this belongs at the end of every
        #                          step; with :linit_only it is applied once,
        #                          after the initial condition, and reports the
        #                          drift it removed.
        #---------------------------------------------------------------------------
        :llagrange_projection   => true,
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 5,
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #
        # A CLOSED all-quadrilateral surface grid of the sphere. The one that
        # ships with this case is a 10x10-per-panel cubed sphere (600 quads,
        # 602 vertices, V - E + F = 2). Regenerate it with gmsh from the .geo
        # shipped next to it:
        #
        #   gmsh -2 problems/ShallowWater/SWsphere/cubed_sphere.geo \
        #        -o problems/ShallowWater/SWsphere/cubed_sphere.msh
        #
        # Any gmsh format works (MSH 2.2 or 4.1, parametric nodes or not): the
        # file is read by GridapGmsh, i.e. by gmsh itself.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/SWsphere/cubed_sphere.msh",
        #---------------------------------------------------------------------------
        # Time integration.
        #
        # src/kernel/solvers/sphere_time_loop.jl builds an ODEProblem and solves
        # it with :ode_solver, like every other case. The Lagrange projection has
        # to run after EVERY RK stage, which is exactly what an integrator's
        # stage_limiter! hook is for; the modal filter is the step_limiter!.
        # :ode_solver must therefore name an integrator that takes limiters.
        #
        # :lcfl_dt => true (the default) takes the step from the CFL condition,
        # Δt = :cfl * Δmin / max(|u| + sqrt(φ)), with Δmin the smallest LGL node
        # spacing. Set it false to use a :Δt of your own.
        #
        # It is an explicit switch rather than "omit :Δt and you get the CFL
        # step" because mod_inputs_user_inputs! fills a missing :Δt with 0.1 s,
        # so an absent :Δt cannot be told from a deliberate one — and taking
        # 0.1 s literally turns this one-day run into 864 000 steps.
        #
        # The instability needs TIME: the perturbation is linearly unstable but
        # grows slowly, and the jet only rolls up into the published train of
        # vortices around day 4-6. A one-day run shows a laminar jet and looks
        # like nothing is happening. :tend is therefore the full 144 h of the
        # test (Galewsky et al. 2004; Marras et al. section 5).
        #---------------------------------------------------------------------------
        # The integrator is now HONOURED: sphere_time_loop.jl builds an
        # ODEProblem and hands it to OrdinaryDiffEq, with the Lagrange
        # projection as the stage_limiter! and the modal filter as the
        # step_limiter!. Any SSPRK integrator works — it has to be one that
        # takes limiters, because the projection must run after every stage.
        # SSPRK33 is the Shu-Osher SSP-RK3 of Marras et al. (2015) section 4.2,
        # i.e. the published scheme for this test; SSPRK54 costs 5 RHS
        # evaluations per step instead of 3 and permits a larger CFL.
        :ode_solver           => SSPRK33(),
        :lcfl_dt              => true,           # take Δt from the CFL condition
        :cfl                  => 0.35,
        :tinit                => 0.0,
        :tend                 => 518400.0,       # 144 h = 6 days, as in the test
        :ndiagnostics_outputs => 24,             # a VTK dump every 6 h
        :ndiagnostics_prints  => 200,            # steps between diagnostic lines
        :case                 => "swsphere",
        :SOL_VARS_TYPE        => TOTAL(),
        :lsource              => true,
        #---------------------------------------------------------------------------
        # Artificial viscosity, the δν∇²(φu) term of Eq. (8b), with ν the paper's
        # 1e5 m²/s. On the shell this is the SURFACE (Laplace-Beltrami)
        # Laplacian, assembled in weak form from the manifold metrics in
        # src/kernel/operators/sphere_rhs.jl — the flat viscous path of rhs.jl is
        # built on the 2×2 inverse Jacobian and has no third metric direction, so
        # it cannot be used here.
        #
        # :ivisc_equations selects which equations are diffused. [2,3,4] is the
        # MOMENTUM only, which is where the paper puts it: the continuity
        # equation ∂φ/∂t + ∇·(φu) = 0 carries no diffusion. Adding 1 diffuses
        # the geopotential too; it makes almost no difference to stability here
        # (measured — see problems/ShallowWater/SWsphere_visc/README.md).
        #
        # IT IS AN ADDITION TO :lfilter, NOT A SUBSTITUTE FOR IT. At this
        # resolution ν = 1e5 with the filter off blows up at 2 days: the coarse
        # elements sit at ~219 km, i.e. a grid Reynolds number of ~175, and
        # viscosity alone does not hold that until about ν = 5e5 — which by then
        # has erased half of max|ζ|, the field the test is judged on. The paper
        # runs the filter AND ν = 1e5, and so should you; ShallowWater/
        # SWsphere_visc is that configuration, verified to 3 days.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :ivisc_equations      => [2, 3, 4],
        :μ                    => 2.0e5,      # set to 1.0e5 together with :lvisc => true
        #---------------------------------------------------------------------------
        # Stabilization: the modal filter, the stand-in for the Boyd-Vandeven
        # filter of the paper's section 4.2 ("chosen to reduce by 5% the highest
        # modes only" -- hence :filter_alpha => 0.05). Applied after every step
        # and followed by a mass-weighted DSS average, so it conserves ∫φ.
        #
        # The paper shows the inviscid solution is badly resolution-sensitive on
        # the cubed sphere, so LEAVE THIS ON. Switching on the artificial
        # viscosity above is not a licence to switch this off: at this
        # resolution that combination fails at 2 days (see the note there). With
        # both off the run blows up sooner still, and the time loop warns at
        # startup that it will.
        #---------------------------------------------------------------------------
        :lfilter              => false,
        :filter_alpha         => 0.05,     # damping of the HIGHEST mode
        :filter_order         => 8,        # exponent: larger = sharper cutoff
        :filter_kcut          => 2/3,      # first damped mode, as a fraction of nop
        #---------------------------------------------------------------------------
        # Plotting parameters
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
