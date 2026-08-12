function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # SWsphere_visc — the Galewsky jet on the shell, stabilised by the
        # ARTIFICIAL DIFFUSION of Eq. (8b) instead of by the modal filter.
        #
        # Same equations, same grid, same initial condition as
        # ShallowWater/SWsphere — the five other user_*.jl files in this
        # directory are one-line includes of that case's. What differs is
        # exactly two switches at the bottom:
        #
        #     SWsphere        :lfilter => true    :lvisc => false
        #     SWsphere_visc   :lfilter => false   :lvisc => true, :μ => 1e5
        #
        # WHY IT IS A SEPARATE CASE. The run needs at least ONE of the two: the
        # inviscid unfiltered high-order solution on the cubed sphere grows
        # grid-scale modes and blows up (Marras, Kopera & Giraldo 2015, QJRMS
        # 141: 1727-1739, section 4.2). SWsphere keeps the filter because that
        # is the paper's own choice for the published test; this case exercises
        # the other branch, so both paths are covered by something runnable
        # rather than by a comment saying "set these two flags".
        #
        # Physics, provenance and the discretization: see the README next to
        # problems/ShallowWater/SWsphere/user_inputs.jl. Everything down to the
        # STABILIZATION block below is that deck, unchanged.
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
        # Initial condition: Galewsky, Scott & Polvani (2004), Tellus 56A,
        # 429-440 — barotropically unstable mid-latitude jet.
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
        # Mesh: the cubed sphere that ships with SWsphere (600 quads, 602
        # vertices, 10 elements per panel edge). Pointed at rather than copied —
        # the two cases are the same grid by construction.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/SWsphere/cubed_sphere.msh",
        #---------------------------------------------------------------------------
        # Time integration.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK33(),
        :lcfl_dt              => true,           # take Δt from the CFL condition
        :cfl                  => 0.35,
        :tinit                => 0.0,
        :tend                 => 518400.0,       # 144 h = 6 days, as in the test
        :ndiagnostics_outputs => 24,             # a VTK dump every 6 h
        :ndiagnostics_prints  => 200,            # steps between diagnostic lines
        :case                 => "swsphere_visc",
        :SOL_VARS_TYPE        => TOTAL(),
        :lsource              => true,
        #---------------------------------------------------------------------------
        # STABILIZATION — this is the whole point of the case.
        #
        # The artificial diffusion δν∇²(φu) of Eq. (8b), with the paper's
        # ν = 1e5 m²/s. On the shell this is the SURFACE (Laplace-Beltrami)
        # Laplacian, assembled in weak form from the manifold metrics in
        # src/kernel/operators/sphere_rhs.jl.
        #
        # :ivisc_equations => [2,3,4] is the MOMENTUM only, where the paper puts
        # it: the continuity equation ∂φ/∂t + ∇·(φu) = 0 carries no diffusion.
        #
        # ν costs no time step here. Diffusion is a second derivative, so it is
        # explicit-stable only for Δt ~ Δ²/ν, and sphere_cfl_dt takes the min of
        # that and the wave-speed limit — but at this resolution the two are
        # ~1e5 s and ~1e2 s, so the gravity waves still set Δt by three orders
        # of magnitude. Raising ν far enough WILL start to cut the step, which
        # is the intended behaviour: the alternative is a silent blow-up.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :ivisc_equations      => [2, 3, 4],
        :μ                    => 1.0e5,
        #---------------------------------------------------------------------------
        # ... and the modal filter is OFF, which is what makes this a test of the
        # diffusion rather than a test of the two together. Switching it back on
        # is legitimate — the two mechanisms compose, they are simply both
        # dissipative — but then the case no longer isolates anything.
        #---------------------------------------------------------------------------
        :lfilter              => false,
        :filter_alpha         => 0.05,
        :filter_order         => 8,
        :filter_kcut          => 2/3,
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
