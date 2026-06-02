function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # Stratified-atmosphere rising thermal bubble (Robert, 1993; Giraldo &
        # Restelli, 2008) integrated with the conservation-form Euler-θ
        # equations and stabilised by the Marras et al. Dynamic SGS model.
        #
        # The case is identical in setup to problems/CompEuler/theta — same
        # mesh, source, fluxes, BCs — but uses :SOL_VARS_TYPE => TOTAL() so
        # the residual-based DSGS coefficient is built from the full
        # conservative variables (the perturbation branch in compute_dsgs_
        # viscosity! would have to add back qe to recover them).
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 0.5,
        :tinit                => 0.0,
        :tend                 => 1000.0,
        :diagnostics_at_times => (0:100:1000),
        :lrestart             => false,
        :case                 => "rtb",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => DSGS(),     # Marras et al. Dynamic SGS
        :energy_equation      => "theta",
        # Per-equation multiplier on the DSGS coefficient. The Marras
        # paper value is 1.0; turn an equation off with 0.0; throttle a
        # too-aggressive coefficient with a value in (0, 1). The ρ
        # entry stays at 0.0 to keep the mass equation conservative
        # (Marras eq. 10).
        :μ                    => [0.0, 1.0, 1.0, 1.0],
        :Pr                   => 0.1,        # artificial Prandtl number (Marras eq. 10b)
        #---------------------------------------------------------------------------
        # Mesh
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_10x10.msh",
        #---------------------------------------------------------------------------
        # Plotting
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,      # plot total state, not perturbation
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :init_refine_lvl      => 1,
        #---------------------------------------------------------------------------
        # AMR
        #---------------------------------------------------------------------------
        :ladapt               => false,
        :amr_freq             => 200,
        :amr_max_level        => 2,
    ) #Dict

    return inputs

end
