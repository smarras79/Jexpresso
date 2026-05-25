function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # SoliWaveIsland: 2D non-linear shallow water equations.
        #
        # Section 5.5 of Marras, Kopera, Constantinescu, Suckale, Giraldo,
        # "A residual-based shock capturing scheme for the continuous/
        #  discontinuous spectral element solution of the 2D shallow water
        #  equations", Advances in Water Resources 114 (2018) 45-63.
        #
        # A right-propagating solitary wave runs up and runs down a circular
        # conical island that sits in the middle of a closed rectangular basin.
        # Solved with CG SEM (no Riemann flux).
        #
        # Conservative variables:
        #   q = [H, Hu, Hv]
        # where H is the local water depth, Hb(x,y) is the (time-independent)
        # bathymetry, and the still-water depth far from the island is h0.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 0.02,
        :tinit                => 0.0,
        :tend                 => 25.0,
        :diagnostics_at_times => (0:1.0:25.0),
        :case                 => "soliwave_island",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Physical parameters / artificial viscosity.
        # We use a small constant Laplacian viscosity on all three equations
        # to keep the CG solution stable around the wave front. This is a
        # placeholder for the Dyn-SGS scheme of Marras et al. (2018), which
        # would replace the constant μ with a residual-based coefficient.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => AV(),
        :μ                    => [0.05, 0.05, 0.05],
        #---------------------------------------------------------------------------
        # CG filter for additional stabilisation.
        #---------------------------------------------------------------------------
        :lfilter              => true,
        :mu_x                 => 0.05,
        :mu_y                 => 0.05,
        :filter_type          => "erf",
        #---------------------------------------------------------------------------
        # Mesh
        # Generate with:
        #   gmsh -2 problems/ShallowWater/SoliWaveIsland/SoliWaveIsland.geo \
        #        -o meshes/gmsh_grids/SoliWaveIsland.msh
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/SoliWaveIsland.msh",
        #---------------------------------------------------------------------------
        # Plotting / output
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
