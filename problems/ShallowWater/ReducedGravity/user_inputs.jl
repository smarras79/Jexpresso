function user_inputs()

    day = 86400.0 #seconds in a day
    year = 360.0 * day #seconds in a year (model year of 360 days)

    inputs = Dict(
        #---------------------------------------------------------------------------
        # ReducedGravity: 2D non-linear shallow water equations with reduced gravity.
        #
        # Based on Sec 3.1 of Choi et al. (2004), "A numerical study 
        # of the wind-driven circulation in a rectangular basin with flat bottom", 
        # J. Geophys. Res., 109, C03023, doi:10.1029/2003JC002007.
        #
        # A widnd driven double gyre circulation in a square basin with flat bottom, 
        # assuming that the simulated layer of fluid rests on a deep, motionless layer.
        # Solved with CG SEM (no Riemann flux).
        #
        # Conservative variables:
        #   q = [H, Hu, Hv]
        # where H is the local water depth.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 1200.0, 
        :tinit                => 0.0,
        :tend                 => 6.0 * year, #6 years in seconds
        :diagnostics_at_times => (0: 30*day : 6.0*year),
        :case                 => "double_gyre_reduced_gravity",
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
        # Note: the continuity equation diffuses the depth perturbation
        # H - He (see user_primitives.jl), not the full cone-shaped depth,
        # so the lake at rest stays an exact equilibrium.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => AV(),
        :μ                    => [0.05, 0.05, 0.05],
        #---------------------------------------------------------------------------
        # CG filter: OFF. The filter acts on the full depth H (it only
        # subtracts qe from the momentum components), and at rest H is
        # cone-shaped with a kink at the wet/dry ring: re-projecting it
        # every step perturbs the lake-at-rest equilibrium that the
        # well-balanced flux/source split preserves. The constant
        # artificial viscosity above is enough to stabilise the fronts.
        #---------------------------------------------------------------------------
        :lfilter              => false,
        #---------------------------------------------------------------------------
        # Mesh
        # Generate with:
        #   gmsh -2 problems/ShallowWater/ReducedGravity/DoubleGyre.geo \
        #        -o meshes/gmsh_grids/DoubleGyre.msh
        #---------------------------------------------------------------------------

        #NOTE - need to create that geo file and mesh for the double gyre case, and check if it is appropriate for the time step I want to use
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/ReducedGravity/DoubleGyre.msh",
        #---------------------------------------------------------------------------
        # Plotting / output: one PNG per variable at every diagnostic time
        # (H-it<n>.png, Hu-it<n>.png, Hv-it<n>.png). Set :lplot_surf3d to
        # true for the Spline2D surface rendering instead of the nodal map,
        # or switch back to "vtk" for ParaView output.
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => true,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs

end
