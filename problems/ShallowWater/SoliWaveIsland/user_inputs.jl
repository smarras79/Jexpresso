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
        :Δt                   => 0.01,
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
        # JACC.jl portable RHS.
        #
        # Flip to true to evaluate the right-hand side through
        # src/kernel/operators/rhs_jacc.jl — ONE kernel body that runs on
        # multi-core CPU, NVIDIA, AMD, Intel or Apple hardware. Everything else
        # about the run is unchanged: same mesh, same Float64 setup, same
        # SSPRK54, same output, and :backend stays at its CPU() default (that
        # switch selects the separate KernelAbstractions path and the two must
        # not both be on).
        #
        # Which device JACC uses is NOT a deck setting — it is a Preferences
        # value written once per project:
        #
        #     julia> using JACC; JACC.set_backend("cuda")   # amdgpu|oneapi|metal|threads
        #     (restart Julia)
        #
        # On the default "threads" backend, ask for threads on the command line:
        #
        #     julia --project -t 8 src/Jexpresso.jl ShallowWater SoliWaveIsland
        #
        # Measured on this deck (5 steps of SSPRK54, HDF5 output diffed against
        # the serial run): H and Hu agree to 4e-16 and 8e-16 RELATIVE, i.e. to
        # floating-point round-off. Not bit-for-bit — the direct stiffness
        # summation adds the same numbers in a different order — but two JACC
        # runs of this deck do agree to the last bit, on CPU and GPU alike.
        # problems/JACCtestz/{README.md, JACC_2d_swe.jl} have the details.
        #---------------------------------------------------------------------------
        :ljacc                => false,
        #---------------------------------------------------------------------------
        # Mesh
        # Generate with:
        #   gmsh -2 problems/ShallowWater/SoliWaveIsland/SoliWaveIsland.geo \
        #        -o meshes/gmsh_grids/SoliWaveIsland.msh
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/SoliWaveIsland/SoliWaveIsland.msh",
        #---------------------------------------------------------------------------
        # Plotting / output: one PNG per variable at every diagnostic time
        # (H-it<n>.png, Hu-it<n>.png, Hv-it<n>.png). Set :lplot_surf3d to
        # true for the Spline2D surface rendering instead of the nodal map,
        # or switch back to "vtk" for ParaView output.
        #---------------------------------------------------------------------------
        :outformat            => "png",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs

end
