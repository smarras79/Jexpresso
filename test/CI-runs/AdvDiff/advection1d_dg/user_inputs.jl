function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # 1D scalar linear advection, discontinuous Galerkin.
        #
        #     ∂q/∂t + ∂(U q)/∂x = 0,   U = 2,   periodic on [-1, 1]
        #
        # This case is the Jexpresso-side twin of the validated Python 1D DG
        # reference. The configuration below
        # deliberately mirrors that reference so the two discretizations are
        # the *same discrete problem*: identical domain, wave speed, initial
        # condition, polynomial order, flux, and boundary treatment. That is
        # what makes an operator-level comparison meaningful rather than
        # merely an order-of-accuracy coincidence.
        #
        #   domain length 2, U = 2  =>  traversal period is exactly 1.0
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK33(),
        :tend                 => 1.0,      # exactly one traversal
        :Δt                   => 4.0e-4,
        :diagnostics_at_times => (0:0.1:1.0),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 4,         # matches the Python reference default
        :lexact_integration  => false,     # collocated LGL (lumped mass) -- the
                                           # convention the Python study calls
                                           # "collocated" (nq = ngl)
        :lsource             => false,
        :lperiodic_1d        => true,
        #---------------------------------------------------------------------------
        # Discontinuous Galerkin (DG) discretization
        #   :AD             => DiscGal()      selects the DG path (default ContGal)
        #   :numerical_flux => upwind_flux()  interface flux (also rusanov_flux())
        #---------------------------------------------------------------------------
        :AD                  => DiscGal(),
        :numerical_flux      => upwind_flux(),
        #---------------------------------------------------------------------------
        # Physical parameters/constants
        #   Viscous DG is not implemented (no _expansion_visc!(::DiscGal)) --
        #   :lvisc must stay false on the DG path.
        #---------------------------------------------------------------------------
        :lvisc               => false,
        #---------------------------------------------------------------------------
        # Mesh parameters and files
        #---------------------------------------------------------------------------
        :lread_gmsh          => false,     # 1D native mesh built by Jexpresso
        #---------------------------------------------------------------------------
        # Output formats: "png" -> plots, "ascii" -> data per npoin, "hdf5" -> .h5
        #   Switch to "hdf5" for machine-readable capture and for the CI ref.
        #---------------------------------------------------------------------------
        :outformat           => "png",
        :loverwrite_output   => true,
        :output_dir          => "./output",
        #---------------------------------------------------------------------------
        # 1D grid built by Jexpresso
        #   Python reference sweeps nelem = 8, 16, 32, 64 on this same domain.
        #---------------------------------------------------------------------------
        :xmin          =>  -1.0,
        :xmax          =>   1.0,
        :nelx          =>   50,
    ) #Dict

    return inputs

end
