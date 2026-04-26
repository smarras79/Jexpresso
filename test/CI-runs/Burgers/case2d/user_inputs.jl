function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # CI version: reduced tend for fast execution (5 time steps)
        #
        # 2D Burgers Riemann problem (flux form)
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :tend                 => 5.0e-3,
        :Δt                   => 1.0e-3,
        :diagnostics_at_times => (5.0e-3,),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 7,
        :lexact_integration  => false,
        :lsource             => false,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc               => true,
        :μ                   => [1.0e-2, 1.0e-2],
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10_burgers2d.msh",
        #---------------------------------------------------------------------------
        # Output formats
        #---------------------------------------------------------------------------
        :outformat           => "hdf5",
        :loverwrite_output   => true,
        :output_dir          => "none",
    ) #Dict

    return inputs
end
