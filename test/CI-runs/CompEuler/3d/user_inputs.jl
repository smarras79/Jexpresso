function user_inputs()
    inputs = Dict(
    #---------------------------------------------------------------------------
    # CI version: reduced tend for fast execution (2 time steps)
    #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 0.5,
        :tinit                => 0.0,
        :tend                 => 1.0,
        :diagnostics_at_times => (1.0,),
        :lsource              => true,
        :lrestart             => false,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => VREM(),
        :μ                    => [0.0, 1.0, 1.0, 1.0, 2.0],
        :energy_equation     => "theta",
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lwarp => false,
        :lstretch => false,
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x1x10.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "none",
        :loverwrite_output   => true,
        :loutput_pert        => false,
        #---------------------------------------------------------------------------
        # AMR (disabled for CI)
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 1,
        :ladapt              => false,
        :lamr                => false,
        :amr_freq            => 100,
        :amr_max_level       => 1,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
