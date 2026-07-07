function user_inputs()
    inputs = Dict(
    #---------------------------------------------------------------------------
    # User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.2,
        :tinit                => 0.0,
        :tend                 => 500.0,
        :diagnostics_at_times => (0:100:1000.0),
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # restart options
        #---------------------------------------------------------------------------
        :lrestart             => false,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:lwall_model          => true,
        :lvisc                => true, #false by default
        :visc_model           => AV(), #VREM(), #SMAG(),
        # :visc_model           => SMAG(),
        # :visc_model           => VREM(),
        :μ                    => [0.0, 60.0, 60.0, 60.0, 60.0],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        # Warping:
        :lwarp => false,
        :mount_type => "agnesi",
        :a_mount => 4000.0,
        :h_mount => 1000.0,
        :c_mount => 5000.0,

        # Stretching factors:
        :lstretch => false,
        :stretch_factor => 1.5,
        :stretch_type => "fixed_first_twoblocks_strong", #strong means that the top is constrained
        :first_zelement_size => 250.0,
        :zlevel_transition => 5000.0,

        # GMSH files:
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x1x10.msh",
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        :lfilter             => false,
        :mu_x                => 0.5,
        :mu_y                => 0.5,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "./output/",
        :loverwrite_output   => true,
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 1,
        #---------------------------------------------------------------------------
        # preadapt: refine the middle of the domain (a box around the mesh
        # center in x,z) up to :preadapt_max_level before t=0, via
        # user_get_preadapt_flags! in initialize.jl. See docs/amr_setup.md
        # section 3.
        #---------------------------------------------------------------------------
        :lpreadapt           => true,
        :preadapt_max_level  => 1,
        #---------------------------------------------------------------------------
        # AMR
        #---------------------------------------------------------------------------
        :lamr                 => true,
        #---------------------------------------------------------------------------
        # AMR parameters
        #---------------------------------------------------------------------------
        :amr_freq            => 250,
        :amr_max_level       => 2,
        #---------------------------------------------------------------------------
        # AMR restart: resume from iter_6 (t=250, given diagnostics_at_times
        # = 50:50:500 and iter_1 = the t=0 IC write). Loads the p4est forest
        # saved at iter_6/iter_6.p4est and the matching VTK solution data.
        # See docs/amr_setup.md section 4.
        #---------------------------------------------------------------------------
        :lrestart_amr        => false,
        :restart_vtk_iout    => 6,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
end
