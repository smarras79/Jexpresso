function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(),
        :Δt                   => 0.5,
        :tinit                => 0.0,
        :tend                 => 1000.0,
        :diagnostics_at_times => (0:100:1000),
        #:restart_input_file_path => "problems/equations/CompEuler/theta",
        :case                 => "rtb",
        :lsource              => true,
        #:SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc          => true, #false by default NOTICE: works only for Inexact
        :visc_model     => AV(),
        # :visc_model     => VREM(),
        #:visc_model     => SMAG(),
        :energy_equation => "theta",
        # :μ              => [0.0, 1.0, 1.0, 2.0], #horizontal viscosity constant for momentum
        :μ              => [0.0, 125.0, 125.0, 125.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10.msh", #for nop=4
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter             => true,
        #:mu_x                => 0.01,
        #:mu_y                => 0.01,
        #:filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => true,
        :lwrite_initial      => true,
        :output_dir          => "./output",
        #:output_dir          => "./test/CI-run",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 1,
        #---------------------------------------------------------------------------
        # preadapt: refine the middle of the domain (a box around the mesh
        # center) up to :preadapt_max_level before t=0, via
        # user_get_preadapt_flags! in initialize.jl. See docs/amr_setup.md
        # section 3.
        #---------------------------------------------------------------------------
        :lpreadapt           => true,
        :preadapt_max_level  => 1,
        #---------------------------------------------------------------------------
        # AMR
        #---------------------------------------------------------------------------
        :lamr              => true,
        #---------------------------------------------------------------------------
        # AMR parameters
        #---------------------------------------------------------------------------
        :amr_freq            => 200,
        :amr_max_level       => 2,
        #---------------------------------------------------------------------------
        # AMR restart: resume from iter_5 (t=400, given diagnostics_at_times
        # = 0:100:1000 and iter_1 = the t=0 IC write). Loads the p4est forest
        # saved at iter_5/iter_5.p4est and the matching VTK solution data.
        # See docs/amr_setup.md section 4.
        #---------------------------------------------------------------------------
        :lrestart_amr        => false,
        :restart_vtk_iout    => 5,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs

end
