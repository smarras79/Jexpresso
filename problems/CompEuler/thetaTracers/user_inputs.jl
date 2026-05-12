function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(),
        :Δt                   => 0.05,
        :tinit                => 0.0,
        :tend                 => 1000.0,
        :diagnostics_at_times => (0:100:1000),
        :restart_time         => 500,
        :lrestart             => false,
        :case                 => "rtb",
        :lsource              => true, 
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc          => true, #false by default NOTICE: works only for Inexact
        #:visc_model     => AV(),
        #:visc_model     => VREM(),
        :visc_model     => SMAG(),
        :μ              => [0.0, 1.0, 1.0, 2.0, 3.0, 3.0], #horizontal viscosity constant for momentum
        #:μ             => [0.0, 40.0, 40.0, 60.0, 60.0, 60.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10.msh", #for nop=4
        :gmsh_filename       => "./meshes/gmsh_grids/square_UNSTR_20el.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/square_UNSTR_5el.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk", #"hdf5",
        :loverwrite_output   => true,
        :output_dir          => "./output",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
