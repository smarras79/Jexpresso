function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 1.0,
        :ode_solver           => SSPRK33(), #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :Δt                   => 0.2,
        :ndiagnostics_outputs => 20,
        :lsource              => false,
        #:CL                   => NCL(), #CL() is defaults
        :SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :νx                   => 30.0, #horizontal viscosity constant for momentum
        :νy                   => 30.0, #vertical   viscosity constant for momentum
        :κ                    => 60.0, #kinematic viscosity constant for θ equation
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x20_periodic.msh", #for nop=4
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "png",
        :output_dir          => "./output/",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
