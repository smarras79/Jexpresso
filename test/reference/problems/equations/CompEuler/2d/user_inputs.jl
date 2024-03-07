function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 300,
        :ode_solver           => SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :Δt                   => 0.2,
        :ndiagnostics_outputs => 2,
        :case                 => "rtb",
        :equationSet          => NCL(),
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        #:lexact_integration  => true,
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #:luser_bc            => true,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:lvisc                => true, #false by default NOTICE: works only for Inexact
        :νx                   => 30.0, #horizontal viscosity constant for momentum
        :νy                   => 30.0, #vertical   viscosity constant for momentum
        :κ                    => 60.0, #kinematic viscosity constant for θ equation
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB20x20.msh", #for nop=4
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB27x27.msh", #for nop=3
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB40x40.msh", #for nop=2
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "hdf5",
        :output_dir          => "./output/",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
