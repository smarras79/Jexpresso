function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 150,

        #:ode_solver           => Tsit5(),   #28s

        :ode_solver           => SSPRK54(), #13.8s dt=0.35
        :Δt                   => 0.35,
        
        #:ode_solver           => SSPRK33(),  #14s dt=0.2
        #:Δt                   => 0.2,
        
        #:ode_solver           => ExplicitRK(tableau = constructDormandPrince()), #WORKING (but slow)
        #:ode_solver           => ExplicitRK(tableau = constructEuler()), #MEMORY BROKEN
        #:ode_solver           => ExplicitRK(tableau = constructSSPRK33()), #MEMORY BROKEN
        
        :ndiagnostics_outputs => 2,
        :case                 => "rtb",
        :lsource              => true,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        :luser_bc            => true,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default
        :νx                   => 30.0, #horizontal viscosity constant for momentum
        :νy                   => 30.0, #vertical   viscosity constant for momentum
        :κ                    => 60.0, #kinematic viscosity constant for θ equation
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
