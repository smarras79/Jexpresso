function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 15000.0, #2π,
        :Δt                   => 0.1,#8.75e-4,
        :ode_solver           => SSPRK54(),
        :ndiagnostics_outputs => 2,
        :case                 => "rtb",
        #:CL                   => NCL(),
        :SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",   # Choice: lgl, cgl 
        :nop                 => 4,      # Polynomial order
        :luser_bc            => true,
        :lsource             => true,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => false, #false by default
        #:visc_model           => "dsgs", #"none", "dsgs"
        :νx                   => 0.000015, #kinematic viscosity constant
        :νy                   => 0.000015, #kinematic viscosity constant
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/agnesi240kmX30km_coarse.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_10x10_laguerre_top.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/agnesi-120kmx30km-hm5000.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_RTB.msh",
        :gmsh_filename         => "./meshes/gmsh_grids/hexa_TFI_120x31_top_lag.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/agnesi240kmX30km_coarse_laguerreTopLateral.msh",
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
