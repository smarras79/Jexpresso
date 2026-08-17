function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 4.0,
        :ode_solver           => SSPRK33(), #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :Δt                   => 0.025,
        :ndiagnostics_outputs => 10,
        :diagnostics_at_times => collect(0.0:0.25:4.0),  # 17 frames: iter_1 (IC) … iter_17 (t=4); drives VTK cadence
        :lsource              => false,
        #:backend              => MetalBackend(),
        #:CL                   => NCL(), #CL() is defaults
        #:SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        :lexact_integration  => false,  # collocated LGL — required for the entrywise operator verification
        :AD                  => DiscGal(),
        :numerical_flux      => upwind_flux(),
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => false, # DG Phase 5: inviscid; the viscous path is flag-gated (rhs.jl), not μ-gated. 
        # NOTE: :lkep must stay unset/false under DiscGal — the KEP volume path
        # (_expansion_inviscid_KEP!) dispatches on AD and has no DiscGal method.
        :ivisc_equations      => [1],
        :μ                    => [0.0], #kinematic viscosity constant for θ equation
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x20_periodic.msh", #for nop=4
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "./output/",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
