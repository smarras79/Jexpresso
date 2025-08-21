
function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 2,
        :ode_solver           => SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :Δt                   => 1e-3,
        :ndiagnostics_outputs => 10,
        :case                 => "rtb",
        :lsource              => false,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 3,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :luser_bc             => true,
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :ivisc_equations      => [1],
        :μ                   => [0.01], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_2x2_domain01.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk", #"hdf5",
        :loverwrite_output   => true,
        :output_dir          => "./output",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        :lwrite_initial      => true
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
