function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                => 3.0/π,
        :Δt                  => 5e-4,
        :diagnostics_interval=> 10, #these are steps, not seconds
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",   # Choice: lgl, cgl 
        :lexact_integration  => false,
        :nop                 => 4,      # Polynomila order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :νx                   => 0.0, #kinematic viscosity constant
        :νy                   => 0.0, #kinematic viscosity constant
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_20x1.msh",
        #---------------------------------------------------------------------------
        # Boundary conditions:
        #---------------------------------------------------------------------------
        :xmin_bc       => "dirichlet", #Use either dirichlet or periodic
        :ymin_bc       => "dirichlet", #Use either dirichlet or periodic
        :zmin_bc       => "dirichlet", #Use either dirichlet or periodic
        :xmax_bc       => "dirichlet", #Use either dirichlet or periodic
        :ymax_bc       => "dirichlet", #Use either dirichlet or periodic
        :zmax_bc       => "dirichlet", #Use either dirichlet or periodic
        :bc_exact_xmin => [0.0 0.0 0.0],
        :bc_exact_xmax => [0.0 0.0 0.0],
        :bc_exact_ymin => [0.0 0.0 0.0],
        :bc_exact_ymax => [0.0 0.0 0.0],
        :bc_exact_zmin => [0.0 0.0 0.0],
        :bc_exact_zmax => [0.0 0.0 0.0],
        
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
