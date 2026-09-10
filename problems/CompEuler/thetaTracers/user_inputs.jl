function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 1000.0,
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #MSRK5(), #SSPRK54(),
        :Δt                   => 0.2,
        :diagnostics_at_times => (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
        :case                 => "rtb",
        :lsource              => true,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        #:visc_model     => AV(),
        #:visc_model     => VREM(),
        #:visc_model     => SMAG(),
        :visc_model     => DSGS(),  # residual-based DynSGS (DSGS.md §3): :μ are multipliers of
                                     # its ν, which sits at the cap 0.5·Δ·c ≈ 1e4 m²/s at the
                                     # bubble's cone edge from the first step. On this mesh
                                     # (Δx_min = 49 m) that needs :Δt => 0.1 (ν·Δt/Δx² ≈ 0.4;
                                     # 0.2 blows up at the first step, the run warns about it),
                                     # or :μ[5] => 1.0 at Δt = 0.2. :Pr (artificial Prandtl
                                     # number of the θ slot) defaults to 0.1.
        :energy_equation => "theta",
        #:μ                   => [0.0, 1.0, 1.0, 2.0, 3.0, 1.0], #horizontal viscosity constant for momentum
        :μ                   => [0.0, 1.0, 1.0, 1.0, 1.0, 1.0], #horizontal viscosity constant for momentum
        #:μ                   => [0.0, 40.0, 40.0, 60.0, 60.0, 60.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/square_UNSTR_20el.msh",
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
