# FV version of ../advection2d_dg: same Gaussian, c = (0.5, 1), [-5,5]×[0,20] periodic, t = 4,
# on 40 × 80 cells. For linear advection all fluxes reduce to exact upwind.
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
        :AD                  => FV(),         # finite volumes = order-zero DG (no :nop)
        :numerical_flux      => hll_flux(),   # upwind_flux() | rusanov_flux() | hll_flux() | hllc_flux() | roe_flux()
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => false, # inviscid: FV has no viscous term yet
        :ivisc_equations      => [1],
        :μ                    => [0.0], #kinematic viscosity constant for θ equation
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./problems/AdvDiff/advection2d_fv/quad_40x80_periodic.msh",
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
