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
        :visc_model     => DSGS(),   # Marras et al. Dynamic SGS; :μ below are multipliers of its ν
        #
        # WHAT CHANGES WHEN THIS DECK RUNS DynSGS INSTEAD OF SMAG(), and why
        # the :μ line below is not the one the SMAG version of this deck used.
        #
        # Measured here, 3 ranks, Δt = 0.2, to t = 1000, max ν over the run:
        #
        #     SMAG()                                3.0 m²/s, flat
        #     DSGS(), :dsgs_sensor => "legacy"      1.0e2 - 1.8e3 m²/s
        #     DSGS(), :dsgs_sensor => "residual"    2.8e2 - 2.0e3 m²/s
        #
        # The :μ entries multiply that ν, so the [0, 1, 1, 2, 3, 1] this case
        # carried under SMAG() (6 and 9 m²/s there) becomes 2-6e3 m²/s under
        # the residual sensor. That run DIES before t = 50 with a DomainError
        # on ρθ < 0 in p = C0(ρθ)^γ -- measured. At [0, 1, 1, 1, 1, 1] it
        # completes 1000 s on either sensor; at Δt = 0.1 the 2-3 vector
        # completes too.
        #
        # "legacy" is the sensor every θ case of this code was validated with
        # (problems/CompEuler/theta_dsgs sets the same pair): the assembled
        # rate against a fixed BDF2, which at hydrostatic rest reads ~0. It is
        # the least dissipative of the three over the first half of this run
        # and it carries the 2-3 :μ vector at Δt = 0.2 where "residual" does
        # not. :dsgs_reference => true is inert under it and is kept only so
        # that flipping to "residual" is a one-line experiment; on THIS mesh
        # it does not lower ν (measured: 4.0e2 - 9.6e2 m²/s with it against
        # 2.8e2 - 2.0e3 without), the rising bubble of theta_dsgs on 1 km
        # elements is where it earns its keep.
        #
        :dsgs_sensor     => "legacy",
        :dsgs_reference  => true,
        :Pr              => 0.1,     # artificial Prandtl number of the θ slot (the default)
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
