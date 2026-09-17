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
        # THE TWO KEYS BELOW ARE NOT OPTIONAL FOR THIS CASE, and without them
        # the run is slow and dies at t ≈ 200 s whatever :Δt is.
        #
        # This deck advances the TOTAL variables (:SOL_VARS_TYPE defaults to
        # TOTAL()) on top of the hydrostatic qe built in initialize.jl, with
        # the full flux and the full gravity source. The element-wise strong
        # residual that DynSGS uses by default (:dsgs_sensor => "residual",
        # the default since September 2026) therefore reads, at rest, the
        # interpolation error of the hydrostatic balance at every element
        # interface — about 30 % of ρg, steady, and nothing to do with the
        # flow. It holds the normalized ratio above 1 for the whole run, so ν
        # sits at its cap Cmax·Δ·(|u|+c) ≈ 10³-10⁴ m²/s from the first step:
        # three orders of magnitude above the ≈ 2 m²/s that SMAG() gives this
        # same case. The bubble is then destroyed by its own stabilization
        # (diffusion time over an element Δx²/ν ≈ 1 s), and the failure time
        # does not move when :Δt is reduced, because the cap does not contain
        # Δt. DSGS.md §1.2, "The reference state"; rhs.jl warns about it now.
        #
        # "legacy" is the sensor every θ case of this code was validated with
        # (problems/CompEuler/theta_dsgs sets the same pair): the assembled
        # rate against a fixed BDF2, which at hydrostatic rest reads ≈ 0.
        # :dsgs_reference => true is what the "residual" sensor needs instead
        # — the element RHS of qe, time-independent, evaluated once and
        # subtracted so the residual measures the DEPARTURE from qe. It is
        # kept here so that flipping :dsgs_sensor to "residual" is a one-line
        # experiment rather than a blow-up.
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
