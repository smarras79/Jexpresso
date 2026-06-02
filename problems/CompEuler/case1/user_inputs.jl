function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 1D CompEuler: Marras et al. DSGS (Dynamic SGS) stabilization test case.
        #
        # Initial condition is the smooth Gaussian acoustic pulse from
        # Kopriva, Implementing Spectral Methods for PDEs, Sec. 7.4.3:
        #
        #     ρ(x,0) = 1
        #     u(x,0) = 0
        #     p(x,0) = 1 + exp(-log(2) * (x - xs)^2 / ω²),  xs = 1.5, ω = 0.125
        #
        # The DSGS viscosity (Marras et al.) is computed from the BDF2 inviscid
        # residual and added to the inviscid RHS to stabilize the Euler system.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK53(),
        :tend                 => 3.0,
        :Δt                   => 5.0e-4,
        :diagnostics_at_times => (0:0.15:3.0),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lsource              => false,
        :lperiodic_1d         => true,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => DSGS(),     # Marras et al. Dynamic SGS
        :μ                    => [0.0, 0.0, 0.0],
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh           => false,
        :xmin                 => -5.0,
        :xmax                 =>  5.0,
        :nelx                 => 50,
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "png",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    ) #Dict

    return inputs

end
