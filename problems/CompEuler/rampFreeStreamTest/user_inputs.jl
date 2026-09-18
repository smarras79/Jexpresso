function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # FREE-STREAM PRESERVATION TEST on the Cao et al. (2021) ramp mesh.
        # See initialize.jl for why this case exists and user_bc.jl for the
        # two constants that ARE the experiment.
        #
        # EVERYTHING THAT COULD CONFOUND THE MEASUREMENT IS OFF.  That is
        # the entire design of this deck, and each "off" below is load
        # bearing:
        #
        #   :lvisc        => false   no molecular viscosity, no DynSGS.  An
        #                            artificial viscosity would DAMP the very
        #                            error being measured and turn a positive
        #                            result into a smaller positive result.
        #   :lpositivity  => false   MUST stay false.  The repair exists to
        #                            hide negative pressure; here negative
        #                            pressure is the signal.
        #   :lfilter      => false   same argument.
        #   :lsource      => false   no source terms to manufacture anything.
        #   :lkep         => false   the BASELINE is the production inviscid
        #                            operator.  Set it true (with
        #                            :volume_flux => ranocha()) for the second
        #                            run -- flux differencing has its own
        #                            free-stream-preservation requirement on
        #                            the metric terms, and user_flux.jl
        #                            already carries the methods.
        #
        # The mesh, the polynomial order, the state and :Δt are the
        # production ones, unchanged, so a result here transfers directly to
        # rampCaoEtAl2021_M7.
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        # 500 steps at the production Δt.  The first positivity repair in the
        # production run happened within the first 1000 RHS calls, i.e. inside
        # step 200, so 500 steps is twice the window in which the real failure
        # first shows -- and it is seconds of wall clock, not hours.
        :tend                 => 5.0e-7,
        :Δt                   => 1.0e-9,
        :diagnostics_at_times => (0:5.0e-8:5.0e-7),   # every 50 steps, 10 frames
        :lrestart             => false,
        :restart_time         => 0.0,
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        :energy_equation      => "energy",
        :lvisc                => false,
        :lfilter              => false,
        :lpositivity          => false,
        :lkep                 => false,
        #---------------------------------------------------------------------------
        # The PRODUCTION mesh, read from the production directory.  Swap for
        # ramp15.msh to ask the same question of the wall-stretched grid.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/CompEuler/rampCaoEtAl2021/ramp15_uniform.msh",
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        :lschlieren           => false,
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :ladapt               => false,
        :lxy_partition        => false,
    ) #Dict

    return inputs

end
