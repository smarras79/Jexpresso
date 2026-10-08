#---------------------------------------------------------------------------------
# PODbenchmark — THE POD REFERENCE BENCHMARK.
#
# Linear advection of a multi-harmonic wave on a periodic line,
#
#     ∂u/∂t + c ∂u/∂x = 0 ,   u(x,0) = Σ_{j=1}^3 A_j cos(2πj x/L + ϕ_j) ,
#
# run for exactly ONE revolution. This is the standard test problem of the
# transport-dominated model-reduction literature, and it is here because its POD
# can be written down in closed form: the code's answer is checked against
# arithmetic rather than against another run. The derivation, the expected
# numbers and how to read the output are in the README next to this file; the
# same closed form is asserted, to 1e-10, by test/test_pod_benchmark.jl.
#
# THE POD BLOCK AT THE BOTTOM IS THREE LINES, and two of them are about the
# benchmark rather than about POD. `:lpod => true` alone gives a decomposition
# of every solution variable over the whole run, sampled at the output cadence,
# written as modes, spectrum, coefficients, CSV and a .jld2 basis — for this
# case or any other.
#---------------------------------------------------------------------------------
function user_inputs()

    #-----------------------------------------------------------------------------
    # The benchmark's own geometry and timing. Kept as locals because the POD
    # window below is DERIVED from them — see the note on :pod_tend.
    #-----------------------------------------------------------------------------
    L     = 2.0          # domain length
    c     = 1.0          # advection speed — must match user_flux.jl
    T     = L/c          # one full revolution: the exact solution returns to u(x,0)
    nsnap = 40           # POD sampling INTERVALS (so 41 snapshots)

    inputs = Dict(
        #---------------------------------------------------------------------------
        # Time integration. One revolution, 2000 steps of an SSP-RK3: the
        # discretization error over the run is far below the agreement the
        # README quotes, so what the POD sees is the exact solution sampled.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK33(),
        :tinit                => 0.0,
        :tend                 => T,
        :Δt                   => 1.0e-3,
        :diagnostics_at_times => (0:T/8:T),
        #---------------------------------------------------------------------------
        # Space. 25 elements at nop = 4 is 101 nodes, i.e. ~33 per wavelength of
        # the shortest harmonic — spectrally converged for this problem.
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lexact_integration   => false,
        :lsource              => false,
        :lperiodic_1d         => true,
        :lvisc                => false,
        :lread_gmsh           => false,     # 1-D: the grid is built by Jexpresso
        :xmin                 => 0.0,
        :xmax                 => L,
        :nelx                 => 25,
        #---------------------------------------------------------------------------
        # POD — see docs/POD.md. In a case that is not a benchmark, the whole
        # block is the first line.
        #---------------------------------------------------------------------------
        :lpod                 => true,
        :pod_nsnapshots       => nsnap,
        #
        # WHY :pod_tend IS ONE INTERVAL SHORT OF :tend, which is the one subtlety
        # of this benchmark and the reason it is written out here.
        #
        # A travelling wave has DEGENERATE pairs of modes: λ₁ = λ₂ exactly, and
        # the pair is one structure in quadrature rather than two. That equality
        # is a statement about averaging over a WHOLE period. Sampling [0,T] at
        # both ends repeats the zero phase, which counts one phase twice and
        # splits every pair by (K/2+1)/(K/2) — 5 % at K = 41 — for a reason that
        # has nothing to do with the decomposition.
        #
        # Ending the POD window one interval early makes the 41 snapshots tile
        # exactly one period with no repeat, and the pairs come out degenerate to
        # round-off. The RUN still goes to T; only the sampling window stops
        # short. (The sampling times are added to the integrator's tstops, so the
        # snapshots are taken AT them and not at the first step after.)
        #
        :pod_tend             => T*nsnap/(nsnap + 1),
        :pod_nmodes_plot      => 6,
        # Keep the raw snapshots too, so the decomposition can be redone offline
        # over a different window without re-running the case — and so that this
        # benchmark's data can be handed to another POD implementation.
        :pod_write_snapshots  => true,
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "png",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    ) #Dict

    return inputs

end
