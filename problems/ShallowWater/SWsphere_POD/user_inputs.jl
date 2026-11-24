#---------------------------------------------------------------------------------
# SWsphere_POD — the Galewsky jet, configured so that its POD is worth looking at.
#
# THE SAME PHYSICS, GRID AND DISCRETISATION as problems/ShallowWater/SWsphere:
# the Galewsky, Scott & Polvani (2004) barotropically unstable jet on the cubed
# sphere shipped next to this file, the conservative Cartesian form of Marras,
# Kopera & Giraldo (2015), continuous-Galerkin spectral elements, SSP-RK3 with
# the Lagrange projection after every stage. The five sibling user_*.jl are
# verbatim copies of SWsphere's (a case directory owns its six files; see
# test/test_case_includes.jl for what happens when one tries to share them).
#
# WHAT IS DIFFERENT, and why each difference matters TO THE DECOMPOSITION:
#
#   :tend => 6 days    SWsphere runs 20. Six days is where the barotropic
#                      instability has rolled the jet up (it is the day
#                      Galewsky et al. plot) — and, for POD, it is where the
#                      FLUCTUATION about the temporal mean is the instability
#                      rather than numerical noise. At one day the jet has
#                      barely moved, the mean absorbs it, and the modes come out
#                      as grid-scale noise: measured on this case, sub-element
#                      content 0.94 at one day against 0.55 at six.
#
#   :lfilter => true   SWsphere ships :lfilter => false, which contradicts its
#                      own comments (line 15: "+ diffusion :lfilter => true …
#                      ← as shipped") and its README, which records that the
#                      filter-off configuration goes non-finite at 2.005 days.
#                      With the filter off the grid-scale content of the
#                      vorticity fluctuation rises to 0.90+ and the POD, quite
#                      correctly, returns modes of that. This deck turns it on.
#
#   POD on three       :vorticity is what the test is judged on. :h and
#   fields             :velocity are PRIMITIVE fields, whose fluctuation is the
#                      flow rather than the curl of the flow's noise; comparing
#                      the three "sub-element content" lines in the report is
#                      the quickest way to see the difference.
#
# WHAT IT COSTS: 2420 steps at Δt = 214.2 s, a few minutes on one core, and
# 25 snapshots × 15 002 nodes × 4 components ≈ 12 MB held for the run.
#
# WHAT IT SHOULD PRODUCE — measured, on this deck, serially:
#
#   vorticity   Σλ = 2.127e+04   E = 60.8 / 18.1 / 10.1 / 4.8 %
#               90 % of the energy in 4 modes, 99 % in 8
#               sub-element content: fluctuation 0.581, mode 1 0.594
#   h           sub-element content 0.366 / 0.376 — the primitive field is
#   velocity    markedly smoother than the derivative: 0.426 / 0.404
#   the modes   the unstable wave train in the 30–60°N band, wavenumber rising
#               with mode index; the (a₁,a₂) phase portrait a growing spiral,
#               which is the instability
#
# If a run reproduces those numbers, the POD machinery is behaving. If it does
# not, the report's own lines — the snapshot spread, the field rms, the
# sub-element content — say which half of the chain to look at first.
#
#   julia> Jexpresso.run_case("ShallowWater", "SWsphere_POD")
#
# See docs/POD.md, and the README next to this file.
#---------------------------------------------------------------------------------
function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # The shell, the grid and its checks — SWsphere's, unchanged.
        #---------------------------------------------------------------------------
        :lspherical_shell     => true,
        :lgrid_only           => false,
        :linit_only           => false,
        :lcheck_grid          => true,
        :lstop_on_bad_grid    => true,
        :lproject_to_sphere   => true,
        :sphere_radius        => 6.37122e6,
        :cubed_sphere_map     => :none,        # this grid is already equiangular
        :sphere_metrics       => :curl_invariant,
        #---------------------------------------------------------------------------
        # Initial condition: Galewsky, Scott & Polvani (2004), Tellus 56A, 429-440.
        #---------------------------------------------------------------------------
        :lgalewsky_perturbation => true,
        :galewsky_umax          => 80.0,
        :galewsky_phi0          => π/7,
        :galewsky_phi1          => π/2 - π/7,
        :galewsky_Omega         => 7.292e-5,
        :galewsky_gravity       => 9.80616,
        :galewsky_hmean         => 10000.0,
        :galewsky_hhat          => 120.0,
        :galewsky_alpha         => 1.0/3.0,
        :galewsky_beta          => 1.0/15.0,
        :galewsky_phi2          => π/4,
        :galewsky_nquad         => 512,
        :llagrange_projection   => true,
        #---------------------------------------------------------------------------
        # Space
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 5,
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/SWsphere_POD/cubed_sphere.msh",
        #---------------------------------------------------------------------------
        # Time. Six days: the day the instability is judged on, and the window
        # over which the fluctuation is the instability rather than the noise.
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :lcfl_dt              => true,
        :cfl                  => 1.0,
        :tinit                => 0.0,
        :tend                 => 6*24*3600,
        :ndiagnostics_outputs => 6,              # a VTK dump a day
        :ndiagnostics_prints  => 400,
        :case                 => "swsphere",
        :SOL_VARS_TYPE        => TOTAL(),
        :lsource              => true,
        #---------------------------------------------------------------------------
        # Stabilisation: the artificial diffusion of Eq. (8b) on the momentum,
        # AND the modal filter. See the header on why the filter is on here.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :ivisc_equations      => [2, 3, 4],
        :μ                    => 1.0e5,
        :lfilter              => true,
        :filter_alpha         => 0.05,
        :filter_order         => 8,
        :filter_kcut          => 2/3,
        #---------------------------------------------------------------------------
        # POD — docs/POD.md. `:lpod => true` alone would do; the rest is this
        # case being a POD case.
        #---------------------------------------------------------------------------
        :lpod                 => true,
        # A DERIVATIVE field and two PRIMITIVE ones, so the report's
        # "sub-element content" lines can be compared directly.
        :pod_fields           => [:vorticity, "h", :velocity],
        :pod_nsnapshots       => 24,             # 25 snapshots, one every 5.76 h
        :pod_tstart           => 0.0,
        :pod_nmodes           => 12,
        :pod_nmodes_plot      => 4,
        :pod_subtract_mean    => true,
        :pod_method           => :auto,
        :pod_nlon             => 720,
        :pod_nlat             => 360,
        :pod_write_vtk        => true,
        :pod_write_png        => true,
        :pod_write_data       => true,
        # The raw snapshots too, so the decomposition can be redone offline over
        # a shorter window or a different rank without re-running the case.
        :pod_write_snapshots  => true,
        :pod_time_scale       => 1.0/86400.0,
        :pod_time_label       => "t [days]",
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
    ) #Dict

    return inputs

end
