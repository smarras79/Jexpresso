function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 2D CompEuler: Mach-3 supersonic flow over a forward-facing step.
        #
        # Wind tunnel 3 m x 1 m with a 0.2 m step 0.6 m from the inflow; the
        # tunnel is filled with, and fed from the left by, a uniform Mach-3
        # stream of air at p = 101325 Pa, T = 293 K (|u| ≈ 1029 m/s), per the
        # Loci/STREAM "2D Supersonic Forward Step" tutorial. This is the
        # dimensional form of Emery (1968) / Woodward & Colella (1984), also
        # Section 5.1 of Nazarov & Hoffman (IJNMF 71:339-357, 2013).
        #
        # A bow shock stands off the step, reflects from the roof, and the
        # reflections merge into a Mach stem. Everything interesting in this
        # case is a discontinuity, which drives the two decisions below.
        #
        # (1) :energy_equation => "energy".  Slot 4 is ρE, not ρθ. ρθ is an
        #     entropy variable: it is conserved across a contact but NOT
        #     across a shock, so the Euler-θ system carries the wrong shock
        #     speed no matter how it is stabilized.
        #
        # (2) :visc_model => DSGS().  Residual-based artificial viscosity as
        #     shock capturing, Nazarov & Hoffman eq. (3.4)-(3.7): the
        #     viscosity is proportional to the local residual of the
        #     conservation laws, so it appears at the shocks and stays near
        #     zero in the smooth 90% of the field. A constant AV() coefficient
        #     large enough to hold the Mach-3 shocks would smear the whole
        #     domain. The total-energy branch of compute_dsgs_viscosity! (2D)
        #     is selected by :energy_equation above.
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => 8.0e-3,          # ≈ 2.7 tunnel flow-throughs
        :lrestart             => false,
        :restart_time         => 0.0,
        # CFL. The grid is h = 0.025 m with :nop => 4, so the tightest LGL
        # node spacing is ≈ 0.0043 m; against the free-stream wave speed
        # |u| + c ≈ 1372 m/s that puts the ADVECTIVE limit near 1.4e-6 s.
        # With :init_refine_lvl => 2 below, every one of those numbers
        # shrinks by 4 (element size) and the viscous one by 16 (Δx²).
        #
        # The binding constraint is NOT advective, it is VISCOUS. DynSGS
        # saturates its own μ_max bound at the step corner (measured: μ =
        # 14.4 Pa·s against a local cap of 14.95), and μΔt/(ρΔx²) at
        # Δt = 5e-7 is already ≈ 0.22 with :μ => 1.0. Scaling :μ up without
        # scaling Δt down therefore blows the viscous limit — see the note
        # on :μ below. The pair (Δt, :μ) has to move together.
        :Δt                   => 0.5e-7,
        :diagnostics_at_times => (0:5.0e-5:8.0e-3),
        # Wall-clock note, not a setting: at Δt = 1.25e-7 the diagnostics
        # above are 3200 steps apart, so the CFL/VTK lines are ~35-40 min
        # apart and the whole run is 64000 steps, order 12 h on one core.
        # A long silence after "Integrator warm-up with real callbacks" is
        # the run working, not a hang. If you ever want to watch it step,
        # `JEXPRESSO_STEP_HEARTBEAT=1` turns on a per-step trace without
        # editing this deck.
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,               # polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #---------------------------------------------------------------------------
        :energy_equation      => "energy",        # slot 4 is ρE — see note (1)
        :lvisc                => true,
        :visc_model           => DSGS(),          # residual-based shock capturing     
        :μ                    => [1.0, 4.0, 4.0, 4.0],
        # Artificial Prandtl number P of eq. (3.7): κ = P/(γ-1)·μ. Nazarov &
        # Hoffman use P ≈ 0.1.
        :Pr                   => 0.1,
        #:ldsgs_global_norms   => true,
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/ffs_step_transfinite.msh",
        :gmsh_filename        => "./problems/CompEuler/shock_circle/plate_hole_circle_unit.msh",
        #---------------------------------------------------------------------------
        # Plotting
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        #:output_dir           => "/scratch/smarras/smarras/output/shock_circle/",
        :output_dir           => "./output",
        :loutput_pert         => false,           # plot the total state
        # Numerical schlieren from ρ, computed at output times only
        # (kernel/physics/schlieren.jl). Adds two point-data fields to the
        # VTU on top of :outvars —
        #   schlieren_grad_rho  |∇ρ| [kg/m⁴], quantitative
        #   schlieren           exp(-k|∇ρ|/max|∇ρ|), the picture
        # For the familiar dark-shock look, colour "schlieren" with a
        # REVERSED greyscale in ParaView. This is the field to look at on
        # this case: the bow shock, its reflection off the roof, the Mach
        # stem and the slip line downstream of the triple point are all
        # density features, and the exponential map keeps the weak ones
        # visible next to the strong bow shock.
        :lschlieren           => true,
        :schlieren_k          => 20.0,            # contrast; Hadjadj uses 10-100
        #---------------------------------------------------------------------------
        # AMR off: the mesh already resolves the shocks at h/nop = 1/80, and
        # DynSGS is what handles what is left under-resolved.
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :init_refine_lvl      => 1,
        :ladapt               => false,
    ) #Dict

    return inputs

end
