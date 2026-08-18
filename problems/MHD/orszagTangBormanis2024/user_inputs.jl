function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # 2D Orszag-Tang vortex, ideal GLM-MHD.
        #
        # A. Bormanis, C. A. Leon, A. Scheinker,
        # "Solving the Orszag-Tang vortex magnetohydrodynamics problem with
        #  physics-constrained convolutional neural networks",
        # Phys. Plasmas 31, 012101 (2024), Sec. I A (Eqs. 6-9) and Sec. III.
        #
        # Run with:
        #   julia --project=. src/Jexpresso.jl MHD orszagTangBormanis2024
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        # Δt = 5e-4 is CFL ≈ 0.24 for this grid: the initial maximum wave
        # speed (= the GLM cleaning speed) is c_h ≈ 2.603 and the smallest
        # LGL spacing at :nop => 4 on a 1/32 element is ≈ 5.4e-3. The margin
        # is deliberate — the vortex steepens into shocks by t ≈ 0.5 and the
        # local wave speeds grow. The reference simulation of the paper used
        # Δt = 8e-4 on its 128² finite-volume grid.
        :Δt                   => 5.0e-4,
        :tinit                => 0.0,
        :tend                 => 1.0,   # the paper's t ∈ [0, 1] interval
        :diagnostics_at_times => (0.0:0.05:1.0),
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,   # GLM ψ-damping source (Dedner mixed cleaning; see user_source.jl)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #
        # :nop => 4 on the 32×32 element grid below gives 32*4 = 128 unique
        # points per direction, i.e. exactly the 128 × 128 resolution of the
        # reference data set of the paper.
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 4,
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #
        # The paper solves IDEAL MHD (no viscosity, no resistivity) with a
        # finite-volume scheme whose Riemann solver supplies the numerical
        # dissipation. A collocated continuous-Galerkin SEM has none, so the
        # Orszag-Tang shocks must be regularized explicitly.
        #
        # Marras-Nazarov DynSGS: the eddy viscosity is set by the LOCAL
        # RESIDUAL of the governing equations rather than by a tuned
        # constant, so it appears at the shocks and stays near zero in the
        # smooth 90% of the domain. That is exactly the failure mode of the
        # Smagorinsky alternative here: ρ Cs² Δ² |S| cannot tell "resolved"
        # from "unresolved", so it had to be scaled up 8x globally to survive
        # the shocks, which over-damped everything else.
        #
        # :μ are per-equation multipliers on the DynSGS coefficient for
        # (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ): no mass diffusion, everything
        # else at full strength. The B entries act as a turbulent
        # resistivity. C1/C2 are Marras's residual and wave-speed-cap
        # coefficients; dsgs_gamma must match γ_mhd in user_flux.jl.
        #---------------------------------------------------------------------------
        :lvisc            => true,
        :μ                => [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model       => DSGS_MHD(),
        :dsgs_C1          => 1.0,
        :dsgs_C2          => 0.5,
        :dsgs_gamma       => 5.0/3.0,
        :dsgs_Prt         => 0.7,
        :lrichardson      => false,      # no gravity/stratification in this problem
        # Slot 4 carries the TOTAL ENERGY ρE: "energy" keeps the kernel's τ·u
        # viscous-work augmentation of the energy equation active.
        :energy_equation  => "energy",
        #---------------------------------------------------------------------------
        # No entropy-stable / kinetic-energy-preserving machinery:
        #---------------------------------------------------------------------------
        :lkep              => false,
        :entropy_variables => false,
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #
        # Doubly periodic UNIT SQUARE [0,1]², 32×32 quads. The mesh ships with
        # the case so that it runs out of the box; regenerate (or relocate it
        # to the JexpressoMeshes tree) with
        #
        #   gmsh -2 problems/MHD/orszagTangBormanis2024/OT_32x32_periodic.geo \
        #        -o meshes/gmsh_grids/OT_32x32_periodic.msh
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./problems/MHD/orszagTangBormanis2024/OT_32x32_periodic.msh",
        #---------------------------------------------------------------------------
        # Filter parameters.
        #
        # OFF on purpose. The Boyd-Vandeven "erf" filter is Jexpresso's other
        # stabilization mechanism, but it filters the CONSERVATIVE variables
        # independently, and filtering ρ, ρu and ρE separately can drive
        # ρE - ½ρ|v|² - ½|B|² negative even where the unfiltered state was
        # fine. Measured here: :mu_x = :mu_y = 0.1 (with :μ = 1) does reach
        # t = 1, but only just — min p = 5.1e-4, a 0.4% margin — and turning
        # it UP to 0.2 aborts at t ≈ 0.55 with sqrt(negative). Scaling the
        # Smagorinsky coefficient above is both safer and better targeted:
        # it adds dissipation in proportion to the local strain rate, i.e.
        # at the shocks and nowhere else.
        #---------------------------------------------------------------------------
        :lfilter             => false,
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => false,
        :lwrite_initial      => true,
        :output_dir          => "./output/",
        :loutput_pert        => false,
        #---------------------------------------------------------------------------
        # AMR (off)
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :ladapt              => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
