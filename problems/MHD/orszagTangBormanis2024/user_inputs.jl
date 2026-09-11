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
        #:Δt                   => 1.5e-4,
        :Δt                   => 0.7e-5,
        :tinit                => 0.0,
        :tend                 => 1.0,   # the paper's t ∈ [0, 1] interval
        :diagnostics_at_times => (0.0:0.5:1.0),
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
        # (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ), all at full strength. The
        # coefficients by equation are those of Dao & Nazarov (2022, JSC
        # 92:77, §4.4), who apply the residual viscosity to these same
        # equations with continuous elements: one kinematic ν from the max
        # of the normalized residuals (their eq. 4.8; :dsgs_CR and :dsgs_Cmax
        # are their C_R and C_max),
        # then ν on ∇ρ (:μ[1] = 1, their eq. 4.4, the term that keeps ρ
        # positive), the dynamic ρν in the momentum stress, κ = ρν/Pr on
        # ∇T with T = p/ρ (:dsgs_nazarov_energy; the default is the
        # Fourier-law c_p ρν/Pr, 2.5× larger at γ = 5/3), η = ν on B. Pr = 1
        # as in their runs. dsgs_gamma must match γ_mhd in user_flux.jl.
        # The VTK output carries one mu_dsgs_<slots> field per distinct
        # coefficient; with Pr = 1 that is mu_dsgs_ρ_Bx_By_Bz_ψ (ν) and
        # mu_dsgs_ρu_ρv_ρE_ρw (ρ̄ν = κ), a third field mu_dsgs_ρE (ρ̄ν/Pr)
        # appearing for Pr ≠ 1.
        #---------------------------------------------------------------------------
        :lvisc            => true,
        :μ                => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        # DynSGS sensor: "legacy" = the sensor this case was validated with
        # (the assembled RHS against a fixed BDF2 of the stage state, in
        # effect a |∂ₜq| sensor); "residual" (the default) = the element-wise
        # strong residual with the stage-consistent stencil, DSGS.md §1.2.
        :visc_model       => DSGS_MHD(),
        :dsgs_sensor      => "legacy",
        :dsgs_CR          => 1.0,
        :dsgs_Cmax        => 0.5,
        :dsgs_gamma       => 5.0/3.0,
        :dsgs_Prt         => 1.0,
        :dsgs_nazarov_energy => true,
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
        :outformat           => "vtk",   # ParaView: the output variables plus the DynSGS coefficient fields
        :loverwrite_output   => false,
        :lwrite_initial      => true,
        :output_dir          => "/scratch/smarras/smarras/MHD/",
        #:output_dir          => "./output/",
        :loutput_pert        => false,
        #---------------------------------------------------------------------------
        # AMR (off)
        #---------------------------------------------------------------------------
        :linitial_refine     => true,
        :init_refine_lvl     => 2,
        :ladapt              => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
