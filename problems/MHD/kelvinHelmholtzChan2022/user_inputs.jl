function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # 2D magnetized Kelvin-Helmholtz instability, ideal GLM-MHD.
        # Chan et al., Frontiers in Physics 10:898028 (2022), Section 3.2.1.
        #
        # Run with:
        #   julia --project=. src/Jexpresso.jl MHD kelvinHelmholtzChan2022
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 2.0e-3,
        :tinit                => 0.0,
        :tend                 => 15.0,  # T_final of the paper
        :diagnostics_at_times => (0.0:0.5:15.0),
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,   # GLM ψ-damping source (Dedner mixed cleaning; see user_source.jl)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 7,       # paper: N = 7 on a 32×32 mesh (Fig. 4)
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #
        # No entropy-stable/KEP fluxes for now: plain weak form stabilized with
        # Smagorinsky SGS viscosity. :μ are per-equation multipliers on the SGS
        # coefficient for (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ): no mass diffusion,
        # everything else on. The B/ψ entries act as a turbulent resistivity.
        #---------------------------------------------------------------------------
        :lvisc            => true,
        :μ                => [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model       => SMAG(),
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
        # Same doubly periodic [-1,1]² grid used by the companion Euler case
        # problems/CompEuler/kelvinHelmholtzChan2022 (from JexpressoMeshes).
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_32x32_unitsquare.msh",
        #---------------------------------------------------------------------------
        # Filter parameters (off: Smagorinsky provides the dissipation)
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
