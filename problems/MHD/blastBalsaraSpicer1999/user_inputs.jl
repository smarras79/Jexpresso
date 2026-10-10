#---------------------------------------------------------------------------------
# MHD blast wave (Balsara & Spicer 1999b; Balsara 2004 §7.3): (0, 1)², periodic, t = 0.01, γ = 1.4.
# 64×64 elements of degree 3 (193² nodes, the paper's 200² zones); positivity as in rotorDaoNazarov2022.
#   julia --project=. src/Jexpresso.jl MHD blastBalsaraSpicer1999
#---------------------------------------------------------------------------------
function user_inputs()

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 1.0e-5,   # CFL ≈ 0.15 at |v| + c_f ≈ 65, smallest LGL spacing ≈ 4.3e-3
        :tinit                => 0.0,
        :tend                 => 0.01,
        :diagnostics_at_times => (0.0:0.001:0.01),
        :conservation_every   => 1,        # mass and energy totals every step: conservation.dat
        :conservation_slots   => [1, 4],
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,     # GLM ψ damping (user_source.jl)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        :interpolation_nodes  => "lgl",
        :nop                  => 3,
        #---------------------------------------------------------------------------
        # DynSGS in its conserved form: one ν per element and ∇·(ν∇q) on every conserved variable, with the
        # legacy residual normalized per element; ρ or p that still undershoot are fixed by the conservative
        # element limiter. As problems/MHD/rotorDaoNazarov2022, plus fast-speed floors and no startup hold.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model           => DSGS_MHD(),
        :dsgs_sensor          => "legacy",
        :dsgs_rel             => 1.0,
        :dsgs_fast_floors     => true,     # floors on c_f = √((γp + |B|²)/ρ): ambient β = 2.5e-4, c_f/c = 75
        :dsgs_hold_steps      => 0,        # ν from the first step: the initial jump is the most violent instant
        :dsgs_norms           => "element",
        :dsgs_CR              => 1.0,
        :dsgs_Cmax            => 0.5,
        :dsgs_gamma           => 1.4,
        :dsgs_Prt             => 1.0,
        :dsgs_conserved       => true,
        :dsgs_nazarov_energy  => false,
        :lpositivity          => true,
        :positivity_method    => "conservative",   # elements scaled toward their mean: ρ, ρv, E conserved
        :positivity_rho_min   => 1.0e-6,           # 1e-6 of the ambient ρ = 1
        :positivity_p_min     => 1.0e-7,           # 1e-6 of the ambient p = 0.1
        :lrichardson          => false,
        :energy_equation      => "energy",
        :lkep                 => false,
        :entropy_variables    => false,
        #---------------------------------------------------------------------------
        # Mesh: the Orszag-Tang 32×32 periodic square refined once by p4est: 64×64.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/MHD/orszagTangBormanis2024/OT_32x32_periodic.msh",
        :linitial_refine      => true,
        :init_refine_lvl      => 1,
        :ladapt               => false,
        :lfilter              => false,
        :outformat            => "vtk",
        :loverwrite_output    => false,
        :lwrite_initial       => true,
        :output_dir           => "./output/",
        :loutput_pert         => false,
    ) #Dict

    return inputs
end
