#---------------------------------------------------------------------------------
# MHD rotor (Balsara & Spicer 1999) as in Dao & Nazarov (2022, JSC 92:77, §5.5): P3 on 300×300 nodes, t = 0.15.
# Here :nop => 3 on the Orszag-Tang 32×32 mesh refined twice by p4est: 128×128 elements, 385×385 nodes.
#   julia --project=. src/Jexpresso.jl MHD rotorDaoNazarov2022
#---------------------------------------------------------------------------------
function user_inputs()

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 2.0e-4,   # CFL ≈ 0.24: max wave speed ≈ 2.6, smallest LGL spacing ≈ 2.2e-3
        :tinit                => 0.0,
        :tend                 => 0.15,
        :diagnostics_at_times => (0.0, 0.05, 0.1, 0.15),
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,     # GLM ψ damping (user_source.jl)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        :interpolation_nodes  => "lgl",
        :nop                  => 3,
        #---------------------------------------------------------------------------
        # DynSGS in its conserved form: one ν per element and ∇·(ν∇q) on every conserved variable, with the
        # legacy residual normalized per element (the low-pressure interior is judged on its own scales).
        # ρ or p that still undershoot are fixed by the conservative element limiter (:positivity_method).
        # :dsgs_gamma must equal γ_mhd = 1.4 (user_flux.jl).
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model           => DSGS_MHD(),
        :dsgs_sensor          => "legacy",
        :dsgs_rel             => 1.0,
        :dsgs_hold_steps      => 2,
        :dsgs_norms           => "element",
        :dsgs_CR              => 1.0,
        :dsgs_Cmax            => 0.5,
        :dsgs_gamma           => 1.4,
        :dsgs_Prt             => 1.0,
        :dsgs_conserved       => true,
        :dsgs_nazarov_energy  => false,   # true: κ = ρν/Pr on the thermal part of E (p dipped to −0.02 at 64² elements)
        :lpositivity          => true,
        :positivity_method    => "conservative",   # elements scaled toward their mean: ρ, ρv, E conserved
        :positivity_rho_min   => 1.0e-6,
        :positivity_p_min     => 1.0e-6,
        :lrichardson          => false,
        :energy_equation      => "energy",
        :lkep                 => false,
        :entropy_variables    => false,
        #---------------------------------------------------------------------------
        # Mesh: the doubly periodic unit square of the Orszag-Tang case, refined by p4est.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/MHD/orszagTangBormanis2024/OT_32x32_periodic.msh",
        :linitial_refine      => true,
        :init_refine_lvl      => 2,
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
