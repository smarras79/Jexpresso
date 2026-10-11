#---------------------------------------------------------------------------------
# Magnetized Mach 800 jet (Wu & Shu, SISC 40 (2018) B1302, Example 5.6) on the right half [0, 0.5] × [0, 1.5]
# with a symmetry axis at x = 0, as in the paper; γ = 1.4, B = (0, B_a, 0), B_a² = 200 (user_flux.jl).
#   julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
#---------------------------------------------------------------------------------
function user_inputs()

    lvl = 0     # p4est refinement of the mesh; Δt halves with each level (CFL ≈ 0.1 at c_h ≈ 812)

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 5.0e-7/2^lvl,   # h = 0.025 at lvl = 0 (2.0e-7 on AJ_half_50x150, h = 0.01)
        :tinit                => 0.0,
        :tend                 => 2.0e-3,
        :diagnostics_at_times => sort(unique(vcat(collect(0.0:2.0e-5:1.0e-4), collect(0.0:1.0e-4:2.0e-3)))),
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,     # GLM ψ damping (user_source.jl)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # DynSGS in its conserved form (∇·(ν∇q) on every conserved variable, one ν per element); nodes with
        # ρ < ρ_min or p < p_min are fixed by the conservative element limiter (stage limiter of the RK).
        # :dsgs_gamma must equal γ_mhd = 1.4 (user_flux.jl).
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model           => DSGS_MHD(),
        :dsgs_sensor          => "residual",   # carries the Dirichlet-node treatment needed at the nozzle
        :dsgs_CR              => 1.0,
        :dsgs_Cmax            => 0.5,
        :dsgs_Cmin            => 0.0,
        :dsgs_cutoff          => 0.0,
        :dsgs_hold_steps      => 0,
        :dsgs_norms           => "element",   # the beam's p is 5.6e-6 of E: only an element scale sees its noise
        :dsgs_rel             => 1.0,
        :dsgs_fast_floors     => true,     # floors on c_f = √((γp + |B|²)/ρ): β_a = 1e-2, c_f/c = 12
        :dsgs_gamma           => 1.4,
        :dsgs_Prt             => 1.0,
        :dsgs_conserved       => true,
        :dsgs_ref_weight      => false,
        :dsgs_nazarov_energy  => false,
        :ldsgs_nodal          => false,
        :lpositivity          => true,
        :positivity_method    => "conservative",   # elements scaled toward their mean: ρ, ρv, E conserved
        :positivity_rho_min   => 1.4e-7,           # 1e-6 of the ambient ρ = 0.14
        :positivity_p_min     => 1.0e-6,           # 1e-6 of the ambient p = 1
        :positivity_report    => true,
        :positivity_report_every => 250,           # every 50 steps
        :lrichardson          => false,
        :energy_equation      => "energy",
        :lkep                 => false,
        :entropy_variables    => false,
        #---------------------------------------------------------------------------
        # Mesh: AJ_half_20x60.msh (h = 0.025) or AJ_half_50x150.msh (h = 0.01); AJ_40x60.msh and
        # AJ_100x150.msh are the full width [-0.5, 0.5] × [0, 1.5].
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/MHD/astroJetWuShu2018/AJ_half_20x60.msh",
        :linitial_refine      => lvl > 0,
        :init_refine_lvl      => lvl,
        :ladapt               => false,
        :lfilter              => false,
        :outformat            => "vtk",
        :loverwrite_output    => false,
        :lwrite_initial       => true,
        :output_dir           => "./output/",
        :loutput_pert         => false,
        :lschlieren           => true,
        :schlieren_k          => 20.0,
    ) #Dict

    return inputs
end
