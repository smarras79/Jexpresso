#---------------------------------------------------------------------------------
# Mach 15 MHD fast shock hitting a round cloud 10 times denser (Dai & Woodward 1998, JCP 142:331, p. 365,
# Figs. 23-25): [0, 1.5] × [0, 1], open boundaries, γ = 5/3. 96×64 elements of degree 4 carry 385 × 257
# nodes, DW's 384 × 256 cells; SC15_48x32.msh (with :Δt doubled, :sc_smooth => 1/32) is the quick look.
# The physics and the boundaries are those of problems/MHD/shockCloudDaiWoodward1998 with :sc_problem => "cloud".
#   julia --project=. src/Jexpresso.jl MHD shockCloudMach15DaiWoodward1998
#---------------------------------------------------------------------------------
function user_inputs()

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :nop                  => 4,
        :Δt                   => 2.0e-5,   # CFL ≈ 0.12 at |v| + c_f = 16.7, smallest LGL spacing 2.7e-3 (96×64)
        :tinit                => 0.0,
        :tend                 => 0.16,     # Fig. 23: the shock is at x ≈ 1.43 (README)
        :diagnostics_at_times => (0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16),
        :conservation_every   => 10,       # mass and energy totals: conservation.dat (they change through the open sides)
        :conservation_slots   => [1, 4],
        :sc_problem           => "cloud",     # Fig. 23 (the letters of Fig. 20: problems/MHD/shockCloudDaiWoodward1998)
        :sc_shock_x           => 0.6,      # initial shock position (Tóth 2000 §6.5, see initialize.jl)
        :sc_smooth            => 1.0/64,   # shock and cloud edge spread over one element (1/32 on SC15_48x32); 0 = DW's sharp data
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,     # GLM ψ damping (user_source.jl)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        :interpolation_nodes  => "lgl",
        #---------------------------------------------------------------------------
        # DynSGS as in problems/MHD/blastBalsaraSpicer1999: conserved form, legacy residual normalized per
        # element, fast-speed floors, thermal-energy residual, ν from the first step; ρ or p that still
        # undershoot are fixed by the conservative element limiter. :dsgs_gamma must equal γ_mhd = 5/3.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model           => DSGS_MHD(),
        :dsgs_sensor          => "legacy",
        :dsgs_rel             => 1.0,
        :dsgs_fast_floors     => true,
        :dsgs_thermal_residual => true,
        :dsgs_hold_steps      => 0,        # ν from the first step
        :dsgs_norms           => "element",
        :dsgs_CR              => 1.0,
        :dsgs_Cmax            => 0.5,
        :dsgs_gamma           => 5.0/3.0,
        :dsgs_Prt             => 1.0,
        :dsgs_conserved       => true,
        :dsgs_nazarov_energy  => false,
        :lpositivity          => true,
        :positivity_method    => "conservative",   # elements scaled toward their mean: ρ, ρv, E conserved
        :positivity_rho_min   => 1.0e-6,           # 1e-6 of the pre-shock ρ = 1
        :positivity_p_min     => 1.0e-6,           # 1e-6 of the pre-shock p = 1
        :lrichardson          => false,
        :energy_equation      => "energy",
        :lkep                 => false,
        :entropy_variables    => false,
        #---------------------------------------------------------------------------
        # Mesh: SC15.geo, open on all four sides (tags left/right/bottom/top).
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/MHD/shockCloudMach15DaiWoodward1998/SC15_96x64.msh",
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
