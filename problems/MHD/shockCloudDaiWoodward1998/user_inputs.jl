#---------------------------------------------------------------------------------
# Mach 12 MHD fast shock hitting the dense letters "MHD" (Dai & Woodward 1998, JCP 142:331, Fig. 20):
# [0, 2] × [0, 1.5], open boundaries, γ = 5/3, DW's output times. 128×96 elements of degree 4 carry
# 513 × 385 nodes, DW's 512 × 384 cells; SC_64x48.msh (with :Δt doubled) is the quick look.
#   julia --project=. src/Jexpresso.jl MHD shockCloudDaiWoodward1998
#---------------------------------------------------------------------------------
function user_inputs()

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :nop                  => 4,
        :Δt                   => 1.6e-5,   # CFL ≈ 0.15 at |v| + c_f = 24.4, smallest LGL spacing 2.7e-3 (128×96)
        :tinit                => 0.0,
        :tend                 => 0.10411,
        :diagnostics_at_times => (0.0, 0.02507, 0.05004, 0.08003, 0.10411),   # the four panels of Fig. 20
        :conservation_every   => 10,       # mass and energy totals: conservation.dat (they change through the open sides)
        :conservation_slots   => [1, 4],
        :sc_shock_x           => 0.2,      # initial shock position (read off Fig. 20, see initialize.jl)
        :sc_smooth            => 1.0/64,   # shock and letters' edges spread over one element (1/32 on SC_64x48); 0 = DW's sharp data
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
        :dsgs_hold_steps      => 0,        # ν from the first step: the initial shock is sharp
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
        # Mesh: SC.geo, open on all four sides (tags left/right/bottom/top).
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/MHD/shockCloudDaiWoodward1998/SC_128x96.msh",
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
