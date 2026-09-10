#---------------------------------------------------------------------------------
# Brio-Wu MHD shock tube, 1D ideal MHD (Brio & Wu 1988), as set up in Dao &
# Nazarov, J. Sci. Comput. 92:77 (2022), Sec. 5.2: domain (0, 1), γ = 2,
# t ∈ [0, 0.2], left/right states in initialize.jl. The MHD counterpart of
# problems/CompEuler/sod1d: stabilized by DynSGS (kernel/physics/SGS.jl,
# the 1D DSGS_MHD kernel) in its conserved form, integrated with SSPRK53.
#
#   julia --project=. src/Jexpresso.jl MHD brioWu1d
#
# 150 elements at N = 4 are 600 LGL points, the "600 DOFs" of the paper's
# Fig. 2(a); :nelx => 300 gives its 1200-DOF case.
#---------------------------------------------------------------------------------
function user_inputs()
    inputs = Dict(
        :ode_solver           => SSPRK53(),
        # Fastest signal: the fast magnetosonic speed of the right state,
        # √((γp + |B|²)/ρ) = 3.75, against the smallest LGL spacing
        # 0.146·(1/150) ≈ 1e-3: Δt = 5e-5 is a Courant number of 0.19.
        :Δt                   => 5.0e-5,
        :tinit                => 0.0,
        :tend                 => 0.2,
        :diagnostics_at_times => (0:0.05:0.2),
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        :lperiodic_1d         => false,
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # DynSGS-MHD, conserved form: one residual-based kinematic ν per
        # element (max over the equations of the normalized residual, Dao &
        # Nazarov eq. 4.8; C1 = their C_R, C2 = C_max), capped at
        # C2·Δ·(|u| + c_f), applied as ∇·(ν∇q) to every conserved variable
        # (user_primitives.jl). :μ are the per-slot multipliers for
        # (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz); :dsgs_gamma must equal γ_mhd = 2
        # of user_flux.jl.
        #---------------------------------------------------------------------------
        :lvisc            => true,
        :μ                => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model       => DSGS_MHD(),
        :dsgs_C1          => 1.0,
        :dsgs_C2          => 0.5,
        :dsgs_C0          => 0.0,
        :dsgs_gamma       => 2.0,
        :dsgs_Prt         => 1.0,
        :dsgs_conserved   => true,
        :dsgs_local_norms => false,   # the whole tube is the reference scale
        :ldsgs_global_norms => true,  # (only matters under MPI)
        :energy_equation  => "energy",
        :lkep              => false,
        :entropy_variables => false,
        :lfilter           => false,
        #---------------------------------------------------------------------------
        # Mesh: uniform 1D grid built from xmin/xmax/nelx
        #---------------------------------------------------------------------------
        :lread_gmsh           => false,
        :xmin                 => 0.0,
        :xmax                 => 1.0,
        :nelx                 => 150,
        #---------------------------------------------------------------------------
        # Output, as for CompEuler/sod1d: one figure fields-it<n>.png per
        # output time with a panel per output variable (ρ, u, v, p, By;
        # Jexpresso in blue, the reference solution of reference_hll.dat
        # dashed at t = 0.2) and a last panel with the DynSGS coefficient
        # per element. :plot_matrix => false writes one PNG per panel instead.
        #---------------------------------------------------------------------------
        :outformat            => "png",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    )
    return inputs
end
