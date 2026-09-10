#---------------------------------------------------------------------------------
# Brio-Wu MHD shock tube, 1D ideal MHD (Brio & Wu 1988), as set up in Dao &
# Nazarov, J. Sci. Comput. 92:77 (2022), Sec. 5.2: domain (0, 1), γ = 2,
# t ∈ [0, 0.1] (see below), left/right states in initialize.jl. The MHD counterpart of
# problems/CompEuler/sod1d: stabilized by DynSGS (kernel/physics/SGS.jl,
# the 1D DSGS_MHD kernel) in its conserved form, integrated with SSPRK53.
#
#   julia --project=. src/Jexpresso.jl MHD brioWu1d
#
# 150 elements at N = 4 are 600 LGL points, the "600 DOFs" of the paper's
# Fig. 2(a); :nelx => 300 gives its 1200-DOF case. Final time 0.1 on (0, 1)
# = Brio & Wu's 0.2 on (−1, 1), the state the paper's Fig. 2 shows.
#---------------------------------------------------------------------------------
function user_inputs()
    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        # Fastest signal: the fast magnetosonic speed of the right state,
        # √((γp + |B|²)/ρ) = 3.75, against the smallest LGL spacing
        # 0.146·(1/150) ≈ 1e-3: Δt = 5e-5 is a Courant number of 0.19.
        :Δt                   => 5.0e-5,
        :tinit                => 0.0,
        # t = 0.1 on (0, 1) is Brio & Wu's t = 0.2 on (−1, 1), the state of
        # the paper's Fig. 2 (its "t̂ = 0.2"): see README.md.
        :tend                 => 0.1,
        :diagnostics_at_times => (0:0.025:0.1),
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
        # Nodal form of the kernel (Dao & Nazarov's own): ν at every node from
        # the assembled residual, so the coefficient is a continuous field and
        # the diffusive flux has no jump at element interfaces; the element
        # form (one ν per element) left one wiggle per element in the plateau
        # behind the compound wave. :dsgs_Cl is their C_l = 0.4 (eq. 4.7).
        :dsgs_nodal       => true,
        :dsgs_Cl          => 0.4,
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
        # Output. Default: the density figure of the paper's Fig. 2
        # (user_plot.jl), density-it<n>.png, with the reference solution of
        # reference_hll.dat and its three zoom boxes at the final time.
        # :plot_user => false gives instead the format of CompEuler/sod1d: one
        # figure fields-it<n>.png per output time with a panel per output
        # variable (ρ, u, v, p, By; the reference dashed at the final time)
        # and a last panel with the DynSGS coefficient per element
        # (:plot_matrix => false writes those panels as separate files).
        #---------------------------------------------------------------------------
        :outformat            => "png",
        :plot_user            => true,
        :loverwrite_output    => true,
        :output_dir           => "./output",
    )
    return inputs
end
