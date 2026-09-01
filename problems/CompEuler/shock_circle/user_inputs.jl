function user_inputs()

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => 2.0e-2,          # ≈ 2.7 tunnel flow-throughs
        :lrestart             => false,
        :restart_time         => 0.0,
        :Δt                   => 0.5e-7,
        :diagnostics_at_times => (0:5.0e-5:2.0e-2),
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,               # polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #---------------------------------------------------------------------------
        :energy_equation      => "energy",        # slot 4 is ρE — see note (1)
        :lvisc                => true,
        :visc_model           => DSGS(),          # residual-based shock capturing     
        :μ                    => [1.0, 4.0, 4.0, 4.0],
        # Artificial Prandtl number P of eq. (3.7): κ = P/(γ-1)·μ. Nazarov &
        # Hoffman use P ≈ 0.1.
        :Pr                   => 0.1,
        #:ldsgs_global_norms   => true,
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        #:gmsh_filename        => "./meshes/gmsh_grids/ffs_step_transfinite.msh",
        #:gmsh_filename        => "./problems/CompEuler/shock_circle/plate_hole_circle_unit.msh",
        :gmsh_filename        => "./problems/CompEuler/shock_circle/cylinder_symmetric.msh",
        :exact_geometry       => Dict("circle_boundary" => (:circle, 1.0, 0.0, 0.2)),
        #---------------------------------------------------------------------------
        # Plotting
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        #:output_dir           => "/scratch/smarras/smarras/output/shock_circle/",
        :output_dir           => "./output",
        :loutput_pert         => false,
        :lschlieren           => true,
        :schlieren_k          => 20.0,            # contrast; Hadjadj uses 10-100
        #---------------------------------------------------------------------------
        # AMR off: the mesh already resolves the shocks at h/nop = 1/80, and
        # DynSGS is what handles what is left under-resolved.
        #---------------------------------------------------------------------------
        :linitial_refine      => true,
        :init_refine_lvl      => 1,
        :ladapt               => false,
    ) #Dict

    return inputs

end
