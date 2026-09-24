function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # Dirichlet Poisson problem  -∇²u = f  on [0,2π]², u = u_ex on ∂Ω
        # (see user_source.jl). One deck, four solves of the SAME SEM system:
        #
        #   default                          full SEM system, sparse direct solve
        #   :linsolve_amg => true            full SEM system, AMG-preconditioned CG
        #                                    on the free (non-Dirichlet) nodes
        #   :lelementLearning     => true,
        #   :lstatic_condensation => true    STATIC CONDENSATION of element learning
        #                                    (elementLearning_Axb!), T^ie computed
        #                                    from the SEM matrix (no network): an
        #                                    exact reduction to the skeleton unknowns
        #     :EL_skeleton_solver => "direct"  sparse direct on the condensed system
        #     :EL_skeleton_solver => "amg"     AMG-preconditioned CG on it
        #
        # AMG options (both AMG solves): :amg_method => "sa" (smoothed
        # aggregation, default) | "rs" (Ruge–Stüben), :amg_rtol => 1e-12.
        #---------------------------------------------------------------------------
        :llinsolve            => true,
        :linsolve_amg         => false,
        :lelementLearning     => false,
        :lstatic_condensation => false,
        :EL_skeleton_solver   => "direct",
        :amg_method           => "sa",
        :amg_rtol             => 1e-12,
        #--- SEM -------------------------------------------------------------------
        :ode_solver           => "BICGSTABLE",
        :ndiagnostics_outputs => 1,
        :lsource              => true,
        :lsparse              => true,
        :ldss_laplace         => true,
        :ldss_differentiation => false,
        :rconst               => [0.0],
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,          # polynomial order
        #---------------------------------------------------------------------------
        # Mesh: the 16×16 mesh of Elliptic/poisson_periodic_sem with its boundary
        # tags renamed (no periodicity), mapped to [0,2π]² (xscale = 4π).
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/Elliptic/poisson_dirichlet_sc/square_16x16.msh",
        :xscale               => 4π,
        :yscale               => 4π,
        :xdisp                => 0.0,
        :ydisp                => 0.0,
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :output_dir           => "./output/",
        :loverwrite_output    => true,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
