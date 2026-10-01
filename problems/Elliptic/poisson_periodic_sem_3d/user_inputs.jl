function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # Triply periodic Poisson problem  -∇²u = f  on [0,2π]³  (see user_source.jl)
        #
        # ONE deck, three solvers of the SAME problem:
        #   :lfft => false   SEM: Jexpresso's native spectral-element
        #                    discretisation on the periodic mesh below, solved
        #                    directly (standard_linsolve!). The periodic SEM
        #                    Laplacian is singular (constants), so the solve
        #                    fixes the zero-mean solution — see
        #                    solve_periodic_sem_system in src/kernel/solvers/Axb.jl.
        #   :lfft => true    FFT: the Fourier spectral solver on a uniform
        #                    :fft_N³ grid of the same box.
        #   :lpseudospectral => true
        #                    pseudo-spectral: Fourier collocation with Kopriva's
        #                    derivative matrix (dense, physical space) on the
        #                    same uniform grid (:ps_N, default :fft_N; even).
        # (:lfft wins if both spectral flags are set.)
        #
        # SEM solves of the periodic system (with :lfft/:lpseudospectral false):
        #   default                        sparse direct (factorise + solve)
        #   :linsolve_amg => true          AMG-preconditioned CG, full system
        #   :lstatic_condensation => true  element-learning static condensation
        #                                  (elementLearning_Axb!, T^ie from the
        #                                  SEM matrix, no network);
        #                                  skeleton system solved by
        #                                  :EL_skeleton_solver => "direct" | "amg"
        # AMG options: :amg_method => "sa" | "rs", :amg_rtol (CG stopping tol).
        #---------------------------------------------------------------------------
        :llinsolve            => true,
        :lfft                 => false,
        :lpseudospectral      => false,
        :linsolve_amg         => false,
        :lstatic_condensation => false,
        :EL_skeleton_solver   => "direct",
        :amg_method           => "sa",
        :amg_rtol             => 1e-12,
        #--- FFT / pseudo-spectral (used only with :lfft / :lpseudospectral) -------
        # 3D: an :fft_N^3 grid of the box [0,2π]^3 (the deck's :nsd = 3)
        :fft_N                => 64,
        :fft_Lx               => 2π, :fft_Ly => 2π, :fft_Lz => 2π,
        :fft_x0               => 0.0, :fft_y0 => 0.0, :fft_z0 => 0.0,
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
        # Mesh: Jexpresso's built-in Cartesian grid, 8×8×8 hexahedra on
        # [0,2π]³, periodic in x, y and z. :nelz makes the deck 3D (:nsd = 3).
        #---------------------------------------------------------------------------
        :lcartesian_grid      => true,
        :nelx                 => 8, :xmin => 0.0, :xmax => 2π,
        :nely                 => 8, :ymin => 0.0, :ymax => 2π,
        :nelz                 => 8, :zmin => 0.0, :zmax => 2π,
        :cartesian_bdy        => Dict(:xmin => "periodicx", :xmax => "periodicx",
                                      :ymin => "periodicy", :ymax => "periodicy",
                                      :zmin => "periodicz", :zmax => "periodicz"),
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat            => "none",
        :output_dir           => "./output/",
        :loverwrite_output    => true,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
