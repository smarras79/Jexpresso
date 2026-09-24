function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # Doubly periodic Poisson problem  -∇²u = f  on [0,2π]²  (see user_source.jl)
        #
        # ONE deck, three solvers of the SAME problem:
        #   :lfft => false   SEM: Jexpresso's native spectral-element
        #                    discretisation on the periodic mesh below, solved
        #                    directly (standard_linsolve!). The periodic SEM
        #                    Laplacian is singular (constants), so the solve
        #                    fixes the zero-mean solution — see
        #                    solve_periodic_sem_system in src/kernel/solvers/Axb.jl.
        #   :lfft => true    FFT: the Fourier spectral solver on a uniform
        #                    :fft_N × :fft_N grid of the same box.
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
        #--- FFT (used only with :lfft => true) -----------------------------------
        :fft_N                => 64,
        :fft_Lx               => 2π,
        :fft_Ly               => 2π,
        :fft_x0               => 0.0,
        :fft_y0               => 0.0,
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
        # Mesh: 16×16 quads on the unit square, periodic in x and y, mapped to
        # [0,2π]² (x ← (x + xdisp)·xscale/2, so xscale = 4π).
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/Elliptic/poisson_periodic_sem/square_periodic_16x16.msh",
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
