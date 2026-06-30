function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # Classical FFT (Fourier spectral) solver for the Laplace/Poisson equation
        #
        #   -∇²u = f   on a PERIODIC rectangle, solved by diagonalizing the
        #   Laplacian with the FFT (Kopriva's FFT routines, see
        #   src/kernel/infrastructure/Kopriva_functions.jl and
        #   src/kernel/solvers/fft_laplace.jl).
        #
        # :lfft => true diverts the linear-solve branch in problems/drivers.jl to
        # fft_linsolve!. That solver builds its OWN uniform periodic grid (an FFT
        # cannot run on the unstructured LGL/GMSH mesh), so the mesh below only
        # feeds the generic setup pipeline and is otherwise unused by the solve.
        #---------------------------------------------------------------------------
        :llinsolve            => true,
        :lfft                 => true,
        #--- FFT grid / geometry (periodic) --------------------------------------
        :fft_N                => 64,        # points per direction (power of 2)
        :fft_Lx               => 2π,        # domain lengths
        :fft_Ly               => 2π,
        :fft_x0               => 0.0,       # lower-left corner
        :fft_y0               => 0.0,
        #--- generic setup flags (kept consistent with the other Elliptic cases) -
        :ode_solver           => "BICGSTABLE",
        :ndiagnostics_outputs => 1,
        :lsource              => true,
        :lsparse              => true,
        :ldss_laplace         => true,
        :ldss_differentiation => false,
        :rconst               => [0.0],
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_50x50.msh",
        #---------------------------------------------------------------------------
        # grid modification parameters
        #---------------------------------------------------------------------------
        :xscale              => 5.0,
        :yscale              => 3.14,
        :xdisp               => 1.0,
        :ydisp               => 0.0,
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "./output/",
        :loverwrite_output   => true,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs

end
