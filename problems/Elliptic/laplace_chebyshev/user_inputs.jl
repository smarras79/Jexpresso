function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # NON-PERIODIC Laplace / Poisson problem on a box, solved by a global
        # CHEBYSHEV SPECTRAL (collocation) method.
        #
        #   -∇²u = f   on [cheb_xmin,cheb_xmax] × [cheb_ymin,cheb_ymax],
        #    u = g     on the boundary  (Dirichlet).
        #
        # This is the NON-periodic counterpart of the Fourier/FFT solver: for a
        # non-periodic problem the natural spectral basis is Chebyshev, whose CGL
        # nodes x_j = cos(jπ/N) cluster toward the boundaries (which defeats the
        # Runge phenomenon and resolves boundary layers). The 1-D CGL derivative
        # matrix is Jexpresso's existing Kopriva `CGLDerivativeMatrix`; the 2-D
        # Laplacian is D²ₓ⊗I + I⊗D²ᵧ with Dirichlet rows, solved directly.
        #
        # The method order is N = :cheb_N (spectral accuracy) and is INDEPENDENT of
        # the SEM :nop. The solver builds its own CGL grid (the mesh is read only
        # to satisfy setup; :nop is the mesh-read order, not the spectral order).
        #
        # As shipped: the harmonic field u = exp(x)cos(y) (∇²u = 0 ⇒ a true
        # LAPLACE problem, f = 0) with the matching non-homogeneous Dirichlet BC.
        #---------------------------------------------------------------------------
        :llinsolve            => true,
        :lcheb                => true,
        :cheb_N               => 24,        # Chebyshev resolution (nodes/dir)
        :cheb_xmin            => -1.0,
        :cheb_xmax            =>  1.0,
        :cheb_ymin            => -1.0,
        :cheb_ymax            =>  1.0,
        #--- generic setup flags --------------------------------------------------
        :ode_solver           => "BICGSTABLE",
        :ndiagnostics_outputs => 1,
        :lsource              => false,     # the RHS comes from user_cheb_rhs
        :lsparse              => true,
        :ldss_laplace         => false,     # the SEM Laplace matrix is not needed
        :ldss_differentiation => false,
        :rconst               => [0.0],
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 4,          # mesh-read order only; not the spectral order
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        # Incidental: read only to satisfy the setup pipeline. Any small 2D mesh
        # works — the Chebyshev solve uses its own CGL grid over [cheb_xmin..] etc.
        :gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_50x50.msh",
        #---------------------------------------------------------------------------
        # grid modification parameters
        #---------------------------------------------------------------------------
        :xscale              => 1.0,
        :yscale              => 1.0,
        :xdisp               => 0.0,
        :ydisp               => 0.0,
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "./output/",
        :loverwrite_output   => true,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs

end
