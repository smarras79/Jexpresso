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
        # CHEBYSHEV member of the 3-method comparison trio: all three solve the SAME
        # manufactured problem on the SAME domain [-π,π]² (only the grid differs):
        #     -∇²u = 2 sin(x)cos(y) ,   u_ex = sin(x)cos(y) ,  Dirichlet g = u_ex.
        #     • laplace_periodic     (Fourier/FFT, periodic)
        #     • laplace_chebyshev    (this case, Chebyshev/Dirichlet)
        #     • elementLearning_2pi  (element learning / SEM, Dirichlet)
        #---------------------------------------------------------------------------
        :llinsolve            => true,
        :lcheb                => true,
        # SAME DOF as the FFT and EL cases: 64 points/dir = 4096 nodes.
        # Chebyshev grid is (cheb_N+1) points/dir ⇒ cheb_N = 63 → 64.
        :cheb_N               => 63,
        :cheb_xmin            => -π,         # SAME domain [-π,π]² as the FFT and EL cases
        :cheb_xmax            =>  π,
        :cheb_ymin            => -π,
        :cheb_ymax            =>  π,
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
