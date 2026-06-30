function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # PERIODIC Laplace / Poisson problem solved by the classical FFT solver
        # directly on a structured periodic GMSH grid over [-π, π]ᵈ (2D or 3D).
        #
        #   -∇²u = f,   u periodic.   (d = 2 ⇒ [-π,π]²,  d = 3 ⇒ [-π,π]³)
        #
        # The solver is dimension-agnostic: it reads the mesh's space dimension
        # (mesh.nsd) and runs a 2-D or 3-D FFT accordingly.
        #
        # :lfft => true        → use fft_linsolve! (problems/drivers.jl dispatch)
        # :fft_use_mesh => true → solve ON the mesh nodes (NOT a synthetic grid):
        #     fft_linsolve! folds the mesh nodes onto the uniform periodic lattice,
        #     FFT-solves, and scatters the result back to the mesh for the normal
        #     Jexpresso (VTK) output.
        #
        # REQUIREMENTS for the mesh path (a classical FFT needs an equispaced grid):
        #   • :nop => 1   so the global SEM nodes are the (uniform) element corners.
        #                 With nop>1 the in-element LGL nodes are non-uniform and
        #                 the solver will stop with an explanatory error.
        #   • A structured grid with a POWER-OF-TWO number of cells per direction
        #     (e.g. 32×32 or 32×32×32) so each axis count is a power of two.
        #   • The period per axis is taken from :fft_Lx/:fft_Ly[/:fft_Lz] (2π here).
        #
        # Point :gmsh_filename at YOUR periodic mesh. (Both a GMSH-periodic mesh and
        # a plain closed structured box work — the solver folds the +π seam back.)
        #---------------------------------------------------------------------------
        :llinsolve            => true,
        :lfft                 => true,
        :fft_use_mesh         => true,
        :fft_Lx               => 2π,        # period in x  (domain [-π, π])
        :fft_Ly               => 2π,        # period in y
        :fft_Lz               => 2π,        # period in z  (ignored for a 2D mesh)
        #--- generic setup flags --------------------------------------------------
        :ode_solver           => "BICGSTABLE",
        :ndiagnostics_outputs => 1,
        :lsource              => false,     # the FFT RHS comes from user_fft_rhs
        :lsparse              => true,
        :ldss_laplace         => false,     # the SEM Laplace matrix is not needed
        :ldss_differentiation => false,
        :rconst               => [0.0],
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 1,          # REQUIRED: linear elements ⇒ uniform grid
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        # User-supplied periodic mesh spanning [-π,π] per direction.
        # NOTE: the repo's other meshes live under "gmsh_grids" (no 'e'); change
        # "gmesh_grids" → "gmsh_grids" if the run reports the file is not found.
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_2d_2pi.msh",
        #---------------------------------------------------------------------------
        # grid modification parameters (identity: the mesh already spans [-π,π]²)
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
