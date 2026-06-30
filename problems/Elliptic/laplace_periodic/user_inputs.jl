function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # PERIODIC Laplace / Poisson problem  -∇²u = f  on [-π,π]²
        # solved by a GLOBAL FOURIER SPECTRAL method (FFT).
        #
        # The order of THIS method is N = :fft_N, the number of Fourier modes /
        # grid points — i.e. spectral accuracy. It is INDEPENDENT of the SEM
        # polynomial order :nop. The FFT uses its own uniform N×N grid over the
        # domain below (the natural collocation grid for the trigonometric basis
        # e^{i k·x}); the SEM/LGL nodal basis is not used by the solve. So pick
        # :fft_N for the Fourier resolution you want (power of 2); set it to your
        # grid's resolution and the FFT grid coincides with your grid points.
        #
        # (:nop only affects how Jexpresso reads the mesh during setup — kept at a
        #  normal value so a GMSH-periodic mesh reads cleanly; it does NOT change
        #  the FFT accuracy. Output is a STRUCTURED_POINTS VTK.)
        #
        # NOTE: Fourier is the right basis here because the problem is PERIODIC.
        # For NON-periodic BCs the analogous spectral choice is a Chebyshev basis
        # (Kopriva's FastChebyshevTransform) — not used in this periodic case.
        #---------------------------------------------------------------------------
        :llinsolve            => true,
        :lfft                 => true,
        :fft_use_mesh         => false,     # FFT on its own spectral grid (mesh-independent)
        :fft_N                => 64,        # Fourier resolution N (modes) — power of 2
        :fft_x0               => -π,        # domain corner  (x ∈ [-π, π])
        :fft_y0               => -π,        # domain corner  (y ∈ [-π, π])
        :fft_Lx               => 2π,        # period in x
        :fft_Ly               => 2π,        # period in y
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
        :nop                 => 4,          # mesh-read order only; does NOT affect the FFT
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        # Your periodic [-π,π]² grid. (Used to set up the run; the FFT in mode (A)
        # solves on its own uniform :fft_N grid over the domain given above.)
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

        #=== OPTIONAL: write the FFT result on the SEM mesh nodes ==================
        #   Only useful if you specifically want the output sampled at the mesh
        #   nodes rather than on the spectral grid. It does NOT make the method
        #   higher/lower order — the FFT order is still N. Because it samples the
        #   solution at the SEM nodes, those nodes must be equispaced, i.e. :nop=>1,
        #   and (since the nop=1 2D periodic-mesh reader crashes) the mesh must NOT
        #   carry periodic tags — the FFT folds the seam itself. To use it:
        #       :fft_use_mesh => true,  :nop => 1,
        #       :fft_Lx => 2π, :fft_Ly => 2π,
        #   with :gmsh_filename a NON-periodic structured square mesh.
        #=========================================================================#
    ) #Dict

    return inputs

end
