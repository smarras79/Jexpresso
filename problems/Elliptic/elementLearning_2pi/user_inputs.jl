function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # ELEMENT-LEARNING member of the 3-method comparison trio. All three solve
        # the SAME manufactured Laplace problem on the SAME domain [-π,π]²
        # (only the discretization / grid differs):
        #
        #     -∇²u = f ,   u_ex(x,y) = sin(x)·cos(y) ,   f = 2·sin(x)·cos(y)
        #
        #   • Elliptic/laplace_periodic   — Fourier/FFT  (periodic BCs, uniform grid)
        #   • Elliptic/laplace_chebyshev  — Chebyshev    (Dirichlet BCs, CGL grid)
        #   • Elliptic/elementLearning_2pi (this case)   — element learning / SEM
        #                                                  (Dirichlet BCs, LGL grid)
        #
        # u_ex = sin(x)cos(y) is the manufactured solution already used by the
        # base `elementLearning` case (MMS_KX = MMS_KY = 1 in user_source.jl); here
        # the only change is the DOMAIN, mapped to [-π,π]² via the grid scaling
        # below. The mesh is the [-1,1]² `square_dirichletT` grid, so
        #   x_phys = (x_mesh + xdisp)·(xscale·0.5) = x_mesh·π ∈ [-π,π]   (xdisp=0)
        # and likewise in y. (The other elementLearning_* cases are left as-is.)
        #
        # NOTE on the grid (the three cases share the DOMAIN, not the grid):
        # element learning solves a DIRICHLET problem, so it needs a Dirichlet-
        # tagged mesh — it cannot use the PERIODIC hexa_TFI_2d_2pi.msh that the FFT
        # case uses. A SMALL square_dirichletT mesh is used here: a 15×15 grid at
        # nop=6 is ~91×91 ≈ 8k nodes, which builds in moments. (The base case's
        # 100×100 mesh is ~360k nodes at nop=6 and makes the matrix wrapper crawl —
        # use it only for production EL, not this comparison.)
        #
        # As shipped this runs the DIRECT SEM solve (Ax=b). To run the element-
        # learning inference instead — and get its "# SOLVER TIMING" line for the
        # comparison — uncomment :lelementLearning and point :NNfile at a trained
        # model (â ↦ Tie), exactly as for the base elementLearning case.
        #---------------------------------------------------------------------------
        :tend                 => 1.0,
        :ode_solver           => "BICGSTABLE",
        :ndiagnostics_outputs => 1,
        :lsource              => true,
        :llinsolve            => true,
        :ldss_laplace         => true,
        :lsparse              => true,
        #:lelementLearning     => true,   # ← uncomment to time the EL inference
        #:lEL_Sample           => true,
        :Nsamp                => 50000,
        :rconst               => [0.0],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        # Small DIRICHLET-tagged square (base extent [-1,1]²) — fast to build.
        :gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_15x15.msh",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 6,      # Polynomial order
        #---------------------------------------------------------------------------
        # grid modification parameters → map [-1,1]² mesh to [-π,π]²
        #---------------------------------------------------------------------------
        :xscale              => 2π,     # x_phys = x_mesh·(2π·0.5) = x_mesh·π
        :yscale              => 2π,
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
