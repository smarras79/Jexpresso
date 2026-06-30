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
        # SAME NUMBER OF DOF as the FFT and Chebyshev cases: 64 points/dir = 4096
        # nodes, so the three solve timings are comparable. The global SEM node
        # count per direction is  nelem·nop + 1 ; with the 3×3 square_dirichletT
        # mesh and nop = 21 that is  3·21 + 1 = 64.
        #
        # Why nop = 21 (and not the usual 6): the FFT needs a power-of-two grid
        # (64), and 6·nelem+1 is ALWAYS odd, so nop=6 can never match a power of
        # two. With the meshes that ship with the base case (3×3, 15×15, …), 3×3
        # at nop=21 is the clean way to hit exactly 64 nodes/dir. (If you have a
        # 9×9 Dirichlet mesh, nop=7 → 64 is gentler; a 21×21 mesh gives nop=3.)
        #
        # NOTE: element learning solves a DIRICHLET problem, so it needs a
        # Dirichlet-tagged mesh — it cannot use the PERIODIC hexa_TFI_2d_2pi.msh
        # the FFT case uses (the three share the DOMAIN, not the grid).
        #
        # As shipped this runs the DIRECT SEM solve (Ax=b). To run the element-
        # learning INFERENCE instead — and get its "# SOLVER TIMING" line for the
        # comparison — uncomment :lelementLearning and point :NNfile at a trained
        # model (â ↦ Tie). That model must be trained at the SAME nop used here;
        # the base model is nop=6, which (see above) cannot share the FFT's DOF —
        # so matched-DOF timing applies to the direct-SEM solve, while EL-inference
        # timing should be compared at the model's own nop.
        #---------------------------------------------------------------------------
        :tend                 => 1.0,
        :ode_solver           => "GMRES", #"BICGSTABLE", #ORK256()
        :ndiagnostics_outputs => 1,
        :lsource              => true, 
        :llinsolve            => true,
        :ldss_laplace         => true,
        :lsparse              => true,
        :lelementLearning     => true,
        :lEL_Sample           => true,        
        :NNfile               => "JX_NN_model_9x9_nop7.onnx",
        #:NNfile               => "JX_RFRC_model.onnx",
        #:NNfile               => "JX_RFRC_final.jld2",
        :Nsamp                => 50000,
        :rconst               => [0.0],        
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        # 3×3 DIRICHLET-tagged square (base extent [-1,1]²). With nop=21 this gives
        # 3·21+1 = 64 nodes/dir = 4096 DOF, matching the FFT and Chebyshev cases.
        :gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_1x1.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_3x3.msh",
 #       :gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_9x9.msh",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        #:nop                 => 21,    # 3 elements × nop 21 + 1 = 64 nodes/dir (= 4096 DOF)
        :nop                 => 7,     # 9 elements × nop 7 + 1 = 64 nodes/dir (= 4096 DOF)
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
