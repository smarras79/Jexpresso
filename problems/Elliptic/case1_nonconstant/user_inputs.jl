function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # NON-CONSTANT DIFFUSIVITY test (derived from Elliptic/case1).
        # Standard direct solver for  -∇·(a(x,y)∇u) = f , verified against the
        # manufactured solution in user_source.jl.
        #---------------------------------------------------------------------------
        :tend                 => 1000.0,
        :ode_solver           => "BICGSTABLE",
        :ndiagnostics_outputs => 1,
        :lsource              => true,
        :llinsolve            => true,
        :ldss_laplace         => true,
        :ldss_differentiation => false,
        #:lsparse              => true,
        :rconst               => [0.0],
        #---------------------------------------------------------------------------
        # Non-constant diffusivity a(x,y) handed to the operator assembly.
        # Must be the SAME a used to build the manufactured source f.
        #---------------------------------------------------------------------------
        :diffusivity          => nonconstant_diffusivity,   # (x,y) -> A0 + A1 x
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT.msh",
        #---------------------------------------------------------------------------
        # grid modification parameters  →  domain [1,6] x [0, 3.14]
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

    return inputs

end
