function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # GEOMETRY-INDUCED NON-CONSTANT DIFFUSIVITY test (from elementLearning_hole).
        #
        # The non-constant diffusion is dictated by the GEOMETRY, exactly as in the
        # element-learning notes: the physical diffusivity is constant (a = 1, the
        # standard Laplacian), but the REFERENCE diffusivity
        #
        #     â(ξ) = (|K|/|[-1,1]^d|) · a · J_K⁻¹ J_K⁻ᵀ
        #
        # is non-constant because the element Jacobian J_K varies over the
        # plate-with-hole mesh (curved/unstructured element shapes — Option 1 in
        # the notes). No artificial a(x,y) gradient is imposed.
        #
        # Workflow:
        #   • Physical problem  -Δu = f  is verified now via the manufactured
        #     solution (user_source.jl, el_source_mode()=:mms): run with the direct
        #     solver (comment :lelementLearning) or via element learning.
        #   • :lEL_nonconstant builds the per-element â feature from the real mesh
        #     metrics for the EL surrogate inference (needs a model trained on the
        #     3·(k+1)² â feature — see tools/EL_training/train_rfrc.py).
        #   • The geometry-induced â is written to the VTU (ahat_eff, ahat_aniso
        #     cell fields) so the non-constant diffusion can be visualised.
        #---------------------------------------------------------------------------
        :tend                 => 1.0,
        :ode_solver           => "GMRES",
        :ndiagnostics_outputs => 1,
        :lsource              => true,
        :llinsolve            => true,
        :ldss_laplace         => true,
        :lsparse              => true,
#        :lelementLearning     => true,
        :lEL_nonconstant      => true,   # geometry-induced â feature (a = 1)
#        :lEL_Sample           => true,   # uncomment to (re)generate training data
        :NNfile               => "JX_NN_model.onnx",
        #:NNfile               => "JX_RFRC_model.onnx",
        :Nsamp                => 50000,
        :rconst               => [0.0],
        #---------------------------------------------------------------------------
        # NOTE: physical diffusivity is constant (a = 1) → no :diffusivity key.
        #       The non-constancy is geometric (J_K), not a physical a(x,y).
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "./output-nonconstant/",
        :loverwrite_output   => true,
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/plate_hole_circle_unit.msh",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 6,      # Polynomial order
        #---------------------------------------------------------------------------
        # grid modification parameters
        #---------------------------------------------------------------------------
        :xscale              => 5.0,
        :yscale              => 3.14,
        :xdisp               => 1.0,
        :ydisp               => 0.0,
    ) #Dict

    return inputs

end
