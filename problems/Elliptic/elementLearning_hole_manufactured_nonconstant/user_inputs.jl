function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # MANUFACTURED-SOLUTION test that ALSO exercises the GEOMETRY-INDUCED â
        # non-constant path, on the plate-with-hole geometry.
        #
        # Physical problem (a = 1, standard Laplacian): user_source.jl is hardwired
        # to el_source_mode() = :mms, so this solves  -Δu = f  with f = -Δu_ex,
        # g = u_ex on ∂Ω, u_ex = sin(x)cos(y) stored in qe. The post-solve L2 error
        # ("# MMS verification: ... ‖e‖_L2 = ...") checks the result against u_ex.
        #
        # Element-learning twist (:lEL_nonconstant => true): the surrogate consumes
        # the per-element 2×2 SPD reference diffusivity â = (|K|/|ref|) J_K⁻¹ J_K⁻ᵀ
        # built from the REAL hole-mesh metrics — non-constant purely because the
        # element shapes (Jacobians) vary. So the manufactured solution verifies
        # the surrogate ON the geometry-induced â feature it is meant to learn.
        # The â field is written to the VTU (see :ahat_output) for inspection.
        #
        #   • As shipped (:lelementLearning + :lEL_nonconstant): verifies the EL
        #     surrogate; needs a model trained on the 3·(k+1)² â feature
        #     (sample with :lEL_Sample, train via tools/EL_training/train_rfrc.py).
        #   • Comment :lelementLearning to verify the physical MMS with the direct
        #     Ax=b solver (runnable now); the â VTU fields are still written.
        #---------------------------------------------------------------------------
        :tend                 => 1.0,
        :ode_solver           => "GMRES", #"BICGSTABLE", #ORK256()
        :ndiagnostics_outputs => 1,
        :lsource              => true,
        :llinsolve            => true,
        :ldss_laplace         => true,
        :lsparse              => true,
        :lelementLearning     => true,
        :lEL_nonconstant      => true,   # geometry-induced â feature (a = 1)
        :lEL_xidependent      => true,   # (sampling) within-element-varying â for bilinear quads
        :EL_sample_shape      => :quad,  # (sampling) :affine | :quad (bilinear, recommended) | :warp
        # Physical-amplitude sampling (match a non-trivial el_diffusivity at inference;
        # defaults below = a≡1, geometry-only). See EL_nonconstant_diffusivity.jl header.
        #:EL_amin              => 1.0,    # (sampling) min of a ~ U(EL_amin,EL_amax) per element
        #:EL_amax              => 1.0,    # (sampling) max of a ~ U(EL_amin,EL_amax) per element
        #:EL_avar              => 0.0,    # (sampling) >0 ⇒ smooth within-element variation of a
        :ahat_output          => :cell, # VTU â format: :cell | :nodal | :tensor
        :lEL_Sample           => true,   # uncomment to (re)generate training data
        :NNfile               => "JX_NN_model.onnx",
        #:NNfile               => "JX_RFRC_model.onnx",
        #:NNfile               => "JX_RFRC_final.jld2",
        :Nsamp                => 100000,
        :rconst               => [0.0],
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        #:output_dir          => "/project/smarras/smarras/Jexpresso/Jexpresso/EL/",
        #:output_dir          => "./output-RNN/",
        :output_dir          => "./output-manufactured-nonconstant/",
        #:output_dir          => "./output-RFRC/",
        #:output_dir          => "./output-RFRC-JLD2/",
        #:output_dir          => "./output-Axb/",
        :loverwrite_output   => true,
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        :gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_1x1.msh",
#        :gmsh_filename       => "./meshes/gmsh_grids/plate_hole_circle_unit.msh",

        
        #:gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_3x3.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_15x15.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_50x50.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/square_dirichletT_100x100.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/jexpresso_domain_unique_bcs.msh",
        #:gmsh_filename        => "./meshes/gmsh_grids/jex-el_domain_unique_bcs.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/plate_word_unit.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/square_UNSTR_unit_square_10x10el.msh",
        #---------------------------------------------------------------------------
        # static adaptivity
        #---------------------------------------------------------------------------
        #:lpreadapt       => true,
        #:amr_max_level   => 1,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 6,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:lvisc                => true, #false by default NOTICE: works only for Inexact
        #:ivisc_equations      => (1, 2, 3, 4),
        #:μ                   => (0.0, 75.0, 75.0, 75.0), #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # grid modification parameters
        #---------------------------------------------------------------------------
        :xscale              => 1.0,
        :yscale              => 1.0,
        :xdisp               => 0.0,
        :ydisp               => 0.0,
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
