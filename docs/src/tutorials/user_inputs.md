# `user_inputs.jl`
Every user-defined problem has its own namelist/input file called `user_inputs.jl`.
A sample user_inputs.jl file with all possible entries (notice that not all entries are necessary).
Use `#` to comment a line. When a line is commented or not explicitaly give, default values will be used. 
The entries (with default values) are defined in `src/io/mod_inputs.jl`.


```
function user_inputs()
    
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver              => SSPRK54(), #From the suite of solvers of DifferentialEquations.jl
        :Δt                      => 0.4,
        :tinit                   => 0.0,
        :tend                    => 1000.0,
        #:tinit                  => 100.0,
        #:tend                   => 1000.0,
        #:lrestart               => true,
        #:restart_input_file_path => "./output/CompEuler/theta/",
        :diagnostics_at_times    => (0:100:1000),
        :lsource                 => true, 
        #:backend                => MetalBackend(),
        #:SOL_VARS_TYPE          => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes     =>"lgl",
        :nop                     => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                   => true,                   #false by default
        :ivisc_equations         => [1, 2, 3, 4],           # use [], not ()!
        :μ                       => [0.0, 20.0, 20.0, 60.0], # use [], not ()!
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter                => true,
        #:mu_x                   => 0.01,
        #:mu_y                   => 0.01,
        #:filter_type            => "erf",
	#---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh              => true, #If false, a 1D problem will be enforced
        :gmsh_filename           => "./meshes/gmsh_grids/hexa_TFI_RTB20x20.msh",
	#---------------------------------------------------------------------------
	# init_refinement
	#---------------------------------------------------------------------------
	# When true this takes the input mesh (which can be as coarse as you'd like)
	# and statically refine it as many times as the value of :init_refine_lvl
	# E.g. Given a grid of 20x10x3 elements, 
	#  :init_refine_lvl => 1 will generate a 40x20x6 grid (octree)
	#---------------------------------------------------------------------------
	:linitial_refine         => false,
	:init_refine_lvl         => 1,
	#---------------------------------------------------------------------------
        # AMR: Adaptive Mesh Refinement (different from initial_refinement)
        #---------------------------------------------------------------------------
        :ladapt                  => false,
        :amr_freq                => 200, #AMR is triggered every :amr_freq steps
        :amr_max_level           => 2,   #levels of oct/quadtree subdivisions
	#---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat               => "vtk",
        :loverwrite_output       => true,
        :lwrite_initial          => false,
        :output_dir              => "./output",
        :loutput_pert            => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
```
