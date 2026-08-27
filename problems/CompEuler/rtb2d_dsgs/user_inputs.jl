#=============================================================================
 rtb2d_dsgs -- the 2D rung of the DynSGS ladder.

 A rising thermal bubble (Robert 1993; Giraldo & Restelli 2008) stabilised by
 the Marras-Nazarov Dynamic SGS model. 10000 x 10000 m, 10 x 10 elements at
 p = 4: 1681 gridpoints, seconds.

 WHY THIS EXISTS ALONGSIDE CompEuler/theta_dsgs
 ----------------------------------------------
 Two reasons, and the first is the practical one:

 1. theta_dsgs points at ./meshes/gmsh_grids/hexa_TFI_10x10.msh, and `meshes`
    is a symlink to the separate JexpressoMeshes repository -- gitignored,
    absent from a fresh clone and from CI. So the one 2D DynSGS case in the
    repository could not be run by someone who has only this repository, which
    is the wrong property for the case you reach for when the 3D one
    misbehaves. This one's mesh is committed next to it.

 2. It is DELIBERATELY the same problem as CompEuler/rtb3d_dsgs: same extents,
    same 10 x 10 elements in the x-z plane, same nop, same bubble (centre
    (0, 2500) m, radius 2000 m, 2 K, linear taper), same dt, same tend. A 3D
    slab of a 2D problem IS the 2D problem, so the pair is a ladder -- if 2D
    is right and 3D is not, the difference is in the 3D path and nowhere else.

 THE ONE MODEL DIFFERENCE BETWEEN THE TWO RUNGS, and it is not a bug
 -------------------------------------------------------------------
 The theta equation's coefficient:

     2D DSGS path   kappa = Pr/(gamma-1) * mu  = 0.25 mu   with :Pr => 0.1
     3D DSGS path   kappa = mu / Pr_t          = 1.43 mu   with Pr_t = 0.7

 a factor of 5.7. The 2D form is Nazarov & Hoffman's ARTIFICIAL heat
 conduction, written for the internal energy of a shock-capturing scheme; the
 3D path routes theta through SGS_diffusion like every other 3D closure, and
 an atmospheric potential temperature wants a turbulent Prandtl number. To put
 the two rungs on exactly equal footing, run this one with

     :Pr => 0.7/(gamma-1) ... no. Simpler, on the 3D side:
     :mu => [0.0, 1.0, 1.0, 1.0, 0.175]        # 0.25 * Pr_t

 which gives the 3D theta equation the 2D path's effective coefficient.

 The other differences are geometric and known: the 2D path's element scale is
 Delta_elem/ngl = 200 m against the 3D path's Delta_elem_filter/nop = 250 m
 (a factor 1.56 in Delta^2), and the 3D viscous kernel is the full stress
 tensor where the 2D one is a scalar Laplacian on the primitives.

 WHAT TO LOOK AT
 ---------------
   theta        the classic mushroom, symmetric about x = 0.
   mu_rhou      this case's initialize.jl puts the DynSGS coefficients in
                qoutvars, so they come out as mu_rhou / mu_rhov / kappa_theta
                alongside the solution. Expect them CONCENTRATED on the
                bubble's edge -- that is the whole claim of a residual sensor.
                Uniform, or pinned at the C2*Delta*(|v|+c) cap, means the
                normalisation or the residual is wrong, not that the case is
                hard.

     DBG_DT=...  DBG_TEND=...  DBG_VISC=0 (closure off)  DBG_MESH=...
=============================================================================#

function user_inputs()

    tend = parse(Float64, get(ENV, "DBG_TEND", "1000.0"))

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => parse(Float64, get(ENV, "DBG_DT", "0.5")),
        :tinit                => 0.0,
        :tend                 => tend,
        :diagnostics_at_times => (0:100:tend),
        :lrestart             => false,
        :case                 => "rtb",
        :lsource              => true,
        # TOTAL, not PERT: the DynSGS coefficient is built from the
        # conservative variables, and theta_dsgs runs TOTAL for the same
        # reason. rtb3d_dsgs defaults to TOTAL too, so the rungs match.
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # The closure
        #---------------------------------------------------------------------------
        :lvisc                => parse(Bool, get(ENV, "DBG_VISC", "true")),
        :visc_model           => DSGS(),
        :energy_equation      => "theta",
        # Per-equation MASK on the coefficient the model computed. Slot 1 stays
        # 0 so mass is strictly conservative (Marras eq. 10).
        :μ                    => [0.0, 1.0, 1.0, 1.0],
        # Artificial Prandtl number of Marras eq. (10b) / Nazarov & Hoffman
        # eq. (3.7). This is the 2D path's theta coefficient -- see the header
        # for how it relates to the 3D path's Pr_t.
        :Pr                   => 0.1,
        #---------------------------------------------------------------------------
        # Mesh -- committed next to this deck, see rtb2d_dsgs_10x10.geo
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => get(ENV, "DBG_MESH",
                                     "./problems/CompEuler/rtb2d_dsgs/rtb2d_dsgs_10x10.msh"),
        #---------------------------------------------------------------------------
        # Plotting
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,      # plot the total state
        #---------------------------------------------------------------------------
        # AMR off: the point is the closure, not the grid.
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :init_refine_lvl      => 1,
        :ladapt               => false,
        :amr_freq             => 200,
        :amr_max_level        => 2,
    ) #Dict

    return inputs

end
