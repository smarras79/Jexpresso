#=============================================================================
 rtb3d_dsgs -- the smallest honest test of DynSGS in 3D.

 A rising thermal bubble (Robert 1993; Giraldo & Restelli 2008) that is a
 CYLINDER along y, on a mesh one element thick in y. 10000 x 1000 x 10000 m,
 10 x 1 x 10 elements at p = 4: 41 x 5 x 41 = 8405 gridpoints, one core,
 minutes. It exists so that the 3D DynSGS path -- new, and structurally
 different from the 1D and 2D ones -- can be exercised somewhere a wrong
 answer is recognisable, instead of only on LESICP2, where 15.9 M nodes and
 three hours stand between a mistake and noticing it.

 WHAT MAKES IT A TEST AND NOT JUST A SMALL RUN
 ---------------------------------------------
 Three things are checkable here and are not checkable on a real LES:

 1. THE SOLUTION MUST BE y-INVARIANT.  initialize.jl measures the bubble
    radius in the x-z plane only, the y faces are free-slip, and nothing else
    in the problem has a y dependence. So the exact solution has v == 0 for
    all time, and `v` in the VTU is a free correctness check on the ENTIRE 3D
    path -- the stress tensor tau_ij including its -2/3 div(u) term, the
    metric terms, the DSS, and the DynSGS coefficient itself. Anything above
    round-off in v means something is wrong, and it is much easier to see than
    a slightly-too-diffuse bubble.

 2. IT SHOULD REPRODUCE THE 2D CASE.  The geometry, the bubble and the
    resolution are those of CompEuler/theta and CompEuler/theta_dsgs -- centre
    (0, 2500) m, radius 2000 m, amplitude 2 K, linear taper theta_c(1 - r/r0),
    1000 m cells at nop = 4. A 3D slab of a 2D problem is the 2D problem, so
    the theta field at t = 1000 s should match theta_dsgs's.

    NOT BIT-FOR-BIT, and the reason is worth knowing before comparing: the 2D
    DSGS path gives slot 4 the Nazarov & Hoffman artificial conduction
    Pr/(gamma-1)*mu with :Pr => 0.1, i.e. 0.25*mu, while the 3D path routes
    theta through SGS_diffusion like every other 3D closure and gives it
    mu/Pr_t = 1.43*mu. That is a factor 5.7 on the THETA diffusivity for the
    same coefficient -- deliberate (see DSGS.md sec. 5.1: theta is not
    internal energy and an atmospheric scalar wants a turbulent Prandtl
    number), and it means the 3D bubble will be slightly more diffuse in
    theta. To compare the two paths on equal footing set

        :mu => [0.0, 1.0, 1.0, 1.0, 0.175]      # 0.25 * Pr_t

    which reproduces the 2D path's effective theta coefficient exactly.

 3. THE PERTURBATION NORMALISATION IS VISIBLE.  This is the one formulation
    change the 3D path makes (DSGS.md sec. 5.2), and this case has the
    stratification that motivates it in miniature. Over the 10 km column the
    reference state runs rho 1.17 -> 0.4 and rho*theta 350 -> 120, while the
    bubble's departure from it is O(2). Normalising eq. (9) on the TOTAL field
    would divide the residual by ~230 instead of ~2 -- a factor 100 -- and
    switch DynSGS off. If the bubble here comes out visibly UNDER-stabilised
    (Gibbs ringing on the theta front, mu_dsgs_* near zero everywhere), that
    is the first thing to check.

 WHAT TO LOOK AT
 ---------------
   v                   must be round-off.  This is the primary check.
   theta               the classic mushroom, symmetric about x = 0, tops out
                       around z = 5-6 km at t = 1000 s.
   mu_dsgs_ρu          the DynSGS momentum coefficient, per element (the VTK
                       fields are named after qvars, so mu_dsgs_ρ, mu_dsgs_ρu,
                       ..., mu_dsgs_ρθ). Expect it CONCENTRATED on the
                       bubble's edge -- that is the whole claim of a residual
                       sensor. If it is uniform, the normalisation or the
                       residual is wrong. It is piecewise constant per element
                       by construction, so it will look blocky.
   mu_dsgs_ρθ          should be mu_dsgs_ρu / Pr_t = 1.43x it, everywhere.
   mu_dsgs_ρ           identically zero, because :μ[1] = 0.

 MEASURED (this deck has been run)
 ---------------------------------
 Both rungs of the ladder run to t = 1000 s. On this one the DynSGS
 coefficient settles at

     nu  max ~90-120 m^2/s, mean ~60          (kinematic)

 against a C2*Delta*(|v|+c) cap of ~4.3e4 -- i.e. the cap is inert and mu_res
 governs, which is the regime the model is designed for. rtb2d_dsgs on the
 same bubble gives max ~190-215; the difference is Delta (250 m here against
 the 2D path's 200 m, so 1.56x in Delta^2) and the theta split noted above.

 The worst residual in the domain sits on the BUBBLE EDGE, where the initial
 condition's linear taper theta_c(1 - r/r0) has its slope discontinuity. That
 is the whole claim of a residual sensor and it is worth checking rather than
 assuming:

     JEXPRESSO_DSGS_MONITOR=1

 prints one line per step with nu max/mean, the normalising denominators, and
 the coordinates of the node that produced the worst ratio -- flagging it if
 that node is ON A BOUNDARY, which is the signature of the residual picking up
 a boundary condition rather than a discretisation error.

 It is not registered in test/ci_cases.jl because that needs a committed
 reference solution in test/CI-ref/ and a copy of the deck under
 test/CI-runs/CompEuler/rtb3d_dsgs/ (see test/CIdescription.md). It is small
 enough and fast enough to belong there; generating the reference is the only
 missing step.

     DBG_PERT=1      run with :SOL_VARS_TYPE => PERT() instead of TOTAL()
     DBG_SMAG=1      :dsgs_add_smagorinsky => true
     DBG_VISC=0      closure off entirely, to see what it is holding
     DBG_DT=...  DBG_TEND=...  DBG_MESH=...
=============================================================================#

function user_inputs()

    # TOTAL is the default because it is what CompEuler/theta_dsgs runs and
    # what LESICP2 runs, so it is the path a mistake would actually reach.
    # DBG_PERT=1 exercises the other branch of compute_dsgs_viscosity!'s
    # `lpert` handling, which is otherwise not covered by any deck: under PERT
    # the stored state IS the departure from qe, so the norms are taken on it
    # directly and the total is recovered by ADDING qe -- the opposite of the
    # TOTAL branch, and silently swapping the two costs a factor of 100 in the
    # denominator on this problem (see note 3 in the header).
    lpert = parse(Bool, get(ENV, "DBG_PERT", "false"))

    inputs = Dict(
        #---------------------------------------------------------------------------
        # Time integration
        #
        # Fully explicit. There is no reason to put an IMEX split in front of a
        # test whose purpose is the CLOSURE: the acoustic stiffness here is
        # mild (isotropic 1000 m cells, so no anisotropy to exploit) and a
        # stage solve between the residual and the viscosity is one more thing
        # to be wrong about. dt = 0.5 s is CompEuler/theta's, on the same mesh.
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => parse(Float64, get(ENV, "DBG_DT", "0.5")),
        :tinit                => 0.0,
        :tend                 => parse(Float64, get(ENV, "DBG_TEND", "1000.0")),
        :diagnostics_at_times => (0:100:parse(Float64, get(ENV, "DBG_TEND", "1000.0"))),
        :lrestart             => false,
        :case                 => "rtb",
        :lsource              => true,     # -rho*g on the w equation
        :SOL_VARS_TYPE        => lpert ? PERT() : TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # THE CLOSURE -- the reason this case exists
        #---------------------------------------------------------------------------
        :lvisc                => parse(Bool, get(ENV, "DBG_VISC", "true")),
        :visc_model           => DSGS(),
        :energy_equation      => "theta",
        # Marras eq. (9). C1 is DELIBERATELY NOT THE PAPER'S 1.0 -- see below.
        #
        # The sensor this code applies is a TENDENCY sensor, not the paper's
        # residual (see :dsgs_residual, and the long comment on DSGS_STRICT in
        # kernel/physics/SGS.jl), so R ~ 1.5|dq/dt| and C1 = 1 has no authority
        # here. On this bubble it over-fires during the spin-up: for the first
        # ~10 s the flow accelerates from rest on a time scale of a few
        # seconds, and 1.5|dq/dt| over the field's own still-small norm comes
        # out at O(0.5) 1/s. At Delta^2 = 62500 m^2 that is
        #
        #     nu ~ 2-5e4 m^2/s
        #
        # i.e. an eddy frequency of 0.4 1/s for a 250 m eddy -- a velocity
        # scale of 100 m/s for a flow that reaches 10. About 4x too much, and
        # it happens in 2D too: rtb2d_dsgs peaks at 2.5e4 at t = 1 s and decays
        # to ~250 by t = 1000. 2D SURVIVES it because its theta slot carries
        # Pr/(gamma-1)*mu = 0.25 mu; the 3D path routes theta through
        # SGS_diffusion at mu/Pr_t = 1.43 mu and applies the full stress tensor
        # to the momenta, and it does not.
        #
        # MEASURED on this deck, dt = 0.5:
        #
        #     C1 = 1.0    dies at t ~ 7.5 s, nu pinned at the 4.3e4 cap
        #     C1 = 0.25   runs to t = 200.  nu max ~400, mean ~175
        #     C1 = 0.1    runs to t = 200.  nu max ~160, mean ~70
        #
        # And it is NOT the theta split alone -- C1 = 1 with :mu[5] = 0.175
        # (the 2D path's effective theta coefficient) still dies at t ~ 9 s --
        # nor the time step: dt = 0.25 dies at t ~ 8.
        :dsgs_C1              => parse(Float64, get(ENV, "DBG_C1", "0.25")),
        :dsgs_C2              => 0.5,
        # Add Smagorinsky to the residual viscosity instead of replacing it.
        # OFF: on this case there is no wall and no surface layer, so the
        # argument that motivates it in a PBL (DSGS.md sec. 5.7) does not
        # apply, and "does the residual term alone hold a CG-SEM bubble" is
        # exactly what this case is asking.
        :dsgs_add_smagorinsky => parse(Bool, get(ENV, "DBG_SMAG", "false")),
        :C_s                  => 0.16,       # read only when the above is true
        :lrichardson          => false,      # multiplies the Smagorinsky part only
        :lwall_damping        => false,      # no wall in this problem
        # Per-equation MASK on the coefficient the closure computed (not a
        # viscosity). Slot 1 stays 0 so the mass equation is strictly
        # conservative, Marras eq. (10). See note 2 in the header for the one
        # value to change if you are comparing against theta_dsgs.
        :μ                    => [0.0, 1.0, 1.0, 1.0,
                                   parse(Float64, get(ENV, "DBG_MU5", "1.0"))],
        # Domain rather than rank-local <q'> and ‖q' - <q'>‖. Costs two
        # Allreduce per RHS call and is a no-op in serial, which is how this
        # case is meant to be run; set it if you take it to many ranks and want
        # mu independent of the partitioning.
        # :ldsgs_global_norms   => true,
        #---------------------------------------------------------------------------
        # Mesh
        #
        # :les_filter_width is DELIBERATELY NOT SET. The mesh is isotropic --
        # dx = dy = dz = 1000 m -- so :max (the default), :geometric and :min
        # all give the same 1000 m, and no choice here can be wrong. That is a
        # property of the .geo, and its header explains why the y extent is one
        # element wide rather than a slab.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => get(ENV, "DBG_MESH",
                                     "./problems/CompEuler/rtb3d_dsgs/rtb3d_dsgs_10x1x10.msh"),
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
