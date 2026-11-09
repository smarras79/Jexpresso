#---------------------------------------------------------------------------------
# 3D emergence of a twisted convection-zone flux tube into a coronal dipole:
# Leake, Linton & Antiochos, ApJ 787:46 (2014), simulations SD, MD, WD, with the
# setup of Leake, Linton & Török, ApJ 778:99 (2013). See README.md.
#
#   julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2014
#   mpiexec -n 128 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2014
#
# LK_REF picks the mesh (LK.geo) and the time step:
#   1: 16³ elements  (smoke test; 19 of the paper's cells per element)
#   2: 38³ elements  (half the paper's resolution at :nop => 4)
#   3: 76³ elements  (the paper's 304³ cells at :nop => 4; generate it first:
#                     gmsh -3 LK.geo -setnumber ref 3 -o LK_76x76x76.msh)
#---------------------------------------------------------------------------------

# Positivity floors as the Runge-Kutta stage limiter (as problems/MHD/fluxEmergenceSon2025):
# ρ ≥ f ρ_e(z) and p ≥ f p_e(z) of the reference atmosphere, p raised through E at fixed
# velocity and field. Not conservative where it fires; every hit is counted and reported.
if !@isdefined(lk_floor_frac)
    const lk_floor_frac = Ref{Float64}(0.2)
    const lk_floor_hits = Ref{Int}(0)
    const lk_floor_next = Ref{Int}(1)
end

function lk_positivity_limiter!(u, integrator, p, t)
    npoin = p.mesh.npoin
    qe    = p.qp.qe
    frac  = lk_floor_frac[]
    γm1   = γ_mhd - 1.0
    nhit  = 0
    @inbounds for ip = 1:npoin
        iρ = ip; iρu = npoin + ip; iρv = 2npoin + ip; iE = 3npoin + ip; iρw = 4npoin + ip
        iBx = 5npoin + ip; iBy = 6npoin + ip; iBz = 7npoin + ip; iψ = 8npoin + ip
        ρmin = frac*qe[ip,1]
        pmin = frac*qe[ip,end]
        ρ = u[iρ]
        if ρ < ρmin
            u[iρ] = ρmin; ρ = ρmin; nhit += 1
        end
        ke = 0.5*(u[iρu]^2 + u[iρv]^2 + u[iρw]^2)/ρ
        me = 0.5*(u[iBx]^2 + u[iBy]^2 + u[iBz]^2) + 0.5*u[iψ]^2
        if γm1*(u[iE] - ke - me) < pmin
            u[iE] = pmin/γm1 + ke + me; nhit += 1
        end
    end
    lk_floor_hits[] += nhit
    if nhit > 0 && lk_floor_hits[] >= lk_floor_next[]
        println(" # positivity floors, rank ", MPI.Comm_rank(MPI.COMM_WORLD), ": ", lk_floor_hits[],
                " stage-node hits so far (", nhit, " this stage) at t = ", round(t, digits = 3))
        lk_floor_next[] = 10*lk_floor_next[]
    end
    return nothing
end

# Explicit dissipation of the papers (Re = Rm = 100, i.e. 0.01 in units of L₀v₀),
# per element centroid, added to the DynSGS coefficient of the momentum (ν) and
# field (η) slots by compute_dsgs_viscosity!(::DSGS_MHD, ::NSD_3D):
#   η = 0.01 [tanh(−(r_η − 100)/5) + 1]/2,  r_η = √(x² + y²)   (Paper I Eq. 9)
#   ν = 0.01, kinematic (README: the papers' dynamic ν would give Re ~ 10⁻⁶ in the corona)
lk_dissipation(x, y, z) = (0.01, 0.01*0.5*(tanh(-(sqrt(x*x + y*y) - 100.0)/5.0) + 1.0))

function user_inputs()

    LK_REF = 2
    mesh   = ("LK_16x16x16.msh", "LK_38x38x38.msh", "LK_76x76x76.msh")[LK_REF]
    # Δt for c_h Δt/Δh ≈ 0.28 at c_h ≈ 200 (the fast speed just above the transition
    # region for SD; 155 for MD, 106 for WD) on the smallest LGL gap Δh of each mesh (README).
    Δt     = (2.0e-3, 8.0e-4, 4.0e-4)[LK_REF]

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(stage_limiter! = lk_positivity_limiter!,
                                                      williamson_condition = false),
        :Δt                   => Δt,
        :tinit                => 0.0,
        :tend                 => 160.0,      # Paper II Figs. 10-11: MD erupts at 120-160 t₀
        :diagnostics_at_times => (0.0:5.0:160.0),
        :conservation_every   => 10,
        :conservation_slots   => [1, 4],
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,       # gravity, ψ damping, velocity damping
        :lwell_balanced       => true,       # subtract the discrete RHS of qe (atmosphere + dipole)
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # The papers' problem (README):
        #---------------------------------------------------------------------------
        :lk_dipole            => "MD",       # SD, MD, WD (Φ_dip/Φ_tube = 0.13, 0.1, 0.067) or ND
        :lk_tube              => true,
        :lk_cz_gradient       => 1.0,        # a of Paper I Eq. 11 (1.25 = exactly adiabatic)
        :lk_glm_alpha         => 1.0,        # ψ damping rate α_p c_h/Δh (README: GLM transient)
        :lk_damp_xd           => 96.0,
        :lk_damp_zd           => 180.0,      # Paper II (Paper I: 170)
        :mhd_dissipation      => lk_dissipation,
        #---------------------------------------------------------------------------
        # DynSGS in conserved form on the departure from the reference state, as the
        # 2D flux-emergence case: legacy sensor, element norms floored at the local
        # scales, fast-speed floors (β ~ 10⁻² in the corona), thermal-energy residual.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model           => DSGS_MHD(),
        :dsgs_sensor          => "legacy",
        :dsgs_CR              => 1.0,
        :dsgs_Cmax            => 0.5,
        :dsgs_Cmin            => 0.0,
        :dsgs_gamma           => 5.0/3.0,
        :dsgs_Prt             => 1.0,
        :dsgs_norms           => "element",
        :dsgs_local_rel       => 1.0,
        :dsgs_conserved       => true,
        :dsgs_fast_floors     => true,
        :dsgs_thermal_residual => true,
        :lrichardson          => false,
        :energy_equation      => "energy",
        :lkep                 => false,
        :entropy_variables    => false,
        #---------------------------------------------------------------------------
        # Mesh: LK.geo (the papers' stretching, element boundaries every 4 cells at ref 3)
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/MHD/fluxEmergenceLeake2014/" * mesh,
        :lfilter              => false,
        :outformat            => "vtk",
        :loverwrite_output    => false,
        :lwrite_initial       => true,
        :output_dir           => "./output/",
        :loutput_pert         => false,
        :linitial_refine      => false,
        :ladapt               => false,
    ) #Dict

    return inputs
end
