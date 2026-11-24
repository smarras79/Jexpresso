#---------------------------------------------------------------------------------
# Mach-7 LAMINAR flow over a circular cylinder.
#
# CompEuler/shock_circle raised from Mach 3 to Mach 7 and made VISCOUS, so
# that it carries a resolved boundary layer and a surface heat flux. It is
# the curved-wall rung on the way to rampCaoEtAl2021.
#
# Geometry (grid cylinder_M7.msh): box [0,3] x [-1,1] with a cylinder of
# radius R = 0.2 m centred at (1,0). The cylinder is a NO-SLIP ISOTHERMAL
# wall; the box's own top and bottom carry the free stream (the bow shock
# never reaches them — see user_bc.jl).
#
# FREE STREAM. Not sea level. A Mach-7 cylinder at p = 101325 Pa would have
# Re_D = 6e7: the boundary layer would be micrometres thick, turbulent, and
# unresolvable, and there would be no point asking for a heat flux. These
# are shock-tunnel conditions of the kind Cao et al. (2021) use for the
# ramp — low density, T = 125 K — with the pressure chosen to put the
# cylinder Reynolds number at 1.0e4:
#
#     M∞ = 7,  T∞ = 125 K,  p∞ = 5 Pa
#     ρ∞ = 1.3937e-4 kg/m³,  c∞ = 223.98 m/s,  u∞ = 1567.8 m/s
#     μ(T∞) = 8.656e-6 Pa·s (Sutherland),  Re_D = 1.01e4
#     T₀ = T∞(1 + (γ-1)/2 M²) = 1345 K
#
# so the laminar boundary layer is δ ~ R/√(Re_R) = 2.8 mm on a 200 mm
# radius, which the 2.2 mm wall cells of the grid resolve with about six
# LGL nodes at :nop => 4. Knudsen number is 1.0e-3, comfortably continuum.
# Stagnation temperature 1345 K keeps calorically-perfect air defensible;
# at sea-level T∞ the same Mach number would demand 3150 K and real-gas
# chemistry.
#
# The gas is still ideal and calorically perfect: this is a NUMERICAL test
# at Mach 7, not a physical model of Mach-7 air.
#---------------------------------------------------------------------------------

#
# Free-stream state, shared by initialize() and by user_bc.jl.
# Returns (ρ, u, v, p, T, ρE) in SI units.
#
function cyl_freestream()

    PhysConst = PhysicalConst{Float64}()

    M∞ = 7.0                    # free-stream Mach number
    p∞ = 5.0                    # Pa
    T∞ = 125.0                  # K

    ρ∞ = p∞/(PhysConst.Rair*T∞)          # kg/m³
    c∞ = sqrt(PhysConst.γ*p∞/ρ∞)         # m/s
    u∞ = M∞*c∞
    v∞ = 0.0

    ρE∞ = p∞/PhysConst.γm1 + 0.5*ρ∞*(u∞*u∞ + v∞*v∞)

    return ρ∞, u∞, v∞, p∞, T∞, ρE∞
end

#
# Isothermal wall temperature. A COLD wall: T_w/T₀ = 300/1345 = 0.22, which
# is the regime hypersonic heating measurements are made in, and it is what
# makes the surface heat flux a meaningful output rather than a number that
# drifts with the interior solution.
#
cyl_Twall() = 300.0   # K


function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        println(" Initialize fields for 2D CompEuler (shock_circle_M7: Mach-7 laminar cylinder) ... ")
    end

    qvars    = ["ρ", "ρu", "ρv", "ρE"]
    qoutvars = ["ρ", "u", "v", "p", "T"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    ρ∞, u∞, v∞, p∞, T∞, ρE∞ = cyl_freestream()

    if rank == 0
        PhysConst = PhysicalConst{Float64}()
        M∞  = u∞/sqrt(PhysConst.γ*p∞/ρ∞)
        μ∞  = inputs[:sutherland_muref]*(T∞/inputs[:sutherland_Tref])^1.5 *
              (inputs[:sutherland_Tref] + inputs[:sutherland_S])/(T∞ + inputs[:sutherland_S])
        ReD = ρ∞*u∞*0.4/μ∞
        @printf("    free stream: ρ = %.5e kg/m³, u = %.1f m/s, p = %.2f Pa, T = %.1f K, M = %.3f\n",
                ρ∞, u∞, p∞, T∞, M∞)
        @printf("    Re_D = %.3e,  T0 = %.0f K,  T_wall = %.0f K (isothermal)\n",
                ReD, T∞*(1.0 + 0.5*PhysConst.γm1*M∞*M∞), cyl_Twall())
    end

    #
    # STARTING FIELD: free stream everywhere EXCEPT a thin layer on the
    # cylinder, where it is blended to the wall condition.
    #
    # A uniform free stream is NOT a legal starting field for a no-slip
    # isothermal wall, and this case found that out the hard way. At the wall
    # nodes the boundary condition sets u = v = 0 and ρE = ρ cv T_w, so ρE
    # goes from 183.85 to 30.02 — a factor of 6 — while the node one cell
    # away is still at 183.85. The jump is dominated by the KINETIC energy
    # being removed, not by the temperature: of the -153.8 change, -171.3 is
    # the kinetic energy going and +17.5 is the wall being hotter (300 K)
    # than the free stream (125 K). Across one 2.2 mm cell that is a violent
    # startup transient, and it drove the pressure to -0.296 Pa against
    # p∞ = 5 Pa within the first steps.
    #
    # rampCaoEtAl2021 already learned this and starts from a compressible
    # laminar boundary-layer profile rather than a uniform stream (its
    # initialize.jl, "STARTING FIELD: a compressible laminar boundary layer,
    # not a uniform stream"). The same idea, adapted: a cylinder has no
    # similarity profile to lay down before the bow shock even exists, so
    # this is not a boundary layer — it is a BOUNDARY-CONDITION-CONSISTENT
    # field, nothing more. The velocity is taken to zero and the temperature
    # to T_w over a layer of thickness δ₀ using the same Pohlhausen blend the
    # ramp uses,
    #
    #     ζ  = n/δ₀,  n = r - R the wall distance
    #     su = 2ζ - 2ζ³ + ζ⁴          (su = 0 at the wall, 1 and smooth at δ₀)
    #     u  = su·u∞,  v = su·v∞,  T = T_w + (T∞ - T_w)·su
    #
    # δ₀ = 5 mm is about two wall cells and about twice the estimated laminar
    # δ ~ 2.8 mm, so the layer is resolved from the first step. The pressure
    # is left at p∞ throughout — a boundary layer has no pressure rise across
    # it — and the density follows from p∞ and the blended T. This is a
    # starting guess, not an answer: the real boundary layer grows out of it.
    #
    PhysConst = PhysicalConst{Float64}()
    cv   = PhysConst.Rair/PhysConst.γm1
    Tw   = cyl_Twall()
    δ₀   = 5.0e-3          # m, blend thickness: ~2 wall cells, ~2δ
    xc, yc, R = 1.0, 0.0, 0.2

    for ip = 1:mesh.npoin

        # mesh.coords[dim, ip] — the per-axis @view(mesh.coords[1,:])/@view(mesh.coords[2,:]) are deprecated.
        n  = sqrt((mesh.coords[1,ip] - xc)^2 + (mesh.coords[2,ip] - yc)^2) - R
        ζ  = clamp(n/δ₀, 0.0, 1.0)
        su = 2.0*ζ - 2.0*ζ^3 + ζ^4        # 0 at the wall, 1 and smooth at δ₀

        uL = su*u∞
        vL = su*v∞
        TL = Tw + (T∞ - Tw)*su
        ρL = p∞/(PhysConst.Rair*TL)          # constant pressure across the layer
        ρEL = ρL*cv*TL + 0.5*ρL*(uL*uL + vL*vL)

        q.qn[ip,1]   = ρL
        q.qn[ip,2]   = ρL*uL
        q.qn[ip,3]   = ρL*vL
        q.qn[ip,4]   = ρEL
        q.qn[ip,end] = p∞

        # Reference state = the undisturbed free stream, NOT the blended
        # layer. qe is what the perturbation output and the PERT() branches
        # read; this case runs TOTAL(), so nothing integrates against it.
        q.qe[ip,1]   = ρ∞
        q.qe[ip,2]   = ρ∞*u∞
        q.qe[ip,3]   = ρ∞*v∞
        q.qe[ip,4]   = ρE∞
        q.qe[ip,end] = p∞
    end

    if rank == 0
        println(" Initialize fields for 2D CompEuler (shock_circle_M7: Mach-7 laminar cylinder) ... DONE ")
    end

    return q
end
