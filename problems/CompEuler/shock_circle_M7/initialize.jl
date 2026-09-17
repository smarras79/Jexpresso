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
    # Uniform free stream everywhere at t = 0. The cylinder is a hole in the
    # mesh, so every node is fluid; the no-slip wall is imposed by the BC and
    # the boundary layer grows out of this impulsive start.
    #
    for ip = 1:mesh.npoin
        q.qn[ip,1]   = ρ∞
        q.qn[ip,2]   = ρ∞*u∞
        q.qn[ip,3]   = ρ∞*v∞
        q.qn[ip,4]   = ρE∞
        q.qn[ip,end] = p∞

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
