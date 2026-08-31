#---------------------------------------------------------------------------------
# Mach-3 flow over a FOUR-ELEMENT high-lift section.
#
# This is problems/CompEuler/naca64A210 with the single aerofoil replaced by a
# slat, a main element and two trailing-edge flaps, deployed as a transport wing
# in the landing configuration — a 747-style layout. It is not a reproduction of
# any published 747 section; those ordinates are not public. Every element is a
# NACA 64A210 of its own chord, deflection and trailing-edge radius, which is
# what lets the whole four-body configuration be stated exactly in
# ./user_exactGeo.jl.
#
# WHAT THIS CASE DEMONSTRATES that the single aerofoil did not: :exact_geometry
# is a Dict of BOUNDARY TAGS, and the kernel curves each named boundary against
# whatever the case says that boundary is. Four elements are therefore four
# entries in the deck and four physical groups in the mesh — no new kernel code,
# no new hook, and user_bc.jl does not change at all.
#
# GEOMETRY (grid b747_highlift.msh, written by ./make_mesh.jl), in the
# shock_circle tunnel, 3 m long and 2 m high:
#
#     slat    chord 0.160 m, LE (0.860, -0.085), drooped 25°
#     main    chord 0.600 m, LE (1.000,  0.000), 0°
#     flap1   chord 0.240 m, LE (1.570, -0.045), lowered 30°
#     flap2   chord 0.150 m, LE (1.795, -0.172), lowered 55°
#
# The flaps sit AFT of the main element's base, which ends at x = 1.548: flap1
# begins at 1.570 and the pair extends the section to x = 1.869. They are behind
# the main, not tucked beneath it.
#
# EVERY EDGE IS ROUND. The 64A210 has a rounded leading edge (radius 0.0056c);
# each trailing edge is truncated where a circle of the stated radius is tangent
# to both surfaces and capped with that circle. So no point on any of the four
# boundaries is a cusp — a cusp is a point at which the boundary has no normal,
# and the free-slip condition is applied edge by edge, so such a node receives
# two projections and neither is the wall it is on. On one element that was
# enough to drive the pressure negative in the wake; on four it would have been
# four times over.
#
# THE SLOTS ARE THE CASE. The gaps are 0.032 m (slat/main), 0.049 m (main/flap1)
# and 0.044 m (flap1/flap2), carrying 8 to 14 elements each — check_mesh.jl
# prints them. A slot two cells wide is not a slot, and the whole point of a
# high-lift section is what goes through them.
#
# Free stream, identical to shock_circle and naca64A210: a uniform Mach-3 stream
# of air at
#
#     p = 101325 Pa,  T = 293 K,  |u| = M*c ≈ 1029 m/s,  v = 0.
#
# A WORD ON THE PHYSICS. High-lift devices are low-speed hardware; running one at
# Mach 3 is not a flight condition anyone flies. This is a geometry and
# boundary-condition demonstration — four curved bodies in close proximity, each
# exactly defined — and the free stream is inherited from the sibling cases so
# the comparison is against something already understood.
#
# Total energy (NOT ρθ) is the prognostic energy variable, because ρθ is an
# entropy variable and is not conserved across a shock.
#---------------------------------------------------------------------------------

#
# Free-stream state, shared by initialize() and by the supersonic-inflow
# boundary condition in user_bc.jl.
#
# Returns (ρ, u, v, p, ρE) in SI units.
#
# The velocity is built from M∞ and the code's own PhysConst (γ = cp/cv) rather
# than hardcoded, so the stream is exactly Mach 3 for the gas the solver
# actually integrates.
#
function naca_freestream()

    PhysConst = PhysicalConst{Float64}()

    M∞ = 3.0                    # inflow Mach number
    p∞ = 101325.0               # Pa
    T∞ = 293.0                  # K

    ρ∞ = p∞/(PhysConst.Rair*T∞)          # kg/m³
    c∞ = sqrt(PhysConst.γ*p∞/ρ∞)         # m/s
    u∞ = M∞*c∞
    v∞ = 0.0

    ρE∞ = p∞/PhysConst.γm1 + 0.5*ρ∞*(u∞*u∞ + v∞*v∞)

    return ρ∞, u∞, v∞, p∞, ρE∞
end


function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        println(" Initialize fields for 2D CompEuler (b747_highlift: Mach-3 flow over a 4-element high-lift section) ... ")
    end

    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q().
    # qoutvars can hold at most neqs+1 entries — that is the width of q.qout.
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρE"]
    qoutvars = ["ρ", "u", "v", "p", "T"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------

    ρ∞, u∞, v∞, p∞, ρE∞ = naca_freestream()

    if rank == 0
        M∞ = sqrt(u∞*u∞ + v∞*v∞)/sqrt(PhysicalConst{Float64}().γ*p∞/ρ∞)
        @printf("    free stream: ρ = %.5f kg/m³, u = %.2f m/s, p = %.1f Pa, M = %.3f\n",
                ρ∞, u∞, p∞, M∞)
    end

    #
    # Uniform free stream everywhere at t = 0. The four elements are not in the
    # mesh — the grid covers the fluid only, each element is a hole — so there is
    # no interior region to mask out: every node is fluid.
    #
    for ip = 1:mesh.npoin
        q.qn[ip,1]   = ρ∞
        q.qn[ip,2]   = ρ∞*u∞
        q.qn[ip,3]   = ρ∞*v∞
        q.qn[ip,4]   = ρE∞
        q.qn[ip,end] = p∞

        # Reference state = the free stream. Nothing in this case runs in
        # PERT() mode, but qe is what the perturbation output and the PERT()
        # branches of the user_* routines read, so it is filled with a
        # meaningful state rather than zeros.
        q.qe[ip,1]   = ρ∞
        q.qe[ip,2]   = ρ∞*u∞
        q.qe[ip,3]   = ρ∞*v∞
        q.qe[ip,4]   = ρE∞
        q.qe[ip,end] = p∞
    end

    if rank == 0
        println(" Initialize fields for 2D CompEuler (b747_highlift: Mach-3 flow over a 4-element high-lift section) ... DONE ")
    end

    return q
end
