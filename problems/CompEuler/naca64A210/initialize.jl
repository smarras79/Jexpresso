#---------------------------------------------------------------------------------
# Mach-3 flow over a NACA 64A210 aerofoil.
#
# This is problems/CompEuler/shock_circle with the circular hole replaced by an
# aerofoil, and it exists to show :exact_geometry on a shape that has no closed
# form. The tunnel, the free stream, the equations, the shock capturing and the
# boundary conditions are shock_circle's; the whole of the difference is the
# geometry, and the whole of the geometry is ./user_exactGeo.jl.
#
# Geometry (grid naca64A210.msh, written by ./make_mesh.jl): the shock_circle
# tunnel, 3 m long and 2 m high, with a NACA 64A210 of chord 0.6 m at zero
# incidence, leading edge at (1, 0). The section is cambered — it is a 10%-thick
# 6A-series aerofoil on the a = 0.8 (modified) mean line, design lift
# coefficient 0.2 — so the flow is asymmetric even at α = 0.
#
# Free stream, identical to shock_circle: a uniform Mach-3 stream of air at
#
#     p = 101325 Pa,  T = 293 K,  |u| = M*c ≈ 1029 m/s,  v = 0.
#
# The TRAILING EDGE IS ROUNDED, not the published cusp: the section is truncated
# at x/c = 0.905 and capped with a circle of radius 0.010c tangent to both
# surfaces (see ./user_exactGeo.jl). A cusp is a point at which the boundary has
# no normal — the two boundary edges meeting there carry normals 168° apart, and
# the free-slip condition is applied edge by edge, so that node gets both
# projections and neither is the wall it is on. Rounding removes it: the largest
# turn of the wall normal between adjacent boundary edges is now 40°, and it is
# at the LEADING edge. The base is 0.012 m across, about 5 elements.
#
# What to expect. At Mach 3 the rounded nose (radius 3.4e-3 m, see below) carries
# a small detached bow shock; the flow then expands over both surfaces, and the
# two streams are turned back parallel through a recompression shock system off
# the base, with a slip line between them because the camber has given the two
# sides different total pressure losses. A rounded base also has its own small
# separated wake, which the published cusp would not. The bow shock reflects
# off the tunnel walls at roughly x = 1.9 and comes back onto the wake. All of
# it is density structure, which is what :lschlieren is on for.
#
# WHY THE CURVED BOUNDARY MATTERS HERE. The grid is straight-sided: nodes on the
# section, high-order nodes in between on the chords. Around the nose a boundary
# edge is 3.4e-3 m long across a 3.4e-3 m radius, so its mid-edge nodes sit
# 3.2e-4 m INSIDE the aerofoil (check_mesh.jl measures exactly this). That is a
# ring of small forward-facing steps through the stagnation region of a Mach-3
# flow, and raising :nop does not remove it — the polygon is in the grid, not in
# the polynomial. :exact_geometry puts those nodes back on the section.
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
        println(" Initialize fields for 2D CompEuler (naca64A210: Mach-3 flow over a NACA 64A210) ... ")
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
    # Uniform free stream everywhere at t = 0. The aerofoil is not in the mesh —
    # the grid covers the fluid only, the section is a hole — so there is no
    # interior region to mask out: every node is fluid.
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
        println(" Initialize fields for 2D CompEuler (naca64A210: Mach-3 flow over a NACA 64A210) ... DONE ")
    end

    return q
end
