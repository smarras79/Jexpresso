#---------------------------------------------------------------------------------
# Mach-3 supersonic flow over a forward-facing step.
#
# Geometry (grid ffs_step_transfinite.msh): a wind tunnel 3 m long and 1 m
# high with a step 0.2 m high starting 0.6 m from the inflow, so the fluid
# region is the L-shape
#
#     ([0,3] x [0,1]) \ ([0.6,3] x [0,0.2]).
#
# Initial and boundary state follow the Loci/STREAM "2D Supersonic Forward
# Step" tutorial: the tunnel is filled with, and continuously fed from the
# left by, a uniform Mach-3 stream of air at
#
#     p = 101325 Pa,  T = 293 K,  |u| = M*c ≈ 1029 m/s,  v = 0.
#
# This is the dimensional form of the classical Emery (1968) / Woodward &
# Colella (1984) problem, which is also Section 5.1 of Nazarov & Hoffman
# (IJNMF 71:339-357, 2013) — the reference for the DynSGS shock capturing
# this case is stabilized with.
#
# The flow is discontinuous from the first instant: a bow shock forms in
# front of the step, reflects off the top wall, and the reflections
# interact to build the familiar Mach stem. Total energy (NOT ρθ) is the
# prognostic energy variable, because ρθ is an entropy variable and is not
# conserved across a shock.
#---------------------------------------------------------------------------------

#
# Free-stream state, shared by initialize() and by the supersonic-inflow
# boundary condition in user_bc.jl.
#
# Returns (ρ, u, v, p, ρE) in SI units.
#
# The tutorial quotes |u| = 1029 m/s, which is Mach 3 for air at 293 K.
# The velocity is rebuilt here from M∞ and the code's own PhysConst
# (γ = cp/cv) instead of being hardcoded, so the stream is exactly Mach 3
# for the gas the solver actually integrates; the number that comes out is
# 1028.9 m/s, i.e. the tutorial value.
#
function ffs_freestream()

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
        println(" Initialize fields for 2D CompEuler (ffs_step: Mach-3 forward-facing step) ... ")
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

    ρ∞, u∞, v∞, p∞, ρE∞ = ffs_freestream()

    if rank == 0
        M∞ = sqrt(u∞*u∞ + v∞*v∞)/sqrt(PhysicalConst{Float64}().γ*p∞/ρ∞)
        @printf("    free stream: ρ = %.5f kg/m³, u = %.2f m/s, p = %.1f Pa, M = %.3f\n",
                ρ∞, u∞, p∞, M∞)
    end

    #
    # Uniform free stream everywhere at t = 0. The step is not in the mesh
    # (the grid covers the fluid L-shape only), so there is no interior
    # region to mask out: every node is fluid.
    #
    for ip = 1:mesh.npoin
        q.qn[ip,1]   = ρ∞
        q.qn[ip,2]   = ρ∞*u∞
        q.qn[ip,3]   = ρ∞*v∞
        q.qn[ip,4]   = ρE∞
        q.qn[ip,end] = p∞

        # Reference state = the free stream. Nothing in this case runs in
        # PERT() mode, but qe is what the perturbation output and the
        # PERT() branches of the user_* routines read, so it is filled
        # with a meaningful state rather than zeros.
        q.qe[ip,1]   = ρ∞
        q.qe[ip,2]   = ρ∞*u∞
        q.qe[ip,3]   = ρ∞*v∞
        q.qe[ip,4]   = ρE∞
        q.qe[ip,end] = p∞
    end

    if rank == 0
        println(" Initialize fields for 2D CompEuler (ffs_step: Mach-3 forward-facing step) ... DONE ")
    end

    return q
end
