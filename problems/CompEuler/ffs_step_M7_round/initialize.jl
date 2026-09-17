#---------------------------------------------------------------------------------
# Mach-7 supersonic flow over a forward-facing step, FILLETED CORNER.
#
# This case is CompEuler/ffs_step_M7 with ONE change: the convex step corner
# at (0.6, 0.2) is filleted with an arc of radius r = 0.05 (two elements)
# instead of being a sharp 90 degrees. Free stream, fluxes, boundary
# conditions, DynSGS settings and time step are identical; only the mesh
# differs, and the corner special case in user_bc.jl is gone because the
# fillet makes it unnecessary.
#
# WHY. At a sharp (0.6, 0.2) the fluid's interior angle is 270 degrees — a
# POINT SINGULARITY of the Euler solution. The expansion there is a
# Prandtl-Meyer fan centred on that one point, and what it sheds downstream
# is an entropy layer, a contact-type feature that convects and never
# self-heals. Woodward & Colella (1984) section IV call it "a singular point
# of the flow" and reset the state in the cells beside it for exactly this
# reason, at Mach 3.
#
# A residual-based viscosity cannot fix that: an expansion fan is a SMOOTH
# solution, so the DynSGS sensor correctly returns almost nothing in it. On
# ffs_step_M7 that shows up as spurious structures streaming downstream of
# the corner, which survived a COMPLETE removal of the top mode by the modal
# filter — so they are not grid-scale ringing, they come from the corner.
# This deck removes the singularity from the geometry instead of patching
# the solution. See README.md.
#
# Geometry (grid ffs_step_round.msh): a wind tunnel 3 m long and 1 m
# high with a step 0.2 m high starting 0.6 m from the inflow, so the fluid
# region is the L-shape
#
#     ([0,3] x [0,1]) \ ([0.6,3] x [0,0.2]),
#
# with the step's top-left corner rounded by an arc of radius r = 0.05
# centred at (0.65, 0.15), tangent to the step face at (0.6, 0.15) and to
# the step top at (0.65, 0.2).
#
# The tunnel is filled with, and continuously fed from the left by, a
# uniform Mach-7 stream of air at
#
#     p = 101325 Pa,  T = 293 K,  |u| = M*c = 2400.4 m/s,  v = 0.
#
# The thermodynamic state is that of the Loci/STREAM "2D Supersonic Forward
# Step" tutorial (the dimensional form of Emery 1968 / Woodward & Colella
# 1984, Section 5.1 of Nazarov & Hoffman, IJNMF 71:339-357, 2013); only the
# Mach number is raised.
#
# WHAT CHANGES AT MACH 7. The bow shock is stronger and stands closer to the
# step, the shock layer is thinner, and the post-shock state is far more
# extreme: the stagnation temperature is T∞(1 + (γ-1)/2 M²) ≈ 3150 K against
# 820 K at Mach 3, and the normal-shock density ratio is 5.5 against 3.9
# (both at γ = cp/cv = 1.398). The gas stays calorically perfect here — this
# is the ideal-gas Euler system, no dissociation, no vibrational excitation —
# so the run is a numerical test at Mach 7, not a physical model of a
# Mach-7 air flow.
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
# The velocity is built from M∞ and the code's own PhysConst (γ = cp/cv)
# rather than hardcoded, so the stream is exactly Mach 7 for the gas the
# solver actually integrates: c∞ = 342.9 m/s gives u∞ = 2400.4 m/s. (The
# same expression at M∞ = 3 returns 1028.7 m/s, the value the Loci/STREAM
# tutorial quotes — this is the one number that differs from ffs_step.)
#
function ffs_freestream()

    PhysConst = PhysicalConst{Float64}()

    M∞ = 7.0                    # inflow Mach number  (ffs_step: 3.0)
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
        println(" Initialize fields for 2D CompEuler (ffs_step_M7_round: Mach-7 step, filleted corner) ... ")
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
        println(" Initialize fields for 2D CompEuler (ffs_step_M7_round: Mach-7 step, filleted corner) ... DONE ")
    end

    return q
end
