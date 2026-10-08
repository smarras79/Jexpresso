#---------------------------------------------------------------------------------
# FREE-STREAM PRESERVATION TEST on the Cao et al. (2021) compression-ramp mesh.
#
# THIS IS NOT A FLOW.  It is a one-property unit test of the discretisation,
# and it exists because eight Mach-7 runs across three geometries all failed
# the same way and none of them ever measured this property.
#
# THE PROPERTY.  A uniform free stream is an exact solution of the Euler
# equations: every flux is constant, so dF/dx + dG/dy = 0 POINTWISE.  A
# correct discretisation must therefore return an RHS of exactly zero
# (machine zero) at every node, for ever.  If it does not, the scheme
# manufactures a source term out of nothing, and at hypersonic Mach that
# source is amplified into the pressure by
#
#     (gamma-1) rhoE/p = 1 + gamma(gamma-1) M^2/2
#
# which is 3.5 at M = 3 and 17.5 at M = 7.7.  The SAME relative error in
# rhoE that a Mach-3 case shrugs off takes a Mach-7.7 case below p = 0.
#
# WHY NOW.  On rampCaoEtAl2021_M7 the positivity repair reported its first
# intervention at
#
#     (x, y) = (0.19659258262890678, 0.08585606812651451)
#
# x = 0.1 + 0.1 cos15 is the OUTFLOW PLANE to sixteen digits, and y is
# 2.5836e-5 below the top-right corner (0.19659, 0.08588) -- exactly one
# wall-normal spacing of ramp15_uniform.msh.  So the first negative pressure
# in that run was one node below the corner where the free-stream "top"
# Dirichlet meets the "outflow" boundary that imposes NOTHING -- in
# undisturbed free stream, at step ~200, 1.7 mm of flow travel, with no
# shock anywhere near it.  Nothing physical happens there.  Either the free
# stream is preserved or it is not, and this case answers that in seconds.
#
# THE CONFIGURATION.  Everything that could confound the measurement is
# off: no viscosity, no DynSGS, no filter, NO POSITIVITY REPAIR (a repair
# would silently hide exactly the quantity being measured), no wall.  Every
# boundary is the free stream except the outflow, which keeps the
# production condition.  What is left is the bare inviscid operator on the
# production mesh with the production state.
#
# READING THE RESULT.  The output carries dp = p - p_inf and its relative
# form.  In ParaView the Information tab's range on "dp_rel" IS the answer:
#
#   ~1e-14, everywhere, not growing   -> the free stream is preserved.  The
#                                        failure is elsewhere and this whole
#                                        line of enquiry is closed.
#   anything larger, growing in time  -> the scheme manufactures a source.
#                                        WHERE it lives is then the bug, and
#                                        the picture shows it directly.
#
# The prediction under test is that it is nonzero and confined to the
# OUTFLOW PLANE, because the outflow is the only boundary of this case
# where nothing is imposed and a missing boundary term would therefore be
# the only place it is not overwritten every stage.  FSP_OUTFLOW below
# turns that into a controlled experiment.
#---------------------------------------------------------------------------------

#
# Exactly rampCaoEtAl2021's free stream -- Table 1 of the paper, M = 7.7,
# p = 760 Pa, T = 125 K.  Duplicated here rather than included so this test
# cannot be perturbed by edits to the production case.
#
function fsp_freestream()

    PhysConst = PhysicalConst{Float64}()

    M∞ = 7.7
    p∞ = 760.0
    T∞ = 125.0

    ρ∞ = p∞/(PhysConst.Rair*T∞)
    c∞ = sqrt(PhysConst.γ*PhysConst.Rair*T∞)
    u∞ = M∞*c∞
    v∞ = 0.0

    ρE∞ = p∞/PhysConst.γm1 + 0.5*ρ∞*(u∞*u∞ + v∞*v∞)

    return ρ∞, u∞, v∞, p∞, T∞, ρE∞
end


function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    qvars    = ["ρ", "ρu", "ρv", "ρE"]
    qoutvars = ["ρ", "u", "v", "dp", "dp_rel"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    ρ∞, u∞, v∞, p∞, T∞, ρE∞ = fsp_freestream()

    if rank == 0
        println(" Initialize fields for 2D CompEuler (rampFreeStreamTest) ... ")
        @printf("    UNIFORM free stream everywhere: rho = %.6f, u = %.1f, p = %.1f, M = 7.7\n",
                ρ∞, u∞, p∞)
        @printf("    outflow condition under test: %s\n", string(FSP_OUTFLOW))
        @printf("    wall condition under test:    %s\n", string(FSP_WALL))
        println("    The exact answer is dp = 0 at every node for all time.")
        println("    Plot dp_rel; its RANGE in ParaView is the measurement.")
    end

    #
    # No profile, no blending, no boundary layer: the uniform state, to the
    # last bit, at every node.  Any departure from it later is the
    # discretisation talking.
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
        println(" Initialize fields for 2D CompEuler (rampFreeStreamTest) ... DONE ")
    end

    return q
end
