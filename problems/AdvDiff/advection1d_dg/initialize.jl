#
# Initial conditions for the 1D scalar advection DG case.
#
# Both ICs below are transcriptions of cases in the Python reference's
# `exact_solution` (drl-amr/backends/python_1d/solvers/utils.py), so the two
# codes can be compared on identical data:
#
#   :sine     -> Python icase 8 : sin(pi*x)
#                Analytic, exactly periodic on [-1,1]. The order-study IC.
#                Integrates to zero, so it gives no relative scale for a
#                conservation check.
#
#   :gaussian -> Python icase 1 : exp(-256 (x - xbar)^2) plus periodic images
#                Nonzero mean, so mass conservation can be reported relative to
#                a nonzero integral. Sharp (half-width ~0.044): expect it to be
#                pre-asymptotic below roughly 32 elements.
#
# Flip IC_KIND to switch. Kept as a constant rather than an input key because
# mod_inputs.jl validates the deck and an unrecognized key is a risk we do not
# need to take for a two-way switch.
#
const IC_KIND = :sine


function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for 1D scalar advection (DG) ........... ")

    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is used to allocate all necessary equation-dependent arrays
    #---------------------------------------------------------------------------------
    qvars = ["q"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat,
                 inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------

    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl

            ip = mesh.connijk[iel_g, i, 1]

            # NOTE: mesh.coords, not mesh.x.
            # St_mesh carries both as separate arrays. AdvDiff/case1 reads
            # mesh.x[ip]; wave1d_dg reads mesh.coords[ip,1] and is known to work
            # under the DG numbering. Until it is confirmed that the DG builder
            # (add_high_order_nodes_1D_native_mesh_dg!) fills `x` as well as
            # `coords`, reading `x` here risks silently wrong node coordinates
            # -- which would look like a broken operator rather than a broken IC.
            x = mesh.coords[ip, 1]

            q.qn[ip, 1] = initial_profile(x, IC_KIND)

            # Background state for plotting / perturbation analysis
            q.qe[ip, 1] = 0.0

        end
    end

    println(" Initialize fields for 1D scalar advection (DG) ........... DONE ")

    return q
end


#
# q(x, 0). Kept separate so the same profile can be evaluated at an arbitrary
# advected position when an exact solution is needed for error measurement.
#
function initial_profile(x, kind::Symbol)

    if kind === :sine
        # Python icase 8. Period 2 = domain length, so exactly periodic.
        return sin(pi * x)

    elseif kind === :gaussian
        # Python icase 1: main pulse plus the two nearest periodic images, so
        # the profile is continuous to machine precision across the wrap.
        beta = 256.0
        L    = 2.0
        return exp(-beta * x^2) +
               exp(-beta * (x - L)^2) +
               exp(-beta * (x + L)^2)

    else
        error("initialize.jl: unknown IC_KIND = $kind (expected :sine or :gaussian)")
    end
end
