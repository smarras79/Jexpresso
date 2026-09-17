#---------------------------------------------------------------------------------
# positivity_driver.jl — the Jexpresso side of the realizability repair.
#
# Everything that has to know about params, inputs, MPI and the solution layout
# lives here; the numerics are in Positivity.jl next door and know none of it.
#
# WHERE IT RUNS. Inside rhs!, immediately after the stage state has been loaded
# into params.uaux and before any flux is evaluated. That placement is the whole
# point: user_flux!, user_fluxaux! (which takes log(p) for ranocha), the DSGS
# sensor and the sound speed then CANNOT be handed a non-realizable state. It
# also makes the FLUXAUX_FLOOR guard in the Mach-7 decks' user_flux.jl
# redundant — that guard was this idea, done badly and in the wrong place.
#
# It repairs the stage state of a low-storage Runge-Kutta, which means it edits
# the integrator's own registers. :lfilter already does exactly that
# (filter.jl), so the pattern is the house one, but it is worth knowing.
#
# SCOPE, KEPT DELIBERATELY SMALL. 2D/3D CompEuler, TOTAL(), ρE in the energy
# slot, CPU, and neqs == nsd + 2 exactly. Anything else is a clear error at the
# first call rather than a silent wrong repair: the θ-form carries ρθ, which is
# positive for a different reason, and the MHD system has a magnetic energy this
# repair knows nothing about. Both deserve their own treatment, not this one.
#
# OFF BY DEFAULT (:lpositivity => false), so no existing case changes at all.
#---------------------------------------------------------------------------------

# Everything below is called through the qualified name Positivity.*, so no
# `using` is needed and no exported name can collide with a Jexpresso one.

# One audit trail per run. A const Ref/instance at module scope is the pattern
# the DSGS kernels already use for this kind of cross-call state.
const POSITIVITY_STATS   = Positivity.PositivityStats()
const POSITIVITY_CHECKED = Ref(false)

#---------------------------------------------------------------------------------
# Configuration check. Runs once, on the first call, then never again.
#---------------------------------------------------------------------------------
function positivity_validate(inputs, params, neqs::Int, ien::Int)

    why = String[]

    get(inputs, :energy_equation, "theta") == "energy" ||
        push!(why, "  :energy_equation must be \"energy\" (slot $(ien) must hold ρE, not ρθ)")

    params.SOL_VARS_TYPE == TOTAL() ||
        push!(why, "  :SOL_VARS_TYPE must be TOTAL() (a perturbation state has no realizable set of its own)")

    neqs == ien ||
        push!(why, "  neqs = $(neqs) but this repair is written for exactly nsd+2 = $(ien) equations")

    get(inputs, :backend, CPU()) == CPU() ||
        push!(why, "  CPU backend only: the GPU path is not wired")

    ρmin = Float64(get(inputs, :positivity_rho_min, 0.0))
    pmin = Float64(get(inputs, :positivity_p_min,   0.0))
    (ρmin > 0.0 && pmin > 0.0) ||
        push!(why, "  :positivity_rho_min and :positivity_p_min must both be set > 0.\n" *
                   "     They are absolute and there is no safe default — pick them from the\n" *
                   "     case's own scales, e.g. 1e-6 of the free-stream ρ and p, so the repair\n" *
                   "     engages only outside the realizable set and never in the solution.")

    if !isempty(why)
        error(" # ERROR positivity_driver.jl: :lpositivity => true is not valid for this case:\n" *
              join(why, "\n") * "\n" *
              " # Set :lpositivity => false, or fix the above.")
    end

    rank = MPI.Comm_rank(get_mpi_comm())
    if rank == 0
    # println_rank, NOT @info: Julia sends @info/@warn to STDERR, and the submit
    # scripts split the streams (--output=%x.%j.out, --error=%x.%j.err), so the
    # audit trail landed in a file nobody reads while the CFL narration it has
    # to be compared against went to the other one.
        println_rank(@sprintf(" # POSITIVITY REPAIR ON: ρ_min = %.3e, p_min = %.3e (absolute). This is a REPAIR, not a preserving scheme — every intervention is counted and reported.",
                       ρmin, pmin); msg_rank = rank)
    end
    return nothing
end

#---------------------------------------------------------------------------------
# apply_positivity!(u, params, SD)
#
# No-op unless :lpositivity. Repairs params.uaux in place, writes it back to the
# integrator state u, and reports on a self-throttling schedule.
#---------------------------------------------------------------------------------
function apply_positivity!(u, params, SD)

    inputs = params.inputs
    get(inputs, :lpositivity, false) || return nothing

    npoin = Int(params.mesh.npoin)
    neqs  = Int(params.neqs)
    nsd   = (SD == NSD_3D()) ? 3 : 2
    ien   = nsd + 2

    if !POSITIVITY_CHECKED[]
        POSITIVITY_CHECKED[] = true
        positivity_validate(inputs, params, neqs, ien)
    end

    T    = eltype(params.uaux)
    γm1  = T(PhysicalConst{Float64}().γm1)
    ρmin = T(inputs[:positivity_rho_min])
    pmin = T(inputs[:positivity_p_min])

    nrep = Positivity.positivity_limit!(@view(params.uaux[:, :]), npoin, ien,
                                        γm1, ρmin, pmin, POSITIVITY_STATS;
                                        coords = params.mesh.coords,   # [dim, ip]
                                        t = NaN)

    # Only write back when something was actually repaired. A case that never
    # needs the repair is then BIT-IDENTICAL to running with :lpositivity off —
    # it pays one read sweep per RHS call and nothing else. That is what makes
    # it safe to leave on in a validated case.
    nrep > 0 && uaux2u!(u, @view(params.uaux[:, :]), neqs, npoin)

    if get(inputs, :positivity_report, true) &&
       Positivity.positivity_should_report(POSITIVITY_STATS)
        if MPI.Comm_rank(get_mpi_comm()) == 0
            println_rank(string(" # POSITIVITY REPAIR ENGAGED — ",
                         Positivity.positivity_summary(POSITIVITY_STATS), "\n",
                         " #   A few node-visits near a shock is the repair doing its job.\n",
                         " #   Engagement growing without bound, or anywhere inside the\n",
                         " #   boundary layer, means the ANSWER is wrong and the repair is\n",
                         " #   only hiding it — check the first-engagement coordinates above\n",
                         " #   against the wall before trusting any heat flux from this run."))
        end
    end

    return nothing
end
