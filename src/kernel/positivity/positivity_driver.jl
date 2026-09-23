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
# SCOPE. Two state layouts, each with its own repair in Positivity.jl:
#
#   * 2D/3D CompEuler, neqs == nsd + 2 exactly  ->  positivity_limit!
#   * 2D ideal GLM-MHD, the nine-field state
#     (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)        ->  positivity_limit_mhd!
#
# both requiring TOTAL(), ρE in the energy slot and the CPU backend. The MHD
# layout is recognised from the case's own qvars, not guessed from neqs: a
# nine-equation system that is not this one must not be repaired as if it were.
# It differs in two structural ways — a non-contiguous momentum slot map and an
# internal energy that owes ½|B|² + ½ψ² to the field — so it is a separate
# function; see its header. The MHD γ comes from :dsgs_gamma, NOT from
# PhysicalConst: the MHD cases here run γ = 1.4, 5/3 and 1.05, and air's 1.4
# would silently be wrong for two of the three.
#
# Anything else is a clear error at the first call rather than a silent wrong
# repair: the θ-form carries ρθ, which is positive for a different reason.
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
#
# Is this the nine-field ideal GLM-MHD state? Decided from the case's OWN
# solution-variable names (initialize.jl's qvars), so it can only be true for
# the state the MHD repair is actually written for. Slot 4 must be the energy
# and 6:9 the field — exactly the ordering every MHD case in problems/MHD uses.
#
function positivity_is_mhd(params, neqs::Int)
    neqs == 9 || return false
    qv = params.qp.qvars
    (qv isa AbstractVector && length(qv) >= 9) || return false
    nm = String.(string.(qv[1:9]))
    return nm == ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
end

function positivity_validate(inputs, params, neqs::Int, ien::Int)

    why = String[]
    lmhd = positivity_is_mhd(params, neqs)

    get(inputs, :energy_equation, "theta") == "energy" ||
        push!(why, "  :energy_equation must be \"energy\" (slot $(ien) must hold ρE, not ρθ)")

    params.SOL_VARS_TYPE == TOTAL() ||
        push!(why, "  :SOL_VARS_TYPE must be TOTAL() (a perturbation state has no realizable set of its own)")

    if lmhd
        # The MHD repair hard-codes the slot map of the problems/MHD state.
        haskey(inputs, :dsgs_gamma) ||
            push!(why, "  the GLM-MHD repair needs :dsgs_gamma, the γ of the MHD equation of state\n" *
                       "     (PhysicalConst's 1.4 is air's and is wrong for an MHD case that is not γ = 1.4)")
    else
        neqs == ien ||
            push!(why, "  neqs = $(neqs) but this repair is written for exactly nsd+2 = $(ien) equations,\n" *
                       "     or for the nine-field GLM-MHD state (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ) —\n" *
                       "     and this case's qvars are not that either")
    end

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
        println_rank(@sprintf(" # POSITIVITY REPAIR ON (%s state): ρ_min = %.3e, p_min = %.3e (absolute), γ = %.6g. This is a REPAIR, not a preserving scheme — every intervention is counted and reported.",
                       lmhd ? "9-field GLM-MHD" : "CompEuler", ρmin, pmin,
                       lmhd ? Float64(inputs[:dsgs_gamma]) : PhysicalConst{Float64}().γ);
                     msg_rank = rank)
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
    lmhd = positivity_is_mhd(params, neqs)
    # The MHD equation of state has its own γ (1.4 on the jet, 5/3 on
    # Orszag-Tang, 1.05 on flux emergence); PhysicalConst's is air's.
    γm1  = lmhd ? T(Float64(inputs[:dsgs_gamma]) - 1.0) : T(PhysicalConst{Float64}().γm1)
    ρmin = T(inputs[:positivity_rho_min])
    pmin = T(inputs[:positivity_p_min])

    nrep = if lmhd
        # Slot map of the problems/MHD state: energy in 4, momentum in (2,3,5),
        # B in 6:8, ψ in 9. See the header of positivity_limit_mhd!.
        Positivity.positivity_limit_mhd!(@view(params.uaux[:, :]), npoin,
                                         γm1, ρmin, pmin, POSITIVITY_STATS;
                                         coords = params.mesh.coords,   # [dim, ip]
                                         t = NaN)
    else
        Positivity.positivity_limit!(@view(params.uaux[:, :]), npoin, ien,
                                     γm1, ρmin, pmin, POSITIVITY_STATS;
                                     coords = params.mesh.coords,   # [dim, ip]
                                     t = NaN)
    end

    # Only write back when something was actually repaired. A case that never
    # needs the repair is then BIT-IDENTICAL to running with :lpositivity off —
    # it pays one read sweep per RHS call and nothing else. That is what makes
    # it safe to leave on in a validated case.
    nrep > 0 && uaux2u!(u, @view(params.uaux[:, :]), neqs, npoin)

    # ---- REPORT. COLLECTIVE, and it has to be. --------------------------------
    #
    # POSITIVITY_STATS is rank-local. On 256 ranks, printing rank 0's copy
    # describes 1/256 of the domain — which is how the first Mach-7 report of
    # this feature came back saying "1 node-visit" before anyone had asked the
    # other 255 ranks. Counts are SUM-reduced and minima MIN-reduced before
    # anything is said.
    #
    # The trigger is `ncalls % every == 0`, identical on every rank because
    # every rank makes the same number of RHS calls, so the reduction is
    # entered in lockstep and cannot deadlock. Triggering on rank 0's own
    # engagement count would hang the run the moment the ranks disagreed.
    every = Int(get(inputs, :positivity_report_every, 1000))
    if get(inputs, :positivity_report, true) && every > 0 &&
       POSITIVITY_STATS.ncalls % every == 0

        comm = get_mpi_comm()
        sums = MPI.Allreduce([Float64(POSITIVITY_STATS.nrho),
                              Float64(POSITIVITY_STATS.nmom),
                              Float64(POSITIVITY_STATS.nenergy),
                              POSITIVITY_STATS.dmass,
                              POSITIVITY_STATS.denergy], MPI.SUM, comm)
        mins = MPI.Allreduce([POSITIVITY_STATS.rho_min,
                              POSITIVITY_STATS.p_min], MPI.MIN, comm)

        #
        # WHERE THE FIRST REPAIR HAPPENED, GLOBALLY.
        #
        # This used to print POSITIVITY_STATS.first_x/y from rank 0 and label
        # it "first on THIS rank". It cost two wrong diagnoses on the Mach-7.7
        # ramp: rank 0's partition holds no near-wall nodes, so it reported the
        # outflow plane while a 1-rank run of the same case put the true first
        # repair at the leading edge, 0.7 mm from the origin. A rank-local
        # coordinate that READS like a global one is worse than no coordinate.
        #
        # So: reduce on the RHS-call index of each rank's first intervention,
        # find the earliest, let the lowest-numbered rank holding it win the
        # tie, and broadcast its coordinate. Every rank reaches all three
        # collectives because the trigger above is the call count, which is
        # identical everywhere.
        #
        myrank    = MPI.Comm_rank(comm)
        mycall    = Float64(POSITIVITY_STATS.first_call)      # typemax if never
        firstcall = MPI.Allreduce(mycall, MPI.MIN, comm)
        owner     = MPI.Allreduce(mycall == firstcall ? Float64(myrank) : Inf,
                                  MPI.MIN, comm)
        where_    = [POSITIVITY_STATS.first_x, POSITIVITY_STATS.first_y,
                     POSITIVITY_STATS.first_t]
        if isfinite(owner)
            MPI.Bcast!(where_, Int(owner), comm)
        end

        total  = sums[1] + sums[2] + sums[3]
        decade = total > 0.0 ? floor(Int, log10(total)) + 1 : 0
        if decade > POSITIVITY_STATS.nreported
            POSITIVITY_STATS.nreported = decade     # global, so identical everywhere
            println_rank(string(" # POSITIVITY REPAIR ENGAGED — ",
                                Positivity.positivity_summary(sums[1], sums[2], sums[3],
                                                              sums[4], sums[5],
                                                              mins[1], mins[2],
                                                              POSITIVITY_STATS.ncalls), "\n",
                                " #   GLOBAL first repair at (x, y) = (",
                                where_[1], ", ", where_[2],
                                ")  on RHS call ", Int(firstcall),
                                ", rank ", Int(owner), "\n",
                                " #   A few node-visits near a shock is the repair doing its job.\n",
                                " #   Engagement growing without bound, or anywhere inside the\n",
                                " #   boundary layer, means the ANSWER is wrong and the repair is\n",
                                " #   only hiding it — check that coordinate against the wall\n",
                                " #   before trusting any heat flux from this run.");
                         msg_rank = MPI.Comm_rank(comm))
        end
    end

    return nothing
end
