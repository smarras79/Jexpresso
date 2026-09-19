#---------------------------------------------------------------------------------
# Positivity.jl — node-wise REALIZABILITY REPAIR for the compressible Euler state.
#
# SELF-CONTAINED ON PURPOSE. This module knows nothing about Jexpresso: no
# types, no MPI, no mesh, no inputs. It takes a (npoin × neqs) matrix of
# conserved variables, the slot index of the energy equation, γ-1 and two
# floors, and repairs the state in place. Everything that has to know about
# Jexpresso lives in positivity_driver.jl next door. That split is deliberate —
# the numerics here can be read, reasoned about and unit-tested without loading
# the solver.
#
# WHAT THIS IS, AND WHAT IT IS NOT
# --------------------------------
# It is NOT a positivity-PRESERVING scheme. A preserving scheme carries a proof
# that the update cannot leave the realizable set under a stated CFL condition.
# For a continuous-Galerkin SEM the theorem that does that is invariant-domain
# preservation with convex limiting (Guermond, Popov & Tomas; the `ryujin` code
# of Maier & Kronbichler): a low-order update with graph viscosity along every
# node-neighbour edge, which is provably positive, then a node-by-node convex
# combination with the high-order update, limited to enforce ρ > 0, internal
# energy > 0 and a local minimum-entropy principle. That is the real answer and
# it is a project, not a file. See README.md.
#
# This is the REPAIR that sits under it until then: the bound enforcement
# without the invariant-domain guarantee. It cannot promise the solution is
# right. It promises three things, which are worth having on their own:
#
#   1. the RHS is never evaluated on a non-realizable state, so a single bad
#      node stops producing NaN fluxes, NaN sound speeds and log(p<0);
#   2. the repair is LOCAL and BOUNDED, so a defect cannot spread by arithmetic;
#   3. it is AUDITED — every intervention is counted and reported, including
#      where the first one happened, so it can never quietly rescue a run whose
#      answer is meaningless.
#
# (3) is the point. The open question on the Mach-7 cases is whether the failure
# is a handful of nodes at the bow shock that then poison the field, or a field
# that is globally garbage. If this repair engages at ten nodes and the run
# continues sensibly, the first is true and convex limiting is worth building.
# If it engages at ten thousand, the second is true and no limiter would have
# helped. Either answer is cheaper to buy here than after writing the limiter.
#
# THE REPAIR
# ----------
# Per node, in order:
#
#   1. ρ < ρ_min          ->  ρ = ρ_min.            Injects mass; recorded.
#
#   2. p = (γ-1)(ρE - ke) < p_min, with ke = |ρu|²/(2ρ):
#
#      a. if ρE > e_min = p_min/(γ-1) and ke > 0, scale the momentum by
#
#             θ = sqrt( (ρE - e_min) / ke )   ∈ [0,1)
#
#         which makes p = p_min exactly. ρE IS UNTOUCHED, so TOTAL ENERGY IS
#         CONSERVED EXACTLY: the repair converts kinetic energy into internal
#         energy and nothing else. That is what a viscous term does, so it is
#         dissipative and entropy-increasing — the right sign for a repair.
#         Nothing is created; a wrong answer here is a locally over-damped one,
#         not a locally energised one.
#
#      b. only if ρE ≤ e_min (momentum scaling cannot reach p_min because the
#         TOTAL energy is already too small): zero the momentum and raise ρE to
#         e_min. This one DOES inject energy, and it is recorded separately for
#         exactly that reason. If this branch is firing, the state is badly
#         broken, not marginally so.
#
# NaN IS LEFT ALONE, DELIBERATELY. `NaN < ρ_min` is false, so a NaN node passes
# through untouched and the solver's own non-finite check still fires. Repairing
# NaN would turn a detectable failure into a silent one; this module exists to
# keep a FINITE state realizable, not to resurrect a dead one.
#---------------------------------------------------------------------------------
module Positivity

export PositivityStats, positivity_limit!, positivity_reset!,
       positivity_touched, positivity_should_report, positivity_summary

#---------------------------------------------------------------------------------
# The audit trail. One instance per run, owned by the driver.
#---------------------------------------------------------------------------------
mutable struct PositivityStats
    ncalls  ::Int        # RHS evaluations the limiter has seen
    nrho    ::Int        # nodes whose density was floored
    nmom    ::Int        # nodes whose momentum was scaled (ρE conserved)
    nenergy ::Int        # nodes whose energy had to be raised (energy INJECTED)
    dmass   ::Float64    # total mass injected by the ρ floor
    denergy ::Float64    # total energy injected by branch 2b
    rho_min ::Float64    # smallest FINITE ρ seen, before repair
    p_min   ::Float64    # smallest FINITE p seen, before repair
    first_x ::Float64    # where the first intervention happened
    first_y ::Float64
    first_t ::Float64
    first_call::Int      # the RHS call the first intervention happened on, so
                         # the report can pick the GLOBALLY earliest one across
                         # ranks instead of announcing rank 0's local view
    nreported::Int       # how many powers of ten have been announced
end

PositivityStats() = PositivityStats(0, 0, 0, 0, 0.0, 0.0,
                                    Inf, Inf, NaN, NaN, NaN, typemax(Int), 0)

function positivity_reset!(s::PositivityStats)
    s.ncalls = 0; s.nrho = 0; s.nmom = 0; s.nenergy = 0
    s.dmass = 0.0; s.denergy = 0.0
    s.rho_min = Inf; s.p_min = Inf
    s.first_x = NaN; s.first_y = NaN; s.first_t = NaN
    s.first_call = typemax(Int)
    s.nreported = 0
    return s
end

positivity_touched(s::PositivityStats) = s.nrho + s.nmom + s.nenergy

#---------------------------------------------------------------------------------
# Kinetic energy at one node: ke = |ρu|²/(2ρ). Momentum occupies slots
# 2 : ien-1, whatever the dimension, because the energy slot is ien = nsd + 2.
#---------------------------------------------------------------------------------
@inline function _kinetic(uaux::AbstractMatrix{T}, ip::Integer, ien::Integer,
                          ρ::T) where {T<:AbstractFloat}
    acc = zero(T)
    @inbounds for k = 2:(ien - 1)
        m = uaux[ip, k]
        acc += m*m
    end
    return T(0.5)*acc/ρ
end

#---------------------------------------------------------------------------------
# Record where the first intervention of the run happened. Coordinates are
# optional: pass empty vectors and it stores NaN.
#---------------------------------------------------------------------------------
# COORDINATES ARE coords[dim, ip], the (3 x npoin) array, NOT the deprecated
# per-axis fields. Indexing matches the rest of the kernel, e.g. rhs.jl:1131
# `x=coords[1, ip], y=coords[2, ip]`. Size-guarded so the module still works
# when handed nothing: coordinates only ever say WHERE the first repair
# happened, and a missing coordinate must never cost a repair.
@inline function _mark_first!(s::PositivityStats, ip::Integer,
                              coords::AbstractArray, t::Real)
    if positivity_touched(s) == 0
        ok = (ndims(coords) == 2 && size(coords, 1) >= 2 && size(coords, 2) >= ip)
        s.first_x = ok ? Float64(coords[1, ip]) : NaN
        s.first_y = ok ? Float64(coords[2, ip]) : NaN
        s.first_t = Float64(t)
        s.first_call = s.ncalls
    end
    return nothing
end

#---------------------------------------------------------------------------------
# positivity_limit!(uaux, npoin, ien, γm1, ρmin, pmin, stats; coords, t)
#
# Repairs `uaux` in place. `ien` is the energy slot (nsd + 2); slots 2:ien-1 are
# momentum.
#
# RETURNS THE NUMBER OF REPAIRS MADE IN THIS CALL, which matters: the caller
# skips the write-back to the integrator state when it is zero, so enabling the
# repair on a case that never needs it is BIT-IDENTICAL to leaving it off, not
# merely equivalent. A healthy case pays one read sweep and nothing else.
#---------------------------------------------------------------------------------
function positivity_limit!(uaux::AbstractMatrix{T},
                           npoin::Integer, ien::Integer,
                           γm1::T, ρmin::T, pmin::T,
                           s::PositivityStats;
                           coords::AbstractArray = zeros(T, 0, 0),
                           t::Real = NaN) where {T<:AbstractFloat}

    s.ncalls += 1
    nrep = 0
    emin = pmin/γm1

    @inbounds for ip = 1:npoin

        ρ = uaux[ip, 1]
        if isfinite(ρ) && Float64(ρ) < s.rho_min
            s.rho_min = Float64(ρ)
        end

        # ---- 1. density floor -------------------------------------------------
        if ρ < ρmin                       # false for NaN: left alone on purpose
            _mark_first!(s, ip, coords, t)
            s.dmass += Float64(ρmin - ρ)
            s.nrho  += 1
            nrep   += 1
            ρ = ρmin
            uaux[ip, 1] = ρ
        end

        # ---- 2. pressure / internal energy ------------------------------------
        ke = _kinetic(uaux, ip, ien, ρ)
        ρE = uaux[ip, ien]
        e  = ρE - ke                      # internal energy per unit volume
        if isfinite(e) && Float64(γm1*e) < s.p_min
            s.p_min = Float64(γm1*e)
        end

        if e < emin                       # false for NaN: left alone on purpose

            if ρE > emin && ke > zero(T)
                # 2a. scale the momentum. ρE untouched -> total energy conserved.
                θ = sqrt(max((ρE - emin)/ke, zero(T)))
                _mark_first!(s, ip, coords, t)
                for k = 2:(ien - 1)
                    uaux[ip, k] *= θ
                end
                s.nmom += 1
                nrep   += 1
            else
                # 2b. the total energy itself is too small. Inject, and say so.
                _mark_first!(s, ip, coords, t)
                for k = 2:(ien - 1)
                    uaux[ip, k] = zero(T)
                end
                s.denergy += Float64(emin - ρE)
                uaux[ip, ien] = emin
                s.nenergy += 1
                nrep      += 1
            end
        end
    end

    return nrep
end

#---------------------------------------------------------------------------------
# Self-throttling report: announce the first intervention, then once per decade
# of the cumulative count. No step counter needed, and it cannot flood a log.
#---------------------------------------------------------------------------------
function positivity_should_report(s::PositivityStats)
    n = positivity_touched(s)
    n == 0 && return false
    decade = floor(Int, log10(n)) + 1
    if decade > s.nreported
        s.nreported = decade
        return true
    end
    return false
end

#---------------------------------------------------------------------------------
# Format from EXPLICIT numbers, so the caller can hand in MPI-reduced totals.
# The struct-taking method below is rank-local and is only good for a serial
# run — on many ranks it describes 1/nranks of the domain, which is how the
# first report of this feature came back saying "1 node-visit" when nobody had
# yet asked the other 255 ranks.
#---------------------------------------------------------------------------------
function positivity_summary(nrho::Real, nmom::Real, nen::Real,
                            dmass::Real, den::Real,
                            rmin::Real, pmin::Real, ncalls::Real)
    return string("repaired ", Int(nrho + nmom + nen), " node-visits in ",
                  Int(ncalls), " RHS calls",
                  "  [ρ-floor ", Int(nrho),
                  ", momentum-scaled ", Int(nmom),
                  ", energy-RAISED ", Int(nen), "]",
                  "  injected: mass ", dmass, ", energy ", den,
                  "  |  GLOBAL min ρ ", rmin, ", GLOBAL min p ", pmin)
end

function positivity_summary(s::PositivityStats)
    return string("repaired ", positivity_touched(s), " node-visits in ",
                  s.ncalls, " RHS calls",
                  "  [ρ-floor ", s.nrho,
                  ", momentum-scaled ", s.nmom,
                  ", energy-RAISED ", s.nenergy, "]",
                  "  injected: mass ", s.dmass, ", energy ", s.denergy,
                  "  |  min ρ seen ", s.rho_min, ", min p seen ", s.p_min,
                  "  |  first at (x, y, t) = (", s.first_x, ", ",
                  s.first_y, ", ", s.first_t, ")")
end

end # module Positivity
