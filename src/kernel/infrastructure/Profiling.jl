# ==============================================================================
# Profiling.jl  --  Jexpresso's optional Extrae (Paraver) instrumentation
#
# A thin, OPT-IN wrapper around the Extrae.jl HPC profiler
# (https://github.com/bsc-quantic/Extrae.jl). It mirrors the standalone shim
# in tools/Extrae/ExtraeShim.jl but is wired into the Jexpresso package and is
# gated by an environment variable so that NORMAL runs are completely
# unaffected:
#
#   * Tracing is OFF unless  JEXPRESSO_EXTRAE  is set to 1/true/yes/on.
#   * When OFF, init() returns immediately WITHOUT loading Extrae, and every
#     other entry point is a single `Ref{Bool}` read that returns `nothing`.
#   * When ON, it loads the real Extrae package (Linux-only binary) and, if
#     that succeeds, emits real Paraver events; if Extrae is not installed or
#     cannot initialise (e.g. on macOS) it degrades silently to no-ops.
#
# So the instrumentation calls sprinkled through the solver are always safe to
# leave in place: they cost essentially nothing on an ordinary run.
#
# To capture a trace, see tools/Extrae/README.md — the same LD_PRELOAD /
# EXTRAE_CONFIG_FILE setup applies, plus `export JEXPRESSO_EXTRAE=1`.
#
# Reference: S. Sanchez-Ramirez & M. Giordano, "Extrae.jl: Julia bindings for
# the Extrae HPC Profiler", Proceedings of JuliaCon (arXiv:2504.12087v1).
# ==============================================================================
module Profiling

# Extrae's General-registry UUID (bsc-quantic/Extrae.jl).
const _EXTRAE_UUID = "8a0c07fa-ade5-4b2a-b81a-b192b2bedf88"

const _EXTRAE = Ref{Union{Module,Nothing}}(nothing)
const _ACTIVE = Ref{Bool}(false)        # true only after a successful init()

# --- event type codes (extrae_type_t / UInt32) ------------------------------
# A single "Jexpresso phase" event type whose VALUE says which coarse phase of
# the run we are in. 6_700_000 is well clear of Extrae's built-in event types.
const EV_PHASE = UInt32(6_700_000)

# phase value codes (extrae_value_t / UInt64)
const PHASE_NONE      = 0     # between phases / idle
const PHASE_SEM_SETUP = 1     # sem_setup (mesh + SEM operators)
const PHASE_INIT      = 2     # initialize() (initial condition)
const PHASE_PARAMS    = 3     # params_setup
const PHASE_TIMELOOP  = 4     # time_loop! (the main workload)
const PHASE_RHS       = 5     # (reserved for a later step) RHS evaluation
const PHASE_HALO      = 6     # (reserved for a later step) MPI halo exchange
# coupling (Jexpresso <-> Alya) phases
const PHASE_CPL_SETUP  = 7    # one-time coupling handshake / data receive
const PHASE_CPL_INTERP = 8    # per-step interpolation of the solution to Alya pts
const PHASE_CPL_COMM   = 9    # per-step MPI send of the packed data to Alya

"""
    Profiling.enabled() -> Bool

Whether the user opted into tracing via the `JEXPRESSO_EXTRAE` environment
variable. Read once-ish; cheap.
"""
enabled() = lowercase(get(ENV, "JEXPRESSO_EXTRAE", "")) in ("1", "true", "yes", "on")

# Flushed per-rank diagnostic print to stderr, used to pinpoint where a traced
# run might stall during Extrae start-up. Active only while tracing is opted in.
function _dbg(rank, msg)
    enabled() || return nothing
    println(stderr, "[extrae rank=$rank] $msg")
    flush(stderr)
    return nothing
end

"""
    Profiling.is_active() -> Bool

`true` only when a real Extrae session is running (opted in AND the native
library initialised). Use it to gate prints, never the instrumentation calls
themselves (those are always safe).
"""
is_active() = _ACTIVE[]

_tcode(x) = UInt32(x)
_vcode(x) = UInt64(x)

function _load_extrae()
    if _EXTRAE[] === nothing
        try
            _EXTRAE[] = Base.require(Base.PkgId(Base.UUID(_EXTRAE_UUID), "Extrae"))
        catch err
            @debug "Extrae not loadable; Jexpresso tracing stays off" exception = err
            _EXTRAE[] = nothing
        end
    end
    return _EXTRAE[]
end

"""
    Profiling.init(rank = 0) -> Bool

Begin a tracing session, but ONLY if `JEXPRESSO_EXTRAE` is set. Registers the
Jexpresso phase labels so Paraver shows names instead of integers. Returns
whether tracing is now active. No-op (returns false) on a normal run.
"""
function init(rank::Integer = 0)
    if !enabled()
        _ACTIVE[] = false
        return false
    end
    _dbg(rank, "init: loading Extrae package (Base.require) ...")
    m = _load_extrae()
    if m === nothing
        _dbg(rank, "init: Extrae package NOT loadable -> tracing off")
        _ACTIVE[] = false
        return false
    end
    _dbg(rank, "init: Extrae package loaded")
    try
        # IMPORTANT: when Extrae is loaded via LD_PRELOAD (the MPI-tracing
        # workflow), the library auto-initialises inside the MPI_Init
        # interception — by the time we get here it is ALREADY initialised.
        # Calling Extrae_init() a second time under an active MPI session can
        # deadlock (it was the cause of a hang right after the Extrae banner).
        # So only initialise if it has not been initialised yet (e.g. a rare
        # run with JEXPRESSO_EXTRAE set but no preload).
        already = false
        try
            _dbg(rank, "init: calling Extrae.isinit() ...")
            already = m.isinit() != 0
            _dbg(rank, "init: isinit -> $already")
        catch
            already = false
            _dbg(rank, "init: isinit threw -> assuming not initialised")
        end
        if !already
            _dbg(rank, "init: calling Extrae.init() ...")
            m.init()
            _dbg(rank, "init: Extrae.init() returned")
        end
        _ACTIVE[] = true
        # Name the phase event + its values for the Paraver timeline.
        _dbg(rank, "init: registering phase event ...")
        m.register(EV_PHASE, "Jexpresso phase",
                   UInt64[PHASE_NONE, PHASE_SEM_SETUP, PHASE_INIT,
                          PHASE_PARAMS, PHASE_TIMELOOP, PHASE_RHS, PHASE_HALO,
                          PHASE_CPL_SETUP, PHASE_CPL_INTERP, PHASE_CPL_COMM],
                   String["idle", "sem_setup", "initialize", "params_setup",
                          "time_loop", "rhs", "halo_exchange",
                          "coupling_setup", "coupling_interp", "coupling_comm"])
        _dbg(rank, "init: register done -> tracing ACTIVE")
        if rank == 0
            @info "Jexpresso: Extrae tracing ACTIVE (JEXPRESSO_EXTRAE set; already_initialised=$already)."
        end
    catch err
        @debug "Extrae.init() failed; Jexpresso tracing stays off" exception = err
        _ACTIVE[] = false
    end
    return _ACTIVE[]
end

"""
    Profiling.finish()

Close the tracing session and flush the Paraver trace. No-op when inactive.
"""
function finish()
    if _ACTIVE[]
        m = _EXTRAE[]
        if m !== nothing
            try
                m.finish()
            catch err
                @debug "Extrae.finish() failed" exception = err
            end
        end
    end
    _ACTIVE[] = false
    return nothing
end

"""
    Profiling.emit(tcode, value)

Emit a punctual `(type, value)` event. No-op when inactive.
"""
function emit(tcode, value)
    _ACTIVE[] || return nothing
    m = _EXTRAE[]
    m === nothing && return nothing
    try
        m.emit(_tcode(tcode), _vcode(value))
    catch err
        @debug "Extrae.emit failed" exception = err
    end
    return nothing
end

"""
    Profiling.user_function(enter::Bool)

Open (`true`) / close (`false`) a user-coded state region. No-op when inactive.
"""
function user_function(enter)
    _ACTIVE[] || return nothing
    m = _EXTRAE[]
    m === nothing && return nothing
    try
        m.user_function(enter)
    catch err
        @debug "Extrae.user_function failed" exception = err
    end
    return nothing
end

"""
    Profiling.region(f, phase)

Run `f()` inside a traced region: emit `phase`, open a user-code region, run
`f`, then close the region and return to `PHASE_NONE` — even if `f` throws.
Returns `f()`'s value. When tracing is inactive this is just `f()` plus two
`Ref{Bool}` reads, so it is safe to wrap hot-ish call sites.

    solution = Profiling.region(Profiling.PHASE_TIMELOOP) do
        time_loop!(...)
    end
"""
function region(f, phase::Integer)
    _ACTIVE[] || return f()
    emit(EV_PHASE, phase)
    user_function(true)
    try
        return f()
    finally
        user_function(false)
        emit(EV_PHASE, PHASE_NONE)
    end
end

"""
    Profiling.region_begin(phase)
    Profiling.region_end()

Open / close a traced region WITHOUT a closure — for bracketing a span inline
when wrapping it in a `do` block would change variable scope (e.g. inside the
coupling exchange). Always pair them, and do not nest one pair inside another
(use sequential, non-overlapping pairs). No-op when inactive.
"""
function region_begin(phase::Integer)
    _ACTIVE[] || return nothing
    emit(EV_PHASE, phase)
    user_function(true)
    return nothing
end

function region_end()
    _ACTIVE[] || return nothing
    user_function(false)
    emit(EV_PHASE, PHASE_NONE)
    return nothing
end

end # module Profiling
