#--------------------------------------------------------
# Extrae tracing integration (opt-in).
#
# Wraps https://github.com/bsc-quantic/Extrae.jl so that Jexpresso can
# emit a Paraver trace of an MPI run without making Extrae a hard
# dependency: when the env var JEXPRESSO_EXTRAE is unset (or Extrae.jl
# is not installed in the active project) every public symbol here is a
# no-op and the code path is byte-for-byte the historical one.
#
# Typical usage at the BSC clusters:
#
#     export JEXPRESSO_EXTRAE=1
#     export EXTRAE_CONFIG_FILE=/path/to/extrae.xml
#     export LD_PRELOAD=$EXTRAE_HOME/lib/libmpitrace.so
#     mpirun -np N julia --project=. ./src/run.jl ...
#
# The libmpitrace.so preload intercepts MPI_Init/MPI_Finalize, so Extrae
# is already tracing MPI calls before Julia gets a chance to run.  The
# extra `Extrae.init` call below is documented as safe (idempotent) and
# is kept so the trace also covers runs where the LD_PRELOAD flavour is
# `libseqtrace.so` (no MPI auto-init).
#
# The user-event palette registered below mirrors the high-level phases
# already named in the code (driver, time loop, RHS, ...) so the
# resulting Paraver trace is readable without manual relabelling.
#--------------------------------------------------------

const _JEXPRESSO_EXTRAE_ENV   = lowercase(get(ENV, "JEXPRESSO_EXTRAE", ""))
const JEXPRESSO_EXTRAE_ENABLED = _JEXPRESSO_EXTRAE_ENV in ("1", "true", "yes", "on")

# Try to load Extrae.jl only when the user actually asked for tracing.
# A failed import is non-fatal: we fall back to no-ops and emit a single
# warning from rank 0 so the run keeps going.
const EXTRAE_LOADED = Ref(false)
if JEXPRESSO_EXTRAE_ENABLED
    try
        @eval using Extrae
        EXTRAE_LOADED[] = true
    catch err
        @warn "JEXPRESSO_EXTRAE is set but `using Extrae` failed; tracing disabled" exception=(err, catch_backtrace())
    end
end

#--------------------------------------------------------
# Event-type codes.  Extrae groups events by `type` (an Int32 channel)
# and `value` (the payload).  Using one type per phase keeps each phase
# on its own Paraver row; value 0 means "leave region", non-zero means
# "enter region" (matches Extrae's user_function convention).
#--------------------------------------------------------
const JE_EVT_DRIVER     = Int32(91_000)
const JE_EVT_MESH       = Int32(91_001)
const JE_EVT_SEM_SETUP  = Int32(91_002)
const JE_EVT_TIME_LOOP  = Int32(91_003)
const JE_EVT_RHS        = Int32(91_004)
const JE_EVT_BC         = Int32(91_005)
const JE_EVT_IO         = Int32(91_006)

const _JE_EVT_NAMES = [
    (JE_EVT_DRIVER,    "Jexpresso: driver"),
    (JE_EVT_MESH,      "Jexpresso: mesh"),
    (JE_EVT_SEM_SETUP, "Jexpresso: SEM setup"),
    (JE_EVT_TIME_LOOP, "Jexpresso: time loop"),
    (JE_EVT_RHS,       "Jexpresso: RHS"),
    (JE_EVT_BC,        "Jexpresso: boundary conditions"),
    (JE_EVT_IO,        "Jexpresso: I/O"),
]

#--------------------------------------------------------
# Public helpers.  Each guards on EXTRAE_LOADED[] so they are cheap
# (single Ref dereference + branch) when tracing is off.
#--------------------------------------------------------
@inline function je_extrae_init()
    EXTRAE_LOADED[] || return nothing
    Base.invokelatest(Extrae.init)
    return nothing
end

@inline function je_extrae_finish()
    EXTRAE_LOADED[] || return nothing
    Base.invokelatest(Extrae.finish)
    return nothing
end

function je_extrae_register_events()
    EXTRAE_LOADED[] || return nothing
    for (code, name) in _JE_EVT_NAMES
        Base.invokelatest(Extrae.register, code, name)
    end
    return nothing
end

@inline function je_trace_enter(evt::Int32)
    EXTRAE_LOADED[] || return nothing
    Base.invokelatest(Extrae.emit, evt, Int64(1))
    return nothing
end

@inline function je_trace_leave(evt::Int32)
    EXTRAE_LOADED[] || return nothing
    Base.invokelatest(Extrae.emit, evt, Int64(0))
    return nothing
end

@inline function je_trace_event(evt::Int32, value::Integer)
    EXTRAE_LOADED[] || return nothing
    Base.invokelatest(Extrae.emit, evt, Int64(value))
    return nothing
end

"""
    @je_trace evt expr

Surround `expr` with Extrae enter/leave events on channel `evt` (an
`Int32` constant from this file).  A no-op when tracing is disabled.
The leave event is emitted from a `finally` block so the trace stays
consistent across exceptions.
"""
macro je_trace(evt, expr)
    quote
        local _je_evt = $(esc(evt))
        je_trace_enter(_je_evt)
        try
            $(esc(expr))
        finally
            je_trace_leave(_je_evt)
        end
    end
end
