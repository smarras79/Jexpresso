# ==============================================================================
# ExtraeShim.jl
#
# A thin, *portable* wrapper around the Extrae.jl HPC profiler
# (https://github.com/bsc-quantic/Extrae.jl).
#
# Why a shim?
# -----------
# The Extrae profiler ships as a binary artifact (Extrae_jll) that is built
# ONLY for Linux (x86_64 / aarch64 / powerpc64le, glibc).  There is **no
# macOS build**.  That means the bindings cannot emit real Paraver traces on
# a MacBook Air.  Rather than make the examples crash on a laptop, this shim:
#
#   * loads the real `Extrae` package and uses it when it is installed AND the
#     underlying native library actually initialises (i.e. on Linux/HPC), and
#   * otherwise degrades to no-op stubs so the *exact same instrumented code*
#     still runs to completion on macOS/Windows (it just produces no trace).
#
# This lets you develop and run the instrumented example on your Mac, then run
# the byte-for-byte identical code on a Linux cluster to get a real
# `*.prv` Paraver trace.
#
# Public surface (mirrors Extrae.jl so call sites are identical):
#   init()                      -> begin a tracing session
#   init(Val(:Distributed))     -> Distributed.jl flavour (see paper, Listing 3)
#   finish()                    -> close the session and flush the trace
#   register(tcode, desc)               -> name an event type
#   register(tcode, desc, vcodes, vdescs) -> name an event type + its values
#   emit(tcode, value)          -> emit a punctual event (type, value)
#   user_function(enter::Bool)  -> open/close a user-code "state" region
#   @user_function expr         -> wrap `expr` in a user-code region
#   is_active()                 -> true when a real Extrae trace is being taken
#
# Reference: S. Sanchez-Ramirez & M. Giordano, "Extrae.jl: Julia bindings for
# the Extrae HPC Profiler", Proceedings of JuliaCon (arXiv:2504.12087v1).
# ==============================================================================
module ExtraeShim

export init, finish, register, emit, user_function, @user_function, is_active

# ---------------------------------------------------------------------------
# Try to bring in the real package.  We keep a handle to the module rather
# than `using` it so that, on platforms without the binary, a load failure is
# caught here and turned into "shim mode" instead of an error at `using` time.
# ---------------------------------------------------------------------------
const _EXTRAE = Ref{Union{Module,Nothing}}(nothing)
const _ACTIVE = Ref{Bool}(false)        # true only after a successful init()

function _load_extrae()
    if _EXTRAE[] === nothing
        try
            _EXTRAE[] = Base.require(Base.PkgId(
                Base.UUID("8a0c07fa-ade5-4b2a-b81a-b192b2bedf88"), "Extrae"))
        catch err
            @debug "Extrae package not loadable; using no-op shim" exception = err
            _EXTRAE[] = nothing
        end
    end
    return _EXTRAE[]
end

"""
    is_active() -> Bool

`true` when a real Extrae session is running (Linux + library initialised),
`false` when the no-op shim is in effect (e.g. on macOS).  Use it to gate
output messages, not the instrumentation itself — the instrumentation calls
are always safe.
"""
is_active() = _ACTIVE[]

# Extrae's C API uses `extrae_type_t == unsigned int` (UInt32) for event type
# codes and `extrae_value_t == unsigned long long` (UInt64) for values. We
# accept any Integer at the call site and convert here so the examples can use
# plain `Int` literals.
_tcode(x) = UInt32(x)
_vcode(x) = UInt64(x)

# ---------------------------------------------------------------------------
# init / finish
# ---------------------------------------------------------------------------
"""
    init()
    init(::Val{:Distributed})

Start a tracing session.  In shim mode this is a no-op that simply records
that we are NOT producing a real trace.  The `Val(:Distributed)` method
mirrors Extrae.jl's helper for `Distributed.jl` workers (paper, Listing 3).
"""
function init(arg = nothing)
    m = _load_extrae()
    if m === nothing
        _ACTIVE[] = false
        return false
    end
    try
        arg === nothing ? m.init() : m.init(arg)
        _ACTIVE[] = true
    catch err
        # Package present but the native library cannot initialise on this
        # platform (the macOS case) — fall back to no-op cleanly.
        @debug "Extrae.init() failed; using no-op shim" exception = err
        _ACTIVE[] = false
    end
    return _ACTIVE[]
end

"""
    finish()

Close the session and flush the trace to disk.  No-op in shim mode.
"""
function finish()
    m = _load_extrae()
    if m !== nothing && _ACTIVE[]
        try
            m.finish()
        catch err
            @debug "Extrae.finish() failed" exception = err
        end
    end
    _ACTIVE[] = false
    return nothing
end

# ---------------------------------------------------------------------------
# register / emit  (events)
# ---------------------------------------------------------------------------
"""
    register(tcode, desc)
    register(tcode, desc, vcodes, vdescs)

Give a human-readable name to an event type (and, optionally, to each of the
discrete values that type can take) so Paraver shows labels instead of raw
integers.  No-op in shim mode.
"""
function register(tcode, desc::AbstractString)
    m = _load_extrae()
    (m === nothing || !_ACTIVE[]) && return nothing
    try
        m.register(_tcode(tcode), String(desc))
    catch err
        @debug "Extrae.register failed" exception = err
    end
    return nothing
end

function register(tcode, desc::AbstractString, vcodes, vdescs)
    m = _load_extrae()
    (m === nothing || !_ACTIVE[]) && return nothing
    try
        m.register(_tcode(tcode), String(desc),
                   UInt64[_vcode(v) for v in vcodes],
                   String[String(d) for d in vdescs])
    catch err
        @debug "Extrae.register(values) failed" exception = err
    end
    return nothing
end

"""
    emit(tcode, value)

Emit a punctual `(type, value)` event on the current task/thread at the
current time.  This is the workhorse used to mark phases, iteration numbers,
problem sizes, etc.  No-op in shim mode.
"""
function emit(tcode, value)
    m = _load_extrae()
    (m === nothing || !_ACTIVE[]) && return nothing
    try
        m.emit(_tcode(tcode), _vcode(value))
    catch err
        @debug "Extrae.emit failed" exception = err
    end
    return nothing
end

# ---------------------------------------------------------------------------
# user_function: open/close a "user code" state region
# ---------------------------------------------------------------------------
"""
    user_function(enter::Bool)

Open (`true`) or close (`false`) a user-coded state region in the trace.
No-op in shim mode.
"""
function user_function(enter)
    m = _load_extrae()
    (m === nothing || !_ACTIVE[]) && return nothing
    try
        m.user_function(enter)
    catch err
        @debug "Extrae.user_function failed" exception = err
    end
    return nothing
end

"""
    @user_function expr

Wrap `expr` in a user-code region: opens the region, evaluates `expr`
(returning its value), then closes the region even if `expr` throws.
Mirrors Extrae.jl's `@user_function` (paper, Listing 1).
"""
macro user_function(expr)
    quote
        ExtraeShim.user_function(true)
        try
            $(esc(expr))
        finally
            ExtraeShim.user_function(false)
        end
    end
end

end # module ExtraeShim
