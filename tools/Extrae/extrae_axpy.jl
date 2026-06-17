# ==============================================================================
# extrae_axpy.jl  --  serial Extrae.jl instrumentation example
#
# Reproduces the canonical examples from the Extrae.jl paper
# (arXiv:2504.12087v1, Listings 1 and 2): instrument an `axpy!` kernel with a
# user-code region (@user_function) and emit a custom event recording the
# vector length.
#
# Run (from the Jexpresso project root):
#
#     julia --project=. tools/Extrae/extrae_axpy.jl
#
# On Linux with Extrae installed this writes a Paraver trace into ./set-0/ and
# a `*.prv` after merge (see README).  On macOS it runs identically but the
# instrumentation calls are no-ops (see ExtraeShim.jl for why).
# ==============================================================================

include(joinpath(@__DIR__, "ExtraeShim.jl"))
using .ExtraeShim

# --- event setup (paper, Listing 2) ----------------------------------------
# Custom event TYPE code (extrae_type_t, a UInt32). Pick a value in the
# "user" range that does not collide with Extrae's built-in event types.
const CODE_VEC_LEN = 84210   # "Vector length"
const CODE_DTYPE   = 84211   # which floating type the kernel ran on

# --- instrumented kernel (paper, Listing 1) --------------------------------
function axpy!(a, x, y)
    @user_function begin
        @inbounds @simd for i in eachindex(x, y)
            y[i] = muladd(a, x[i], y[i])
        end
    end
    return y
end

function run_benchmark(::Type{T}, n::Int; reps::Int = 50) where {T}
    a = T(2)
    x = ones(T, n)
    y = ones(T, n)

    # Mark the problem size and the element type as events so Paraver can
    # colour the timeline by them.
    emit(CODE_VEC_LEN, n)
    emit(CODE_DTYPE, _dtype_code(T))

    local out
    for _ in 1:reps
        out = axpy!(a, x, y)
    end
    return out
end

_dtype_code(::Type{Float16}) = 16
_dtype_code(::Type{Float32}) = 32
_dtype_code(::Type{Float64}) = 64

function main()
    init()                                   # begin tracing (paper, Listing 1)
    try
        # Name the event types so the trace is human-readable in Paraver.
        register(CODE_VEC_LEN, "Vector length")
        register(CODE_DTYPE, "Element type (bits)",
                 (16, 32, 64), ("Float16", "Float32", "Float64"))

        n = 1_000_000
        for T in (Float16, Float32, Float64)
            run_benchmark(T, n)
        end
    finally
        finish()                             # flush trace (paper, Listing 1)
    end

    if is_active()
        println("Extrae trace written. Merge it with: mpi2prv -f TRACE.mpits -o axpy.prv")
    else
        println("Ran in no-op shim mode (Extrae native library not available " *
                "on this platform). The instrumented code executed correctly; " *
                "no Paraver trace was produced.")
    end
end

main()
