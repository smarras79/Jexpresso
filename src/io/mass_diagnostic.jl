#
# Mass-conservation writer (opt-in, :lmass_diagnostic => true).
#
# Every :mass_diagnostic_interval accepted steps, rank 0 appends one line
#
#     step,t_s,t_days,mass,rel_change
#
# to <output_dir>/mass.csv, where mass = Σᵢ Mᵢ q₁ᵢ is the integral of the first
# conserved variable with the lumped mass matrix (ρ for the Euler systems, the
# layer volume ∫H dA in m³ for shallow water) and rel_change = (mass − mass₀)/mass₀,
# the relative change against the initial (or restart) state. Every line is flushed, so the file
# can be followed while the run goes: tools/plot_mass_live.py redraws it live.
#
# One in-place Allreduce! of one number per write, over owned nodes only
# (shared interface nodes carry the full assembled value on every rank), and
# the line is formatted into a preallocated buffer: a write allocates nothing
# (test/test_mass_diagnostic.jl). CPU, continuous Galerkin, lumped mass.
#

mutable struct MassDiagnostic
    active::Bool
    io::Union{IOStream, Nothing}
    dir::String                    # output directory the open file lives in
    mass0::Float64
    step::Int
    interval::Int
    lpert::Bool                    # state holds perturbations: add qe[:, 1]
    red::Vector{Float64}           # Allreduce buffer
    buf::Vector{UInt8}             # line buffer
end

MassDiagnostic(; interval::Int = 1, lpert::Bool = false) =
    MassDiagnostic(false, nothing, "", NaN, 0, max(1, interval), lpert,
                   zeros(Float64, 1), Vector{UInt8}(undef, 256))

# every field has a bounded width, so a line always fits the 256-byte buffer
const _MASS_LINE = Printf.Format("%d,%.10e,%.10g,%.16e,%+.6e\n")

function mass_diagnostic_enabled(inputs)
    return get(inputs, :lmass_diagnostic, false) == true
end

# Σᵢ∈owned Mᵢ q₁ᵢ on this rank (function barrier: concrete argument types).
function _local_mass(u::AbstractVector, M::AbstractVector, owner::AbstractVector,
                     qe::AbstractMatrix, npoin::Int, rank::Int, lpert::Bool)
    s = 0.0
    @inbounds for ip in 1:npoin
        if owner[ip] == rank
            q1 = lpert ? u[ip] + qe[ip, 1] : u[ip]
            s += M[ip] * q1
        end
    end
    return s
end

# Σᵢ∈owned Mᵢ on this rank: the domain area (2D) or volume (3D).
function _local_area(M::AbstractVector, owner::AbstractVector, npoin::Int, rank::Int)
    s = 0.0
    @inbounds for ip in 1:npoin
        owner[ip] == rank && (s += M[ip])
    end
    return s
end

# Name of the first conserved variable (e.g. "H", "ρ"), for the CSV header.
function _first_variable_name(params)
    qp = params.qp
    if hasproperty(qp, :qvars)
        v = getproperty(qp, :qvars)
        if v !== nothing && length(v) >= 1 && v[1] !== nothing
            return string(v[1])
        end
    end
    return "q1"
end

function global_mass(d::MassDiagnostic, u, params, comm)
    rank = MPI.Comm_rank(comm)
    d.red[1] = _local_mass(u, params.M, params.mesh.gip2owner, params.qp.qe,
                           Int(params.mesh.npoin), rank, d.lpert)
    MPI.Allreduce!(d.red, MPI.SUM, comm)
    return d.red[1]
end

function _write_mass_line!(d::MassDiagnostic, step::Int, t::Float64, mass::Float64)
    io = d.io
    io === nothing && return nothing
    rel = mass / d.mass0 - 1.0
    n = Printf.format(d.buf, 1, _MASS_LINE, step, t, t / 86400.0, mass, rel) - 1
    buf = d.buf
    GC.@preserve buf unsafe_write(io, pointer(buf), n)
    flush(io)
    return nothing
end

"""
    open_mass_diagnostic!(d, u, params, inputs, t0)

Activate the writer right before the real solve (after the integrator warm-up,
which runs one step with the output directory swapped to a temporary one):
create `<output_dir>/mass.csv`, take mass₀ from the state `u` at `t0` and write
it as the first line. Collective.
"""
function open_mass_diagnostic!(d::MassDiagnostic, u, params, inputs, t0)
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    if !(params.backend isa CPU && params.AD isa ContGal && params.M isa AbstractVector)
        rank == 0 && @warn ":lmass_diagnostic needs the CPU backend, continuous Galerkin and the lumped mass matrix; disabled."
        return nothing
    end
    d.mass0 = global_mass(d, u, params, comm)
    d.step  = 0
    # domain area (2D) / volume (3D) = Σᵢ Mᵢ over owned nodes: (mass − mass₀)/area
    # is the domain-mean change of q₁ — for shallow water the equivalent change of
    # the layer thickness in m, which the live plotter shows
    d.red[1] = _local_area(params.M, params.mesh.gip2owner, Int(params.mesh.npoin), rank)
    MPI.Allreduce!(d.red, MPI.SUM, comm)
    area = d.red[1]
    q1name = _first_variable_name(params)
    if rank == 0
        d.dir = string(inputs[:output_dir])
        mkpath(d.dir)
        d.io = open(joinpath(d.dir, "mass.csv"), "w")
        println(d.io, "# Jexpresso mass diagnostic: mass = sum_i M_i q1_i, the integral of the first conserved")
        println(d.io, "# variable with the lumped mass matrix (shallow water: the layer volume, m^3);")
        println(d.io, "# rel_change = (mass - mass0)/mass0, the relative change. Written every ", d.interval, " accepted step(s).")
        println(d.io, "# area = ", area)
        println(d.io, "# q1 = ", q1name)
        println(d.io, "step,t_s,t_days,mass,rel_change")
        _write_mass_line!(d, 0, Float64(t0), d.mass0)
        println(" # Mass diagnostic: ", joinpath(d.dir, "mass.csv"), " (every ", d.interval,
                " steps; live plot: python3 tools/plot_mass_live.py ", d.dir, ")")
        flush(stdout)
    end
    d.active = true
    return nothing
end

function mass_diagnostic_condition(d::MassDiagnostic)
    d.active || return false
    d.step += 1
    return d.step % d.interval == 0
end

function mass_diagnostic_affect!(d::MassDiagnostic, integrator)
    mass = global_mass(d, integrator.u, integrator.p, get_mpi_comm())
    _write_mass_line!(d, d.step, Float64(integrator.t), mass)
    return nothing
end

function close_mass_diagnostic!(d::MassDiagnostic)
    d.active = false
    if d.io !== nothing
        close(d.io)
        d.io = nothing
    end
    return nothing
end
