# Lasher-Trapp et al. (2001) CCOPE deep convection: balanced initialization.
#
#   * The reference (qe) is the HYDROSTATIC rest state, dP/dz = -ρ_model(z, P) g, with ρ solved from the
#     model's own equation of state (saturation_adjustment_sam_microphysics! on one point), so
#     p_model(qe) = P_ref exactly and the rest state stays at round-off. The sounding only supplies the
#     target T(z) and qv(z); its pressure only anchors the integration at the lowest level.
#   * The sounding is read directly (all rows are data) and interpolated linearly in z. Below its lowest
#     level, T and qv are extrapolated linearly from the first two levels; above the highest level they are
#     held constant.
#   * Initiation (paper): a Gaussian surface heat flux (user_bc.jl) plus "a small random velocity
#     perturbation (0.5 m/s on the coarse grid, 0.2 m/s on the fine grid) added to the horizontal wind within
#     the lowest kilometer". The perturbation here is zero-mean, uniform in ±amp (amp = :lt_noise_amp,
#     default 0.5 m/s), and a deterministic hash of the periodically wrapped node coordinates, so periodic
#     copies, MPI ranks and restarts agree.
#   * Any sounding wind lives in the perturbation; the reference momentum is zero.

using MPI

# ---------------------------------------------------------------------------------------------------------------
# Sounding: columns z [m], p [Pa], T [°C], qv [g/kg], u [m/s], v [m/s]
# ---------------------------------------------------------------------------------------------------------------
const LT_SND = Ref{Matrix{Float64}}(zeros(0, 0))
# global domain extent (xmin, xmax, ymin, ymax, zmax), set in initialize() (also on restart): the sponges in
# user_source! need it, and the kernel does not pass the extent to user_source! in this path
const LT_DOMAIN = Ref((0.0, 0.0, 0.0, 0.0, 0.0))

function lt_read_sounding(fname)
    rows = Vector{Vector{Float64}}()
    for l in eachline(fname)
        s = strip(l)
        isempty(s) && continue
        push!(rows, parse.(Float64, split(s)))
    end
    return permutedims(reduce(hcat, rows))
end

# T [K], qv [kg/kg], u, v at height z
function lt_profile(z)
    S = LT_SND[]
    n = size(S, 1)
    if z <= S[1, 1]
        k, t = 1, (z - S[1, 1]) / (S[2, 1] - S[1, 1])            # t <= 0: linear extrapolation
    elseif z >= S[n, 1]
        k, t = n - 1, 1.0                                          # held at the top value
    else
        k = searchsortedlast(@view(S[:, 1]), z)
        k = min(k, n - 1)
        t = (z - S[k, 1]) / (S[k+1, 1] - S[k, 1])
    end
    lerp(c) = (1 - t) * S[k, c] + t * S[k+1, c]
    return lerp(3) + 273.15, max(lerp(4), 0.0) / 1000.0, lerp(5), lerp(6)
end

# target thermodynamic state of the sounding: hl = cp T + g z (no condensate in the sounding), qt = qv
function lt_target(z, PhysConst)
    T, qv, _, _ = lt_profile(z)
    return PhysConst.cp * T + PhysConst.g * z, qv
end

# ---------------------------------------------------------------------------------------------------------------
# The model's own equation of state at one point, and the density giving a target pressure
# ---------------------------------------------------------------------------------------------------------------
function lt_model_eos_1pt(ρ, hl, qt, z, MicroConst, PhysConst)
    uaux = zeros(1, 8)
    uaux[1, 1] = ρ; uaux[1, 5] = ρ * hl; uaux[1, 6] = ρ * qt
    qe = zeros(1, 8); qe[1, 1] = 1.0
    Tabs = zeros(1); qn = zeros(1); qi = zeros(1); qc = zeros(1)
    qr = zeros(1); qs = zeros(1); qg = zeros(1); qsatt = zeros(1)
    saturation_adjustment_sam_microphysics!(uaux, qe, Tabs, qn, qi, qc, qr, qs, qg, qsatt, [z], 1,
                                            MicroConst, PhysConst, NSD_3D(), false)
    return uaux[1, end], Tabs[1], qn[1]
end

function lt_model_density(P, hl, qt, z, MicroConst, PhysConst)
    ρ = P / (PhysConst.Rair * 280.0)
    for it in 1:50
        Pm, _, _ = lt_model_eos_1pt(ρ, hl, qt, z, MicroConst, PhysConst)
        r = P / Pm
        ρ *= r
        abs(r - 1) < 1e-14 && break
    end
    return ρ
end

# ---------------------------------------------------------------------------------------------------------------
# Hydrostatic reference table on z = 0, h, 2h, ...
# ---------------------------------------------------------------------------------------------------------------
struct LTRef
    h::Float64
    P::Vector{Float64}
    dP::Vector{Float64}
end

function build_lt_reference(MicroConst, PhysConst, zmax; h = 0.5)
    S  = LT_SND[]
    z0, P0 = S[1, 1], S[1, 2]                       # anchor: lowest sounding level
    n  = ceil(Int, zmax / h) + 2
    f(z, p) = begin
        hl, qt = lt_target(z, PhysConst)
        -lt_model_density(p, hl, qt, z, MicroConst, PhysConst) * PhysConst.g
    end
    rk4(z, p, dz) = begin
        k1 = f(z, p); k2 = f(z + dz/2, p + dz/2*k1); k3 = f(z + dz/2, p + dz/2*k2); k4 = f(z + dz, p + dz*k3)
        p + dz/6*(k1 + 2k2 + 2k3 + k4)
    end
    P  = zeros(n)
    ks = clamp(round(Int, z0 / h) + 1, 1, n)       # grid point nearest the anchor
    P[ks] = rk4(z0, P0, (ks - 1) * h - z0)          # short step from the anchor onto the grid
    for k in ks:n-1                                 # upward
        P[k+1] = rk4((k - 1) * h, P[k], h)
    end
    for k in ks:-1:2                                # downward (below the lowest sounding level)
        P[k-1] = rk4((k - 1) * h, P[k], -h)
    end
    dP = [f((k - 1) * h, P[k]) for k in 1:n]
    return LTRef(h, P, dP)
end

function lt_ref_pressure(R::LTRef, z)
    s = max(z, 0.0) / R.h
    k = min(floor(Int, s) + 1, length(R.P) - 1)
    t = s - (k - 1)
    h00 = 2t^3 - 3t^2 + 1; h10 = t^3 - 2t^2 + t; h01 = -2t^3 + 3t^2; h11 = t^3 - t^2
    return h00 * R.P[k] + h10 * R.h * R.dP[k] + h01 * R.P[k+1] + h11 * R.h * R.dP[k+1]
end

# ---------------------------------------------------------------------------------------------------------------
# Horizontal-wind perturbation: zero-mean, uniform in ±amp below z_top, deterministic in the wrapped coordinates
# ---------------------------------------------------------------------------------------------------------------
function lt_noise(x, y, z, xmin, ymin, Lx, Ly; amp = 0.5, z_top = 1000.0)
    z > z_top && return 0.0, 0.0
    xi = mod(round(x - xmin, digits = 3), round(Lx, digits = 3)) + 0.0
    yi = mod(round(y - ymin, digits = 3), round(Ly, digits = 3)) + 0.0
    zi = round(z, digits = 3) + 0.0
    u1 = Float64(hash((xi, yi, zi, 1)) >> 11) * 2.0^-53
    u2 = Float64(hash((xi, yi, zi, 2)) >> 11) * 2.0^-53
    return amp * (2u1 - 1), amp * (2u2 - 1)
end

# ---------------------------------------------------------------------------------------------------------------
function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 3D Lasher-Trapp 2001 (hydrostatic reference) ............... "
    end

    qvars    = ("ρ", "ρu", "ρv", "ρw", "hl", "ρqt", "ρqp")
    qoutvars = ["ρ", "ρu", "ρv", "ρw", "hl", "ρqt", "ρqp",
                "T", "qn", "qc", "qi", "qr", "qs", "qg", "qsatt"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat,
                 inputs[:backend]; neqs = length(qvars), qoutvars = qoutvars)

    if inputs[:backend] == CPU()
        PhysConst  = PhysicalConst{Float64}()
        MicroConst = MicrophysicalConst{Float64}()
        LT_SND[]   = lt_read_sounding(inputs[:sounding_file])

        xmin = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)
        xmax = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
        ymin = MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm)
        ymax = MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm)
        zmax = MPI.Allreduce(maximum(mesh.z), MPI.MAX, comm)
        LT_DOMAIN[] = (xmin, xmax, ymin, ymax, zmax)

        if inputs[:lrestart] == true
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path],
                                     inputs, mesh.npoin, HDF5(); nvar = length(qvars))
        else
            Lx, Ly = xmax - xmin, ymax - ymin

            R   = build_lt_reference(MicroConst, PhysConst, zmax + 10.0)
            amp = Float64(get(inputs, :lt_noise_amp, 0.5))

            for ip = 1:mesh.npoin
                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]
                P              = lt_ref_pressure(R, z)
                hl_ref, qt_ref = lt_target(z, PhysConst)
                ρref           = lt_model_density(P, hl_ref, qt_ref, z, MicroConst, PhysConst)
                _, _, us, vs   = lt_profile(z)
                Δu, Δv         = amp > 0 ? lt_noise(x, y, z, xmin, ymin, Lx, Ly; amp = amp) : (0.0, 0.0)
                u, v           = us + Δu, vs + Δv

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip, 1] = 0.0                     # thermodynamics at rest: perturbation zero
                    q.qn[ip, 2] = ρref * u                # reference momentum is zero
                    q.qn[ip, 3] = ρref * v
                    q.qn[ip, 4] = 0.0
                    q.qn[ip, 5] = 0.0
                    q.qn[ip, 6] = 0.0
                else
                    q.qn[ip, 1] = ρref
                    q.qn[ip, 2] = ρref * u
                    q.qn[ip, 3] = ρref * v
                    q.qn[ip, 4] = 0.0
                    q.qn[ip, 5] = ρref * hl_ref
                    q.qn[ip, 6] = ρref * qt_ref
                end
                q.qn[ip, 7]   = 0.0
                q.qn[ip, end] = P

                q.qe[ip, 1]   = ρref
                q.qe[ip, 2]   = 0.0
                q.qe[ip, 3]   = 0.0
                q.qe[ip, 4]   = 0.0
                q.qe[ip, 5]   = ρref * hl_ref
                q.qe[ip, 6]   = ρref * qt_ref
                q.qe[ip, 7]   = 0.0
                q.qe[ip, end] = P
            end

            # ---- self-checks with the model's own EOS ----
            δ = 0.01
            maxp = 0.0; maxres = 0.0; maxqn = 0.0
            for ip = 1:mesh.npoin
                z  = mesh.z[ip]
                ρ  = q.qe[ip, 1]
                Pm, _, qn = lt_model_eos_1pt(ρ, q.qe[ip, 5] / ρ, q.qe[ip, 6] / ρ, z, MicroConst, PhysConst)
                maxp  = max(maxp, abs(Pm - q.qe[ip, end]))
                maxqn = max(maxqn, qn)
                dPdz  = (lt_ref_pressure(R, z + δ) - lt_ref_pressure(R, max(z - δ, 0.0))) / (z + δ - max(z - δ, 0.0))
                maxres = max(maxres, abs(dPdz + ρ * PhysConst.g) / (ρ * PhysConst.g))
            end
            maxp   = MPI.Allreduce(maxp,   MPI.MAX, comm)
            maxres = MPI.Allreduce(maxres, MPI.MAX, comm)
            maxqn  = MPI.Allreduce(maxqn,  MPI.MAX, comm)
            if rank == 0
                @info "LT reference: max |p_model(qe) - P_ref| = $(maxp) Pa   (should be round-off)"
                @info "LT reference: max |dP/dz + ρ g| / (ρ g) = $(maxres)"
                @info "LT reference: max condensate in the reference = $(maxqn) kg/kg   (0 if the sounding is subsaturated)"
                @info "LT reference: P(z=0) = $(lt_ref_pressure(R, 0.0)) Pa, P(z_top) = $(lt_ref_pressure(R, zmax)) Pa; wind noise amp = $(amp) m/s"
            end
        end

        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                ρtot = q.qn[:, 1] .+ q.qe[:, 1]
                q.qn[:, 2] ./= ρtot
                q.qn[:, 3] ./= ρtot
                q.qn[:, 4] ./= ρtot
                q.qn[:, 5] ./= ρtot
                q.qe[:, 5] ./= q.qe[:, 1]
            else
                q.qn[:, 2] ./= q.qn[:, 1]
                q.qn[:, 3] ./= q.qn[:, 1]
                q.qn[:, 4] ./= q.qn[:, 1]
                q.qn[:, 5] ./= q.qn[:, 1]
                q.qn[:, 6] ./= q.qn[:, 1]
                q.qn[:, 7] ./= q.qn[:, 1]
                q.qe[:, 5] ./= q.qe[:, 1]
            end
        end
    end

    if rank == 0
        @info " Initialize fields for 3D Lasher-Trapp 2001 ........................ DONE "
    end
    return q
end
