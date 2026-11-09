# Balanced BOMEX initialization (UNTESTED: written without a Julia install; please run the checks at the bottom).
#
# What changes relative to the original:
#   1. The reference pressure is the HYDROSTATIC solution of dP/dz = -rho(z,P) g, with rho solved from the MODEL's
#      own pressure (your saturation_adjustment_sam_microphysics!, called on one point), so p_model(qe) = P_ref exactly
#      even though the model's constants (0.61, Rair, qsat) differ slightly from TD's.  TD only defines the target
#      theta_l/q_t -> T, q_l profile.  Before: isothermal P with TD's rho: neither hydrostatic nor model-consistent.
#   2. The optional noise is applied to theta_l and q_t BEFORE the thermodynamic state is built (at fixed P), so
#      rho, T, q_l and hl are all consistent with the perturbed state and p' = 0 at t = 0.
#   3. Noise is a deterministic hash of the node coordinates (periodic images wrapped), so shared nodes, MPI ranks
#      and restarts agree.  Default off.  Enable with inputs[:bomex_noise] = true.
#   4. Profile branches cover z < 0 (round-off) correctly.
#
# The reference (qe) and the initial total are built from the same P, so the perturbation q is exactly zero when
# noise is off.

using MPI

# ---------------------------------------------------------------------------------------------------------------
# Profiles (Siebesma et al. 2003, as in your original code)
# ---------------------------------------------------------------------------------------------------------------
function bomex_profiles(z::Float64)
    zl1, zl2, zl3, zl4 = 520.0, 1480.0, 2000.0, 3000.0
    z = max(z, 0.0)
    if z <= zl1
        θl = 298.7
        qt = 17.0 + (z / zl1) * (16.3 - 17.0)
    elseif z <= zl2
        θl = 298.7 + (z - zl1) * (302.4 - 298.7) / (zl2 - zl1)
        qt = 16.3 + (z - zl1) * (10.7 - 16.3) / (zl2 - zl1)
    elseif z <= zl3
        θl = 302.4 + (z - zl2) * (308.2 - 302.4) / (zl3 - zl2)
        qt = 10.7 + (z - zl2) * (4.2 - 10.7) / (zl3 - zl2)
    else
        θl = 308.2 + (z - zl3) * (311.85 - 308.2) / (zl4 - zl3)
        qt = 4.2 + (z - zl3) * (3.0 - 4.2) / (zl4 - zl3)
    end
    return θl, qt / 1000.0
end

@inline function bomex_thermo(param_set, P, θl, qt)
    TS  = TD.PhaseEquil_pθq(param_set, P, θl, qt)
    T   = TD.air_temperature(param_set, TS)
    ρ   = TD.air_density(param_set, TS)
    qpt = TD.PhasePartition(param_set, TS)
    return T, ρ, qpt.liq
end

# ---------------------------------------------------------------------------------------------------------------
# The MODEL's own equation of state, evaluated at one point by calling your saturation adjustment.
# Everything that must balance (rho_ref given P_ref) is solved with this, so p' = 0 holds whatever constants
# (0.61 vs R_v/R_d - 1, Rair, cp, Lc, qsat formulas) the model uses.  TD only defines the target T, q_l profile.
# ---------------------------------------------------------------------------------------------------------------
function model_eos_1pt(ρ, hl, qt, z, MicroConst, PhysConst)
    uaux = zeros(1, 8)                  # 7 prognostic variables + pressure in the last column
    uaux[1, 1] = ρ; uaux[1, 5] = ρ * hl; uaux[1, 6] = ρ * qt
    qe = zeros(1, 8); qe[1, 1] = 1.0
    Tabs = zeros(1); qn = zeros(1); qi = zeros(1); qc = zeros(1)
    qr = zeros(1); qs = zeros(1); qg = zeros(1); qsatt = zeros(1)
    saturation_adjustment_sam_microphysics!(uaux, qe, Tabs, qn, qi, qc, qr, qs, qg, qsatt, [z], 1,
                                            MicroConst, PhysConst, NSD_3D(), false)
    return uaux[1, end], Tabs[1], qc[1]
end

# density for which the model's pressure equals P at fixed (hl, qt)
function model_density(P, hl, qt, z, MicroConst, PhysConst)
    ρ = P / (PhysConst.Rair * 290.0)
    for it in 1:40
        Pm, _, _ = model_eos_1pt(ρ, hl, qt, z, MicroConst, PhysConst)
        r = P / Pm
        ρ *= r
        abs(r - 1) < 1e-14 && break
    end
    return ρ
end

# target profile: TD gives T and q_l at (theta_l, q_t, P); we only keep hl and q_t from it
function bomex_target(param_set, PhysConst, z, P, Δθ, Δqt)
    θl, qt = bomex_profiles(z)
    T, _, ql = bomex_thermo(param_set, P, θl + Δθ, qt + Δqt)
    return PhysConst.cp * T + PhysConst.g * z - PhysConst.Lc * ql, qt + Δqt
end

# ---------------------------------------------------------------------------------------------------------------
# Hydrostatic reference table: dP/dz = -rho_model(z,P) g, g = PhysConst.g (what the source term uses)
# ---------------------------------------------------------------------------------------------------------------
struct BomexRef
    h::Float64
    P::Vector{Float64}
    dP::Vector{Float64}
end

function build_bomex_reference(param_set, MicroConst, PhysConst, zmax; h = 0.5, P_sfc = 1.015e5)
    n  = ceil(Int, zmax / h) + 2
    P  = zeros(n)
    dP = zeros(n)
    f(z, p) = begin
        hl, qt = bomex_target(param_set, PhysConst, z, p, 0.0, 0.0)
        -model_density(p, hl, qt, z, MicroConst, PhysConst) * PhysConst.g
    end
    P[1] = P_sfc
    for k in 1:n-1
        z  = (k - 1) * h
        k1 = f(z, P[k])
        k2 = f(z + h / 2, P[k] + h / 2 * k1)
        k3 = f(z + h / 2, P[k] + h / 2 * k2)
        k4 = f(z + h, P[k] + h * k3)
        P[k+1] = P[k] + h / 6 * (k1 + 2k2 + 2k3 + k4)
        dP[k]  = k1
    end
    dP[n] = f((n - 1) * h, P[n])
    return BomexRef(h, P, dP)
end

function ref_pressure(R::BomexRef, z)
    z = max(z, 0.0)
    s = z / R.h
    k = min(floor(Int, s) + 1, length(R.P) - 1)
    t = s - (k - 1)
    h00 = 2t^3 - 3t^2 + 1
    h10 = t^3 - 2t^2 + t
    h01 = -2t^3 + 3t^2
    h11 = t^3 - t^2
    return h00 * R.P[k] + h10 * R.h * R.dP[k] + h01 * R.P[k+1] + h11 * R.h * R.dP[k+1]
end

# ---------------------------------------------------------------------------------------------------------------
# Noise: deterministic in the (periodically wrapped) coordinates
# ---------------------------------------------------------------------------------------------------------------
function bomex_noise(x, y, z, xmin, ymin, Lx, Ly;
                     z_noise = 1600.0, amp_θ = 0.1, amp_q = 2.5e-5)   # verify amplitudes and depth vs the case spec
    z > z_noise && return 0.0, 0.0
    xi = mod(round(x - xmin, digits = 3), round(Lx, digits = 3)) + 0.0
    yi = mod(round(y - ymin, digits = 3), round(Ly, digits = 3)) + 0.0
    zi = round(z, digits = 3) + 0.0
    u1 = Float64(hash((xi, yi, zi, 1)) >> 11) * 2.0^-53
    u2 = Float64(hash((xi, yi, zi, 2)) >> 11) * 2.0^-53
    return amp_θ * (2u1 - 1), amp_q * (2u2 - 1)
end

# ---------------------------------------------------------------------------------------------------------------
# Point state.  Returns reference (ρref, hl_ref, qt_ref, ql_ref), pressure P, and the perturbed total (ρ, hl, qt).
# ---------------------------------------------------------------------------------------------------------------
function initialize_bomex!(z, param_set, MicroConst, PhysConst, R::BomexRef; Δθ = 0.0, Δqt = 0.0)
    P = ref_pressure(R, z)
    hl_ref, qt_ref = bomex_target(param_set, PhysConst, z, P, 0.0, 0.0)
    ρref = model_density(P, hl_ref, qt_ref, z, MicroConst, PhysConst)
    _, _, ql_ref = model_eos_1pt(ρref, hl_ref, qt_ref, z, MicroConst, PhysConst)
    if Δθ == 0.0 && Δqt == 0.0
        return ρref, hl_ref, qt_ref, ql_ref, P, ρref, hl_ref, qt_ref
    end
    hl1, qt1 = bomex_target(param_set, PhysConst, z, P, Δθ, Δqt)
    ρ1 = model_density(P, hl1, qt1, z, MicroConst, PhysConst)      # same P: p' = 0, rho' carries the buoyancy
    return ρref, hl_ref, qt_ref, ql_ref, P, ρ1, hl1, qt1
end

# ---------------------------------------------------------------------------------------------------------------
function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 3D BOMEX (hydrostatic reference) ........................ "
    end

    qvars    = ("ρ", "ρu", "ρv", "ρw", "hl", "ρqt", "ρqp")
    qoutvars = ["ρ", "ρu", "ρv", "ρw", "hl", "ρqt", "ρqp",
                "T", "qn", "qc", "qi", "qr", "qs", "qg", "qsatt", "flux_lw", "flux_sw", "dTdt_lat"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat,
                 inputs[:backend]; neqs = length(qvars), qoutvars = qoutvars)

    if inputs[:backend] == CPU()
        PhysConst = PhysicalConst{Float64}()
        new_param_set = create_updated_TD_Parameters(PhysConst.potential_temperature_reference_pressure)

        # SAME constructor you use where saturation_adjustment_sam_microphysics! is called elsewhere:
        MicroConst = MicrophysicalConst{Float64}()      # <-- adjust to your code
        if MicroConst.Lc != PhysConst.Lc
            @warn "MicroConst.Lc != PhysConst.Lc: hl built with PhysConst.Lc, model inverts with MicroConst.Lc"
        end

        if inputs[:lrestart] == true
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path],
                                     inputs, mesh.npoin, HDF5(); nvar = length(qvars))
        else
            use_noise = get(inputs, :bomex_noise, false)

            # global extents (mesh arrays are per rank)
            xmin = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)
            xmax = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
            ymin = MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm)
            ymax = MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm)
            zmax = MPI.Allreduce(maximum(mesh.z), MPI.MAX, comm)
            Lx, Ly = xmax - xmin, ymax - ymin

            R = build_bomex_reference(new_param_set, MicroConst, PhysConst, max(zmax, 3000.0) + 10.0)

            u = 0.0; v = 0.0; w = 0.0     # rest state; reinstate the wind profile here when needed

            for ip = 1:mesh.npoin
                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]
                Δθ, Δqt = use_noise ? bomex_noise(x, y, z, xmin, ymin, Lx, Ly) : (0.0, 0.0)

                ρref, hl_ref, qt_ref, ql_ref, P, ρ, hl, qt =
                    initialize_bomex!(z, new_param_set, MicroConst, PhysConst, R; Δθ = Δθ, Δqt = Δqt)

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip, 1] = ρ - ρref
                    q.qn[ip, 2] = (ρ - ρref) * u
                    q.qn[ip, 3] = (ρ - ρref) * v
                    q.qn[ip, 4] = (ρ - ρref) * w
                    q.qn[ip, 5] = ρ * hl - ρref * hl_ref
                    q.qn[ip, 6] = ρ * qt - ρref * qt_ref
                else
                    q.qn[ip, 1] = ρ
                    q.qn[ip, 2] = ρ * u
                    q.qn[ip, 3] = ρ * v
                    q.qn[ip, 4] = ρ * w
                    q.qn[ip, 5] = ρ * hl
                    q.qn[ip, 6] = ρ * qt
                end
                q.qn[ip, 7]   = 0.0
                q.qn[ip, end] = P

                q.qe[ip, 1]   = ρref
                q.qe[ip, 2]   = ρref * u
                q.qe[ip, 3]   = ρref * v
                q.qe[ip, 4]   = ρref * w
                q.qe[ip, 5]   = ρref * hl_ref
                q.qe[ip, 6]   = ρref * qt_ref
                q.qe[ip, 7]   = 0.0
                q.qe[ip, end] = P
            end

            # ---- diagnostics (all with the model's own EOS) ----
            δ = 0.01
            maxres = 0.0; maxp = 0.0; maxrho_td = 0.0
            for ip = 1:mesh.npoin
                z  = mesh.z[ip]
                ρ  = q.qe[ip, 1]
                hl = q.qe[ip, 5] / ρ
                qt = q.qe[ip, 6] / ρ
                Pm, _, _ = model_eos_1pt(ρ, hl, qt, z, MicroConst, PhysConst)
                maxp = max(maxp, abs(Pm - q.qe[ip, end]))
                dPdz = (ref_pressure(R, z + δ) - ref_pressure(R, max(z - δ, 0.0))) / (z + δ - max(z - δ, 0.0))
                maxres = max(maxres, abs(dPdz + ρ * PhysConst.g) / (ρ * PhysConst.g))
                θl, qt0 = bomex_profiles(z)
                ρ_td = bomex_thermo(new_param_set, q.qe[ip, end], θl, qt0)[2]
                maxrho_td = max(maxrho_td, abs(ρ - ρ_td) / ρ)
            end
            maxp      = MPI.Allreduce(maxp,      MPI.MAX, comm)
            maxres    = MPI.Allreduce(maxres,    MPI.MAX, comm)
            maxrho_td = MPI.Allreduce(maxrho_td, MPI.MAX, comm)
            if rank == 0
                @info "BOMEX reference: max |p_model(qe) - P_ref| = $(maxp) Pa   (must be ~1e-9 or smaller)"
                @info "BOMEX reference: max |dP/dz + ρ g| / (ρ g) = $(maxres)   (~1e-8 or smaller)"
                @info "BOMEX reference: max |ρ_model - ρ_TD| / ρ  = $(maxrho_td)   (size of the old EOS inconsistency)"
            end
        end
    end

    if rank == 0
        @info " Initialize fields for 3D BOMEX ........................ DONE "
    end
    return q
end
