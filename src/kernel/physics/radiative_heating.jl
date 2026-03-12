Base.@kwdef mutable struct Atmosphere_State{T <: AbstractFloat, dim}
    t_lay = zeros(T, dim)
    p_lay = zeros(T, dim)
    vmr_h2o = zeros(T, dim)
    q_liq = zeros(T, dim)
    q_ice = zeros(T, dim)
    rho = zeros(T, dim)
    vmr_o3 = zeros(T, dim)
end


"""
    compute_radiative_heating(solution_new, atmos_data, mesh, inputs,
        connijk_spa, extra_meshes_connijk, extra_meshes_coords,
        extra_meshes_extra_Je, extra_meshes_extra_nops,
        extra_meshes_extra_nelems, nelem, ngl, κ, σ, regime::Symbol)

Compute the radiative heating rate Q_rad (W/m³) at each spatial node from the
RTE solution field I(x,Ω).

# Physics

The volumetric heating rate is derived from the first moment of the RTE:

    ∇·F_net = κ_abs(x) [G(x) - 4π B(T(x))]

where:
    G(x) = ∫_{4π} I(x,Ω) dΩ        (mean intensity × 4π, i.e. actinic flux)
    B(T)  = (σ_SB / π) T⁴           (Planck blackbody radiance, gray approximation)

The heating rate seen by the atmospheric dynamical core is:

    Q_rad(x) = -∇·F_net / (ρ cₚ)    (K/s, temperature tendency)

or in flux divergence form (W/m³):

    dT/dt|_rad = κ_abs(x) [4π B(T) - G(x)] / (ρ cₚ)

For the shortwave regime B(T) = 0 (no thermal emission at solar wavelengths)
and the heating is purely absorptive:

    dT/dt|_rad = κ_abs(x) G(x) / (ρ cₚ)    [SW]

# Method

Rather than differencing the flux divergence (which amplifies numerical noise
on the spectral element mesh), we use the equivalent absorption form above,
which only requires the angular integral G = ∫ I dΩ at each spatial node.
This is computed from the same angular quadrature used in the forward solve.

# Arguments
- `solution_new`: prolonged solution vector on full mesh, length n_spa
- `atmos_data`: atmospheric state struct (needs t_lay, rho)
- `mesh`: spatial mesh struct
- `κ`: volumetric absorption coefficient (m⁻¹), length npoin — from atmos_to_rad_*
- `σ`: volumetric scattering coefficient (m⁻¹), length npoin
- `regime`: :longwave or :shortwave

# Returns
- `Q_rad`: radiative heating rate (W/m³), length npoin — add to dynamical core
            energy equation as a source term
- `dTdt_rad`: temperature tendency (K/s), length npoin — directly usable as
              forcing in a theta/T prognostic equation
- `F_net`: net flux magnitude (W/m²) at each spatial node, for diagnostics

# Notes
For coupling to an atmospheric model:
- In the energy equation: ∂(ρ e)/∂t = ... + Q_rad
- In the temperature equation: ∂T/∂t = ... + dTdt_rad
- The sign convention: Q_rad > 0 means warming (net absorption > emission)
- For a two-stream or multi-stream solver the G integral naturally handles
  the angular distribution; no assumption about isotropy is made here.
"""
function compute_rt_radiative_heating(
        solution, atmos_data, mesh,
        connijk_spa, extra_connijk, extra_coords,
        extra_Je, extra_nops, extra_nelem, extra_npoin,
        nelem, ngl, κ, σ, ωθ, ωϕ, node_div, ladaptive, swlw::Bool;
        # Optional direct beam quantities for SW direct-diffuse splitting.
        # If G_dir and Q_dir are nothing the function behaves exactly as before.
        G_dir  :: Union{Vector{Float64}, Nothing} = nothing,
        Q_dir  :: Union{Vector{Float64}, Nothing} = nothing,
        sw_μ₀  :: Float64 = 0.0,
        sw_φ₀  :: Float64 = 0.0)

    PhysConst = PhysicalConst{Float64}()

    npoin    = mesh.npoin
    σ_SB     = 5.670374419e-8
    cₚ       = PhysConst.cp

    Q_rad    = zeros(Float64, npoin)
    dTdt_rad = zeros(Float64, npoin)
    F_net    = zeros(Float64, npoin)
    G_accum  = zeros(Float64, npoin)
    Fx_accum = zeros(Float64, npoin)
    Fy_accum = zeros(Float64, npoin)
    Fz_accum = zeros(Float64, npoin)

    # ── Angular integration of diffuse solution ───────────────────────────────
    for iel = 1:nelem
        for i = 1:ngl, j = 1:ngl, k = 1:ngl
            ip = mesh.connijk[iel, i, j, k]
            if ladaptive
                for e_ext = 1:extra_nelem[iel]
                    nop = extra_nops[iel][e_ext]
                    for iθ = 1:nop+1, iϕ = 1:nop+1
                        ip_ext = extra_connijk[iel][e_ext, iθ, iϕ]
                        θ      = extra_coords[iel][1, ip_ext]
                        ϕ      = extra_coords[iel][2, ip_ext]
                        Je_ang = extra_Je[iel][e_ext, iθ, iϕ]
                        ip_g   = connijk_spa[iel][i, j, k, e_ext, iθ, iϕ]
                        I_val  = solution[ip_g]
                        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
                        dΩ = Je_ang * ωθ[iθ] * ωϕ[iϕ]
                        G_accum[ip]  += I_val * dΩ / node_div[ip]
                        Fx_accum[ip] += I_val * Ωx * dΩ / node_div[ip]
                        Fy_accum[ip] += I_val * Ωy * dΩ / node_div[ip]
                        Fz_accum[ip] += I_val * Ωz * dΩ / node_div[ip]
                    end
                end
            else
                for e_ext = 1:extra_nelem
                    nop = extra_nops[e_ext]
                    for iθ = 1:nop+1, iϕ = 1:nop+1
                        ip_ext = extra_connijk[e_ext, iθ, iϕ]
                        θ      = extra_coords[1, ip_ext]
                        ϕ      = extra_coords[2, ip_ext]
                        Je_ang = extra_Je[e_ext, iθ, iϕ]
                        ip_g   = (ip-1) * extra_npoin + ip_ext
                        I_val  = solution[ip_g]
                        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
                        dΩ = Je_ang * ωθ[iθ] * ωϕ[iϕ]
                        G_accum[ip]  += I_val * dΩ / node_div[ip]
                        Fx_accum[ip] += I_val * Ωx * dΩ / node_div[ip]
                        Fy_accum[ip] += I_val * Ωy * dΩ / node_div[ip]
                        Fz_accum[ip] += I_val * Ωz * dΩ / node_div[ip]
                    end
                end
            end
        end
    end

    # ── MPI reduction of diffuse angular integrals ────────────────────────────
    if MPI.Comm_size(MPI.COMM_WORLD) > 1
        gnpoin_spa = mesh.gnpoin
        g_G  = zeros(Float64, gnpoin_spa)
        g_Fx = zeros(Float64, gnpoin_spa)
        g_Fy = zeros(Float64, gnpoin_spa)
        g_Fz = zeros(Float64, gnpoin_spa)

        for ip = 1:npoin
            gip = mesh.ip2gip[ip]
            g_G[gip]  += G_accum[ip]
            g_Fx[gip] += Fx_accum[ip]
            g_Fy[gip] += Fy_accum[ip]
            g_Fz[gip] += Fz_accum[ip]
        end

        MPI.Allreduce!(g_G,  MPI.SUM, MPI.COMM_WORLD)
        MPI.Allreduce!(g_Fx, MPI.SUM, MPI.COMM_WORLD)
        MPI.Allreduce!(g_Fy, MPI.SUM, MPI.COMM_WORLD)
        MPI.Allreduce!(g_Fz, MPI.SUM, MPI.COMM_WORLD)

        for ip = 1:npoin
            gip = mesh.ip2gip[ip]
            G_accum[ip]  = g_G[gip]
            Fx_accum[ip] = g_Fx[gip]
            Fy_accum[ip] = g_Fy[gip]
            Fz_accum[ip] = g_Fz[gip]
        end
    end

    # ── Add direct beam contribution to G and flux if splitting is active ─────
    # G_dir and Q_dir are non-nothing only for the SW case with direct-diffuse
    # splitting. For LW and for SW without splitting they remain nothing and
    # the function behaves identically to the original.
    splitting = !isnothing(G_dir) && !isnothing(Q_dir)

    if splitting
        bdy_nodes = [ip for ip = 1:mesh.npoin 
             if mesh.z[ip] ≈ mesh.zmin || 
                mesh.x[ip] ≈ mesh.xmin || mesh.x[ip] ≈ mesh.xmax ||
                mesh.y[ip] ≈ mesh.ymin || mesh.y[ip] ≈ mesh.ymax]
                
        G_total = G_accum .+ G_dir
        @info "G_dir    on boundaries: $(round.(extrema(G_dir[bdy_nodes]),   sigdigits=4))"
        @info "G_diffuse on boundaries: $(round.(extrema(G_accum[bdy_nodes]), sigdigits=4))"
        @info "G_total  on boundaries: $(round.(extrema(G_total[bdy_nodes]),  sigdigits=4))"

        θ_sun = π - acos(clamp(sw_μ₀, 0.0, 1.0))
        Ωx_sun = sin(θ_sun)*cos(sw_φ₀)
        Ωy_sun = sin(θ_sun)*sin(sw_φ₀)
        Ωz_sun = cos(θ_sun)
        for ip = 1:npoin
            G_accum[ip]  += G_dir[ip]
            F_dir_ip      = G_dir[ip] * sw_μ₀   # F_dir = G_dir × μ₀
            Fx_accum[ip] += F_dir_ip * Ωx_sun
            Fy_accum[ip] += F_dir_ip * Ωy_sun
            Fz_accum[ip] += F_dir_ip * Ωz_sun
        end
    end

    # ── Heating rate at each spatial node ─────────────────────────────────────
    for ip = 1:npoin
        T   = atmos_data.t_lay[ip]
        ρ   = atmos_data.rho[ip]
        κ_a = κ[ip]

        G = G_accum[ip]   # total G: diffuse + direct (if splitting active)

        if swlw
            # Longwave: absorption minus emission
            B4π       = 4.0 * σ_SB * T^4
            Q_rad[ip] = κ_a * (G - B4π)
        else
            if splitting
                # SW with direct-diffuse splitting:
                # Q_total = κ_abs × G_total
                #         = κ_abs × (G_diffuse + G_dir)
                #         = κ_abs × G_accum  (G_dir already added above)
                # This is equivalent to Q_dir + κ_abs × G_diffuse.
                # We use Q_dir directly since it was computed from the exact
                # Beer-Lambert F_dir, which is more accurate than
                # κ_abs × G_dir where G_dir comes from the discrete τ.
                Q_rad[ip] = Q_dir[ip] + κ_a * (G - G_dir[ip])
            else
                # SW without splitting: standard
                Q_rad[ip] = κ_a * G
            end
        end

        F_net[ip]    = sqrt(Fx_accum[ip]^2 + Fy_accum[ip]^2 + Fz_accum[ip]^2)
        dTdt_rad[ip] = Q_rad[ip] / (ρ * cₚ)
    end

    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    if rank == 0
        @info "Radiative heating$(splitting ? " [SW direct-diffuse split]" : "")"
        @info "  Q_rad   extrema (W/m³): $(extrema(Q_rad))"
        @info "  dT/dt   extrema (K/s):  $(extrema(dTdt_rad))"
        @info "  F_net   extrema (W/m²): $(extrema(F_net))"
        @info "  G       extrema (W/m²): $(extrema(G_accum))"
        if splitting
            @info "  G_dir   extrema (W/m²): $(extrema(G_dir))"
            @info "  G_dif   extrema (W/m²): $(extrema(G_accum .- G_dir))"
            @info "  Q_dir   extrema (W/m³): $(extrema(Q_dir))"
        end
    end

    return Q_rad, dTdt_rad, F_net, G_accum
end

function dycore_to_atmos_data!(q, qe, npoin, T, qc, qi, lmoist, atmos_data, PhysConst, ::TOTAL)
    for i=1:npoin
        if (lmoist)
            atmos_data.t_lay[i] =  T[i]
            atmos_data.qv = ((q[i,6]/q[i,1])-qc[i]-qi[i])*PhysConst.Mol_mass_water/ PhysConst.Mol_mass_air
            atmos_data.q_liq[i] = qc[i]
            atmos_data.q_ice[i] = qi[i]
        else
            atmos_data.t_lay[i] = (q[i,5]/q[i,1])/((PhysConst.pref/q[end])^(PhysConst.Rair/PhysConst.cp))
        end
        atmos_data.p_lay[i] = q[i,end]
        atmos_data.rho[i] = q[i,1]
    end

end

function dycore_to_atmos_data!(q, qe, npoin, T, qc, qi, lmoist, atmos_data, PhysConst, ::PERT)

    for i=1:npoin
        ρ = q[i,1]+qe[i,1]
        θ = (q[i,5]+qe[i,5])/ρ
        qt = (q[i,6]+qe[i,6])/ρ
        if (lmoist)
            atmos_data.t_lay[i] =  T[i]
            atmos_data.vmr_h2o[i] = ((qt)-qc[i]-qi[i])*PhysConst.Mol_mass_water/ PhysConst.Mol_mass_air
            atmos_data.q_liq[i] = qc[i]
            atmos_data.q_ice[i] = qi[i]
        else
            atmos_data.t_lay[i] = (θ)/((PhysConst.pref/q[i,end])^(PhysConst.Rair/PhysConst.cp))
        end
        atmos_data.p_lay[i] = q[i,end]
        atmos_data.rho[i] = ρ
    end
end

function get_RT_heat_fluxes!(q, qe, mesh, micro, metrics, atmos_data, params, dψ, ψ, ω, PhysConst, inputs)

    #First get atmospheric states for optics:
    dycore_to_atmos_data!(q,qe, mesh.npoin, micro.Tabs, micro.qc, micro.qi, inputs[:lmoist], atmos_data, PhysConst, inputs[:SOL_VARS_TYPE])
    #Two RT problems need to be solved
    #First do shortwave
    #Get gray shortwave absorption and scattering coefficients
    κ, σ, inputs[:rad_HG_g] = atmos_to_rad_shortwave(atmos_data,mesh.npoin)
    κ_ext_sw         = κ .+ σ
    #compute sample optical depth for lateral boundary conditions
    z_prof, τ_from_TOA = build_sw_lateral_bc_profile(mesh, κ_ext_sw, mesh.ngl)
    #Solve shortwave RT problem and get heating rate
    inputs[:RT_shortwave] = true
    inputs[:RT_longwave] = false
    Q, dTdt, micro.rt_sol_sw = build_radiative_transfer_problem(mesh, inputs, 1, mesh.ngl, dψ, ψ, ω, metrics.Je,
                                     metrics.dξdx, metrics.dξdy, metrics.dξdz, 
                                     metrics.dηdx, metrics.dηdy, metrics.dηdz,
                                     metrics.dζdx, metrics.dζdy, metrics.dζdz,
                                     metrics.nx, metrics.ny, metrics.nz, 
                                     mesh.elem_to_face, mesh.extra_mesh, κ, σ, atmos_data, z_prof, τ_from_TOA, params.QT, NSD_3D(), params.AD;
                                     rt_sol_sw = micro.rt_sol_sw, rt_sol_sw_available = micro.rt_sol_sw_available)

    if !(inputs[:energy_equation] == "theta") || (inputs[:lmoist])
        micro.flux_sw .= Q
    else
        micro.flux_sw .= dTdt
    end

    #Second do longwave
    #Get gray longwave absorption and scattering coefficients
    κ, σ = atmos_to_rad_longwave(atmos_data,mesh.npoin)
    #Solve longwave RT problem and get heating rate
    inputs[:RT_shortwave] = false
    inputs[:RT_longwave] = true
    inputs[:rad_HG_g] = 0.0
    Q, dTdt, micro.rt_sol_lw = build_radiative_transfer_problem(mesh, inputs, 1, mesh.ngl, dψ, ψ, ω, metrics.Je,
                                     metrics.dξdx, metrics.dξdy, metrics.dξdz, 
                                     metrics.dηdx, metrics.dηdy, metrics.dηdz,
                                     metrics.dζdx, metrics.dζdy, metrics.dζdz,
                                     metrics.nx, metrics.ny, metrics.nz, 
                                     mesh.elem_to_face, mesh.extra_mesh, κ, σ, atmos_data, z_prof, τ_from_TOA, params.QT, NSD_3D(), params.AD;
                                     rt_sol_lw = micro.rt_sol_lw, rt_sol_lw_available = micro.rt_sol_lw_available)
                                     
    if !(inputs[:energy_equation] == "theta") || inputs[:lmoist]
        micro.flux_lw .= -Q
    else
        micro.flux_lw .= -dTdt
    end
    micro.rt_sol_sw_available = true
    micro.rt_sol_lw_available = true
end


