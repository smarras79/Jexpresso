# Turbulence Statistics Module
# Jexpresso analogue of dod2d mod_aver.f90
# Computes temporal averages of flow quantities for turbulence analysis

"""
    initialize_turbulence_stats!(turb_stats::St_turbulence_stats, t_start::AbstractFloat)

Initialize turbulence statistics structure and set the starting time for averaging.

# Arguments
- `turb_stats`: Turbulence statistics structure
- `t_start`: Starting time for statistics accumulation
"""
function initialize_turbulence_stats!(turb_stats::St_turbulence_stats, t_start::AbstractFloat)

    # Enable statistics
    turb_stats.l_enabled[] = true

    # Reset accumulated time
    turb_stats.acutim[] = zero(eltype(turb_stats.avvel))
    turb_stats.elapsed_time[] = t_start

    # Zero all accumulated arrays
    fill!(turb_stats.avvel, zero(eltype(turb_stats.avvel)))
    fill!(turb_stats.avrho, zero(eltype(turb_stats.avrho)))
    fill!(turb_stats.avpress, zero(eltype(turb_stats.avpress)))
    fill!(turb_stats.avve2, zero(eltype(turb_stats.avve2)))
    fill!(turb_stats.avvex, zero(eltype(turb_stats.avvex)))
    fill!(turb_stats.avpr2, zero(eltype(turb_stats.avpr2)))
    fill!(turb_stats.avmueff, zero(eltype(turb_stats.avmueff)))

    println(" # --------------------------------------------------------")
    println(" # Turbulence statistics initialized at t = ", t_start)
    println(" # --------------------------------------------------------")
end

"""
    reset_turbulence_stats!(turb_stats::St_turbulence_stats)

Reset all accumulated statistics to zero for starting a new averaging window.
"""
function reset_turbulence_stats!(turb_stats::St_turbulence_stats)

    # Reset accumulated time
    turb_stats.acutim[] = zero(eltype(turb_stats.avvel))

    # Zero all accumulated arrays
    fill!(turb_stats.avvel, zero(eltype(turb_stats.avvel)))
    fill!(turb_stats.avrho, zero(eltype(turb_stats.avrho)))
    fill!(turb_stats.avpress, zero(eltype(turb_stats.avpress)))
    fill!(turb_stats.avve2, zero(eltype(turb_stats.avve2)))
    fill!(turb_stats.avvex, zero(eltype(turb_stats.avvex)))
    fill!(turb_stats.avpr2, zero(eltype(turb_stats.avpr2)))
    fill!(turb_stats.avmueff, zero(eltype(turb_stats.avmueff)))

    println(" # Turbulence statistics reset for new averaging window")
end

"""
    accumulate_statistics!(turb_stats, uaux, press, mueff, dt, npoin, ndim, SD)

Accumulate turbulence statistics over a time step dt.
This implements the accumulation method similar to dod2d's favre_average.

# Arguments
- `turb_stats`: Turbulence statistics structure
- `uaux`: Auxiliary array containing [ρ, ρu, ρv, ρw, ...] at each point (npoin, neqs+1)
- `press`: Pressure field (npoin)
- `mueff`: Effective viscosity field (npoin) - can be nothing if not available
- `dt`: Time step size
- `npoin`: Number of points
- `ndim`: Number of spatial dimensions (1, 2, or 3)
- `SD`: Spatial dimension type (NSD_1D(), NSD_2D(), NSD_3D())
"""
function accumulate_statistics!(turb_stats::St_turbulence_stats, uaux, press, mueff, dt, npoin, ndim, SD)

    # Skip if statistics are disabled
    if !turb_stats.l_enabled[]
        return
    end

    # Update accumulated time
    turb_stats.acutim[] += dt

    # Accumulate statistics based on spatial dimension
    if SD == NSD_1D()
        accumulate_statistics_1D!(turb_stats, uaux, press, mueff, dt, npoin)
    elseif SD == NSD_2D()
        accumulate_statistics_2D!(turb_stats, uaux, press, mueff, dt, npoin)
    elseif SD == NSD_3D()
        accumulate_statistics_3D!(turb_stats, uaux, press, mueff, dt, npoin)
    end
end

"""
    accumulate_statistics_1D!(turb_stats, uaux, press, mueff, dt, npoin)

Accumulate 1D turbulence statistics.
"""
function accumulate_statistics_1D!(turb_stats::St_turbulence_stats, uaux, press, mueff, dt, npoin)

    for ip = 1:npoin
        ρ = uaux[ip, 1]
        ρu = uaux[ip, 2]
        p = press[ip]

        # Compute velocity
        u = ρu / ρ

        # Accumulate mean quantities (weighted by dt)
        turb_stats.avvel[ip, 1] += u * dt
        turb_stats.avrho[ip] += ρ * dt
        turb_stats.avpress[ip] += p * dt

        # Accumulate second-order moments
        turb_stats.avve2[ip, 1] += (u * u) * dt
        turb_stats.avpr2[ip] += (p * p) * dt

        # Accumulate effective viscosity if available
        if !isnothing(mueff)
            turb_stats.avmueff[ip] += mueff[ip] * dt
        end
    end
end

"""
    accumulate_statistics_2D!(turb_stats, uaux, press, mueff, dt, npoin)

Accumulate 2D turbulence statistics including Reynolds stress <u*v>.
"""
function accumulate_statistics_2D!(turb_stats::St_turbulence_stats, uaux, press, mueff, dt, npoin)

    for ip = 1:npoin
        ρ = uaux[ip, 1]
        ρu = uaux[ip, 2]
        ρv = uaux[ip, 3]
        p = press[ip]

        # Compute velocities
        u = ρu / ρ
        v = ρv / ρ

        # Accumulate mean quantities (weighted by dt)
        turb_stats.avvel[ip, 1] += u * dt
        turb_stats.avvel[ip, 2] += v * dt
        turb_stats.avrho[ip] += ρ * dt
        turb_stats.avpress[ip] += p * dt

        # Accumulate second-order moments
        turb_stats.avve2[ip, 1] += (u * u) * dt
        turb_stats.avve2[ip, 2] += (v * v) * dt
        turb_stats.avpr2[ip] += (p * p) * dt

        # Accumulate Reynolds stress components
        turb_stats.avvex[ip, 1] += (u * v) * dt  # <u*v>

        # Accumulate effective viscosity if available
        if !isnothing(mueff)
            turb_stats.avmueff[ip] += mueff[ip] * dt
        end
    end
end

"""
    accumulate_statistics_3D!(turb_stats, uaux, press, mueff, dt, npoin)

Accumulate 3D turbulence statistics including all Reynolds stress components.
"""
function accumulate_statistics_3D!(turb_stats::St_turbulence_stats, uaux, press, mueff, dt, npoin)

    for ip = 1:npoin
        ρ = uaux[ip, 1]
        ρu = uaux[ip, 2]
        ρv = uaux[ip, 3]
        ρw = uaux[ip, 4]
        p = press[ip]

        # Compute velocities
        u = ρu / ρ
        v = ρv / ρ
        w = ρw / ρ

        # Accumulate mean quantities (weighted by dt)
        turb_stats.avvel[ip, 1] += u * dt
        turb_stats.avvel[ip, 2] += v * dt
        turb_stats.avvel[ip, 3] += w * dt
        turb_stats.avrho[ip] += ρ * dt
        turb_stats.avpress[ip] += p * dt

        # Accumulate second-order moments
        turb_stats.avve2[ip, 1] += (u * u) * dt
        turb_stats.avve2[ip, 2] += (v * v) * dt
        turb_stats.avve2[ip, 3] += (w * w) * dt
        turb_stats.avpr2[ip] += (p * p) * dt

        # Accumulate Reynolds stress components
        turb_stats.avvex[ip, 1] += (u * v) * dt  # <u*v>
        turb_stats.avvex[ip, 2] += (u * w) * dt  # <u*w>
        turb_stats.avvex[ip, 3] += (v * w) * dt  # <v*w>

        # Accumulate effective viscosity if available
        if !isnothing(mueff)
            turb_stats.avmueff[ip] += mueff[ip] * dt
        end
    end
end

"""
    normalize_statistics!(turb_stats::St_turbulence_stats)

Normalize accumulated statistics by dividing by total accumulated time.
This implements the normalization step similar to dod2d's eval_average_window.

After normalization:
- avvel contains time-averaged velocity <u_i>
- avrho contains time-averaged density <ρ>
- avpress contains time-averaged pressure <p>
- avve2 contains <u_i * u_i>
- avvex contains Reynolds stresses <u_i * u_j> for i≠j
- avpr2 contains <p * p>
- avmueff contains time-averaged effective viscosity
"""
function normalize_statistics!(turb_stats::St_turbulence_stats)

    if !turb_stats.l_enabled[]
        @warn "Turbulence statistics are not enabled. Skipping normalization."
        return
    end

    acutim = turb_stats.acutim[]

    if acutim ≈ 0.0
        @warn "Accumulated time is zero. Cannot normalize statistics."
        return
    end

    # Normalize all arrays by dividing by accumulated time
    turb_stats.avvel ./= acutim
    turb_stats.avrho ./= acutim
    turb_stats.avpress ./= acutim
    turb_stats.avve2 ./= acutim
    turb_stats.avvex ./= acutim
    turb_stats.avpr2 ./= acutim
    turb_stats.avmueff ./= acutim

    println(" # --------------------------------------------------------")
    println(" # Turbulence statistics normalized over time = ", acutim)
    println(" # --------------------------------------------------------")
end

"""
    compute_rms_fluctuations!(vel_rms, turb_stats, ndim)

Compute RMS velocity fluctuations from the time-averaged statistics.

# Arguments
- `vel_rms`: Output array for RMS fluctuations (npoin, ndim)
- `turb_stats`: Turbulence statistics structure (must be normalized first!)
- `ndim`: Number of spatial dimensions

# Note
This function assumes normalize_statistics! has already been called.
The RMS fluctuation is computed as: u'_rms = sqrt(<u*u> - <u>^2)
"""
function compute_rms_fluctuations!(vel_rms, turb_stats::St_turbulence_stats, ndim)

    npoin = size(turb_stats.avvel, 1)

    for ip = 1:npoin
        for idim = 1:ndim
            mean_vel = turb_stats.avvel[ip, idim]
            mean_vel2 = turb_stats.avve2[ip, idim]

            # Compute variance: Var(u) = <u^2> - <u>^2
            variance = mean_vel2 - mean_vel * mean_vel

            # Ensure non-negative (numerical errors might make it slightly negative)
            variance = max(variance, zero(variance))

            # RMS = sqrt(variance)
            vel_rms[ip, idim] = sqrt(variance)
        end
    end
end

"""
    compute_reynolds_stresses!(reynolds_stress, turb_stats, ndim)

Compute Reynolds stresses from the time-averaged statistics.

# Arguments
- `reynolds_stress`: Output array for Reynolds stresses
- `turb_stats`: Turbulence statistics structure (must be normalized first!)
- `ndim`: Number of spatial dimensions

# Note
This function assumes normalize_statistics! has already been called.
The Reynolds stress is computed as: R_ij = <u_i * u_j> - <u_i> * <u_j>

For 2D: reynolds_stress[ip,1] = R_12 (u'v')
For 3D: reynolds_stress[ip,1] = R_12 (u'v'), [ip,2] = R_13 (u'w'), [ip,3] = R_23 (v'w')
"""
function compute_reynolds_stresses!(reynolds_stress, turb_stats::St_turbulence_stats, ndim)

    npoin = size(turb_stats.avvel, 1)

    if ndim == 2
        # 2D: Only u'v'
        for ip = 1:npoin
            mean_u = turb_stats.avvel[ip, 1]
            mean_v = turb_stats.avvel[ip, 2]
            mean_uv = turb_stats.avvex[ip, 1]

            # Reynolds stress: <u'v'> = <u*v> - <u>*<v>
            reynolds_stress[ip, 1] = mean_uv - mean_u * mean_v
        end

    elseif ndim == 3
        # 3D: u'v', u'w', v'w'
        for ip = 1:npoin
            mean_u = turb_stats.avvel[ip, 1]
            mean_v = turb_stats.avvel[ip, 2]
            mean_w = turb_stats.avvel[ip, 3]
            mean_uv = turb_stats.avvex[ip, 1]
            mean_uw = turb_stats.avvex[ip, 2]
            mean_vw = turb_stats.avvex[ip, 3]

            # Reynolds stresses
            reynolds_stress[ip, 1] = mean_uv - mean_u * mean_v  # <u'v'>
            reynolds_stress[ip, 2] = mean_uw - mean_u * mean_w  # <u'w'>
            reynolds_stress[ip, 3] = mean_vw - mean_v * mean_w  # <v'w'>
        end
    end
end

"""
    write_turbulence_stats_vtk(turb_stats, mesh, OUTPUT_DIR, filename, ndim, SD)

Write turbulence statistics to VTK file for visualization.

# Arguments
- `turb_stats`: Turbulence statistics structure (normalized)
- `mesh`: Mesh structure containing coordinates
- `OUTPUT_DIR`: Output directory path
- `filename`: Base filename for output (without extension)
- `ndim`: Number of spatial dimensions
- `SD`: Spatial dimension type
"""
function write_turbulence_stats_vtk(turb_stats::St_turbulence_stats, mesh, OUTPUT_DIR::String, filename::String, ndim, SD)

    if !turb_stats.l_enabled[]
        @warn "Turbulence statistics are not enabled. Skipping VTK output."
        return
    end

    println(" # Writing turbulence statistics to VTK: ", filename, ".vtu")

    # Note: VTK output would require WriteVTK.jl integration
    # This is a placeholder showing the structure
    # In practice, this would call write_output_vtk with the statistics arrays

    @info "VTK output function placeholder - integrate with existing write_output_vtk infrastructure"
end

"""
    write_turbulence_stats_ascii(turb_stats, mesh, OUTPUT_DIR, filename, ndim, t_current)

Write turbulence statistics to ASCII file for analysis.

# Arguments
- `turb_stats`: Turbulence statistics structure (normalized)
- `mesh`: Mesh structure containing coordinates
- `OUTPUT_DIR`: Output directory path
- `filename`: Base filename for output
- `ndim`: Number of spatial dimensions
- `t_current`: Current simulation time
"""
function write_turbulence_stats_ascii(turb_stats::St_turbulence_stats, mesh, OUTPUT_DIR::String, filename::String, ndim, t_current)

    if !turb_stats.l_enabled[]
        @warn "Turbulence statistics are not enabled. Skipping ASCII output."
        return
    end

    npoin = size(turb_stats.avvel, 1)
    fname = string(filename, "_stats_t", @sprintf("%.6f", t_current), ".dat")

    println(" # Writing turbulence statistics to ASCII: ", fname)

    open(string(OUTPUT_DIR, "/", fname), "w") do f
        # Write header
        @printf(f, "# Turbulence Statistics\n")
        @printf(f, "# Time: %.6f\n", t_current)
        @printf(f, "# Accumulated time: %.6f\n", turb_stats.acutim[])
        @printf(f, "# Number of points: %d\n", npoin)
        @printf(f, "# Spatial dimensions: %d\n", ndim)
        @printf(f, "#\n")

        if ndim == 1
            @printf(f, "# ip | x | <u> | <ρ> | <p> | <u*u> | <p*p> | <μ_eff>\n")
            for ip = 1:npoin
                @printf(f, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                    ip, mesh.x[ip],
                    turb_stats.avvel[ip,1], turb_stats.avrho[ip], turb_stats.avpress[ip],
                    turb_stats.avve2[ip,1], turb_stats.avpr2[ip], turb_stats.avmueff[ip])
            end

        elseif ndim == 2
            @printf(f, "# ip | x | y | <u> | <v> | <ρ> | <p> | <u*u> | <v*v> | <u*v> | <p*p> | <μ_eff>\n")
            for ip = 1:npoin
                @printf(f, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                    ip, mesh.x[ip], mesh.y[ip],
                    turb_stats.avvel[ip,1], turb_stats.avvel[ip,2],
                    turb_stats.avrho[ip], turb_stats.avpress[ip],
                    turb_stats.avve2[ip,1], turb_stats.avve2[ip,2],
                    turb_stats.avvex[ip,1], turb_stats.avpr2[ip], turb_stats.avmueff[ip])
            end

        elseif ndim == 3
            @printf(f, "# ip | x | y | z | <u> | <v> | <w> | <ρ> | <p> | <u*u> | <v*v> | <w*w> | <u*v> | <u*w> | <v*w> | <p*p> | <μ_eff>\n")
            for ip = 1:npoin
                @printf(f, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                    ip, mesh.x[ip], mesh.y[ip], mesh.z[ip],
                    turb_stats.avvel[ip,1], turb_stats.avvel[ip,2], turb_stats.avvel[ip,3],
                    turb_stats.avrho[ip], turb_stats.avpress[ip],
                    turb_stats.avve2[ip,1], turb_stats.avve2[ip,2], turb_stats.avve2[ip,3],
                    turb_stats.avvex[ip,1], turb_stats.avvex[ip,2], turb_stats.avvex[ip,3],
                    turb_stats.avpr2[ip], turb_stats.avmueff[ip])
            end
        end
    end

    println(" # Turbulence statistics written to: ", OUTPUT_DIR, "/", fname)
end

"""
    print_turbulence_stats_summary(turb_stats::St_turbulence_stats, ndim)

Print a summary of the turbulence statistics to the console.
"""
function print_turbulence_stats_summary(turb_stats::St_turbulence_stats, ndim)

    if !turb_stats.l_enabled[]
        println(" # Turbulence statistics: DISABLED")
        return
    end

    println(" # ========================================================")
    println(" # Turbulence Statistics Summary")
    println(" # ========================================================")
    println(" # Accumulated time: ", turb_stats.acutim[])
    println(" # Spatial dimensions: ", ndim)

    # Compute global statistics
    mean_rho = sum(turb_stats.avrho) / length(turb_stats.avrho)
    mean_p = sum(turb_stats.avpress) / length(turb_stats.avpress)

    println(" # Domain-averaged <ρ>: ", mean_rho)
    println(" # Domain-averaged <p>: ", mean_p)

    if ndim >= 1
        mean_u = sum(turb_stats.avvel[:, 1]) / size(turb_stats.avvel, 1)
        println(" # Domain-averaged <u>: ", mean_u)
    end
    if ndim >= 2
        mean_v = sum(turb_stats.avvel[:, 2]) / size(turb_stats.avvel, 1)
        println(" # Domain-averaged <v>: ", mean_v)
    end
    if ndim >= 3
        mean_w = sum(turb_stats.avvel[:, 3]) / size(turb_stats.avvel, 1)
        println(" # Domain-averaged <w>: ", mean_w)
    end

    println(" # ========================================================")
end
