"""
    ERA5 Data Structures for Jexpresso LES

    This module defines data structures and allocation functions for ERA5 reanalysis
    data integration with Jexpresso's physics solvers.

    Author: Jexpresso Development Team
    Date: 2025-11-30
"""

"""
    St_ERA5Fields{T, dims1, backend}

Structure to store interpolated ERA5 fields on the Jexpresso mesh.

Fields:
- `temperature`: Temperature field [K]
- `u_wind`: Zonal wind component [m/s]
- `v_wind`: Meridional wind component [m/s]
- `w_wind`: Vertical wind component [m/s]
- `specific_humidity`: Specific humidity [kg/kg]
- `pressure`: Pressure field [Pa]
- `theta`: Potential temperature [K] (derived)
- `rho`: Density [kg/m³] (derived)
"""
Base.@kwdef mutable struct St_ERA5Fields{T <: AbstractFloat, dims1, backend}
    # Primary ERA5 fields
    temperature       = KernelAbstractions.zeros(backend, T, dims1)
    u_wind            = KernelAbstractions.zeros(backend, T, dims1)
    v_wind            = KernelAbstractions.zeros(backend, T, dims1)
    w_wind            = KernelAbstractions.zeros(backend, T, dims1)
    specific_humidity = KernelAbstractions.zeros(backend, T, dims1)
    pressure          = KernelAbstractions.zeros(backend, T, dims1)

    # Derived thermodynamic fields
    theta             = KernelAbstractions.zeros(backend, T, dims1)
    rho               = KernelAbstractions.zeros(backend, T, dims1)
end

"""
    St_ERA5Forcing{T, dims1, backend}

Structure to store ERA5-derived forcing terms for LES.

Fields:
- `u_nudging`: Zonal wind nudging tendency [m/s²]
- `v_nudging`: Meridional wind nudging tendency [m/s²]
- `T_nudging`: Temperature nudging tendency [K/s]
- `q_nudging`: Humidity nudging tendency [kg/(kg·s)]
- `nudging_timescale`: Relaxation timescale [s]
"""
Base.@kwdef mutable struct St_ERA5Forcing{T <: AbstractFloat, dims1, backend}
    u_nudging         = KernelAbstractions.zeros(backend, T, dims1)
    v_nudging         = KernelAbstractions.zeros(backend, T, dims1)
    T_nudging         = KernelAbstractions.zeros(backend, T, dims1)
    q_nudging         = KernelAbstractions.zeros(backend, T, dims1)
    nudging_timescale::T = 3600.0  # Default: 1 hour
end

"""
    allocate_ERA5Fields(npoin, mesh, inputs, T, backend; lERA5=false)

Allocate ERA5 field structure and optionally load data.

# Arguments
- `npoin`: Number of mesh points
- `mesh`: Jexpresso mesh structure
- `inputs`: Input parameters dictionary
- `T`: Float type (Float32 or Float64)
- `backend`: Computation backend (CPU() or CUDABackend())
- `lERA5`: Flag to enable ERA5 data loading (default: false)

# Returns
- `St_ERA5Fields`: Allocated ERA5 fields structure
"""
function allocate_ERA5Fields(npoin, mesh, inputs, T, backend; lERA5=false)

    if lERA5
        dims1 = (Int64(npoin), 1)
    else
        dims1 = (Int64(1))
    end

    era5_fields = St_ERA5Fields{T, dims1, backend}()

    if lERA5
        println("Loading ERA5 initial conditions...")
        load_era5_initial_conditions!(backend, inputs, era5_fields, mesh)
    end

    return era5_fields
end

"""
    allocate_ERA5Forcing(npoin, inputs, T, backend; lERA5_forcing=false)

Allocate ERA5 forcing structure.

# Arguments
- `npoin`: Number of mesh points
- `inputs`: Input parameters dictionary
- `T`: Float type (Float32 or Float64)
- `backend`: Computation backend (CPU() or CUDABackend())
- `lERA5_forcing`: Flag to enable ERA5 forcing (default: false)

# Returns
- `St_ERA5Forcing`: Allocated ERA5 forcing structure
"""
function allocate_ERA5Forcing(npoin, inputs, T, backend; lERA5_forcing=false)

    if lERA5_forcing
        dims1 = (Int64(npoin), 1)
    else
        dims1 = (Int64(1))
    end

    era5_forcing = St_ERA5Forcing{T, dims1, backend}()

    if lERA5_forcing && haskey(inputs, :era5_nudging_timescale)
        era5_forcing.nudging_timescale = T(inputs[:era5_nudging_timescale])
    end

    return era5_forcing
end

"""
    load_era5_initial_conditions!(backend, inputs, era5_fields, mesh)

Load ERA5 data and interpolate to mesh for initial conditions.

# Arguments
- `backend`: Computation backend
- `inputs`: Input parameters dictionary with ERA5 configuration
- `era5_fields`: ERA5 fields structure to populate
- `mesh`: Jexpresso mesh structure

# Required inputs dictionary keys:
- `:era5_file`: Path to ERA5 NetCDF file
- `:era5_target_lat`: Target latitude for data extraction
- `:era5_target_lon`: Target longitude for data extraction
- `:era5_time_index`: Time index to read (default: 1)
"""
function load_era5_initial_conditions!(backend, inputs, era5_fields, mesh)

    # Check required inputs
    if !haskey(inputs, :era5_file)
        error("ERA5 file path not specified in inputs[:era5_file]")
    end
    if !haskey(inputs, :era5_target_lat)
        error("ERA5 target latitude not specified in inputs[:era5_target_lat]")
    end
    if !haskey(inputs, :era5_target_lon)
        error("ERA5 target longitude not specified in inputs[:era5_target_lon]")
    end

    era5_file = inputs[:era5_file]
    target_lat = inputs[:era5_target_lat]
    target_lon = inputs[:era5_target_lon]
    time_index = haskey(inputs, :era5_time_index) ? inputs[:era5_time_index] : 1

    # Optional: custom variable names
    var_names = haskey(inputs, :era5_var_names) ? inputs[:era5_var_names] : Dict()

    # Optional: pressure level selection
    pressure_levels = haskey(inputs, :era5_pressure_levels) ? inputs[:era5_pressure_levels] : nothing

    # Optional: lat/lon ranges for subsetting
    lat_range = haskey(inputs, :era5_lat_range) ? inputs[:era5_lat_range] : (nothing, nothing)
    lon_range = haskey(inputs, :era5_lon_range) ? inputs[:era5_lon_range] : (nothing, nothing)

    # Read ERA5 data from NetCDF
    println("Reading ERA5 data from: $era5_file")
    era5_data = read_era5_netcdf(era5_file;
                                 var_names=var_names,
                                 time_index=time_index,
                                 lat_range=lat_range,
                                 lon_range=lon_range,
                                 pressure_levels=pressure_levels)

    # Interpolate to mesh
    println("Interpolating ERA5 data to mesh at lat=$target_lat°, lon=$target_lon°")
    era5_mesh = interpolate_era5_to_mesh(backend, era5_data, mesh, target_lat, target_lon)

    # Copy to ERA5 fields structure
    if backend == CPU()
        era5_fields.temperature .= era5_mesh["temperature"]
        era5_fields.u_wind .= era5_mesh["u_wind"]
        era5_fields.v_wind .= era5_mesh["v_wind"]
        era5_fields.w_wind .= era5_mesh["w_wind"]
        era5_fields.specific_humidity .= era5_mesh["specific_humidity"]
        era5_fields.pressure .= era5_mesh["pressure"]
    else
        KernelAbstractions.copyto!(backend, era5_fields.temperature, era5_mesh["temperature"])
        KernelAbstractions.copyto!(backend, era5_fields.u_wind, era5_mesh["u_wind"])
        KernelAbstractions.copyto!(backend, era5_fields.v_wind, era5_mesh["v_wind"])
        KernelAbstractions.copyto!(backend, era5_fields.w_wind, era5_mesh["w_wind"])
        KernelAbstractions.copyto!(backend, era5_fields.specific_humidity, era5_mesh["specific_humidity"])
        KernelAbstractions.copyto!(backend, era5_fields.pressure, era5_mesh["pressure"])
    end

    # Compute derived fields
    compute_era5_derived_fields!(backend, era5_fields, mesh.npoin)

    println("ERA5 initial conditions loaded successfully!")
end

"""
    compute_era5_derived_fields!(backend, era5_fields, npoin)

Compute derived thermodynamic fields from ERA5 data.

# Arguments
- `backend`: Computation backend
- `era5_fields`: ERA5 fields structure
- `npoin`: Number of points
"""
function compute_era5_derived_fields!(backend, era5_fields, npoin)

    # Physical constants
    R_dry = 287.05    # Gas constant for dry air [J/(kg·K)]
    p0 = 100000.0     # Reference pressure [Pa]
    κ = R_dry / 1004.0  # Poisson constant (R/cp)

    if backend == CPU()
        for i in 1:npoin
            # Potential temperature: θ = T * (p0/p)^κ
            era5_fields.theta[i] = era5_fields.temperature[i] *
                                   (p0 / era5_fields.pressure[i])^κ

            # Density from ideal gas law: ρ = p / (R * T)
            era5_fields.rho[i] = era5_fields.pressure[i] /
                                (R_dry * era5_fields.temperature[i])
        end
    else
        # GPU kernel would go here
        # For now, copy to CPU, compute, copy back
        T_cpu = Array(era5_fields.temperature)
        p_cpu = Array(era5_fields.pressure)
        theta_cpu = zeros(npoin)
        rho_cpu = zeros(npoin)

        for i in 1:npoin
            theta_cpu[i] = T_cpu[i] * (p0 / p_cpu[i])^κ
            rho_cpu[i] = p_cpu[i] / (R_dry * T_cpu[i])
        end

        KernelAbstractions.copyto!(backend, era5_fields.theta, theta_cpu)
        KernelAbstractions.copyto!(backend, era5_fields.rho, rho_cpu)
    end
end

"""
    apply_era5_nudging!(q, era5_fields, era5_forcing, dt, ::PERT)

Apply ERA5 nudging/relaxation forcing to the solution.

# Arguments
- `q`: Solution state vector
- `era5_fields`: ERA5 reference fields
- `era5_forcing`: ERA5 forcing structure with nudging timescale
- `dt`: Time step
- `::PERT`: Perturbation formulation type

# Notes
This function applies Newtonian relaxation:
    dq/dt = -(q - q_ERA5) / τ
where τ is the nudging timescale.
"""
function apply_era5_nudging!(q, era5_fields, era5_forcing, dt, ::PERT)

    τ = era5_forcing.nudging_timescale

    # Apply nudging to momentum components
    # q[2] = ρu, q[3] = ρv, q[4] = ρw
    # Nudge toward ERA5 winds: u_ERA5, v_ERA5, w_ERA5

    # This is a simplified version - actual implementation should
    # consider the specific formulation used in Jexpresso

    npoin = length(era5_fields.u_wind)

    for i in 1:npoin
        # Compute target momentum from ERA5 winds
        ρ_target = era5_fields.rho[i]
        u_target = era5_fields.u_wind[i]
        v_target = era5_fields.v_wind[i]
        w_target = era5_fields.w_wind[i]

        # Apply nudging tendency
        era5_forcing.u_nudging[i] = -(q[i,2]/q[i,1] - u_target) / τ
        era5_forcing.v_nudging[i] = -(q[i,3]/q[i,1] - v_target) / τ
        era5_forcing.T_nudging[i] = -(era5_fields.temperature[i]) / τ  # Placeholder
        era5_forcing.q_nudging[i] = -(q[i,6]/q[i,1] - era5_fields.specific_humidity[i]) / τ
    end
end
