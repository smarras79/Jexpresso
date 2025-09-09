# Physical constants
const BETA_H = 5.0         # stability function parameter
const BETA_M = 5.0         # stability function parameter
const GAMMA_H = 16.0       # stability function parameter
const GAMMA_M = 16.0       # stability function parameter
const ZINT = 10.0          # reference height for interpolation (m)

# Add convergence and stability parameters
const MAX_ITER = 15        # Increased from 10
const CONVERGENCE_TOL = 1e-4  # Tighter convergence tolerance
const RELAXATION_FACTOR = 0.5  # More conservative than 0.7
const MIN_MOL = 1.0        # Minimum |L| to avoid division issues
const MAX_MOL = 2000.0     # Maximum |L| for numerical stability

# Wind averaging parameters
const WIND_AVERAGING_TIME = 200.0  # seconds

# Structure to store wind velocity time averaging data
mutable struct WindAverageData
    u_sum::Float64
    v_sum::Float64
    w_sum::Float64
    count::Int64
    last_update_time::Float64
    
    # Constructor
    WindAverageData() = new(0.0, 0.0, 0.0, 0, 0.0)
end

# Structure to store MOL and wind averaging data
mutable struct SurfaceLayerStorage
    previous_mol::Dict{Tuple{Int,Int,Int}, Float64}  # Store MOL by grid location
    wind_averages::Dict{Tuple{Int,Int,Int}, WindAverageData}  # Store wind averages by grid location
    time_series::Dict{Tuple{Int,Int,Int}, Vector{Tuple{Float64, Float64, Float64, Float64}}}  # (time, u, v, w)
end

# Global instance (you might want to pass this as a parameter instead)
const surface_storage = SurfaceLayerStorage(
    Dict{Tuple{Int,Int,Int}, Float64}(),
    Dict{Tuple{Int,Int,Int}, WindAverageData}(),
    Dict{Tuple{Int,Int,Int}, Vector{Tuple{Float64, Float64, Float64, Float64}}}()
)

"""
Update running wind average using exponential moving average or sliding window
"""
function update_wind_average!(grid_key::Tuple{Int,Int,Int}, u_vel::Real, v_vel::Real, w_vel::Real, 
                              current_time::Real, dt::Real)
    
    # Get or create wind average data for this grid point
    if !haskey(surface_storage.wind_averages, grid_key)
        surface_storage.wind_averages[grid_key] = WindAverageData()
        surface_storage.time_series[grid_key] = Vector{Tuple{Float64, Float64, Float64, Float64}}()
    end
    
    wind_data = surface_storage.wind_averages[grid_key]
    time_series = surface_storage.time_series[grid_key]
    
    # Add current data point
    push!(time_series, (current_time, u_vel, v_vel, w_vel))
    
    # Remove data older than WIND_AVERAGING_TIME
    cutoff_time = current_time - WIND_AVERAGING_TIME
    filter!(x -> x[1] >= cutoff_time, time_series)
    
    # Calculate average from remaining data
    if !isempty(time_series)
        n_points = length(time_series)
        u_avg = sum(x[2] for x in time_series) / n_points
        v_avg = sum(x[3] for x in time_series) / n_points
        w_avg = sum(x[4] for x in time_series) / n_points
        
        # Update the average data structure
        wind_data.u_sum = u_avg
        wind_data.v_sum = v_avg
        wind_data.w_sum = w_avg
        wind_data.count = n_points
        wind_data.last_update_time = current_time
        
        return u_avg, v_avg, w_avg
    else
        # If no data available, return current values
        return u_vel, v_vel, w_vel
    end
end

"""
Alternative: Exponential moving average implementation (more memory efficient)
"""
function update_wind_average_exponential!(grid_key::Tuple{Int,Int,Int}, u_vel::Real, v_vel::Real, w_vel::Real, 
                                          current_time::Real, dt::Real)
    
    # Time constant for exponential averaging (tau = WIND_AVERAGING_TIME / 3 for ~95% contribution)
    tau = WIND_AVERAGING_TIME / 3.0
    alpha = dt / (tau + dt)  # Smoothing factor
    
    # Get or create wind average data for this grid point
    if !haskey(surface_storage.wind_averages, grid_key)
        surface_storage.wind_averages[grid_key] = WindAverageData()
        # Initialize with current values
        wind_data = surface_storage.wind_averages[grid_key]
        wind_data.u_sum = u_vel
        wind_data.v_sum = v_vel
        wind_data.w_sum = w_vel
        wind_data.count = 1
        wind_data.last_update_time = current_time
        return u_vel, v_vel, w_vel
    end
    
    wind_data = surface_storage.wind_averages[grid_key]
    
    # Exponential moving average update
    wind_data.u_sum = (1.0 - alpha) * wind_data.u_sum + alpha * u_vel
    wind_data.v_sum = (1.0 - alpha) * wind_data.v_sum + alpha * v_vel
    wind_data.w_sum = (1.0 - alpha) * wind_data.w_sum + alpha * w_vel
    wind_data.count += 1
    wind_data.last_update_time = current_time
    
    return wind_data.u_sum, wind_data.v_sum, wind_data.w_sum
end

function MOST!(τ_f, wθ,
               lwall_model,
               uprimitiveieq,
               inputs,
               PhysConst, MPConst,
               dist, 
               iel, ieq,
               connijk,
               coords,
               poin_in_bdy_face, elem_to_face, bdy_face_type,
               k, l, m, iface_bdy, idx1, idx2;
               mol_prev::Union{Float64, Nothing} = nothing,  # Add optional previous MOL
               current_time::Real = 0.0,                    # Current simulation time
               dt::Real = 1.0,                              # Time step
               use_exponential_avg::Bool = true)            # Choose averaging method

    # Get physical constants    
    Rair   = PhysConst.Rair # dry air gas constant
    cp     = PhysConst.cp    # specific heat at constant pressure
    GRAV   = PhysConst.g
    KARMAN = PhysConst.karman
    LV     = MPConst.Lc      # latent heat of vaporization
    
    # Get boundary point information
    ip = connijk[iel, k, l, m]
    iface_bdy = elem_to_face[iel, k, l, m, 1]

    # For surface boundary, use second point above surface (iwallnode = 2)
    wall_z    = coords[ip, 3]  # z-coordinate of surface point
    iwallnode = 2              # second point above surface
    ip2       = connijk[iel, k, l, iwallnode]
    z2        = coords[ip2, 3]
    z_level   = z2 - wall_z    # height above surface
    
    # Extract atmospheric variables at second grid point
    ρ_atm = uprimitiveieq[k, l, iwallnode, 1]
    u_vel_instant = uprimitiveieq[k, l, iwallnode, 2]  # Instantaneous values
    v_vel_instant = uprimitiveieq[k, l, iwallnode, 3]
    w_vel_instant = uprimitiveieq[k, l, iwallnode, 4]
    θ_atm = uprimitiveieq[k, l, iwallnode, 5]
    θ_skn = uprimitiveieq[k, l, 1, 5]    
    q_atm = 0.0 #uprimitiveieq[k, l, iwallnode, 6]  # specific humidity (assuming index 6)

    # Update time-averaged wind velocities
    grid_key = (k, l, m)
    
    if use_exponential_avg
        u_vel, v_vel, w_vel = update_wind_average_exponential!(grid_key, u_vel_instant, v_vel_instant, w_vel_instant, 
                                                               current_time, dt)
    else
        u_vel, v_vel, w_vel = update_wind_average!(grid_key, u_vel_instant, v_vel_instant, w_vel_instant, 
                                                   current_time, dt)
    end
    
    # Debug output (optional)
    if current_time > 0.0 && mod(current_time, 10.0) < dt  # Print every 10 seconds
        instant_wspd = sqrt(u_vel_instant^2 + v_vel_instant^2)
        avg_wspd = sqrt(u_vel^2 + v_vel^2)
        println("Grid ($k,$l,$m): Instant wspd = $(round(instant_wspd, digits=2)) m/s, " *
               "Avg wspd = $(round(avg_wspd, digits=2)) m/s")
    end

    p_atm = 100000.0 #initialize to be modified by perfectGasLaw_ρθtoP!(p_atm, ...)
    T_atm = θ_atm    #initialize to be modified by perfectGasLaw_θPtoT!(T_atm, ...)
    perfectGasLaw_ρθtoP!(p_atm, PhysConst; ρ=ρ_atm, θ=θ_atm)
    perfectGasLaw_θPtoT!(T_atm, PhysConst; θ=θ_atm, Press=p_atm)

    # Surface (skin) quantities:
    z0_local   = 0.1 #clamp(rand(dist), 0.0, 0.1) # momentum roughness
    z0t_local  = z0_local*0.1                # thermal roughness
    tsk_local  = θ_skn                       # skin temperature
    qsfc_local =   0.0                       # skin specific humidity
    
    # Calculate wind speed magnitude using time-averaged velocities
    wspd = sqrt(u_vel^2 + v_vel^2)
    
    # Calculate air density
    ρ_air = p_atm / (Rair * T_atm)

    # Get previous MOL value for this grid location
    mol_previous = mol_prev
    if mol_previous === nothing
        # Try to get from storage, or use default for first time step
        mol_previous = get(surface_storage.previous_mol, grid_key, 100.0)  # Default neutral value
    end

    # Call similarity theory routine (non-iterative version using previous MOL)
    # Use time-averaged wind velocities
    result = surface_layer_similarity_prev_mol(z_level, z0_local, z0t_local,
                                               T_atm, q_atm,
                                               u_vel, v_vel, p_atm,  # Time-averaged velocities
                                               tsk_local, qsfc_local,
                                               mol_previous,  # Use previous MOL
                                               PhysConst)

    # Store current MOL for next time step
    surface_storage.previous_mol[grid_key] = result.mol

    # Validate MOST output before proceeding
    #if !validate_most_output(result.ustar, result.tstar, result.mol, z_level, z0_local)
    #    println("WARNING: MOST output failed validation checks")
    #    # Could implement fallback to neutral conditions here
    #end

    # Calculate exchange coefficients
    coeffs = calculate_exchange_coefficients!(result.ustar,
                                              result.tstar,
                                              result.qstar, 
                                              result.wspd, result.psim, result.psih, 
                                              z_level, z0_local, z0t_local, PhysConst)
    
    ### MATIAS USA: L del paso anterior - IMPLEMENTED! ###
    ### Wind velocity time averaging - IMPLEMENTED! ###
    #### u* = ln(1+z/z0) - perturbation
    
    # Calculate surface fluxes using time-averaged wind velocities
    # Sensible heat flux (W/m²)
    hfx_local = ρ_air * cp * coeffs.ch * result.wspd * (tsk_local - T_atm)
    
    # Momentum fluxes (N/m²) - use time-averaged velocities
    taux_local = ρ_air * coeffs.cd * result.wspd * u_vel  # Time-averaged u_vel
    tauy_local = ρ_air * coeffs.cd * result.wspd * v_vel  # Time-averaged v_vel
    
    # Apply fluxes with appropriate signs
    # Momentum: negative because surface drag opposes wind
    τ_f[iface_bdy, idx1, idx2, 1] = -taux_local  # ✓ Correct for drag
    τ_f[iface_bdy, idx1, idx2, 2] = -tauy_local  # ✓ Correct for drag
    
    # Heat: sign depends on your flux convention
    # If wθ represents "flux TO atmosphere FROM surface", use positive when surface is warmer
    # If wθ represents "flux TO surface FROM atmosphere", use as calculated
    wθ[iface_bdy, idx1, idx2, 1] = hfx_local  # Check this based on your convention
    
end

"""
Non-iterative surface_layer_similarity using previous time step MOL
"""
function surface_layer_similarity_prev_mol(z::Real, z0::Real, z0t_in::Real, ta::Real, qa::Real, 
                                           ua::Real, va::Real, ps::Real, tsk::Real, qsfc::Real,
                                           mol_prev::Real,  # Previous time step MOL
                                           PhysConst)

    GRAV   = PhysConst.g
    KARMAN = PhysConst.karman
    cp     = PhysConst.cp
    Rair   = PhysConst.Rair
    
    # Calculate basic variables
    wspd = sqrt(ua^2 + va^2)
    wspd = max(wspd, 0.1)  # minimum wind speed
    
    # Calculate potential temperature
    th = ta * (100000.0/ps)^(Rair/cp)
    
    # Convert specific humidity to mixing ratio
    qv = qa / (1.0 - qa)
    
    # Virtual potential temperature
    tv   = th * (1.0 + 0.61*qv)
    tvs  = tsk * (1.0 + 0.61*qsfc) * (100000.0/ps)^(Rair/cp)
    dthv = tv - tvs
    
    # Enhanced z0t calculation
    z0t = z0t_in
    if z0t <= 0.0
        czil = 0.1
        ustar_init = wspd * KARMAN / log(z/z0)  # Better initial guess
        z0t = z0 * exp(-KARMAN * czil * sqrt(ustar_init * z0 / 1.5e-5))
        z0t = clamp(z0t, 1.0e-7, z0)  # Ensure reasonable bounds
    end
    
    # Use previous MOL (with bounds checking)
    mol = clamp(mol_prev, -MAX_MOL, MAX_MOL)
    
    # Ensure minimum MOL magnitude
    if abs(mol) < MIN_MOL
        mol = sign(mol) * MIN_MOL
    end
    
    # Calculate stability parameters using previous MOL
    zeta   = z / mol
    zeta0  = z0 / mol
    zeta0t = z0t / mol

    # Apply stability parameter limits
    zeta   = clamp(zeta, -5.0, 5.0)
    zeta0  = clamp(zeta0, -5.0, 5.0)
    zeta0t = clamp(zeta0t, -5.0, 5.0)
    
    # Calculate stability functions
    if zeta >= 0.0
        # Stable conditions
        psim_stable = -min(BETA_M * zeta, 10.0)
        psih_stable = -min(BETA_H * zeta, 10.0)
        psim        = psim_stable - (-min(BETA_M * zeta0, 10.0))
        psih        = psih_stable - (-min(BETA_H * zeta0t, 10.0))
        psiq        = psih
    else
        # Unstable conditions
        zeta_bounded   = max(zeta, -2.0)
        zeta0_bounded  = max(zeta0, -2.0)
        zeta0t_bounded = max(zeta0t, -2.0)
        
        x  = (1.0 - GAMMA_M * zeta_bounded)^0.25
        x0 = (1.0 - GAMMA_M * zeta0_bounded)^0.25
        
        # Check for numerical issues
        if x <= 0.0 || x0 <= 0.0
            # Fallback to near-neutral
            psim = -0.1 * zeta
            psih = -0.1 * zeta
        else
            psim_unstable  = 2.0*log((1.0+x)/2.0) + log((1.0+x^2)/2.0) - 2.0*atan(x) + π/2
            psih_unstable  = 2.0*log((1.0+sqrt(1.0-GAMMA_H*zeta_bounded))/2.0)
            
            psim0_unstable = 2.0*log((1.0+x0)/2.0) + log((1.0+x0^2)/2.0) - 2.0*atan(x0) + π/2
            psih0_unstable = 2.0*log((1.0+sqrt(1.0-GAMMA_H*zeta0t_bounded))/2.0)
            
            psim = psim_unstable - psim0_unstable
            psih = psih_unstable - psih0_unstable
        end
        psiq = psih
    end
    
    # Calculate friction velocity
    denom = log(z/z0) - psim
    if abs(denom) < 0.1
        ustar = 0.1  # Fallback value
    else
        ustar = wspd * KARMAN / denom
        ustar = clamp(ustar, 0.01, 5.0)  # Reasonable bounds for u*
    end
    
    # Calculate temperature scale
    denom_h = log(z/z0t) - psih
    if abs(denom_h) < 0.1
        tstar = 0.01 * sign(th - tsk)  # Small fallback value
    else
        tstar = KARMAN * (th - tsk) / denom_h
        tstar = clamp(tstar, -2.0, 2.0)  # Reasonable bounds for T*
    end
    
    # Calculate moisture scale
    qstar = KARMAN * (qv - qsfc) / (log(z/z0t) - psiq)
    qstar = clamp(qstar, -0.01, 0.01)  # Reasonable bounds for q*
    
    # Optional: Update MOL for consistency check (but don't iterate)
    # This gives you a "predicted" MOL for the next time step if desired
    if abs(tstar) > 1.0e-6
        mol_updated = -tv * ustar^3 / (KARMAN * GRAV * tstar)
        mol_updated = clamp(mol_updated, -MAX_MOL, MAX_MOL)
        
        # You could use mol_updated as the MOL for the next time step
        # or blend it with the current value for stability
        mol = 0.9 * mol + 0.1 * mol_updated  # Light blending for stability
    end
    
    return (ustar=ustar, tstar=tstar, qstar=qstar, mol=mol, 
            psim=psim, psih=psih, psiq=psiq, wspd=wspd, th=th, qv=qv)
end

"""
Get current wind averages for a grid point (for debugging/monitoring)
"""
function get_wind_averages(k::Int, l::Int, m::Int)
    grid_key = (k, l, m)
    if haskey(surface_storage.wind_averages, grid_key)
        wind_data = surface_storage.wind_averages[grid_key]
        return wind_data.u_sum, wind_data.v_sum, wind_data.w_sum, wind_data.count
    else
        return 0.0, 0.0, 0.0, 0
    end
end

"""
Reset wind averages (useful for restart simulations)
"""
function reset_wind_averages!()
    empty!(surface_storage.wind_averages)
    empty!(surface_storage.time_series)
    empty!(surface_storage.previous_mol)
end

"""
Validate MOST output for physical reasonableness
"""
function validate_most_output(ustar::Real, tstar::Real, mol::Real, z::Real, z0::Real)
    checks = []
    
    # Check friction velocity bounds
    push!(checks, 0.01 < ustar < 5.0)
    
    # Check temperature scale bounds  
    push!(checks, abs(tstar) < 5.0)
    
    # Check MOL bounds and finite
    push!(checks, isfinite(mol) && abs(mol) > MIN_MOL && abs(mol) < MAX_MOL)
    
    # Check stability parameter bounds
    zeta = z / mol
    push!(checks, abs(zeta) < 10.0)
    
    # Check height ratio reasonableness
    push!(checks, z > 2.0 * z0)
    
    return all(checks)
end

"""
Calculate bulk Richardson number with enhanced error handling
"""
function bulk_richardson(z::Real, ta::Real, tsk::Real, ua::Real, va::Real, ps::Real, PhysConst)
    wspd = sqrt(ua^2 + va^2)
    wspd = max(wspd, 0.1)
    
    th  = ta * (100000.0/ps)^(PhysConst.Rair/PhysConst.cp)
    ths = tsk * (100000.0/ps)^(PhysConst.Rair/PhysConst.cp)
    dth = th - ths
    
    rib = (PhysConst.g * z * dth) / (th * wspd^2)
    
    # Bound Richardson number to reasonable range
    return clamp(rib, -10.0, 10.0)
end

"""
Calculate exchange coefficients with enhanced bounds checking
"""
function calculate_exchange_coefficients!(ustar::Real, tstar::Real, qstar::Real, wspd::Real,
                                          psim::Real, psih::Real, z::Real, z0::Real, z0t::Real, PhysConst)
    
    # Drag coefficient for momentum with bounds checking
    denom_m = log(z/z0) - psim
    if abs(denom_m) < 0.1
        cd = 0.001  # Fallback minimum value
    else
        cd = (PhysConst.karman / denom_m)^2
        cd = clamp(cd, 0.0001, 0.01)  # Reasonable bounds
    end
    
    # Exchange coefficient for heat with bounds checking
    denom_h = log(z/z0t) - psih
    if abs(denom_h) < 0.1 || wspd < 0.1
        ch = 0.001  # Fallback minimum value
    else
        ch = PhysConst.karman * ustar / (wspd * denom_h)
        ch = clamp(ch, 0.0001, 0.01)  # Reasonable bounds
    end
    
    cq = ch  # assuming same roughness length for heat and moisture
    
    return (cd=cd, ch=ch, cq=cq)
end
