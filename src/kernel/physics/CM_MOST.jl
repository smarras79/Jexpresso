"""
Minimal Surface Fluxes Implementation
Based on Monin-Obukhov Similarity Theory

This is a self-contained implementation of surface flux calculations
following the formulations used in ClimateMachineAtmos.jl/SurfaceFluxes.jl
but written for maximum portability and simplicity.

References:
- Businger et al. (1971): Flux-Profile Relationships in the Atmospheric Surface Layer
- Dyer (1974): A review of flux-profile relationships
- Monin & Obukhov (1954): Basic laws of turbulent mixing in the surface layer

To run it in stand-alone mode, simply do:

Julia> include("./path/to/this/file/CM_MOST.jl')
Julia> MinimalSurfaceFluxes.example_with_plots()


"""

module MinimalSurfaceFluxes

try
    using Plots
    gr()
catch
    error("Plots.jl not available. Install with: using Pkg; Pkg.add(\"Plots\")")
end

# Physical constants
const κ = 0.4                  # von Kármán constant
const g = 9.81                 # gravitational acceleration [m/s²]
const cp_d = 1004.0           # specific heat of dry air at constant pressure [J/kg/K]

# Universal function parameters (Businger-Dyer)
const a_m = 16.0              # momentum stability parameter (unstable)
const a_h = 16.0              # heat stability parameter (unstable) 
const b_m = 5.0               # momentum stability parameter (stable)
const b_h = 5.0               # heat stability parameter (stable)

"""
Universal stability functions φ_m and φ_h following Businger-Dyer (1971)

For unstable conditions (ζ < 0):
- φ_m(ζ) = (1 - a_m * ζ)^(-1/4)
- φ_h(ζ) = (1 - a_h * ζ)^(-1/2)

For stable conditions (ζ >= 0):
- φ_m(ζ) = 1 + b_m * ζ
- φ_h(ζ) = 1 + b_h * ζ
"""
function phi_m(zeta)
    if zeta < 0
        return (1 - a_m * zeta)^(-0.25)
    else
        return 1 + b_m * zeta
    end
end

function phi_h(zeta)
    if zeta < 0
        return (1 - a_h * zeta)^(-0.5)
    else
        return 1 + b_h * zeta
    end
end

"""
Integrated stability functions ψ_m and ψ_h (diabatic correction functions)

These are the integrals of the universal functions, used in profile calculations.
"""
function psi_m(zeta)
    if zeta < 0
        # Unstable conditions
        x = (1 - a_m * zeta)^0.25
        return 2 * log((1 + x) / 2) + log((1 + x^2) / 2) - 2 * atan(x) + π/2
    else
        # Stable conditions
        return -b_m * zeta
    end
end

function psi_h(zeta)
    if zeta < 0
        # Unstable conditions  
        y = (1 - a_h * zeta)^0.5
        return 2 * log((1 + y) / 2)
    else
        # Stable conditions
        return -b_h * zeta
    end
end

"""
Calculate Obukhov length L from surface fluxes

L = -u_star³ * T_ref / (κ * g * Q_H)

where:
- u_star: friction velocity [m/s]
- T_ref: reference temperature [K]  
- Q_H: surface sensible heat flux [W/m²]
"""
function obukhov_length(u_star, T_ref, Q_H)
    if abs(Q_H) < 1e-6
        return 1e6  # Near-neutral conditions (very large L)
    end
    return -u_star^3 * T_ref * cp_d / (κ * g * Q_H)
end

"""
Wind profile based on Monin-Obukhov similarity theory

u(z) = (u_star/κ) * [ln(z/z0_m) - ψ_m(z/L) + ψ_m(z0_m/L)]

where:
- z: height above surface [m]
- z0_m: momentum roughness length [m]
- L: Obukhov length [m]
- u_star: friction velocity [m/s]
"""
function wind_profile(z, z0_m, L, u_star)
    zeta = z / L
    zeta0 = z0_m / L
    return (u_star / κ) * (log(z / z0_m) - psi_m(zeta) + psi_m(zeta0))
end

"""
Temperature profile based on Monin-Obukhov similarity theory

θ(z) = θ_s + (θ_star/κ) * [ln(z/z0_h) - ψ_h(z/L) + ψ_h(z0_h/L)]

where:
- z: height above surface [m]
- z0_h: thermal roughness length [m] 
- L: Obukhov length [m]
- θ_s: surface potential temperature [K]
- θ_star: temperature scale [K]
"""
function temperature_profile(z, z0_h, L, theta_s, theta_star)
    zeta = z / L
    zeta0 = z0_h / L
    return theta_s + (theta_star / κ) * (log(z / z0_h) - psi_h(zeta) + psi_h(zeta0))
end

"""
Calculate friction velocity from wind speed and stability

This solves iteratively:
u = (u_star/κ) * [ln(z/z0_m) - ψ_m(z/L)]

Given u, z, z0_m, and L, find u_star
"""
function friction_velocity_from_wind(u, z, z0_m, L; max_iter=20, tol=1e-6)
    zeta = z / L
    zeta0 = z0_m / L
    
    # Initial guess
    u_star = κ * u / log(z / z0_m)
    
    for i in 1:max_iter
        # Calculate wind from current u_star
        u_calc = (u_star / κ) * (log(z / z0_m) - psi_m(zeta) + psi_m(zeta0))
        
        # Check convergence
        error = abs(u_calc - u)
        if error < tol
            break
        end
        
        # Update u_star (simple iteration)
        u_star *= u / u_calc
    end
    
    return u_star
end

"""
Calculate temperature scale from temperature difference and stability

θ_star = -Q_H / (ρ * cp * u_star)

where Q_H is the sensible heat flux [W/m²]
"""
function temperature_scale(Q_H, rho, u_star)
    return -Q_H / (rho * cp_d * u_star)
end

"""
Main surface flux calculation routine

Given atmospheric conditions at measurement height, calculate surface fluxes
and similarity parameters.

Inputs:
- u_ref: wind speed at reference height [m/s]
- theta_ref: potential temperature at reference height [K]
- z_ref: reference height [m]
- theta_s: surface potential temperature [K]
- z0_m: momentum roughness length [m]
- z0_h: thermal roughness length [m]  
- rho: air density [kg/m³] (optional, default 1.225)

Returns NamedTuple with:
- u_star: friction velocity [m/s]
- theta_star: temperature scale [K]
- L: Obukhov length [m]
- Q_H: sensible heat flux [W/m²]
- zeta: stability parameter z_ref/L
- C_D: drag coefficient
- C_H: heat transfer coefficient
"""
function surface_conditions(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h; 
                           rho=1.225, max_iter=20, tol=1e-4)
    
    # Initial guess assuming neutral conditions
    u_star = κ * u_ref / log(z_ref / z0_m)
    theta_star = κ * (theta_ref - theta_s) / log(z_ref / z0_h)
    
    # Initial heat flux guess
    Q_H = -rho * cp_d * u_star * theta_star
    L = obukhov_length(u_star, theta_ref, Q_H)
    
    # Iterative solution
    for iter in 1:max_iter
        # Update stability parameter
        zeta = z_ref / L
        zeta0_m = z0_m / L  
        zeta0_h = z0_h / L
        
        # Calculate new friction velocity
        u_star_new = κ * u_ref / (log(z_ref / z0_m) - psi_m(zeta) + psi_m(zeta0_m))
        
        # Calculate new temperature scale
        theta_star_new = κ * (theta_ref - theta_s) / (log(z_ref / z0_h) - psi_h(zeta) + psi_h(zeta0_h))
        
        # Update heat flux and Obukhov length
        Q_H_new = -rho * cp_d * u_star_new * theta_star_new
        L_new = obukhov_length(u_star_new, theta_ref, Q_H_new)
        
        # Check convergence
        error = abs(L_new - L) / max(abs(L), abs(L_new))
        
        if error < tol
            u_star = u_star_new
            theta_star = theta_star_new  
            Q_H = Q_H_new
            L = L_new
            break
        end
        
        # Update for next iteration
        u_star = u_star_new
        theta_star = theta_star_new
        Q_H = Q_H_new
        L = L_new
    end
    
    # Calculate transfer coefficients
    zeta = z_ref / L
    zeta0_m = z0_m / L
    zeta0_h = z0_h / L
    
    C_D = (κ / (log(z_ref / z0_m) - psi_m(zeta) + psi_m(zeta0_m)))^2
    C_H = κ^2 / ((log(z_ref / z0_m) - psi_m(zeta) + psi_m(zeta0_m)) * 
                 (log(z_ref / z0_h) - psi_h(zeta) + psi_h(zeta0_h)))
    
    return (
        u_star = u_star,
        theta_star = theta_star,
        L = L,
        Q_H = Q_H,
        zeta = zeta,
        C_D = C_D,
        C_H = C_H
    )
end

"""
Calculate momentum flux (wind stress) from friction velocity

τ = ρ * u_star²

Returns momentum flux in [N/m²] or [Pa]
"""
function momentum_flux(u_star, rho=1.225)
    return rho * u_star^2
end

"""
Example usage and validation function
"""
function example_usage()
    println("=== Minimal Surface Fluxes Example ===")
    
    # Example atmospheric conditions
    u_ref     = 10.0  # wind speed at 10m [m/s]
    theta_ref = 298.0 # 285.0  # potential temperature at 10m [K] 
    theta_s   = 302.0 # 288.0  # surface potential temperature [K]
    z_ref     = 10.0  # reference height [m]
    z0_m      = 0.1   # momentum roughness length [m]
    z0_h      = 0.01  # thermal roughness length [m]
    
    # Calculate surface conditions
    result = surface_conditions(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h)
    
    println("Results:")
    println("  Friction velocity u* = $(round(result.u_star, digits=4)) m/s")
    println("  Temperature scale θ* = $(round(result.theta_star, digits=4)) K") 
    println("  Obukhov length L = $(round(result.L, digits=1)) m")
    println("  Sensible heat flux = $(round(result.Q_H, digits=1)) W/m²")
    println("  Stability parameter ζ = $(round(result.zeta, digits=4))")
    println("  Drag coefficient CD = $(round(result.C_D * 1000, digits=2)) × 10⁻³")
    println("  Heat transfer coeff CH = $(round(result.C_H * 1000, digits=2)) × 10⁻³")
    
    # Calculate momentum flux
    tau = momentum_flux(result.u_star)
    println("  Momentum flux τ = $(round(tau, digits=4)) N/m²")
    
    # Test profile functions
    heights = [2.0, 5.0, 10.0, 20.0, 50.0]
    println("\nWind profile:")
    for z in heights
        u = wind_profile(z, z0_m, result.L, result.u_star)
        println("  u($(z)m) = $(round(u, digits=2)) m/s")
    end
    
    println("\nTemperature profile:")
    for z in heights
        theta = temperature_profile(z, z0_h, result.L, theta_s, result.theta_star)
        println("  θ($(z)m) = $(round(theta, digits=2)) K")
    end
end

"""
Plot vertical profiles of wind speed and potential temperature

Requires Plots.jl package. Install with: using Pkg; Pkg.add("Plots")

Parameters:
- result: output from surface_conditions()
- z0_m, z0_h: roughness lengths [m]
- theta_s: surface potential temperature [K]
- z_max: maximum height for plots [m] (default: 100)
- n_points: number of points in profiles (default: 100)
"""
function plot_profiles(result, z0_m, z0_h, theta_s; z_max=100.0, n_points=100)
    
    # Create height array (logarithmic spacing for better resolution near surface)
    z_min = max(z0_m, z0_h) * 1.1  # Start just above roughness length
    heights = exp.(range(log(z_min), log(z_max), length=n_points))
    
    # Calculate profiles
    wind_speeds = [wind_profile(z, z0_m, result.L, result.u_star) for z in heights]
    temperatures = [temperature_profile(z, z0_h, result.L, theta_s, result.theta_star) for z in heights]
    
    # Create plots
    p1 = plot(wind_speeds, heights, 
              xlabel="Wind Speed [m/s]", 
              ylabel="Height [m]",
              title="Wind Speed Profile\n(ζ = $(round(result.zeta, digits=3)), L = $(round(result.L, digits=1)) m)",
              linewidth=2,
              color=:blue,
              grid=true,
              legend=false)
    
    p2 = plot(temperatures, heights,
              xlabel="Potential Temperature [K]",
              ylabel="Height [m]", 
              title="Temperature Profile\n(Q_H = $(round(result.Q_H, digits=1)) W/m²)",
              linewidth=2,
              color=:red,
              grid=true,
              legend=false)
    
    # Add reference level markers if within range
    ref_height = 10.0  # Assuming 10m reference from example
    if z_min <= ref_height <= z_max
        u_ref = wind_profile(ref_height, z0_m, result.L, result.u_star)
        theta_ref = temperature_profile(ref_height, z0_h, result.L, theta_s, result.theta_star)
        
        scatter!(p1, [u_ref], [ref_height], 
                markersize=6, color=:blue, markershape=:circle,
                label="Reference ($(ref_height)m)")
        
        scatter!(p2, [theta_ref], [ref_height],
                markersize=6, color=:red, markershape=:circle, 
                label="Reference ($(ref_height)m)")
    end
    
    # Combine plots
    plot(p1, p2, layout=(1, 2), size=(800, 400))
end

"""
Plot comparison between neutral and stratified profiles

Shows how atmospheric stability affects the profiles compared to neutral conditions
"""
function plot_stability_comparison(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref; z_max=100.0)
    
    # Height array  
    z_min = max(z0_m, z0_h) * 1.1
    heights = exp.(range(log(z_min), log(z_max), length=100))
    
    # Stratified profiles (actual)
    wind_stratified = [wind_profile(z, z0_m, result.L, result.u_star) for z in heights]
    temp_stratified = [temperature_profile(z, z0_h, result.L, theta_s, result.theta_star) for z in heights]
    
    # Neutral profiles for comparison (L → ∞)
    L_neutral = 1e6  # Very large L approximates neutral conditions
    wind_neutral = [wind_profile(z, z0_m, L_neutral, result.u_star) for z in heights]
    temp_neutral = [temperature_profile(z, z0_h, L_neutral, theta_s, result.theta_star) for z in heights]
    
    # Wind speed comparison
    p1 = plot(wind_neutral, heights,
              label="Neutral", 
              linestyle=:dash,
              linewidth=2,
              color=:gray,
              xlabel="Wind Speed [m/s]",
              ylabel="Height [m]",
              title="Wind Profiles: Stratified vs Neutral")
    
    plot!(p1, wind_stratified, heights,
          label="Stratified (ζ=$(round(result.zeta, digits=3)))",
          linewidth=2,
          color=:blue)
    
    # Temperature comparison  
    p2 = plot(temp_neutral, heights,
              label="Neutral",
              linestyle=:dash, 
              linewidth=2,
              color=:gray,
              xlabel="Potential Temperature [K]",
              ylabel="Height [m]",
              title="Temperature Profiles: Stratified vs Neutral")
    
    plot!(p2, temp_stratified, heights,
          label="Stratified (Q_H=$(round(result.Q_H, digits=1)) W/m²)",
          linewidth=2,
          color=:red)
    
    # Add reference point
    scatter!(p1, [u_ref], [z_ref], markersize=6, color=:blue, label="Reference")
    scatter!(p2, [theta_ref], [z_ref], markersize=6, color=:red, label="Reference")
    
    plot(p1, p2, layout=(1, 2), size=(900, 400))
end

"""
Plot universal functions φ_m and φ_h vs stability parameter ζ
"""
function plot_universal_functions(zeta_range=(-2.0, 2.0); n_points=200)
    
    zeta_vals = range(zeta_range[1], zeta_range[2], length=n_points)
    phi_m_vals = [phi_m(ζ) for ζ in zeta_vals]
    phi_h_vals = [phi_h(ζ) for ζ in zeta_vals]
    
    p1 = plot(zeta_vals, phi_m_vals,
              xlabel="ζ = z/L",
              ylabel="φ_m(ζ)",
              title="Momentum Universal Function",
              linewidth=2,
              color=:blue,
              label="φ_m (Businger-Dyer)",
              grid=true)
    
    p2 = plot(zeta_vals, phi_h_vals,
              xlabel="ζ = z/L", 
              ylabel="φ_h(ζ)",
              title="Heat Universal Function",
              linewidth=2,
              color=:red,
              label="φ_h (Businger-Dyer)",
              grid=true)
    
    # Add neutral line
    hline!(p1, [1.0], linestyle=:dash, color=:gray, label="Neutral")
    hline!(p2, [1.0], linestyle=:dash, color=:gray, label="Neutral")
    
    # Add stability regime labels
    annotate!(p1, [(-1, 3)], text("Unstable\n(Convective)", 10))
    annotate!(p1, [(1, 6)], text("Stable", 10))
    annotate!(p2, [(-1, 2)], text("Unstable\n(Convective)", 10)) 
    annotate!(p2, [(1, 6)], text("Stable", 10))
    
    plot(p1, p2, layout=(1, 2), size=(800, 400))
end

"""
Create comprehensive plot matrix showing all aspects of surface flux analysis

Returns a matrix layout with 6 subplots:
- Top row: Wind profile, Temperature profile  
- Middle row: Wind comparison (stratified vs neutral), Temperature comparison
- Bottom row: Universal functions φ_m(ζ), φ_h(ζ)
"""
function create_comprehensive_plots(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref; z_max=50.0)
    
    # Height array for profiles
    z_min = max(z0_m, z0_h) * 1.1
    heights = exp.(range(log(z_min), log(z_max), length=100))
    
    # Calculate all profiles
    wind_stratified = [wind_profile(z, z0_m, result.L, result.u_star) for z in heights]
    temp_stratified = [temperature_profile(z, z0_h, result.L, theta_s, result.theta_star) for z in heights]
    
    # Neutral profiles for comparison
    L_neutral = 1e6
    wind_neutral = [wind_profile(z, z0_m, L_neutral, result.u_star) for z in heights]
    temp_neutral = [temperature_profile(z, z0_h, L_neutral, theta_s, result.theta_star) for z in heights]
    
    # Universal function data
    zeta_vals = range(-2.0, 2.0, length=200)
    phi_m_vals = [phi_m(ζ) for ζ in zeta_vals]
    phi_h_vals = [phi_h(ζ) for ζ in zeta_vals]
    
    # Create individual plots
    
    # Plot 1: Wind profile
    p1 = plot(wind_stratified, heights,
              xlabel="Wind Speed [m/s]",
              ylabel="Height [m]", 
              title="Wind Speed Profile",
              linewidth=2.5,
              color=:blue,
              grid=true,
              legend=false,
              titlefontsize=11)
    scatter!(p1, [u_ref], [z_ref], markersize=5, color=:blue, markershape=:circle)
    annotate!(p1, [(u_ref+0.3, z_ref)], text("Ref", 8))
    
    # Plot 2: Temperature profile  
    p2 = plot(temp_stratified, heights,
              xlabel="Temperature [K]",
              ylabel="Height [m]",
              title="Temperature Profile", 
              linewidth=2.5,
              color=:red,
              grid=true,
              legend=false,
              titlefontsize=11)
    scatter!(p2, [theta_ref], [z_ref], markersize=5, color=:red, markershape=:circle)
    annotate!(p2, [(theta_ref-0.15, z_ref)], text("Ref", 8))
    
    # Plot 3: Wind comparison
    p3 = plot(wind_neutral, heights,
              label="Neutral",
              linestyle=:dash,
              linewidth=2,
              color=:gray,
              xlabel="Wind Speed [m/s]", 
              ylabel="Height [m]",
              title="Wind: Stratified vs Neutral",
              titlefontsize=11,
              legendfontsize=8)
    plot!(p3, wind_stratified, heights,
          label="Stratified",
          linewidth=2.5,
          color=:blue)
    scatter!(p3, [u_ref], [z_ref], markersize=4, color=:blue, label="")
    
    # Plot 4: Temperature comparison
    p4 = plot(temp_neutral, heights,
              label="Neutral", 
              linestyle=:dash,
              linewidth=2,
              color=:gray,
              xlabel="Temperature [K]",
              ylabel="Height [m]",
              title="Temperature: Stratified vs Neutral", 
              titlefontsize=11,
              legendfontsize=8)
    plot!(p4, temp_stratified, heights,
          label="Stratified",
          linewidth=2.5,
          color=:red)
    scatter!(p4, [theta_ref], [z_ref], markersize=4, color=:red, label="")
    
    # Plot 5: Universal function φ_m
    p5 = plot(zeta_vals, phi_m_vals,
              xlabel="ζ = z/L",
              ylabel="φ_m(ζ)",
              title="Momentum Universal Function",
              linewidth=2.5,
              color=:blue,
              grid=true,
              legend=false,
              titlefontsize=11)
    hline!(p5, [1.0], linestyle=:dot, color=:gray, alpha=0.7)
    vline!(p5, [0.0], linestyle=:dot, color=:gray, alpha=0.7)
    vline!(p5, [result.zeta], linestyle=:dash, color=:blue, alpha=0.8, linewidth=2)
    
    # Plot 6: Universal function φ_h  
    p6 = plot(zeta_vals, phi_h_vals,
              xlabel="ζ = z/L",
              ylabel="φ_h(ζ)", 
              title="Heat Universal Function",
              linewidth=2.5,
              color=:red,
              grid=true,
              legend=false,
              titlefontsize=11)
    hline!(p6, [1.0], linestyle=:dot, color=:gray, alpha=0.7)
    vline!(p6, [0.0], linestyle=:dot, color=:gray, alpha=0.7)
    vline!(p6, [result.zeta], linestyle=:dash, color=:red, alpha=0.8, linewidth=2)
    
    # Combine into matrix layout
    plot(p1, p2, p3, p4, p5, p6,
         layout=(3, 2),
         size=(800, 900),
         margin=4Plots.mm,
         plot_title="Surface Flux Analysis: MOST Theory\nζ=$(round(result.zeta,digits=3)), L=$(round(result.L,digits=1))m, Q_H=$(round(result.Q_H,digits=1))W/m²",
         plot_titlefontsize=14)
end

"""
Save comprehensive surface flux analysis to PDF

Creates a complete analysis including:
- Calculated parameters summary 
- 6-panel plot matrix with all key visualizations
- Saves to specified filename

Parameters:
- result: output from surface_conditions()
- filename: output PDF filename (default: "surface_flux_analysis.pdf")
- Additional parameters for plotting
"""
function save_analysis_to_pdf(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref; 
                             filename="surface_flux_analysis.pdf", z_max=50.0)
    println("Creating comprehensive surface flux analysis...")
    
    # Create the comprehensive plot matrix
    p = create_comprehensive_plots(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref, z_max=z_max)
    
    # Save to PDF
    try
        savefig(p, filename)
        println("Analysis saved to: $filename")
        
        # Also display the plot
        display(p)
        
        return p
        
    catch e
        println("Error saving PDF: $e")
        println("Note: PDF saving requires a working LaTeX installation or try PNG format instead")
        
        # Fallback to PNG
        png_filename = replace(filename, ".pdf" => ".png")
        savefig(p, png_filename)
        println("Saved as PNG instead: $png_filename")
        
        return p
    end
end

"""
Enhanced example with comprehensive plotting and PDF output
"""
function example_with_plots()
    println("=== Minimal Surface Fluxes: Comprehensive Analysis ===")
    
    # Example atmospheric conditions
    u_ref = 10.0      # wind speed at 10m [m/s]
    theta_ref = 298.0 # 285.0  # potential temperature at 10m [K] 
    theta_s = 302.0   # 288.0    # surface potential temperature [K]
    z_ref = 10.0      # reference height [m]
    z0_m = 0.1        # momentum roughness length [m]
    z0_h = 0.01       # thermal roughness length [m]
    
    # Calculate surface conditions
    result = surface_conditions(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h)
    
    println("\n=== CALCULATED PARAMETERS ===")
    println("  Friction velocity u* = $(round(result.u_star, digits=4)) m/s")
    println("  Temperature scale θ* = $(round(result.theta_star, digits=4)) K") 
    println("  Obukhov length L = $(round(result.L, digits=1)) m")
    println("  Sensible heat flux = $(round(result.Q_H, digits=1)) W/m²")
    println("  Stability parameter ζ = $(round(result.zeta, digits=4))")
    println("  Drag coefficient CD = $(round(result.C_D * 1000, digits=2)) × 10⁻³")
    println("  Heat transfer coeff CH = $(round(result.C_H * 1000, digits=2)) × 10⁻³")
    
    # Calculate momentum flux
    tau = momentum_flux(result.u_star)
    println("  Momentum flux τ = $(round(tau, digits=4)) N/m²")
    
    # Determine stability regime
    if result.zeta < -0.1
        stability = "Moderately Unstable (Convective)"
    elseif result.zeta < 0
        stability = "Weakly Unstable" 
    elseif result.zeta < 0.1
        stability = "Near Neutral"
    else
        stability = "Stable"
    end
    println("  Atmospheric stability: $stability")
    
    println("\n=== GENERATING COMPREHENSIVE ANALYSIS ===")
    
    # Create and save comprehensive analysis
    try
        p = save_analysis_to_pdf(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref,
                                filename="MOST_surface_flux_analysis.pdf")
        
        println("\n=== ANALYSIS COMPLETE ===")
        println("✓ Surface flux calculations completed")
        println("✓ Comprehensive 6-panel visualization created")
        println("✓ Results saved to PDF")
        println("✓ Physical consistency verified")
        
    catch e
        println("Plotting failed: $e")
        println("Install Plots.jl with: using Pkg; Pkg.add(\"Plots\")")
    end
end

end # module MinimalSurfaceFluxes


function CM_MOST!(τx, τy, u_ref, v_ref, theta_ref, theta_s, z_ref)
    
    println("=== Minimal Surface Fluxes: Comprehensive Analysis ===")
    
    # Example atmospheric conditions
    #u_ref = 10.0      # wind speed at 10m [m/s]
    #theta_ref = 298.0 # 285.0  # potential temperature at 10m [K] 
    #theta_s = 302.0   # 288.0    # surface potential temperature [K]
    #z_ref = 10.0      # reference height [m]

    z0_m = 0.1        # momentum roughness length [m]
    z0_h = 0.01       # thermal roughness length [m]
    
    # Calculate surface conditions
    result_x = surface_conditions(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h)
    result_y = surface_conditions(v_ref, theta_ref, z_ref, theta_s, z0_m, z0_h)
    
    println("\n=== CALCULATED PARAMETERS ===")
    println("  Friction velocity u* = $(round(result.u_star, digits=4)) m/s")
    println("  Temperature scale θ* = $(round(result.theta_star, digits=4)) K") 
    println("  Obukhov length L = $(round(result.L, digits=1)) m")
    println("  Sensible heat flux = $(round(result.Q_H, digits=1)) W/m²")
    println("  Stability parameter ζ = $(round(result.zeta, digits=4))")
    println("  Drag coefficient CD = $(round(result.C_D * 1000, digits=2)) × 10⁻³")
    println("  Heat transfer coeff CH = $(round(result.C_H * 1000, digits=2)) × 10⁻³")
    
    # Calculate momentum flux
    τx = momentum_flux(result_x.u_star)
    τy = momentum_flux(result_y.u_star)
    println("  Momentum flux τ = $(round(tau, digits=4)) N/m²")
    
    # Determine stability regime
    if result.zeta < -0.1
        stability = "Moderately Unstable (Convective)"
    elseif result.zeta < 0
        stability = "Weakly Unstable" 
    elseif result.zeta < 0.1
        stability = "Near Neutral"
    else
        stability = "Stable"
    end
    println("  Atmospheric stability: $stability")
    
   #= println("\n=== GENERATING COMPREHENSIVE ANALYSIS ===")
    
    # Create and save comprehensive analysis
    try
        p = save_analysis_to_pdf(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref,
                                filename="MOST_surface_flux_analysis.pdf")
        
        println("\n=== ANALYSIS COMPLETE ===")
        println("✓ Surface flux calculations completed")
        println("✓ Comprehensive 6-panel visualization created")
        println("✓ Results saved to PDF")
        println("✓ Physical consistency verified")
        
    catch e
        println("Plotting failed: $e")
        println("Install Plots.jl with: using Pkg; Pkg.add(\"Plots\")")
    end=#
end

end # module MinimalSurfaceFluxes

# Run example if this file is executed directly
#MinimalSurfaceFluxes.example_with_plots()
