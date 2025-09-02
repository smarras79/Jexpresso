using Roots
using Printf

# Main function to find uτ given u2 and y2
function find_uτ(u2, y2)
    # Handle the sign of u2 properly
    u2_abs = abs(u2)
    u2_sign = sign(u2)
    
    # Check for edge cases
    if u2_abs < 1e-12
        return 0.0
    end
    
    if y2 <= 0
        error("Wall-normal distance y2 must be positive, got y2 = $y2")
    end
    
    # Define the nonlinear function to find zeros
    # f(uτ) = uτ*(κinv * ln(y2 * uτ / ν) + C) - |u2| = 0
    function wall_model_residual(uτ)
        
        # Wall model parameters (constants)
        κinv = 2.5    # Inverse of von Karman constant (1/κ)
        C = 5.5       # Additive constant in log law
        ν = 1.0e-5    # Kinematic viscosity
                
        return uτ * (κinv * log(y2 * uτ / ν) + C) - u2_abs
    end
    
    # Solve for uτ using Brent's method
    uτ_low = 1e-8
    uτ_high = 10.0 * u2_abs  # Conservative upper bound
    
    try
        uτ_solution = find_zero(wall_model_residual, (uτ_low, uτ_high), Roots.Brent())
        return uτ_solution * u2_sign  # Apply original sign
    catch e
        @warn "Brent method failed for u2=$u2, y2=$y2: $e"
        # Fallback to Newton method
        return find_uτ_newton_fallback(u2_abs, y2) * u2_sign
    end
end

# Fallback Newton method
function find_uτ_newton_fallback(u2_abs, y2)
    function residual_and_derivative(uτ)
        residual = uτ * (κinv * log(y2 * uτ / ν) + C) - u2_abs
        derivative = κinv * log(y2 * uτ / ν) + C + κinv
        return (residual, derivative)
    end
    
    # Initial guess
    uτ_init = u2_abs / (κinv * log(y2 * u2_abs / (ν * 10)) + C)
    uτ_init = max(uτ_init, 1e-8)  # Ensure positive initial guess
    
    try
        uτ_solution = find_zero(residual_and_derivative, uτ_init, Roots.Newton())
        return uτ_solution
    catch e
        @warn "Newton fallback also failed: $e"
        return NaN
    end
end

# Test function to demonstrate usage
function test_wall_model()
    # Test cases with different signs and magnitudes
    test_cases = [
        (0.1, 17.0),     # Positive velocity
        (0.0, 17.0),     # zero velocity
        (10.0, 950.0),     # Positive velocity
        (-10.0, 950.0),    # Negative velocity
        (5.0, 500.0),      # Smaller velocity, closer to wall
        (-15.0, 1200.0),   # Negative velocity, further from wall
        (0.1, 100.0),      # Small positive velocity
        (-0.1, 100.0),     # Small negative velocity
    ]
            
    for (u2_test, y2_test) in test_cases
        uτ = find_uτ(u2_test, y2_test)
        
        if !isnan(uτ)
            # Calculate dimensionless parameters
            y_plus = y2_test * abs(uτ) / ν
            u_plus = u2_test / uτ
            
            @printf("u₂ = %8.3f, y₂ = %8.1f → uτ = %8.5f, y⁺ = %8.1f, u⁺ = %8.3f\n", 
                   u2_test, y2_test, uτ, y_plus, u_plus)
        else
            @printf("u₂ = %8.3f, y₂ = %8.1f → FAILED\n", u2_test, y2_test)
        end
    end
end

#
# MOST
#
function MOST!(τ_f, wθ,
               lwall_model,
               uprimitiveieq,
               inputs,
               PhysConst,
               iel, ieq,
               connijk,
               coords,
               poin_in_bdy_face, elem_to_face, bdy_face_type,
               k, l, m, iface_bdy, idx1, idx2)
    
    ip        = connijk[iel, k, l, m]
    iface_bdy = elem_to_face[iel, k, l, m, 1]
    
    # For top wall, use second point below (m=N-1, where N is max in that direction)
    wall_y = coords[ip, 3]  # y-coordinate of wall point
    ieq    = 2
    u2     = uprimitiveieq[k, l, 2, ieq]
    ip2    = connijk[iel, k, l, 2]
    y2     = coords[ip2, 3]
    
    # Get velocity components and temperature at second grid point
    u_vel   = uprimitiveieq[k, l, 2, 2]  # u-component
    v_vel   = uprimitiveieq[k, l, 2, 3]  # v-component
    vel_mag = sqrt(u_vel^2 + v_vel^2)
    θ_2     = uprimitiveieq[k, l, 2, 5]    # potential temperature (assuming index 5)
    
    # Height of first grid point above surface
    z1 = abs(y2 - wall_y)
    
    # Surface parameters
    z0     = 0.1         # momentum roughness length (m)
    z0h    = z0/10.0     # thermal roughness length (m)
    κ      = 0.4         # von Kármán constant
    g      = PhysConst.g # gravitational acceleration
    θ_surf = 300.0       # surface temperature (K) - should be prescribed
    
    # MOST iterative solution for surface fluxes
    L_old    = 1e6       # initial guess for Obukhov length (neutral conditions)
    max_iter = 20
    tol      = 1e-3
    
    uτ       = 0.0
    θ_star   = 0.0
    for iter = 1:max_iter
        # Stability parameter
        ζ = z1 / L_old
        
        # Limit stability parameter to prevent numerical issues
        ζ = max(-2.0, min(1.0, ζ))
        
        # Universal functions and integrated stability functions
        if ζ < 0  # Unstable conditions
            x = (1 - 16*ζ)^0.25
            φ_m = (1 - 16*ζ)^(-0.25)
            φ_h = (1 - 16*ζ)^(-0.5)
            
            # Integrated stability functions
            Ψ_m = 2*log((1 + x)/2) + log((1 + x^2)/2) - 2*atan(x) + π/2
            Ψ_h = 2*log((1 + x^2)/2)
            
        elseif ζ > 0  # Stable conditions
            φ_m = 1 + 5*ζ
            φ_h = 1 + 5*ζ
            
            # Integrated stability functions
            Ψ_m = -5*ζ
            Ψ_h = -5*ζ
            
        else  # Neutral conditions
            φ_m = 1.0
            φ_h = 1.0
            Ψ_m = 0.0
            Ψ_h = 0.0
        end
        
        # Calculate friction velocity
        if vel_mag > 1e-12
            log_term_m = log(z1/z0) - Ψ_m
            log_term_m = max(log_term_m, 0.1)  # Prevent division by small numbers
            uτ = κ * vel_mag / log_term_m
        else
            uτ = 0.0
        end
        
        # Calculate temperature scale
        if abs(θ_2 - θ_surf) > 1e-6 && uτ > 1e-12
            log_term_h = log(z1/z0h) - Ψ_h
            log_term_h = max(log_term_h, 0.1)
            θ_star = κ * (θ_2 - θ_surf) / log_term_h
        else
            θ_star = 0.0
        end
        
        # Update Obukhov length
        if abs(θ_star) > 1e-12 && uτ > 1e-12
            L_new = -uτ^3 * θ_2 / (κ * g * θ_star)
        else
            L_new = 1e6  # Neutral conditions
        end
        
        # Check convergence
        if abs(L_new - L_old) / max(abs(L_old), 1.0) < tol
            break
        end
        
        # Update with relaxation to improve stability
        L_old = 0.7 * L_old + 0.3 * L_new
    end
    
    # Calculate surface fluxes and apply boundary conditions
    if !isnan(uτ) && uτ > 1e-12
        τw_mag = uprimitiveieq[k, l, m, 1] * uτ^2  # ρ * u_τ^2
        
        if vel_mag > 1e-12
            # Compute wall shear stress components
            # The shear stress always opposes the flow direction
            τw_x = -τw_mag * u_vel / vel_mag
            τw_y = -τw_mag * v_vel / vel_mag
            
            # Apply wall normal sign for gradient direction
            # For bottom wall: gradient = (u_wall - u_interior) / dy < 0 (if u_interior > 0)
            # For top wall: gradient = (u_interior - u_wall) / dy < 0 (if u_interior > 0)
            τ_f[iface_bdy, idx1, idx2, 1] = -τw_x
            τ_f[iface_bdy, idx1, idx2, 2] = -τw_y
            
            # Apply heat flux boundary condition (if temperature equation is being solved)
            # Heat flux: q = -ρ * c_p * u_τ * θ_star
            # This would typically be applied to the temperature equation boundary condition
            #q_wall = -uprimitiveieq[k, l, m, 1] * PhysConst.cp * uτ * θ_star
            q_wall = 0.12 #[Km/s] constant imposition
            wθ[iface_bdy, idx1, idx2, 1] = q_wall
            # τ_f[iface_bdy, idx1, idx2, 5] = q_wall  # assuming temperature is equation 5
        end
    else
        # Fallback to no-slip condition if MOST fails
        τ_f[iface_bdy, idx1, idx2, 1] = 0.0
        τ_f[iface_bdy, idx1, idx2, 2] = 0.0

        # Even in fallback, we might still want to impose the heat flux
        wθ[iface_bdy, idx1, idx2, 1] = 0.0
    end
end
