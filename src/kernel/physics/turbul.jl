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

#test_wall_model()
