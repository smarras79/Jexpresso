using Roots

"""
    nonlinear_function(uτ, κinv, u2, ν, y2, C)

Define the nonlinear function: uτ*(κinv*log(y2*uτ/ν) + C) - u2

# Arguments
- `uτ`: Variable to solve for
- `κinv`: Inverse of κ (1/0.4 = 2.5)
- `u2`: Constant (10)
- `ν`: Constant (1e-4)
- `y2`: Constant (1)
"""
function nonlinear_function(uτ, κinv, u2, ν, y2, C)
    return uτ*(κinv * log(y2 * uτ / ν) + C) - u2
end

function find_zeros_nonlinear(; κinv=2.5, u2=10, ν=1e-4, y2=1, C=5.0,
                              search_range=(1.0, 1.0), method=:brent)
    
    # Define the function with fixed parameters
    f(uτ) = nonlinear_function(uτ, κinv, u2, ν, y2, C)
    
    try
        if method == :newton
            # For Newton's method, use the midpoint as initial guess
            x0 = sum(search_range) / 2
            root = find_zero(f, x0, Roots.Newton())
        elseif method == :bisection
            root = find_zero(f, search_range, Roots.Bisection())
        else
            # Default to Brent's method (robust and fast)
            root = find_zero(f, search_range, Roots.Brent())
        end
        
        return root
    catch e
        println("Error finding root: ", e)
        println("Try adjusting the search_range or checking if a root exists in the interval")
        return nothing
    end
end

"""
    analyze_function(; κinv=2.5, u2=10, ν=1e-4, y2=1, C=5.0, 
                     uτ_range=(0.01, 10.0), n_points=1000)

Analyze the behavior of the nonlinear function over a range of uτ values.
This helps visualize the function and identify potential zero locations.

# Returns
- Tuple of (uτ_values, function_values)
"""
function analyze_function(; κinv=2.5, u2=10, ν=1e-4, y2=1, C=5.0,
                         uτ_range=(0.01, 10.0), n_points=1000)
    
    uτ_vals = range(uτ_range[1], uτ_range[2], length=n_points)
    f_vals = [nonlinear_function(uτ, κinv, u2, ν, y2, C) for uτ in uτ_vals]
    
    return uτ_vals, f_vals
end

# Example usage and testing
function jeFind_uτ(u2, y2, κ, ν, C)

    println("=== Nonlinear Zero Finding Example ===")
    println("Function: uτ*(κinv*log(y2*uτ/ν) + C) - u2 = 0")
    println("Parameters: κ=$κ, u2=$u2, ν=$ν, y2=$y2")
    println()
    
    # Find the zero
    method = :bisection #Newton doesn't seem to converge
    κinv = 1.0/κ
    #uτ_range = (-10.0, 10.0)
    #
    # uτ as zeros of log-law: uτ((1/κ)log(y*uτ/ν) + C) - u = 0.0 where y,u are at some point ABOVE the surface
    #

    f(uτ, p) = uτ*(κinv * log(y2 * uτ / ν) + C) - u2 + p
    u0 = [-1.0, 1.0]
    p = 0.0
    prob = NonlinearProblem(f, u0, p)
    sol = solve(prob)
    println(" Uτ ===== ", sol.uτ)
    #uτ = find_zeros_nonlinear(; κinv=κinv, u2=u2, ν=ν, y2=y2, C=C, search_range=uτ_range, method=method)
    
    #if uτ !== nothing
    #    println("Root found: uτ = ", uτ)
    #    # Verify the solution
    #    verification = nonlinear_function(uτ, κinv, u2, ν, y2, C)
    #    println("Verification f(uτ) = ", verification)
    #    println("(Should be close to zero)")
    #end

    return uτ
end

#----------------------------------------------
# flow/fluid parameters (these should come from Jexpresso)
#----------------------------------------------
#=
u2 = 10.0
y2 = 1.0
κ  = 0.4
ν  = 1.0e-4
ρ  = 1.0
C  = 5.0 # for smooth wall
=#
#----------------------------------------------
# No user changes beyond this poin
#----------------------------------------------
#=
uτ = jeFind_uτ(u2, y2, κ, ν, C)

#verify: 
τw = ρ*uτ^uτ
yp = uτ*y2/ν
up = 1.0/κ*log(yp) + C
upverify = u2/uτ

println(" τw = ", τw)
println(" uτ = ", uτ)
println(" yp = ", yp)
println(" up = ", up)
println(" upverify = ", upverify)
println(" up - upverify = ", up - upverify)
=#
