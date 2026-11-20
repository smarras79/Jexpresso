# Steady state the run converges to, overlaid on the output plots (src/io/plotting/jeplots.jl).
# The transient has no closed form, so only the final plot should match it.
function user_analytic_solution(x, t, outvar, inputs)
    return reshape([rad_steady(xi - inputs[:xmin]) for xi in x], :, 1)
end
