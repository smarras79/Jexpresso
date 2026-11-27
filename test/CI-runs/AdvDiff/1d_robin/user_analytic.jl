# Exact solution, overlaid on the output plots (src/io/plotting/jeplots.jl).
function user_analytic_solution(x, t, outvar, inputs)
    return reshape([r1g_exact(xi - inputs[:xmin], t) for xi in x], :, 1)
end
