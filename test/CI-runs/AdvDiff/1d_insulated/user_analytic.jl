# Exact solution, overlaid on the output plots (src/io/plotting/jeplots.jl).
function user_analytic_solution(x, t, outvar, inputs)
    return reshape([i1_exact(xi - inputs[:xmin], t) for xi in x], :, 1)
end
