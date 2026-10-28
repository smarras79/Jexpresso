#---------------------------------------------------------------------------------
# Exact KdV soliton, overlaid on the 1D output plots by src/io/plotting/jeplots.jl
#
#     user_analytic_solution(x, t, outvar, inputs) -> npoin × nvar matrix
#
# u is compared with η, v with ∂ₓη and w with ∂ₓ²η (kdv_soliton in user_flux.jl).
# This is the exact solution of KdV, the τ → 0 limit, not of KdVH itself, so
# the gap between the curves is discretisation error plus the O(cτ) model
# error of the relaxation.
#---------------------------------------------------------------------------------
function user_analytic_solution(x, t, outvar, inputs)

    nvar = length(outvar)
    q    = fill(NaN, length(x), nvar)

    for i in eachindex(x)
        η, ηx, ηxx = kdv_soliton(x[i], t)

        for ivar = 1:nvar
            name = string(outvar[ivar])
            q[i, ivar] = name == "u" ? η   :
                         name == "v" ? ηx  :
                         name == "w" ? ηxx : NaN
        end
    end

    return q
end
