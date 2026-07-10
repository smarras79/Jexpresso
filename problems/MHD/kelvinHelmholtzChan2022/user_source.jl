#---------------------------------------------------------------------------------
# The GLM-MHD system of Chan et al. (2022), Eq. (6), is solved WITHOUT source
# terms (inputs[:lsource] = false), so this is a no-op kept only because the
# harness includes user_source.jl for every case.
#
# NOTE: the non-conservative term Υ of the paper's Eq. (7) is NOT a pointwise
# source (it involves ∇·B and v·∇ψ) and cannot be expressed through this
# hook; it is intentionally omitted in this first, non-entropy-stable setup.
# See the header of user_flux.jl.
#---------------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=9, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    for ieq = 1:neqs
        S[ieq] = 0.0
    end

end
