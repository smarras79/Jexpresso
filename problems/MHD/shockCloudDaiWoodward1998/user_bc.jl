#---------------------------------------------------------------------------------
# Boundary conditions. Dai & Woodward use an open boundary on all four sides.
#
#   "left"    The post-shock gas enters at u = 13.84 > c_f = 10.52 (super-fast): every characteristic
#             comes in, so the whole post-shock state is imposed (Tóth 2000 §6.5 fixes its inflow the same way).
#   others    Open (do-nothing): kernel/boundaryconditions/BCs.jl copies back only the slots written here,
#             so leaving qbdy alone lets the interior solution leave through the strong-form CG boundary,
#             as for the outflow sides of problems/MHD/astroJetWuShu2018. ψ = 0 is the one incoming
#             condition the GLM pair (Bₙ, ψ) needs.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    if tag == "left"
        s = sc_conserved(SC_POST)
        for ieq = 1:9
            qbdy[ieq] = s[ieq]
        end
    else
        qbdy[9] = 0.0
    end
    return nothing
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::PERT)
    error(" problems/MHD/shockCloudDaiWoodward1998: PERT() solution variables are not supported.")
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, inputs)
    flux = zeros(size(q,2),1)
    return flux
end
