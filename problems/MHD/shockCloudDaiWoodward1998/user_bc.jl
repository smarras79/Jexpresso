#---------------------------------------------------------------------------------
# Boundary conditions. Dai & Woodward use open boundaries on all four sides.
#
# THE DO-NOTHING CONDITION IS UNSTABLE HERE. Leaving qbdy alone, as the outflow sides of
# problems/MHD/astroJetWuShu2018 do, imposes nothing. That suits a super-fast outflow, where every
# characteristic leaves. It goes wrong where gas at rest touches the boundary, as in the post-shock gas
# of the "cloud" test at x = 0. Measured on SC15_48x32: the boundary started to empty uniformly in y,
# u went from −0.007 at t = 0.01 to −2.9 at t = 0.015, and by t = 0.03 the gas within 0.25 of the
# boundary had lost 30% of its density and half its pressure, before any wave from the cloud arrived.
#
# CHARACTERISTIC FAR FIELD. At each boundary node the state q is split into waves about a far-field
# state q∞, using the eigenvectors of the Jacobian A = ∂(F nx + G ny)/∂q of the GLM-MHD flux at q∞
# (outward normal n):
#
#     q_b = q∞ + Σ_{λ_k > 0} r_k l_k (q − q∞)       (A r_k = λ_k r_k,  l_k r_j = δ_kj)
#
# The waves with λ > 0 leave, so they keep their interior value; the waves with λ ≤ 0 enter, so they
# take the far-field value. At a super-fast inflow no MHD wave leaves, so q_b = q∞ apart from the
# outgoing GLM wave (speed c_h > |v|): the whole state is imposed, as Tóth (2000, §6.5) does for DW's
# single-cloud test. At a super-fast outflow every MHD wave leaves and q_b = q. The condition is linear
# about q∞: a strong wave, such as the bow shock crossing the boundary, is let out only approximately.
#
# THE FAR FIELD is the undisturbed solution: DW's planar shock, moving at s = (ρ₂u₂ − ρ₁u₁)/(ρ₂ − ρ₁)
# (18.59 for "letters", 5.18 for "cloud"),
#
#     q∞(x, t) = post-shock state  for x < x_s + s t,   pre-shock state otherwise.
#
# It is exact on every side until the waves from the clouds reach it. The eigenvectors are computed once
# per state and per side, by initialize.jl (sc_build_far_field!), because the Jacobian depends on c_h.
#
# EXCEPT ACROSS THE SHOCK. Where the shock meets the top and bottom, the computed shock is spread over
# about two elements (:sc_smooth and DynSGS), while q∞ jumps at x_s + s t. A node inside that spread was
# pulled toward the far-field state of the wrong side: on SC15_48x32, ρ there was off by 12 after five
# steps, and the run blew up at t = 0.0016. Within |x − x_s − s t| < δ = 2 h + :sc_smooth (h the element
# size) the far field is not known, so those nodes only get ψ = 0 and are otherwise left open.
#---------------------------------------------------------------------------------

if !@isdefined(SC_FAR)
    # (state 1 = pre / 2 = post, side 1..4 = left, right, bottom, top) -> projector onto the outgoing waves
    const SC_FAR    = Dict{Tuple{Int,Int},Matrix{Float64}}()
    const SC_FAR_Q  = Vector{NTuple{9,Float64}}(undef, 2)   # q∞: pre, post (conserved, HL)
    const SC_SHOCK  = Ref((0.2, 18.59, 0.05))                # (x_s, s, δ)
    const SC_SIDES  = ("left" => (-1.0, 0.0), "right" => (1.0, 0.0), "bottom" => (0.0, -1.0), "top" => (0.0, 1.0))
end

# Normal GLM-MHD flux F nx + G ny, the same as user_flux! of problems/MHD/orszagTangBormanis2024.
function sc_normal_flux(q, nx, ny)
    ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ = q
    u, v, w = ρu/ρ, ρv/ρ, ρw/ρ
    c_h  = c_h_mhd[]
    p    = pressure_mhd(ρ, ρu, ρv, ρw, ρE, Bx, By, Bz, ψ)
    magp = 0.5*(Bx*Bx + By*By + Bz*Bz)
    vdB  = u*Bx + v*By + w*Bz
    enfl = 0.5*(ρu*u + ρv*v + ρw*w) + γ_mhd*p/(γ_mhd - 1.0) + 2.0*magp
    un, Bn = u*nx + v*ny, Bx*nx + By*ny
    return [ρ*un,
            ρu*un + (p + magp)*nx - Bx*Bn,
            ρv*un + (p + magp)*ny - By*Bn,
            un*enfl - Bn*vdB + c_h*ψ*Bn,
            ρw*un - Bz*Bn,
            un*Bx - u*Bn + c_h*ψ*nx,
            un*By - v*Bn + c_h*ψ*ny,
            un*Bz - w*Bn,
            c_h*Bn]
end

# Projector R diag(λ > 0) R⁻¹ onto the waves leaving through the side of outward normal n, at state q.
function sc_outgoing_projector(q::NTuple{9,Float64}, nx, ny)
    A = zeros(9, 9)
    for k = 1:9
        h  = 1.0e-6*max(1.0, abs(q[k]))
        qp = collect(q); qp[k] += h
        qm = collect(q); qm[k] -= h
        A[:, k] = (sc_normal_flux(qp, nx, ny) - sc_normal_flux(qm, nx, ny))/(2.0*h)
    end
    E  = eigen(A)
    λ  = real.(E.values)
    R  = real.(E.vectors)
    tol = 1.0e-8*maximum(abs, λ)
    return R*Diagonal(Float64.(λ .> tol))*inv(R)
end

function sc_build_far_field!(st, xs, δ)
    SC_FAR_Q[1] = sc_conserved(st.pre)
    SC_FAR_Q[2] = sc_conserved(st.post)
    s = (st.post[1]*st.post[3] - st.pre[1]*st.pre[3])/(st.post[1] - st.pre[1])
    SC_SHOCK[] = (xs, s, δ)
    for (i, q) in enumerate(SC_FAR_Q), (j, (_, n)) in enumerate(SC_SIDES)
        SC_FAR[(i, j)] = sc_outgoing_projector(q, n...)
    end
    return s
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    j = tag == "left" ? 1 : tag == "right" ? 2 : tag == "bottom" ? 3 : tag == "top" ? 4 : 0
    j == 0 && return nothing
    xs, s, δ = SC_SHOCK[]
    d  = coords[1] - xs - s*t
    if abs(d) < δ                    # inside the shock's spread: open, ψ = 0
        qbdy[9] = 0.0
        return nothing
    end
    i  = d < 0.0 ? 2 : 1
    q∞ = SC_FAR_Q[i]
    P  = SC_FAR[(i, j)]
    for a = 1:9
        acc = q∞[a]
        for b = 1:9
            acc += P[a,b]*(q[b] - q∞[b])
        end
        qbdy[a] = acc
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
