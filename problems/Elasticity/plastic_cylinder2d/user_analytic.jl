#---------------------------------------------------------------------------------
# References for the elastic-perfectly plastic thick cylinder.
#
# There is no closed form for the whole elastic-plastic history of a von Mises
# cylinder in plane strain with ν < 1/2 (σ_zz is not the mean stress while
# part of the wall is still elastic), so user_solid_exact returns `nothing`
# and the solver compares nothing at run time. What IS exact:
#
#   1. FIRST YIELD. Up to it the solution is Lamé's, and yielding starts at the
#      bore, where the von Mises stress of
#          σrr = −p,  σθθ = P(1 + b²/a²),  σzz = 2νP,   P = p a²/(b² − a²)
#      reaches σy0. Below δ_e the computed pressure must BE the Lamé pressure.
#
#   2. THE LIMIT PRESSURE. With the whole wall flowing, plastic incompressibility
#      and ε_zz = 0 give s_zz = 0, so the von Mises condition reads
#      σθθ − σrr = 2σy0/√3, and equilibrium dσrr/dr = (σθθ − σrr)/r integrates
#      to p_L = (2/√3) σy0 ln(b/a) (Hill 1950). Displacement-controlled, the
#      bore pressure rises to p_L and stays there.
#
#   3. IN BETWEEN, the problem stays axisymmetric: u = u(r) e_r. The radial
#      reference below solves exactly that, on a 1D grid of `nel` elements with
#      one quadrature point each (the volumetric constraint then cannot lock),
#      through the SAME constitutive routine and the same load increments as
#      the SEM run — so SEM minus reference is the 2D discretisation's error.
#---------------------------------------------------------------------------------

user_solid_exact(x, y, t) = nothing        # no closed form: see above

"""
    plastic_cylinder_first_yield() -> (; p_e, δ_e)

Pressure and bore displacement at which the bore starts to yield (Lamé).
"""
function plastic_cylinder_first_yield()
    P  = plastic_cylinder_properties()
    c0 = P.a^2/(P.b^2 - P.a^2)                         # P/p
    srr = -1.0; stt = c0*(1.0 + P.b^2/P.a^2); szz = 2.0*P.ν*c0
    svm = sqrt(0.5*((srr - stt)^2 + (stt - szz)^2 + (szz - srr)^2))
    p_e = P.σy0/svm
    M   = 2.0*P.K + 2.0*P.G/3.0
    δ_e = p_e*c0*(P.a/M + P.b^2/(2.0*P.G*P.a))
    return (; p_e, δ_e)
end

"""
    plastic_cylinder_limit_pressure() -> p_L = (2/√3) σy0 ln(b/a)
"""
plastic_cylinder_limit_pressure() = (P = plastic_cylinder_properties(); 2.0*P.σy0*log(P.b/P.a)/sqrt(3.0))

"""
    plastic_cylinder_radial_reference(δs; nel = 2000) -> (; p, ep_bore)

The axisymmetric plane-strain problem on a 1D radial grid, driven through the
bore displacements `δs` (the same increments as the SEM run): the bore
pressure and the accumulated plastic strain at the bore after every increment.
"""
function plastic_cylinder_radial_reference(δs::AbstractVector; nel::Int = 2000)
    P     = plastic_cylinder_properties()
    model = SolidJ2Plasticity(P.E, P.ν, P.K, P.G, P.σy0, P.Hiso, P.Hkin)
    r     = collect(range(P.a, P.b, length = nel + 1))
    u     = zeros(nel + 1)
    H     = zeros(5, nel); Htr = zeros(5, nel)
    σ     = zeros(4); D = zeros(3, 3)
    fint  = zeros(nel + 1)
    dl = zeros(nel); dd = zeros(nel + 1); du = zeros(nel)
    ps = Float64[]; epb = Float64[]

    for δ in δs
        # previous solution, rescaled to the new bore displacement, as first guess
        u[1] != 0 ? (u .*= δ/u[1]) : (u .= δ .* P.a ./ r)
        u[1] = δ
        for it = 1:50
            fill!(fint, 0.0); fill!(dl, 0.0); fill!(dd, 0.0); fill!(du, 0.0)
            for e = 1:nel
                h  = r[e+1] - r[e]; rm = 0.5*(r[e] + r[e+1])
                εrr = (u[e+1] - u[e])/h
                εθθ = 0.5*(u[e] + u[e+1])/rm
                solid_constitutive_update!(σ, D, @view(Htr[:, e]), model, εrr, εθθ, 0.0, @view(H[:, e]), 0.0)
                # B = [dφ₁/dr dφ₂/dr; φ₁/r φ₂/r] at the midpoint; weight r dr
                w  = rm*h
                b1 = (-1.0/h, 0.5/rm); b2 = (1.0/h, 0.5/rm)
                fint[e]   += w*(b1[1]*σ[1] + b1[2]*σ[2])
                fint[e+1] += w*(b2[1]*σ[1] + b2[2]*σ[2])
                k(bi, bj) = w*(bi[1]*(D[1,1]*bj[1] + D[1,2]*bj[2]) + bi[2]*(D[2,1]*bj[1] + D[2,2]*bj[2]))
                dd[e]   += k(b1, b1); du[e] += k(b1, b2)
                dl[e]   += k(b2, b1); dd[e+1] += k(b2, b2)
            end
            res = -fint[2:end]                              # no load on the free nodes
            # Round-off floor of this residual with 2000 elements: ~1e-10 of |f|
            # (nodal sums of nearly cancelling terms). Converged well above it.
            rel = norm(res)/max(norm(fint), floatmin(Float64))
            rel <= 1.0e-9 && break
            it == 50 && rel > 1.0e-7 &&
                error("plastic_cylinder_radial_reference: Newton did not converge at δ = $δ (|r|/|f| = $rel)")
            it == 50 && break
            Δu = Tridiagonal(dl[2:end], dd[2:end], du[2:end]) \ res
            u[2:end] .+= Δu
        end
        H .= Htr
        push!(ps, fint[1]/P.a)                              # bore reaction per radian / a
        push!(epb, H[5, 1])
    end
    return (; p = ps, ep_bore = epb)
end
