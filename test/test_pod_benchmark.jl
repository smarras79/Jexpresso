#---------------------------------------------------------------------------------
# test/test_pod_benchmark.jl — THE REFERENCE BENCHMARK.
#
#   julia test/test_pod_benchmark.jl
#
# test/test_pod.jl checks that the code does what the algorithm says. This file
# checks something stronger: that on a problem whose POD is known in CLOSED FORM
# — the one the model-reduction literature uses to say what POD can and cannot
# do — Jexpresso returns that closed form, to round-off, including the numbers a
# reduced-order model would be built from.
#
# THE BENCHMARK: LINEAR ADVECTION OF A MULTI-HARMONIC WAVE.
#
#     ∂u/∂t + c ∂u/∂x = 0 ,   x ∈ [0,L) periodic,
#     u(x,0) = Σ_{j=1}^J A_j cos(2πj x/L + ϕ_j)   ⟹   u(x,t) = u(x - ct, 0).
#
# This is the standard test problem of the transport-dominated model-reduction
# literature (Ohlberger & Rave 2016; Greif & Urban 2019; and the travelling-wave
# discussion of Holmes, Lumley, Berkooz & Rowley, 2nd ed., §3.3). It is used
# because the whole decomposition can be written down, and because what it says
# is the sharpest known limitation of POD (item 5 below).
#
# THE CLOSED FORM. Average over one period of the translation:
#
#     R(x,x') = ⟨u'(x,t) u'(x',t)⟩ = Σ_j (A_j²/2) cos(2πj(x-x')/L) ,          (B1)
#
# a convolution kernel, so its eigenfunctions are the Fourier modes and, since
# cos(2πj(x-x')/L) = cos·cos + sin·sin, each wavenumber contributes a TWO-
# dimensional eigenspace:
#
#     λ_{2j-1} = λ_{2j} = A_j² L / 4 ,   span{ cos(2πjx/L), sin(2πjx/L) } ,   (B2)
#     E_j = λ_j / Σλ = A_j² / (2 Σ_i A_i²) ,                                  (B3)
#     ⟨‖u'‖²⟩ = Σλ = (L/2) Σ_j A_j² ,                                         (B4)
#     ε(r = 2m)² = Σ_{j>m} A_j² / Σ_j A_j²        (the truncation error).     (B5)
#
# Everything below is a comparison against (B1)-(B5).
#
# WHAT EACH ITEM CATCHES, because a benchmark that only checks the eigenvalues
# is weaker than it looks:
#
#   1. THE SPECTRUM, against (B2)-(B4) to 1e-10 relative. Catches a wrong
#      normalisation of the correlation matrix, a missing 1/K, a mean that was
#      not removed.
#   2. THE EIGENSPACES, not the modes. λ_{2j-1} = λ_{2j} EXACTLY, so the two
#      members of a pair are defined only up to a rotation between them — any
#      code that claims a particular pair of modes there is claiming something
#      the problem does not determine. What IS determined is the plane they
#      span, and that is what is checked: cos(2πjx/L) and sin(2πjx/L) must lie
#      in it to round-off.
#   3. THE INNER PRODUCT, by running the whole thing on a deliberately NON-
#      UNIFORM grid (Chebyshev-Lobatto, clustered towards element edges like the
#      LGL nodes of the real solver) and asserting BOTH that the mass-weighted
#      answer is right AND that the unweighted one is measurably wrong. The
#      second half is the point: on a uniform grid the two agree, so a uniform
#      benchmark cannot tell a correct implementation from one that silently
#      dropped the mass matrix.
#   4. THE ROM QUANTITIES: the truncation error curve (B5) — which is the a
#      priori error bound a reduced model is chosen with — and the phase
#      portrait of each pair, which must be a circle of constant radius because
#      the structure travels.
#   5. THE KOLMOGOROV n-WIDTH. With every A_j equal, (B3) makes the spectrum
#      FLAT: 2J modes each carrying 1/(2J) of the energy, and no truncation
#      captures more than its share. That is the known, and famous, failure of
#      linear model reduction for transport, and a POD implementation that
#      reported a decaying spectrum here would be wrong in the most consequential
#      way possible — it would promise a reduced model that cannot exist.
#
# Needs `Test` and nothing else: the two files under test are free of Jexpresso
# types by design. `problems/AdvDiff/PODbenchmark` is the same problem run
# through the solver, deck and all, with this file's closed form in its README.
#---------------------------------------------------------------------------------

using Test
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "src", "kernel", "rom", "pod_core.jl"))


#
# A 1-D spectral-element grid, open at both ends, with Chebyshev-Lobatto nodes
# inside each element — NON-UNIFORM, clustered at element edges, like the LGL
# nodes of the solver. `w` is the trapezoid mass matrix of that grid, which
# integrates any linear function exactly and sums to L.
#
function sem_grid_1d(nelem::Int, ngl::Int, L::Float64)
    npoin = nelem*(ngl-1) + 1
    x     = zeros(npoin)
    h     = L/nelem
    ξ     = [cos(π*(ngl-1-k)/(ngl-1)) for k = 0:ngl-1]
    for iel = 1:nelem, i = 1:ngl
        x[(iel-1)*(ngl-1) + i] = (iel-1)*h + h*(ξ[i] + 1)/2
    end
    w = zeros(npoin)
    for ip = 1:npoin-1
        d = x[ip+1] - x[ip]
        w[ip] += d/2; w[ip+1] += d/2
    end
    return x, w
end

#
# The exact solution, sampled at K times that tile exactly ONE period of the
# translation with no repeated phase: t_k = k·T/K, k = 0…K-1.
#
# Both ends of one period would REPEAT the zero phase, and that single repeat
# splits every degenerate pair by (K/2+1)/(K/2) — 5 % at K = 40. It is a
# property of the sampling and not of the decomposition, and the benchmark is
# written to avoid it rather than to tolerate it. (The same care is why
# problems/AdvDiff/PODbenchmark sets :pod_tend one interval short of :tend.)
#
function advected_wave_snapshots(x, A, ϕ, L, c, K)
    npoin = length(x)
    T     = L/c
    X     = Array{Float64}(undef, npoin, 1, K)
    t     = Vector{Float64}(undef, K)
    for k = 1:K
        t[k] = (k-1)*T/K
        for ip = 1:npoin
            X[ip,1,k] = sum(A[j]*cos(2π*j*(x[ip] - c*t[k])/L + ϕ[j]) for j in eachindex(A))
        end
    end
    return X, t
end

minner(w, u, v) = sum(w[ip]*u[ip]*v[ip] for ip in eachindex(w))


@testset verbose = true "POD reference benchmark: advection of a multi-harmonic wave" begin

L = 2.0
c = 0.5
A = [1.0, 0.5, 0.25]                # a spectrum that decays by 4× per pair
ϕ = [0.0, 0.7, -1.3]
K = 40

x, w  = sem_grid_1d(25, 5, L)
X, t  = advected_wave_snapshots(x, A, ϕ, L, c, K)

# the grid itself, before anything is decomposed on it
@test sum(w) ≈ L rtol = 1e-12
@test maximum(diff(x))/minimum(diff(x)) > 2      # genuinely non-uniform

P = pod_from_snapshots(X, w, t; name = "u")

#---------------------------------------------------------------------------------
@testset "1. the spectrum, against the closed form (B2)-(B4)" begin

    @test length(P.λ) == 2*length(A)             # two modes per wavenumber, and no more

    for j in eachindex(A)
        λexact = A[j]^2*L/4
        @test P.λ[2j-1] ≈ λexact rtol = 1e-10
        @test P.λ[2j]   ≈ λexact rtol = 1e-10
        Eexact = A[j]^2/(2*sum(abs2, A))
        @test P.energy[2j-1] ≈ Eexact rtol = 1e-10
        @test P.energy[2j]   ≈ Eexact rtol = 1e-10
    end

    @test P.total ≈ (L/2)*sum(abs2, A) rtol = 1e-10
    @test P.cumenergy[end] ≈ 1.0 atol = 1e-12
    # the mean of a travelling wave over a period is zero, and POD must find that
    @test maximum(abs, P.q̄) < 1e-12
    @test P.orthoerr < 1e-10
end

#---------------------------------------------------------------------------------
@testset "2. the eigenSPACES, which is what the problem determines" begin

    for j in eachindex(A)
        # λ_{2j-1} = λ_{2j}, so the pair is a plane and not two modes
        @test P.λ[2j-1] ≈ P.λ[2j] rtol = 1e-10

        φ1 = view(P.Φ, :, 1, 2j-1)
        φ2 = view(P.Φ, :, 1, 2j)
        for f in (Float64[cos(2π*j*xi/L) for xi in x], Float64[sin(2π*j*xi/L) for xi in x])
            nf  = sqrt(minner(w, f, f))
            res = f .- minner(w, f, φ1).*φ1 .- minner(w, f, φ2).*φ2
            @test sqrt(minner(w, res, res))/nf < 1e-8          # lies in the plane
        end
        # …and nothing of a DIFFERENT wavenumber does
        if j < length(A)
            g  = Float64[cos(2π*(j+1)*xi/L) for xi in x]
            ng = sqrt(minner(w, g, g))
            @test abs(minner(w, g, φ1))/ng < 1e-8
            @test abs(minner(w, g, φ2))/ng < 1e-8
        end
    end
end

#---------------------------------------------------------------------------------
@testset "3. the inner product has to be the mass matrix" begin

    # the weighted answer is the closed form …
    @test P.energy[1] ≈ A[1]^2/(2*sum(abs2, A)) rtol = 1e-10

    # … and the unweighted one, on this non-uniform grid, is not. Same snapshots,
    # same code, Euclidean weights: the clustered nodes count for several times
    # their share of the domain and the energies come out wrong. If this ever
    # starts agreeing, the benchmark has stopped being able to see the
    # difference — check the grid is still non-uniform before believing it.
    Pe   = pod_from_snapshots(X, ones(length(w)), t; name = "euclidean")
    worst = maximum(abs(Pe.energy[i] - P.energy[i]) for i = 1:length(P.λ))
    @test worst > 1e-3
    @test Pe.orthoerr < 1e-10                    # …it is still an orthonormal basis,
    @test !isapprox(Pe.energy[1], P.energy[1]; rtol = 1e-6)  # just of the wrong inner product
end

#---------------------------------------------------------------------------------
@testset "4. what a reduced-order model is built from" begin

    #--- the a-priori truncation error, Eq. (B5)
    e   = pod_truncation_error(P)
    tot = sum(abs2, A)
    for m = 0:length(A)
        eexact = sqrt(sum(abs2, A[m+1:end])/tot)
        @test e[2m+1] ≈ eexact atol = 1e-9
    end

    #--- and it is the error actually made, not just the one predicted
    for m = 0:length(A)
        num = 0.0; den = 0.0
        for k = 1:K
            q   = pod_reconstruct(P, P.a[k,1:2m])
            res = X[:,1,k] .- q[:,1]
            num += minner(w, res, res)
            den += minner(w, X[:,1,k], X[:,1,k])
        end
        @test sqrt(num/den) ≈ sqrt(sum(abs2, A[m+1:end])/tot) atol = 1e-8
    end

    #--- each pair is one TRAVELLING structure: constant radius in its own phase
    #    plane, and in quadrature
    for j in eachindex(A)
        rad = [sqrt(P.a[k,2j-1]^2 + P.a[k,2j]^2) for k = 1:K]
        @test (maximum(rad) - minimum(rad)) < 1e-8*maximum(rad)
        @test abs(sum(P.a[:,2j-1].*P.a[:,2j])) < 1e-8*sum(abs2, P.a[:,2j-1])
    end

    #--- the reduced coordinates are the projections
    @test maximum(abs, pod_project(P, w, X[:,:,7]) .- P.a[7,:]) < 1e-9
end

#---------------------------------------------------------------------------------
@testset "5. the Kolmogorov n-width: POD cannot compress a travelling wave" begin

    #
    # Equal amplitudes: (B3) makes every mode carry 1/(2J) of the energy, so the
    # spectrum is FLAT and truncation buys nothing — keeping half the modes
    # leaves exactly half the energy behind. This is the known limitation of
    # linear model reduction for transport, and the benchmark asserts that the
    # code reports it rather than a comfortable decaying spectrum.
    #
    J    = 6
    Aeq  = ones(J)
    Xeq, teq = advected_wave_snapshots(x, Aeq, zeros(J), L, c, 4J)
    Peq  = pod_from_snapshots(Xeq, w, teq; name = "flat")

    @test length(Peq.λ) == 2J
    for i = 1:2J
        @test Peq.energy[i] ≈ 1/(2J) rtol = 1e-8
    end
    @test pod_rank_for_energy(Peq, 0.99) == 2J        # 99 % needs EVERY mode
    eq = pod_truncation_error(Peq)
    @test eq[J+1] ≈ sqrt(0.5) atol = 1e-8             # half the modes, half the energy

    @printf("\n   flat-spectrum case: %d modes, each %.4f %% of the energy; ",
            length(Peq.λ), 100*Peq.energy[1])
    @printf("rank for 99%%: %d of %d\n", pod_rank_for_energy(Peq, 0.99), length(Peq.λ))
end

#---------------------------------------------------------------------------------
# What the benchmark measured, printed so that a run of this file is also a
# report and not only a pass/fail.
#---------------------------------------------------------------------------------
println()
println("   POD of the advected multi-harmonic wave, against the closed form:")
println("   mode     λ (computed)      λ (exact)      E [%]    E exact [%]")
for i = 1:length(P.λ)
    j = (i+1) ÷ 2
    @printf("   %4d   %14.10f  %14.10f  %8.4f  %8.4f\n",
            i, P.λ[i], A[j]^2*L/4, 100*P.energy[i], 100*A[j]^2/(2*sum(abs2, A)))
end
@printf("   Σλ = %.10f   (exact %.10f)   max|ΦᵀMΦ - I| = %.2e\n",
        P.total, (L/2)*sum(abs2, A), P.orthoerr)
println()

end # testset
