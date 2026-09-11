#=
Reference solution of the Brio-Wu MHD shock tube (Brio & Wu 1988) for
problems/MHD/brioWu1d, computed with a first-order finite-volume scheme
(HLL flux, forward Euler, CFL 0.4) on a very fine uniform grid — the same
kind of reference Dao & Nazarov (2022, JSC 92:77, Sec. 5.2) take from Athena
at 10 000 points. Ideal MHD with Bx constant, magnetic pressure B²/2, γ = 2,
domain (0, 1), final time 0.2.

    julia tools/brio_wu_reference.jl [ncells=10000] [nout=2000]

writes problems/MHD/brioWu1d/reference_hll.dat with columns
x ρ u v w p By Bz (nout points, cell averages sampled uniformly).
=#
const γ  = 2.0
const Bx = 0.75

function cons(ρ, u, v, w, p, By, Bz)
    E = p/(γ - 1) + 0.5*ρ*(u*u + v*v + w*w) + 0.5*(Bx*Bx + By*By + Bz*Bz)
    return (ρ, ρ*u, ρ*v, ρ*w, E, By, Bz)
end

@inline function prim(U)
    ρ, mu, mv, mw, E, By, Bz = U
    u = mu/ρ; v = mv/ρ; w = mw/ρ
    p = (γ - 1)*(E - 0.5*ρ*(u*u + v*v + w*w) - 0.5*(Bx*Bx + By*By + Bz*Bz))
    return ρ, u, v, w, p, By, Bz
end

@inline function flux(U)
    ρ, u, v, w, p, By, Bz = prim(U)
    B2 = Bx*Bx + By*By + Bz*Bz
    pt = p + 0.5*B2
    E  = U[5]
    uB = u*Bx + v*By + w*Bz
    return (ρ*u, ρ*u*u + pt - Bx*Bx, ρ*u*v - Bx*By, ρ*u*w - Bx*Bz,
            (E + pt)*u - Bx*uB, By*u - Bx*v, Bz*u - Bx*w)
end

@inline function cfast(U)
    ρ, u, v, w, p, By, Bz = prim(U)
    a2 = γ*max(p, 0.0)/ρ
    b2 = (Bx*Bx + By*By + Bz*Bz)/ρ
    bx2 = Bx*Bx/ρ
    return sqrt(0.5*(a2 + b2 + sqrt(max((a2 + b2)^2 - 4*a2*bx2, 0.0))))
end

@inline function hll(UL, UR)
    uL = UL[2]/UL[1]; uR = UR[2]/UR[1]
    cL = cfast(UL);   cR = cfast(UR)
    SL = min(uL - cL, uR - cR); SR = max(uL + cL, uR + cR)
    FL = flux(UL); FR = flux(UR)
    if SL >= 0
        return FL
    elseif SR <= 0
        return FR
    else
        return ntuple(k -> (SR*FL[k] - SL*FR[k] + SL*SR*(UR[k] - UL[k]))/(SR - SL), 7)
    end
end

function run(ncells, nout, tend)
    dx = 1.0/ncells
    x  = [(i - 0.5)*dx for i = 1:ncells]
    UL = cons(1.0,   0.0, 0.0, 0.0, 1.0,  1.0, 0.0)
    UR = cons(0.125, 0.0, 0.0, 0.0, 0.1, -1.0, 0.0)
    U  = [xi < 0.5 ? UL : UR for xi in x]
    t = 0.0
    F = Vector{NTuple{7,Float64}}(undef, ncells + 1)
    while t < tend - 1e-15
        smax = maximum(abs(Ui[2]/Ui[1]) + cfast(Ui) for Ui in U)
        dt = min(0.4*dx/smax, tend - t)
        F[1] = flux(U[1]); F[end] = flux(U[end])
        @inbounds for i = 2:ncells
            F[i] = hll(U[i-1], U[i])
        end
        @inbounds for i = 1:ncells
            U[i] = ntuple(k -> U[i][k] - dt/dx*(F[i+1][k] - F[i][k]), 7)
        end
        t += dt
    end
    out = joinpath(@__DIR__, "..", "problems", "MHD", "brioWu1d", "reference_hll.dat")
    open(out, "w") do io
        println(io, "# Brio-Wu, gamma = 2, Bx = 0.75, t = $tend: HLL finite volume, $ncells cells, sampled at $nout points")
        println(io, "# x rho u v w p By Bz")
        for j = 1:nout
            i = clamp(round(Int, (j - 0.5)/nout*ncells + 0.5), 1, ncells)
            ρ, u, v, w, p, By, Bz = prim(U[i])
            println(io, join((x[i], ρ, u, v, w, p, By, Bz), " "))
        end
    end
    println("wrote ", out)
end

ncells = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10000
nout   = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2000
tend   = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.1
run(ncells, nout, tend)
