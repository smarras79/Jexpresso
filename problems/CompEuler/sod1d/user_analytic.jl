#---------------------------------------------------------------------------------
# Exact solution of the Sod shock tube, for overlay on the 1D output plots.
#
# The shared 1D plotter (src/io/plotting/jeplots.jl) calls
#
#     user_analytic_solution(x, t, outvar, inputs) -> npoin × nvar matrix
#
# if — and only if — the case being run defines it. Most 1D cases have no
# closed form, so nothing is hardcoded in the plotter; this file supplies the
# solution for sod1d alone and simply does not exist for anything else.
#
# Return NaN for any entry that cannot be supplied: Plots skips non-finite
# points, so a case may provide an exact solution for some variables and not
# others.
#
# Problem (matching initialize.jl):
#
#     ρ,u,p = (1.000, 0, 1.0)    x < 0.5
#     ρ,u,p = (0.125, 0, 0.1)    x > 0.5
#
# Structure at t > 0, left to right: left rarefaction fan, contact
# discontinuity, right shock. The solution is self-similar in ξ = (x-x₀)/t.
#
# Reference: E. F. Toro, "Riemann Solvers and Numerical Methods for Fluid
# Dynamics", 3rd ed., Springer 2009, Ch. 4 (exact Riemann solver) — the
# star-region pressure is found by Newton iteration on the pressure function
# f(p) = f_L(p) + f_R(p) + (u_R - u_L) = 0.
#---------------------------------------------------------------------------------

# --- left/right states and the diaphragm, as in initialize.jl -------------
const SOD_ρL, SOD_uL, SOD_pL = 1.000, 0.0, 1.0
const SOD_ρR, SOD_uR, SOD_pR = 0.125, 0.0, 0.1
const SOD_X0 = 0.5

#
# Toro eq. (4.6)-(4.7): the pressure function and its derivative for one side.
#
@inline function _sod_f(p, ρk, pk, ck, γ)
    if p > pk
        # shock branch
        Ak = 2.0/((γ + 1.0)*ρk)
        Bk = (γ - 1.0)/(γ + 1.0)*pk
        f  = (p - pk)*sqrt(Ak/(p + Bk))
        df = sqrt(Ak/(Bk + p))*(1.0 - 0.5*(p - pk)/(Bk + p))
    else
        # rarefaction branch
        f  = 2.0*ck/(γ - 1.0)*((p/pk)^((γ - 1.0)/(2.0*γ)) - 1.0)
        df = 1.0/(ρk*ck)*(p/pk)^(-(γ + 1.0)/(2.0*γ))
    end
    return f, df
end

#
# Star-region pressure and velocity. Newton with a two-rarefaction initial
# guess, which is positive by construction and converges in a handful of
# iterations for this problem.
#
function _sod_star(γ)
    cL = sqrt(γ*SOD_pL/SOD_ρL)
    cR = sqrt(γ*SOD_pR/SOD_ρR)

    p = ((cL + cR - 0.5*(γ - 1.0)*(SOD_uR - SOD_uL)) /
         (cL/SOD_pL^((γ - 1.0)/(2.0*γ)) + cR/SOD_pR^((γ - 1.0)/(2.0*γ))))^((2.0*γ)/(γ - 1.0))
    p = max(p, 1.0e-12)

    for _ = 1:100
        fL, dfL = _sod_f(p, SOD_ρL, SOD_pL, cL, γ)
        fR, dfR = _sod_f(p, SOD_ρR, SOD_pR, cR, γ)
        Δ = (fL + fR + (SOD_uR - SOD_uL))/(dfL + dfR)
        p_new = max(p - Δ, 1.0e-12)
        if abs(p_new - p)/(0.5*(p_new + p)) < 1.0e-12
            p = p_new
            break
        end
        p = p_new
    end

    fL, _ = _sod_f(p, SOD_ρL, SOD_pL, cL, γ)
    fR, _ = _sod_f(p, SOD_ρR, SOD_pR, cR, γ)
    u = 0.5*(SOD_uL + SOD_uR) + 0.5*(fR - fL)

    return p, u, cL, cR
end

#
# Sample the self-similar solution at ξ = (x - x₀)/t.  Toro Ch. 4.5.
#
function _sod_sample(ξ, pstar, ustar, cL, cR, γ)
    if ξ <= ustar
        # ---- left of the contact -------------------------------------
        if pstar > SOD_pL
            # left shock
            S = SOD_uL - cL*sqrt((γ + 1.0)/(2.0*γ)*pstar/SOD_pL + (γ - 1.0)/(2.0*γ))
            if ξ <= S
                return SOD_ρL, SOD_uL, SOD_pL
            else
                ρ = SOD_ρL*((pstar/SOD_pL + (γ - 1.0)/(γ + 1.0)) /
                            ((γ - 1.0)/(γ + 1.0)*pstar/SOD_pL + 1.0))
                return ρ, ustar, pstar
            end
        else
            # left rarefaction fan
            SH = SOD_uL - cL
            cstar = cL*(pstar/SOD_pL)^((γ - 1.0)/(2.0*γ))
            ST = ustar - cstar
            if ξ <= SH
                return SOD_ρL, SOD_uL, SOD_pL
            elseif ξ >= ST
                ρ = SOD_ρL*(pstar/SOD_pL)^(1.0/γ)
                return ρ, ustar, pstar
            else
                # inside the fan
                u = 2.0/(γ + 1.0)*(cL + 0.5*(γ - 1.0)*SOD_uL + ξ)
                c = 2.0/(γ + 1.0)*(cL + 0.5*(γ - 1.0)*(SOD_uL - ξ))
                ρ = SOD_ρL*(c/cL)^(2.0/(γ - 1.0))
                p = SOD_pL*(c/cL)^((2.0*γ)/(γ - 1.0))
                return ρ, u, p
            end
        end
    else
        # ---- right of the contact ------------------------------------
        if pstar > SOD_pR
            # right shock
            S = SOD_uR + cR*sqrt((γ + 1.0)/(2.0*γ)*pstar/SOD_pR + (γ - 1.0)/(2.0*γ))
            if ξ >= S
                return SOD_ρR, SOD_uR, SOD_pR
            else
                ρ = SOD_ρR*((pstar/SOD_pR + (γ - 1.0)/(γ + 1.0)) /
                            ((γ - 1.0)/(γ + 1.0)*pstar/SOD_pR + 1.0))
                return ρ, ustar, pstar
            end
        else
            # right rarefaction fan
            SH = SOD_uR + cR
            cstar = cR*(pstar/SOD_pR)^((γ - 1.0)/(2.0*γ))
            ST = ustar + cstar
            if ξ >= SH
                return SOD_ρR, SOD_uR, SOD_pR
            elseif ξ <= ST
                ρ = SOD_ρR*(pstar/SOD_pR)^(1.0/γ)
                return ρ, ustar, pstar
            else
                u = 2.0/(γ + 1.0)*(-cR + 0.5*(γ - 1.0)*SOD_uR + ξ)
                c = 2.0/(γ + 1.0)*(cR - 0.5*(γ - 1.0)*(SOD_uR - ξ))
                ρ = SOD_ρR*(c/cR)^(2.0/(γ - 1.0))
                p = SOD_pR*(c/cR)^((2.0*γ)/(γ - 1.0))
                return ρ, u, p
            end
        end
    end
end

#
# Hook called by plot_results(::NSD_1D). outvar carries the names of the
# plotted variables — here the conservative set ("ρ", "ρu", "ρE") — so the
# exact primitive state is mapped onto whichever of them is being drawn.
# Anything unrecognised is left NaN and simply not plotted.
#
function user_analytic_solution(x, t, outvar, inputs)

    γ = PhysicalConst{Float64}().γ
    nvar = length(outvar)
    q = fill(NaN, length(x), nvar)

    # t = 0 is the discontinuous initial state; sampling ξ = ±Inf is exact
    # but pointless to overlay, so decline.
    if t <= 0.0
        return nothing
    end

    pstar, ustar, cL, cR = _sod_star(γ)

    for i in eachindex(x)
        ξ = (x[i] - SOD_X0)/t
        ρ, u, p = _sod_sample(ξ, pstar, ustar, cL, cR, γ)
        E = p/(γ - 1.0) + 0.5*ρ*u*u
        for ivar = 1:nvar
            name = string(outvar[ivar])
            q[i, ivar] = name == "ρ"  ? ρ     :
                         name == "ρu" ? ρ*u   :
                         name == "ρE" ? E     :
                         name == "u"  ? u     :
                         name == "p"  ? p     : NaN
        end
    end

    return q
end
