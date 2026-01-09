function user_source!(S, q, qe, npoin, ::CL, ::TOTAL; 
                     neqs=1, x=0.0, y=0.0, 
                     ymin=0.0, ymax=15000.0,      # ← FIX!
                     ngl=5, nely=10, 
                     xmin=-80000.0, xmax=80000.0)  # ← FIX!
   
    PhysConst = PhysicalConst{Float64}()
    
    ρ = q[1] - qe[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ * PhysConst.g
    S[4] = 0.0
    S[5] = 0.0
    S[6] = 0.0

    # SPONGE
    zs = 10000.0      # Start sponge at 10 km
    xr = 80000.0      # ← FIX!
    xl = -80000.0     # ← FIX!
    
    if (y > zs)
        betay_coe = sinpi(0.5 * (y - zs) / (ymax - zs))^2
    else
        betay_coe = 0.0
    end
    ctop = 0.1 * betay_coe
   
    if (x > xr)
        betaxr_coe = sinpi(0.5 * (x - xr) / (xmax - xr))^2
    else
        betaxr_coe = 0.0
    end
   
    if (x < xl)
        betaxl_coe = sinpi(0.5 * (xl - x) / (xl - xmin))^2
    else
        betaxl_coe = 0.0
    end
   
    cxr = 0.01 * betaxr_coe
    cxl = 0.01 * betaxl_coe
    cs = 1.0 - (1.0 - ctop) * (1.0 - cxr) * (1.0 - cxl)

    if (x >= xmin && x <= xmax)     
      S[1] -= cs * (q[1] - qe[1])
      S[2] -= cs * (q[2] - qe[2])  # ← FIX! (was qe[2]*20.0)
      S[3] -= cs * q[3]
      S[4] -= cs * (q[4] - qe[4])
      S[5] -= cs * (q[5] - qe[5])  # ← ADD!
      S[6] -= cs * (q[6] - qe[6])  # ← ADD!
    end
    
    return S
end


# being modified to do TOTAL
function user_source!(S, q, qe, npoin, ::CL, ::PERT; neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=15000.0, ngl=5, nely=10, xmin=-80000, xmax=80000, ctop_mult=0.1)  # CHANGED: ctop_mult=0.1 (from 0.2), ymax and xmin/xmax

    PhysConst = PhysicalConst{Float64}()

    #
    # S(q(x)) = -ρg
    #
    ρ = q[1]

    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ * PhysConst.g
    S[4] = 0.0
    S[5] = 0.0
    S[6] = 0.0
    #### SPONGE (updated: zs=10km, mult=0.1 for hill2d)

    nsponge_points = 8

    # distance from the boundary. xs in Restelli's thesis
    dsy = (ymax - ymin) / (nely * (ngl - 1))  # equivalent grid spacing
    dbl = ymax - y
    zs = 10000.0  # CHANGED: Start sponge at 10 km (from 16 km)
    dsx = (xmax - xmin) / (nely * (ngl - 1))  # equivalent grid spacing
    dbx = min(xmax - x, x - xmin)
    xr = xmax  # CHANGED: Use domain edge (remove interior ±40km for standard hill2d boundaries)
    xl = xmin  # CHANGED: Use domain edge
    
    if (y >= zs)  # nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe = sinpi(0.5 * (y - zs) / (ymax - zs))^2  # CHANGED: Use ymax, squared sine
    else
        betay_coe = 0.0
    end
    #if (abs(x) <= xmin)
      ctop = betay_coe  # CHANGED: Direct from betay (no extra 0.5)
    #else
     # ctop = 0.0
    #end 

    if (x > xr)  # nsponge_points * dsy) #&& dbl >= 0.0)
        betaxr_coe = sinpi(0.5 * (x - xr) / (xmax - xr))^2  # CHANGED: Squared sine
    else
        betaxr_coe = 0.0
    end

    if (x < xl)  # nsponge_points * dsy) #&& dbl >= 0.0)
        betaxl_coe = sinpi(0.5 * (xl - x) / (xl - xmin))^2  # CHANGED: Squared sine
    else
        betaxl_coe = 0.0
    end
    
    cxr = 0.01 * betaxr_coe  # CHANGED: Minimal lateral (0.01)
    cxl = 0.01 * betaxl_coe  # CHANGED: Minimal lateral
    ctop = ctop_mult * min(ctop, 1.0)  # CHANGED: Use param (now 0.1 max), cap at 1
    cxr = min(cxr, 1.0)
    cxl = min(cxl, 1.0)
    cs = 1.0 - (1.0 - ctop) * (1.0 - cxr) * (1.0 - cxl)
    
    #@info "β x: " ctop,cxr,cxl,cs, zs, y, x, ymin, ymax, dsy, dbl
    S[1] -= (cs) * (q[1])
    S[2] -= (cs) * (q[2])
    S[3] -= (cs) * (q[3])
    S[4] -= (cs) * (q[4])
    S[5] -= (cs) * (q[5] - qe[5])
    S[6] -= (cs) * (q[6] - qe[6])
    return S
end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)

    T = eltype(x)
    # distance from the boundary. xs in Restelli's thesis
    zs = T(10000.0)  # CHANGED: Start sponge at 10 km (from 15 km)
    xr = T(xmax)  # CHANGED: Use domain edge
    xl = T(xmin)  # CHANGED: Use domain edge

    if (y >= zs)  # nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe = T(sinpi(T(0.5) * (y - zs) / (ymax - zs))^2)  # CHANGED: Squared sine, ymax
    else
        betay_coe = T(0.0)
    end
    ctop = T(0.1) * betay_coe  # CHANGED: Explicit 0.1 mult for hill2d
    
    if (x > xr)  # nsponge_points * dsy) #&& dbl >= 0.0)
        betaxr_coe = T(sinpi(T(0.5) * (x - xr) / (xmax - xr))^2)  # CHANGED: Squared sine
    else
        betaxr_coe = T(0.0)
    end

    if (x < xl)  # nsponge_points * dsy) #&& dbl >= 0.0)
        betaxl_coe = T(sinpi(T(0.5) * (xl - x) / (xl - xmin))^2)  # CHANGED: Squared sine
    else
        betaxl_coe = T(0.0)
    end

    cxr = T(0.01) * betaxr_coe  # CHANGED: Minimal lateral (0.01)
    cxl = T(0.01) * betaxl_coe  # CHANGED: Minimal lateral
    ctop = T(0.1) * min(ctop, T(1.0))  # CHANGED: Cap at 0.1 max
    cxr = min(cxr, T(1.0))
    cxl = min(cxl, T(1.0))
    cs = T(1.0) - (T(1.0) - ctop) * (T(1.0) - cxr) * (T(1.0) - cxl)

    #
    # S(q(x)) = -ρg
    #
    ρ = q[1]

    # CHANGED: Align with TOTAL (4 eqs; add if neqs=5/6)
    return T(-cs * q[1]), T(-cs * q[2]), T(-cs * q[3] - ρ * PhysConst.g), T(-cs * q[4])
end