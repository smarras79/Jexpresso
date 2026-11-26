#---------------------------------------------------------------------------------
# TC6: conservative non-linear shallow water fluxes (paper Eq. 41, flat bottom).
# State q = [h, hu, hv]; constants TC6_* are defined in initialize.jl.
#   F = [hu, hu²/h + g h²/2, hu hv/h]
#   G = [hv, hu hv/h,        hv²/h + g h²/2]
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=3, ip=1)
    h  = q[1]
    hu = q[2]
    hv = q[3]
    u  = hu / h
    v  = hv / h
    p  = 0.5 * TC6_G * h * h

    F[1] = hu
    F[2] = hu * u + p
    F[3] = hu * v

    G[1] = hv
    G[2] = hv * u
    G[3] = hv * v + p
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::PERT; neqs=3, ip=1)
    error(" # ERROR TC6: only SOL_VARS_TYPE => TOTAL() is implemented.")
end
