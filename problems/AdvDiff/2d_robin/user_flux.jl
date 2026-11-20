# Advective flux (F, G) = (c u, 0): uniform advection along the periodic x direction, so the walls
# carry no advective flux and the Robin condition is the only thing acting there.
# (TOTAL, PERT and NCL are all spelled out so this file overrides whatever a case run earlier in
# the same session defined; see the note in problems/AdvDiff/kopriva/user_source.jl.)

function user_flux!(F, G, SD::NSD_2D, q, qe, mesh::St_mesh, ::CL, ::TOTAL; neqs = 1, ip = 1)
    F[1] = ROBIN_C[] * q[1]
    G[1] = 0.0
end

function user_flux!(F, G, SD::NSD_2D, q, qe, mesh::St_mesh, ::CL, ::PERT; neqs = 1, ip = 1)
    F[1] = ROBIN_C[] * (q[1] + qe[1])
    G[1] = 0.0
end

function user_flux!(F, G, SD::NSD_2D, q, qe, mesh::St_mesh, ::NCL, ::AbstractPert; neqs = 1, ip = 1)
    F[1] = ROBIN_C[] * q[1]
    G[1] = 0.0
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    return T(q[1] * ROBIN_C[]), T(0.0)
end
