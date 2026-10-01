# The linear solve uses no fluxes; these are the empty hooks the setup expects.
function user_flux!(F, G, H, q, qe, mesh::St_mesh, ::CL, ::TOTAL; neqs=1, ip=1)
    F[1] = 0.0; G[1] = 0.0; H[1] = 0.0
end
function user_flux!(F, G, H, q, qe, mesh::St_mesh, ::CL, ::PERT; neqs=1, ip=1)
    F[1] = 0.0; G[1] = 0.0; H[1] = 0.0
end
