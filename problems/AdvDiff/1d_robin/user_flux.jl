# Pure diffusion: no advective flux. (TOTAL and PERT spelled out so this file overrides methods a
# case run earlier in the same session may have defined.)
function user_flux!(F, G, SD::NSD_1D, q, qe, mesh::St_mesh, ::CL, ::TOTAL; neqs = 1, ip = 1)
    F[1] = 0.0
end
function user_flux!(F, G, SD::NSD_1D, q, qe, mesh::St_mesh, ::CL, ::PERT; neqs = 1, ip = 1)
    F[1] = 0.0
end
