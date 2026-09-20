#
# Primitive helpers for the non-linear shallow water solver (AMR variant).
#
# uprimitive is the quantity differentiated by the viscous Laplacian.
# The continuity equation diffuses the depth PERTURBATION H - He
# (He = qe[1] is the lake-at-rest depth defined in initialize.jl): at
# rest the depth itself is cone-shaped (∇²H ≠ 0 over the island), so
# diffusing the full H would pump mass in/out of the wet/dry ring and
# break the discrete equilibrium that user_flux.jl/user_source.jl
# preserve. The momentum components are diffused in conservation form;
# their reference state is zero, so no subtraction is needed.
#
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[1] - qe[1]   # H - He (free-surface perturbation)
    uprimitive[2] = u[2]           # Hu
    uprimitive[3] = u[3]           # Hv
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    user_primitives!(u, qe, uprimitive, TOTAL())
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    return T(u[1] - qe[1]), T(u[2]), T(u[3])
end

function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    uout[1] = u[1]
    uout[2] = u[2]
    uout[3] = u[3]
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    uout[1] = u[1]
    uout[2] = u[2]
    uout[3] = u[3]
end

#
# AMR restart hook (docs/amr_setup.md section 4). read_vtk_amr_restart! (called
# from initialize.jl when :lrestart_amr is true) reads the VTK point data saved
# for qoutvars = ["H", "Hu", "Hv"] and calls this to rebuild the conserved
# state q.qn on the reloaded, already-adapted mesh. Because user_uout! writes
# the conserved variables straight out (H, Hu, Hv), the reconstruction is a
# direct copy. qe (lake at rest) is regenerated from geometry in initialize().
#
function user_read_vtu_point_data!(q, vars, ip_map, mesh)
    H_arr  = vars["H"]
    Hu_arr = vars["Hu"]
    Hv_arr = vars["Hv"]

    for ip = 1:mesh.npoin
        j = ip_map[ip]
        q.qn[ip, 1] = H_arr[j]
        q.qn[ip, 2] = Hu_arr[j]
        q.qn[ip, 3] = Hv_arr[j]
    end
end
