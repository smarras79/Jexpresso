# Save conservative fields independently of the case's plotted variables.
# DynSGS history is deliberately reinitialized by params_setup on restart.
function write_shallow_water_vtk_restart_data!(vtk, state, reference, mesh, inputs)
    get(inputs, :parsed_equations, "") == "ShallowWater" || return nothing
    size(reference, 2) >= 3 && length(state) == 3mesh.npoin || return nothing
    for (ivar, name) in enumerate(("H", "Hu", "Hv"))
        offset = (ivar - 1) * mesh.npoin
        vtk["restart_" * name, VTKPointData()] = @view state[offset+1:offset+mesh.npoin]
        vtk["restart_reference_" * name, VTKPointData()] = @view reference[:, ivar]
    end
    vtk["restart_nparts", VTKFieldData()] = Float64(mesh.nparts)
    return nothing
end

"""
Restore a 2D shallow-water snapshot from the latest PVD entry or the explicitly
selected :restart_vtk_iout. New snapshots carry conservative and reference
fields; older ReducedGravity snapshots use their H/Hu/Hv arrays and the normal
lake-at-rest reference from initialize. DynSGS history is never loaded.
"""
function read_shallow_water_vtk_restart!(q, mesh, inputs; output_dir="")
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    directory = get(inputs, :restart_vtk_input_dir, output_dir)
    iout = get(inputs, :restart_vtk_iout, -1)
    pvd_time = NaN
    if iout < 0
        last = rank == 0 ? read_pvd_last_entry(joinpath(directory, "simulation.pvd")) : (NaN, -1)
        pvd_time, iout = MPI.bcast(last, 0, comm)
    end
    nparts = MPI.Comm_size(comm)
    part = lpad(rank + 1, ndigits(nparts), '0')
    path = joinpath(directory, "iter_$(iout)", "iter_$(iout)_$(part).vtu")
    # Older files have no rank metadata; their PVTU piece count still records
    # the partitioning, including when local point counts happen to match.
    pvtu = joinpath(directory, "iter_$(iout).pvtu")
    pieces = length(collect(eachmatch(r"<Piece\s+Source=", read(pvtu, String))))
    pieces == nparts || error("Shallow-water VTK restart requires the same MPI rank count (saved $pieces, current $nparts).")
    names = ["H", "Hu", "Hv", "__coords__"]
    append!(names, ["restart_" * name for name in ("H", "Hu", "Hv")])
    append!(names, ["restart_reference_" * name for name in ("H", "Hu", "Hv")])
    data, npoin = read_vtu_point_data(path, names)
    npoin == mesh.npoin || error("Shallow-water VTK restart requires the same mesh and polynomial order (saved $npoin points, current $(mesh.npoin)).")
    metadata, _ = read_vtu_point_data(path, ["TimeValue", "restart_nparts"]; field_data=true)
    if haskey(metadata, "restart_nparts")
        only(metadata["restart_nparts"]) == nparts || error("Shallow-water VTK restart rank count differs.")
    end
    time = haskey(metadata, "TimeValue") ? only(metadata["TimeValue"]) : pvd_time
    isfinite(time) || error("VTK snapshot has no TimeValue; select it through simulation.pvd to recover the restart time.")
    ipmap = vtk_to_mesh_ipmap(data["__coords__"], view(mesh.coords, 1, :),
                             view(mesh.coords, 2, :), zeros(mesh.npoin))
    conservative = all(haskey(data, "restart_" * name) for name in ("H", "Hu", "Hv"))
    references = all(haskey(data, "restart_reference_" * name) for name in ("H", "Hu", "Hv"))
    for (ivar, name) in enumerate(("H", "Hu", "Hv"))
        key = conservative ? "restart_" * name : name
        haskey(data, key) || error("Shallow-water VTK snapshot is missing conservative field $name.")
        copyto!(view(q.qn, :, ivar), data[key][ipmap])
        if references
            copyto!(view(q.qe, :, ivar), data["restart_reference_" * name][ipmap])
        end
    end
    inputs[:tinit] = time
    inputs[:restart_vtk_iout] = iout
    rank == 0 && println(" # Shallow-water VTK restart: iter_$iout at t=$time (DynSGS history reinitialized)")
    return nothing
end
