
function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

            """
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 3D CompEuler with θ equation ........................ "
    end
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ["ρ", "ρu", "ρv", "ρw", "ρθ"]
    qoutvars = ["ρ", "u", "v", "w", "θ"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------
    
    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()
        if inputs[:lrestart] == true
            #
            # READ RESTART HDF5:
            #
           q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars))

            # pvtu_filepath = inputs[:restart_input_file_path]

            # if rank == 0; println("\n--- Starting READ procedure ---"); end
            # read_pvtu_restart!(q, pvtu_filepath; comm=comm, verbose=false)
            # if rank == 0; println("Read successful."); end
            
            
            for ip=1:mesh.npoin
                ρ  = q.qn[ip,1]
                ρθ = q.qn[ip,5]
                θ  = ρθ/ρ
                P = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
                q.qn[ip,end] = P
            
                ρe  = q.qe[ip,1]
                ρθe = q.qe[ip,5]
                θe  = ρθe/ρ
                Pe = perfectGasLaw_ρθtoP(PhysConst, ρ=ρe, θ=θe)
                q.qe[ip,end] = Pe
            end
            
        else
            #
            # INITIAL STATE from scratch:
            #
            comm = MPI.COMM_WORLD
            max_x = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
            min_x = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)
            xc = (max_x + min_x)/2
            zc = 2500.0 #m
            r0 = 2000.0 #m
        
            θref = 300.0 #K
            θc   =   2.0 #K
            for ip = 1:mesh.npoin
            
                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]
            
                r = sqrt( (x - xc)^2 + (z - zc)^2 )
            
                Δθ = 0.0 #K
                if r < r0
                    Δθ = θc*(1.0 - r/r0)
                end
                θ = θref + Δθ
                p    = PhysConst.pref*(1.0 - PhysConst.g*z/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
                pref = PhysConst.pref*(1.0 - PhysConst.g*z/(PhysConst.cp*θref))^(PhysConst.cpoverR)
                ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)    #kg/m³
                ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

                u = 0.0
                v = 0.0
                w = 0.0
            
                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*w - ρref*w
                    q.qn[ip,5] = ρ*θ - ρref*θref
                    q.qn[ip,end] = p
                
                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = u
                    q.qe[ip,3] = v
                    q.qe[ip,4] = w
                    q.qe[ip,5] = ρref*θref
                    q.qe[ip,end] = pref
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*w
                    q.qn[ip,5] = ρ*θ
                    q.qn[ip,end] = p

                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = ρref*u
                    q.qe[ip,3] = ρref*v
                    q.qe[ip,4] = ρref*w
                    q.qe[ip,5] = ρref*θref
                    q.qe[ip,end] = pref
                end
                #end
            end
        end
    
        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[:,2] .= q.qn[:,2]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,3] .= q.qn[:,3]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,4] .= q.qn[:,4]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,5] .= q.qn[:,5]./(q.qn[:,1] + q.qe[:,1])
            
                #Store initial background state for plotting and analysis of pertuebations
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]
                q.qn[:,5] .= q.qn[:,5]./q.qn[:,1]

                #Store initial background state for plotting and analysis of pertuebations
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            end
        end

    else
        if (inputs[:SOL_VARS_TYPE] == PERT())
            lpert = true
        else
            lpert = false
        end
        PhysConst = PhysicalConst{TFloat}()
        xc = TFloat((maximum(mesh.x) + minimum(mesh.x))/2)
        zc = TFloat(2500.0) #m
        rθ = TFloat(2000.0) #m

        θref = TFloat(300.0) #K
        θc   =   TFloat(2.0) #K
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, mesh.z, xc, rθ, zc, θref, θc, PhysConst, lpert; ndrange = (mesh.npoin))
    end
    if rank == 0
        @info " Initialize fields for 3D CompEuler with θ equation ........................ DONE "
    end
    
    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, z, xc, rθ, zc, θref, θc, PhysConst, lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    x = x[ip]
    y = y[ip]
    z = z[ip]
    r = sqrt( (x - xc)^2 + (z - zc)^2 )
    Δθ = T(0.0) #K
    
    if r < rθ
        Δθ = T(θc*(T(1.0) - r/rθ))
    end
    
    θ = θref + Δθ
    p    = PhysConst.pref*(T(1.0) - PhysConst.g*z/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
    pref = PhysConst.pref*(T(1.0) - PhysConst.g*z/(PhysConst.cp*θref))^(PhysConst.cpoverR)
    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)    #kg/m³
    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

    u = T(0.0)
    v = T(0.0)
    w = T(0.0)

    if (lpert)
        qn[ip,1] = ρ - ρref
        qn[ip,2] = ρ*u - ρref*u
        qn[ip,3] = ρ*v - ρref*v
        qn[ip,4] = ρ*w - ρref*w
        qn[ip,5] = ρ*θ - ρref*θref
        qn[ip,end] = p
    else
        qn[ip,1] = ρ
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*w
        qn[ip,5] = ρ*θ
        qn[ip,end] = p
    end

                    #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = ρref
    qe[ip,2] = ρref*u
    qe[ip,3] = ρref*v
    qe[ip,4] = ρref*w
    qe[ip,5] = ρref*θref
    qe[ip,end] = pref

end

function user_get_adapt_flags(inputs, old_ad_lvl, q, qe, connijk, nelem, ngl)
    adapt_flags = KernelAbstractions.zeros(CPU(), TInt, Int64(nelem))
    ips         = KernelAbstractions.zeros(CPU(), TInt, ngl * ngl * ngl)
    tol         = 1.0
    max_level   = inputs[:amr_max_level] 
    
    for iel = 1:nelem
        m = 1
        for i = 1:ngl
            for j = 1:ngl
                for k = 1:ngl
                    ips[m] = connijk[iel, i, j, k]
                    m += 1
                end
            end
        end
        # @info q[ips,4] - qe[ips,4]
        theta      = q[ips, 5] ./ q[ips, 1]
        theta_ref  = qe[ips, 5] ./ qe[ips, 1]
        dtheta     = theta - theta_ref
        # @info dtheta
        if any(dtheta .> tol) && (old_ad_lvl[iel] < max_level)
            adapt_flags[iel] = refine_flag
        end
        if all(dtheta .< tol)
            adapt_flags[iel] = coarsen_flag
        end
    end
    return adapt_flags
end


"""
    read_pvtu_restart!(q::SolutionData, pvtu_filepath::String; comm=MPI.COMM_WORLD, verbose=true)

Reads a parallel VTK file collection (`.pvtu` and its associated `.vtu` files)
and populates the `q.qn` array for a simulation restart.

This function is MPI-aware and must be called by all processes in the communicator.
Rank 0 reads the `.pvtu` file to get the list of partition files and broadcasts
this list. Each rank then reads its corresponding `.vtu` file in parallel.

# Arguments
- `q::SolutionData`: A mutable struct holding the solution arrays. The `q.qn`
  field, a matrix of size `(npoin_local, nvars)`, will be populated with data.
  It is assumed to be pre-allocated with the correct dimensions for the local partition.
- `pvtu_filepath::String`: The path to the main `.pvtu` file.
- `comm::MPI.Comm`: The MPI communicator (defaults to `MPI.COMM_WORLD`).
- `verbose::Bool`: If true, rank 0 will print status messages.
"""
# run_restart.jl
# using MPI
# using ReadVTK
# using XMLDict
# using Printf
# using WriteVTK # Provides MeshCell and VTK writing functions

# # --- Solver-Specific Data Structure ---
# # This is an example of what your solver's data structure might look like.
# # The I/O functions below are generic and will work as long as the
# # object `q` has a `qn` matrix.
# mutable struct SolverState
#     qn::Matrix{Float64}
#     # Your structure could have other fields, which the I/O functions will ignore.
# end

# # --- Generic I/O Functions ---

# """
#     read_pvtu_restart!(q, pvtu_filepath::String; ...)

# Reads parallel VTK data (.pvtu) to populate the `q.qn` array for a restart.
# This function is generic and works on any object `q` with a `q.qn` matrix.
# """
# function read_pvtu_restart!(q, pvtu_filepath::String; comm=MPI.COMM_WORLD, verbose=true)
#     rank = MPI.Comm_rank(comm)
#     nranks = MPI.Comm_size(comm)
    
#     local_vtu_filename = ""

#     # This assumes a fixed filename like "iter_11.pvtu" needs to be appended.
#     # If the full path is already in `pvtu_filepath`, you can remove the joinpath().
#     full_pvtu_path = joinpath(pvtu_filepath, "iter_11.pvtu")
    
#     if rank == 0
#         if !isfile(full_pvtu_path)
#             error("PVTU file not found: $full_pvtu_path")
#         end
#         if verbose
#             println("Rank 0: Reading PVTU manifest from '$full_pvtu_path'...")
#         end

#         xml_string = read(full_pvtu_path, String)
        
#         # --- START OF CORRECTION ---
#         # Add a check to ensure the file is not empty before parsing.
#         if isempty(xml_string)
#             error("The PVTU file at '$full_pvtu_path' was found, but it is empty. Please check the file content.")
#         end
#         # --- END OF CORRECTION ---
        
#         xml_data = xml_dict(xml_string)
        
#         pieces = xml_data["VTKFile"]["PUnstructuredGrid"]["Piece"]
#         if !(pieces isa AbstractVector); pieces = [pieces]; end
#         if length(pieces) != nranks; error("PVTU has $(length(pieces)) pieces, but MPI size is $nranks."); end
        
#         vtu_files = Vector{String}(undef, nranks)
#         for i in 1:nranks
#             local piece = pieces[i]
#             if !(piece isa AbstractDict)
#                 error("Piece $i is not a valid dictionary structure: $piece")
#             end

#             if haskey(piece, :Source)
#                 vtu_files[i] = piece[:Source]
#             elseif haskey(piece, "Source")
#                 vtu_files[i] = piece["Source"]
#             else
#                 error("Could not find Source attribute in piece $i: $piece")
#             end
#         end

#         local_vtu_filename = vtu_files[1]
#         for dest_rank in 1:(nranks-1)
#             filename_bytes = Vector{UInt8}(vtu_files[dest_rank + 1])
#             MPI.Send(filename_bytes, dest_rank, 0, comm)
#         end
#     else
#         status = MPI.Probe(0, 0, comm)
#         count = MPI.Get_count(status, UInt8)
#         filename_bytes = Vector{UInt8}(undef, count)
#         MPI.Recv!(filename_bytes, 0, 0, comm)
#         local_vtu_filename = String(filename_bytes)
#     end

#     MPI.Barrier(comm)

#     pvtu_dir = dirname(full_pvtu_path)
#     local_vtu_filepath = joinpath(pvtu_dir, local_vtu_filename)
#     if !isfile(local_vtu_filepath); error("Rank $rank: Cannot find data file: $local_vtu_filepath"); end
    
#     vtk = VTKFile(local_vtu_filepath)
#     point_data = get_point_data(vtk)
#     var_names = collect(keys(point_data))
    
#     npoin_local, nvars_alloc = size(q.qn)
#     npoin_file = vtk.n_points
#     nvars_file = length(var_names)
    
#     if npoin_local != npoin_file; error("Rank $rank: Point count mismatch. Allocated $npoin_local, file has $npoin_file."); end
    
#     nvars = min(nvars_alloc, nvars_file)
#     for ivar in 1:nvars
#         # Use ReadVTK.get_data to avoid ambiguity with other functions
#         q.qn[:, ivar] .= ReadVTK.get_data(point_data[var_names[ivar]])
#     end
    
#     if verbose && rank == 0
#         println("All ranks successfully loaded restart data.")
#     end
    
#     return nothing
# end

# """
#     write_pvtu_data(q, x, y, z, var_names, output_basename; ...)

# Writes the distributed solution data to a new parallel VTK collection.
# This function is generic and takes coordinates as separate arguments.
# """
# function write_pvtu_data(q, x::Vector, y::Vector, z::Vector, var_names::Vector{String}, output_basename::String; comm=MPI.COMM_WORLD)
#     rank = MPI.Comm_rank(comm)
#     nranks = MPI.Comm_size(comm)
#     if size(q.qn, 2) != length(var_names); error("Variable name count mismatch."); end

#     output_dir = output_basename
#     if rank == 0; mkpath(output_dir); end
#     MPI.Barrier(comm)

#     coords = hcat(x, y, z)'
#     vtu_filename = joinpath(output_dir, "part_$(rank).vtu")

#     vtk_grid(vtu_filename, coords, MeshCell[]) do vtk
#         for (ivar, name) in enumerate(var_names)
#             vtk[name] = q.qn[:, ivar]
#         end
#     end

#     MPI.Barrier(comm)

#     if rank == 0
#         pvtu_filename = joinpath(output_dir, output_basename * ".pvtu")
#         pvtk_grid(pvtu_filename, coords, MeshCell[]; part=rank+1, nparts=nranks) do pvtk
#             for (ivar, name) in enumerate(var_names)
#                 pvtk[name] = q.qn[:, ivar]
#             end
#         end
#         println("Rank 0: Successfully wrote parallel data to '$pvtu_filename'")
#     end

#     return nothing
# end
