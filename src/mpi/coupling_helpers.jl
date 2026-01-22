"""
    coupling_helpers.jl

Helper functions for common coupling operations with external codes.
These functions simplify sending/receiving data between Jexpresso and
external codes that may use different data types and conventions.
"""

using MPI

"""
    send_app_name_to_partner!(ctx::CouplingContext, app_name::String, partner_code_id::Int;
                              max_length::Int=128, tag::Int=100)

Send application name (as character string) to partner code.
This is useful when the external code does an MPI_Gather to collect names from all codes.

# Arguments
- `ctx`: Coupling context
- `app_name`: Application name string (e.g., "Jexpresso")
- `partner_code_id`: ID of the partner code to send to
- `max_length`: Maximum string length (matches external code's buffer size)
- `tag`: MPI tag for the message

# Example
```julia
if !isnothing(Jexpresso.coupling_ctx)
    ctx = Jexpresso.coupling_ctx
    send_app_name_to_partner!(ctx, "Jexpresso", 2; max_length=128)
end
```
"""
function send_app_name_to_partner!(ctx, app_name::String, partner_code_id::Int;
                                   max_length::Int=128, tag::Int=100)
    if !JexpressoCoupling.is_coupling_root(ctx)
        return
    end

    # Pad or truncate to max_length
    padded_name = rpad(app_name, max_length)[1:max_length]

    # Convert to byte array (MPI_CHARACTER equivalent)
    name_bytes = Vector{UInt8}(padded_name)

    partner_rank = JexpressoCoupling.get_partner_root_rank(ctx, partner_code_id)

    MPI.Send(name_bytes, ctx.comm_world; dest=partner_rank, tag=tag)

    if ctx.is_root
        @info "[$(ctx.code_name)] Sent app name '$app_name' to code $partner_code_id (rank $partner_rank)"
    end
end

"""
    recv_app_name_from_partner!(ctx::CouplingContext, partner_code_id::Int;
                                max_length::Int=128, tag::Int=100) -> String

Receive application name from partner code.

# Returns
- String containing the partner's application name (trimmed of trailing spaces)
"""
function recv_app_name_from_partner!(ctx, partner_code_id::Int;
                                     max_length::Int=128, tag::Int=100)
    if !JexpressoCoupling.is_coupling_root(ctx)
        return ""
    end

    name_bytes = Vector{UInt8}(undef, max_length)
    partner_rank = JexpressoCoupling.get_partner_root_rank(ctx, partner_code_id)

    MPI.Recv!(name_bytes, ctx.comm_world; source=partner_rank, tag=tag)

    # Convert bytes to string and trim
    partner_name = String(name_bytes) |> rstrip

    if ctx.is_root
        @info "[$(ctx.code_name)] Received app name '$partner_name' from code $partner_code_id"
    end

    return partner_name
end

"""
    exchange_app_names!(ctx::CouplingContext, my_name::String, partner_code_id::Int;
                        max_length::Int=128) -> String

Exchange application names with partner code (mutual send/receive).

# Returns
- Partner's application name

# Example
```julia
partner_name = exchange_app_names!(ctx, "Jexpresso", 2)
println("Coupled with: ", partner_name)
```
"""
function exchange_app_names!(ctx, my_name::String, partner_code_id::Int;
                            max_length::Int=128)
    if !JexpressoCoupling.is_coupling_root(ctx)
        return ""
    end

    # Pad or truncate to max_length
    my_name_padded = rpad(my_name, max_length)[1:max_length]
    my_bytes = Vector{UInt8}(my_name_padded)

    partner_bytes = Vector{UInt8}(undef, max_length)
    partner_rank = JexpressoCoupling.get_partner_root_rank(ctx, partner_code_id)

    # Use Sendrecv to avoid deadlock
    MPI.Sendrecv!(my_bytes, partner_bytes, ctx.comm_world;
                  dest=partner_rank, sendtag=100,
                  source=partner_rank, recvtag=100)

    partner_name = String(partner_bytes) |> rstrip

    if ctx.is_root
        @info "[$(ctx.code_name)] Exchanged names with code $partner_code_id: '$partner_name'"
    end

    return partner_name
end

"""
    send_field_array!(ctx::CouplingContext, field_data::AbstractArray{T},
                      partner_code_id::Int; tag::Int=200) where T

Send field data array to partner code. Only coupling root participates.

# Arguments
- `ctx`: Coupling context
- `field_data`: Array of field data (e.g., velocity, pressure)
- `partner_code_id`: Partner code ID
- `tag`: MPI tag

# Example
```julia
# Send velocity data at interface
interface_velocity = extract_boundary_velocity(mesh, state)
send_field_array!(ctx, interface_velocity, 2; tag=200)
```
"""
function send_field_array!(ctx, field_data::AbstractArray{T},
                          partner_code_id::Int; tag::Int=200) where T
    if !JexpressoCoupling.is_coupling_root(ctx)
        return
    end

    partner_rank = JexpressoCoupling.get_partner_root_rank(ctx, partner_code_id)
    MPI.Send(field_data, ctx.comm_world; dest=partner_rank, tag=tag)

    if ctx.is_root
        @debug "[$(ctx.code_name)] Sent field array (size=$(length(field_data))) to code $partner_code_id"
    end
end

"""
    recv_field_array!(ctx::CouplingContext, field_data::AbstractArray{T},
                      partner_code_id::Int; tag::Int=200) where T

Receive field data array from partner code. Only coupling root participates.

# Arguments
- `ctx`: Coupling context
- `field_data`: Preallocated buffer to receive into
- `partner_code_id`: Partner code ID
- `tag`: MPI tag
"""
function recv_field_array!(ctx, field_data::AbstractArray{T},
                          partner_code_id::Int; tag::Int=200) where T
    if !JexpressoCoupling.is_coupling_root(ctx)
        return
    end

    partner_rank = JexpressoCoupling.get_partner_root_rank(ctx, partner_code_id)
    MPI.Recv!(field_data, ctx.comm_world; source=partner_rank, tag=tag)

    if ctx.is_root
        @debug "[$(ctx.code_name)] Received field array (size=$(length(field_data))) from code $partner_code_id"
    end
end

"""
    gather_to_coupling_root!(ctx::CouplingContext, local_data::AbstractArray{T},
                             global_data::AbstractArray{T}) where T

Gather data from all local ranks to the coupling root (local rank 0).
This prepares data for sending to partner code.

# Example
```julia
# Each rank has local boundary data
local_boundary = get_local_boundary_data()

# Root gathers all
if JexpressoCoupling.is_coupling_root(ctx)
    global_boundary = zeros(total_size)
else
    global_boundary = zeros(0)  # Non-roots don't need this
end

gather_to_coupling_root!(ctx, local_boundary, global_boundary)

# Now root can send global_boundary to partner
if JexpressoCoupling.is_coupling_root(ctx)
    send_field_array!(ctx, global_boundary, 2)
end
```
"""
function gather_to_coupling_root!(ctx, local_data::AbstractArray{T},
                                  global_data::AbstractArray{T}) where T
    JexpressoCoupling.gather_from_local!(ctx, local_data, global_data; root=0)
end

"""
    scatter_from_coupling_root!(ctx::CouplingContext, global_data::AbstractArray{T},
                                local_data::AbstractArray{T}) where T

Scatter data from coupling root to all local ranks.
This distributes received data from partner code to all Jexpresso ranks.

# Example
```julia
# Root receives from partner
if JexpressoCoupling.is_coupling_root(ctx)
    global_data = zeros(total_size)
    recv_field_array!(ctx, global_data, 2)
else
    global_data = zeros(0)
end

# Scatter to all ranks
local_data = zeros(local_size)
scatter_from_coupling_root!(ctx, global_data, local_data)

# Now all ranks have their portion
apply_boundary_condition!(local_data)
```
"""
function scatter_from_coupling_root!(ctx, global_data::AbstractArray{T},
                                     local_data::AbstractArray{T}) where T
    JexpressoCoupling.scatter_to_local!(ctx, global_data, local_data; root=0)
end

"""
    coupling_exchange_workflow!(ctx::CouplingContext,
                                send_data::AbstractArray{T},
                                recv_data::AbstractArray{T},
                                partner_code_id::Int;
                                tag_send::Int=200,
                                tag_recv::Int=201) where T

Complete workflow for exchanging field data between coupling roots.
Uses non-blocking communication for better performance.

# Example
```julia
# At each time step
send_velocity = extract_boundary_velocity()
recv_pressure = zeros(size(send_velocity))

coupling_exchange_workflow!(ctx, send_velocity, recv_pressure, 2)

# recv_pressure now contains data from partner
apply_pressure_bc!(recv_pressure)
```
"""
function coupling_exchange_workflow!(ctx,
                                     send_data::AbstractArray{T},
                                     recv_data::AbstractArray{T},
                                     partner_code_id::Int;
                                     tag_send::Int=200,
                                     tag_recv::Int=201) where T
    if !JexpressoCoupling.is_coupling_root(ctx)
        return
    end

    partner_rank = JexpressoCoupling.get_partner_root_rank(ctx, partner_code_id)

    # Start non-blocking send and receive
    send_req = MPI.Isend(send_data, ctx.comm_world; dest=partner_rank, tag=tag_send)
    recv_req = MPI.Irecv!(recv_data, ctx.comm_world; source=partner_rank, tag=tag_recv)

    # Wait for both to complete
    MPI.Wait(send_req)
    MPI.Wait(recv_req)
end

"""
    broadcast_from_root_to_all_local!(ctx::CouplingContext, data::AbstractArray)

Broadcast data from coupling root to all local ranks.
Shorthand for common operation after receiving from partner.

# Example
```julia
# Root receives
if JexpressoCoupling.is_coupling_root(ctx)
    data = zeros(100)
    recv_field_array!(ctx, data, 2)
else
    data = zeros(100)
end

# Share with all ranks
broadcast_from_root_to_all_local!(ctx, data)
```
"""
function broadcast_from_root_to_all_local!(ctx, data::AbstractArray)
    JexpressoCoupling.broadcast_to_local!(ctx, data; root=0)
end
