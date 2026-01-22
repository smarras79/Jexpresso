# Coupling Jexpresso with External Codes

This guide shows how to couple Jexpresso with external codes that use MPI for data exchange.

## Quick Reference

### Helper Functions Available

Located in `src/mpi/coupling_helpers.jl`:

**Application Name Exchange:**
- `send_app_name_to_partner!(ctx, name, partner_id; max_length=128)` - Send name string
- `recv_app_name_from_partner!(ctx, partner_id; max_length=128)` - Receive name string
- `exchange_app_names!(ctx, my_name, partner_id)` - Bidirectional name exchange

**Field Data Exchange:**
- `send_field_array!(ctx, data, partner_id; tag=200)` - Send array to partner
- `recv_field_array!(ctx, data, partner_id; tag=200)` - Receive array from partner
- `coupling_exchange_workflow!(ctx, send_data, recv_data, partner_id)` - Bidirectional exchange

**Data Distribution (within Jexpresso):**
- `gather_to_coupling_root!(ctx, local_data, global_data)` - Collect from all Jexpresso ranks
- `scatter_from_coupling_root!(ctx, global_data, local_data)` - Distribute to all Jexpresso ranks
- `broadcast_from_root_to_all_local!(ctx, data)` - Share full array with all ranks

## Common Patterns

### Pattern 1: Send Application Name During Initialization

Your external code does: `MPI_Gather(app_name, 128, MPI_CHARACTER, ...)`

In Jexpresso (e.g., in `initialize.jl` or startup code):

```julia
# Check if coupling is enabled
if !isnothing(Jexpresso.coupling_ctx)
    ctx = Jexpresso.coupling_ctx

    # Send "Jexpresso" name to partner code
    # max_length=128 matches your external code's buffer size
    send_app_name_to_partner!(ctx, "Jexpresso", 2; max_length=128)
end
```

### Pattern 2: Send Boundary Velocity, Receive Pressure (Each Timestep)

This is the most common coupling pattern for fluid-structure or multi-physics coupling.

**In your timestepping loop or boundary condition file:**

```julia
function apply_coupled_boundary!(mesh, state, t, dt)
    ctx = Jexpresso.coupling_ctx

    if isnothing(ctx)
        return  # Not in coupling mode, skip
    end

    # 1. Extract velocity at coupling interface
    #    (This is domain-specific - adjust for your mesh structure)
    local_velocity = extract_boundary_velocity(mesh, state)

    # 2. Gather velocity from all Jexpresso ranks to coupling root
    n_total_interface = get_total_interface_size(mesh)

    if is_coupling_root(ctx)
        global_velocity = zeros(n_total_interface)
    else
        global_velocity = zeros(0)
    end

    gather_to_coupling_root!(ctx, local_velocity, global_velocity)

    # 3. Coupling root exchanges with external code
    if is_coupling_root(ctx)
        global_pressure = zeros(n_total_interface)

        # Send velocity (tag=200), receive pressure (tag=201)
        coupling_exchange_workflow!(ctx, global_velocity, global_pressure, 2;
                                    tag_send=200, tag_recv=201)
    else
        global_pressure = zeros(0)
    end

    # 4. Distribute received pressure to all Jexpresso ranks
    local_pressure = zeros(length(local_velocity))
    scatter_from_coupling_root!(ctx, global_pressure, local_pressure)

    # 5. Apply pressure as boundary condition
    apply_pressure_bc!(mesh, state, local_pressure)
end
```

### Pattern 3: Broadcast Instead of Scatter

When all Jexpresso ranks need the full array (not just their portion):

```julia
# Root receives
if is_coupling_root(ctx)
    global_data = zeros(n_data)
    recv_field_array!(ctx, global_data, 2; tag=300)
else
    global_data = zeros(n_data)
end

# Everyone gets the full array
broadcast_from_root_to_all_local!(ctx, global_data)

# Now all ranks can use global_data
```

### Pattern 4: One-Way Communication (Jexpresso → External)

If you only need to send data (e.g., for monitoring or post-processing):

```julia
# Gather data
gather_to_coupling_root!(ctx, local_data, global_data)

# Send to partner
if is_coupling_root(ctx)
    send_field_array!(ctx, global_data, 2; tag=500)
end
```

## Integration into Jexpresso

### Option 1: Modify User Boundary Conditions

Add coupling to your `user_bc.jl`:

```julia
# In your boundary condition function
function user_bc!(...)
    # ... existing BC code ...

    # Add coupled BC if coupling is active
    if !isnothing(Jexpresso.coupling_ctx)
        apply_coupled_boundary!(mesh, state, t, dt)
    end
end
```

### Option 2: Add to Timestepping Loop

Modify the main driver or timestepping routine:

```julia
for timestep in 1:nsteps
    # ... solve equations ...

    # Exchange data with coupled code
    if !isnothing(Jexpresso.coupling_ctx)
        exchange_coupling_data!(Jexpresso.coupling_ctx, mesh, state, t)
    end

    # ... continue ...
end
```

### Option 3: Initialize at Startup

For one-time exchanges (like mesh info), add to `initialize.jl`:

```julia
# After mesh initialization
if !isnothing(Jexpresso.coupling_ctx)
    ctx = Jexpresso.coupling_ctx

    # Send mesh information to partner
    if is_coupling_root(ctx)
        mesh_coords = get_interface_coordinates(mesh)
        send_field_array!(ctx, mesh_coords, 2; tag=10)
    end
end
```

## MPI Tags Convention

Use distinct tags for different data types to avoid confusion:

```julia
# Recommended tag ranges
TAG_APP_NAME = 100          # Application name exchange
TAG_MESH_INFO = 10-19       # Mesh/geometry data
TAG_VELOCITY = 200          # Velocity field
TAG_PRESSURE = 201          # Pressure field
TAG_TEMPERATURE = 202       # Temperature field
TAG_CONTROL = 300-399       # Control/status messages
```

Example:
```julia
# Send velocity with tag 200
send_field_array!(ctx, velocity, 2; tag=200)

# Receive pressure with tag 201 (different tag!)
recv_field_array!(ctx, pressure, 2; tag=201)
```

## Matching External Code Expectations

### Character String Gathering (MPI_Gather)

Your external code:
```fortran
call MPI_Gather(app_name, 128, MPI_CHARACTER, app_dumm, 128, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
```

If the external code gathers from **all codes' roots** to its rank 0:

```julia
# Jexpresso sends its name to external code's root
if is_coupling_root(ctx)
    # The external code's root will gather this
    send_app_name_to_partner!(ctx, "Jexpresso", 2; max_length=128)
end
```

### Field Data Gathering

If external code expects gathered field data:

```julia
# External code does: MPI_Gather(local_field, ..., global_field, ...)
# Jexpresso sends the full gathered array

# 1. Gather within Jexpresso
gather_to_coupling_root!(ctx, local_field, global_field)

# 2. Send to external code
if is_coupling_root(ctx)
    send_field_array!(ctx, global_field, 2)
end
```

## Synchronization

### When to Synchronize

```julia
# After name exchange (one-time)
exchange_app_names!(ctx, "Jexpresso", 2)
synchronize_coupling(ctx)

# Between timesteps (if needed)
for timestep in 1:nsteps
    solve_timestep!(...)
    exchange_coupling_data!(...)
    synchronize_coupling(ctx)  # Ensure both codes are in sync
end

# Within Jexpresso only (doesn't affect external code)
synchronize_local(ctx)  # Barrier among Jexpresso ranks only
```

## Debugging Tips

### Enable Detailed Logging

```julia
# At the top of your coupling code
ENV["JULIA_DEBUG"] = "Main"  # Enable debug prints

# Or selectively
if is_coupling_root(ctx)
    @info "About to send $(length(data)) elements to code 2"
end
```

### Add Checkpoints

```julia
if is_coupling_root(ctx)
    println("[Jexpresso] Checkpoint 1: Before send")
end
send_field_array!(ctx, data, 2)

if is_coupling_root(ctx)
    println("[Jexpresso] Checkpoint 2: After send")
end

synchronize_coupling(ctx)  # Wait for everyone

if is_coupling_root(ctx)
    println("[Jexpresso] Checkpoint 3: After sync")
end
```

### Verify Data

```julia
if is_coupling_root(ctx)
    println("Sending data:")
    println("  Size: $(length(send_data))")
    println("  Min:  $(minimum(send_data))")
    println("  Max:  $(maximum(send_data))")
    println("  Mean: $(sum(send_data)/length(send_data))")
end
```

## Complete Example

See `src/mpi/example_coupling_workflow.jl` for runnable examples of all patterns.

Run the example:
```bash
mpirun -np 4 julia --project=. -e 'include("src/mpi/example_coupling_workflow.jl"); run_coupling_workflow_demo()'
```

## Performance Considerations

### Use Non-Blocking Communication

Preferred:
```julia
coupling_exchange_workflow!(ctx, send_data, recv_data, 2)  # Non-blocking internally
```

Instead of:
```julia
send_field_array!(ctx, send_data, 2)  # Blocking
recv_field_array!(ctx, recv_data, 2)  # Blocking
```

### Minimize Gathers

If possible, have only boundary ranks participate:
```julia
# Good: Only boundary ranks send
if is_boundary_rank(rank)
    send_to_coupling_root(boundary_data)
end

# Less efficient: All ranks gather, then send
gather_to_coupling_root!(ctx, local_data, global_data)  # All ranks participate
```

### Reuse Buffers

```julia
# Allocate once
global_velocity = zeros(n_interface)
global_pressure = zeros(n_interface)

# Reuse in loop
for timestep in 1:nsteps
    gather_to_coupling_root!(ctx, local_velocity, global_velocity)
    coupling_exchange_workflow!(ctx, global_velocity, global_pressure, 2)
    scatter_from_coupling_root!(ctx, global_pressure, local_pressure)
end
```

## Troubleshooting

**Problem:** Deadlock during exchange
- **Solution:** Ensure both codes use matching send/recv or both use `coupling_exchange_workflow!`
- Check that tags match between sender and receiver

**Problem:** Wrong data received
- **Solution:** Verify buffer sizes match: `length(send_data) == length(recv_data)`
- Check data types match (Float64, Float32, etc.)

**Problem:** "Root ranks not available" error
- **Solution:** This function is only available on coupling roots
- Add check: `if is_coupling_root(ctx)`

**Problem:** Application name truncated or padded strangely
- **Solution:** Ensure `max_length` parameter matches external code's buffer size

## Next Steps

1. **Test basic coupling:** Run `test_simple_coupling.jl`
2. **Identify coupling points:** Where in your simulation do you need to exchange data?
3. **Modify boundary conditions:** Add coupling exchange to your BC routines
4. **Verify data flow:** Use debug prints to confirm data is correct
5. **Performance tuning:** Profile and optimize communication patterns

For more information, see:
- `src/mpi/COUPLING_USAGE.md` - General coupling documentation
- `src/mpi/example_coupling_workflow.jl` - Runnable examples
- `COUPLING_TESTS.md` - Testing guide
