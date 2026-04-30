# Using JexpressoCoupling Module

This guide explains how to use the JexpressoCoupling module to couple Jexpresso with external parallel codes via MPI.

## Overview

The JexpressoCoupling module enables:
- **Multi-code coupling**: Run multiple parallel solvers together with synchronized communication
- **Flexible communication**: Synchronous and asynchronous data exchange between codes
- **Efficient boundaries**: Send/receive boundary conditions, field data, or interface information
- **Collective operations**: Broadcast, reduce, scatter, gather within and across codes

## Quick Start

### 1. Basic Usage with Command Line Options

Jexpresso now supports coupling mode via command-line flags:

```bash
# Run Jexpresso in coupled mode (normal run)
mpirun -np 4 julia --project=. run_jexpresso.jl CompEuler theta

# Run Jexpresso in coupled mode
mpirun -np 4 julia --project=. run_jexpresso.jl \
    --coupling \
    --code-id 1 \
    --n-codes 2 \
    --code-name "Jexpresso" \
    CompEuler theta
```

**Command-line options:**
- `--coupling`: Enable coupling mode (required)
- `--code-id`: Integer ID for this code (1, 2, 3, ...) [default: 1]
- `--n-codes`: Total number of coupled codes [default: 1]
- `--code-name`: Descriptive name for this code [default: "Jexpresso"]

### 2. Launching Multiple Coupled Codes

**Option A: MPMD style (Multiple Program Multiple Data)**
```bash
# Launch Jexpresso and external code together with MPI MPMD mode
mpirun -np 4 julia --project=. run_jexpresso.jl \
    --coupling --code-id 1 --n-codes 2 --code-name "Jexpresso" \
    CompEuler theta : \
    -np 4 ./external_code --code-id 2 --n-codes 2
```

**Option B: Single mpirun with rank splitting**

Create a launcher script (`launch_coupled.jl`):
```julia
#!/usr/bin/env julia
using MPI

MPI.Init()
world_rank = MPI.Comm_rank(MPI.COMM_WORLD)
world_size = MPI.Comm_size(MPI.COMM_WORLD)

# Split ranks: first half = code 1, second half = code 2
jexpresso_ranks = world_size ÷ 2
code_id = (world_rank < jexpresso_ranks) ? 1 : 2

if code_id == 1
    # Launch Jexpresso
    push!(empty!(ARGS), "--coupling", "--code-id", "1", "--n-codes", "2",
          "--code-name", "Jexpresso", "CompEuler", "theta")
    using Jexpresso
else
    # Launch external code
    push!(empty!(ARGS), "--code-id", "2", "--n-codes", "2")
    include("external_code_main.jl")
end

MPI.Finalize()
```

Then launch:
```bash
mpirun -np 8 julia --project=. launch_coupled.jl
```

### 3. Accessing the Coupling Context

When coupling is enabled, a global `coupling_ctx` variable is available:

```julia
# In your user code (e.g., in initialize.jl or boundary conditions)
if !isnothing(Jexpresso.coupling_ctx)
    ctx = Jexpresso.coupling_ctx

    println("This rank belongs to code $(ctx.code_id) named $(ctx.code_name)")
    println("Local rank: $(ctx.local_rank) of $(ctx.local_size)")
    println("World rank: $(ctx.world_rank) of $(ctx.world_size)")
    println("Am I coupling root? $(ctx.is_root)")
end
```

## Common Usage Patterns

### Pattern 1: Exchange Boundary Data at Each Timestep

```julia
using Jexpresso.JexpressoCoupling

function apply_coupled_boundary!(mesh, state, t)
    ctx = Jexpresso.coupling_ctx

    if isnothing(ctx)
        return  # Not in coupling mode
    end

    # Extract interface data
    n_interface = length(mesh.boundary_nodes)
    send_data = extract_velocity_at_boundary(state)

    # Exchange with partner code (e.g., code 2)
    if is_coupling_root(ctx)
        recv_data = zeros(n_interface)
        exchange_field_data!(ctx, send_data, recv_data, 2; tag=0)
    else
        recv_data = zeros(n_interface)
    end

    # Broadcast received data to all local ranks
    broadcast_to_local!(ctx, recv_data)

    # Apply as boundary condition
    apply_pressure_bc!(mesh, state, recv_data)
end
```

### Pattern 2: Asynchronous Exchange for Better Performance

```julia
# During initialization
function setup_coupling_buffers!(ctx)
    if is_coupling_root(ctx)
        n_interface = 100
        partner_code = 2
        register_exchange_buffer!(ctx, :velocity, n_interface, partner_code; tag=10)
        register_exchange_buffer!(ctx, :pressure, n_interface, partner_code; tag=20)
    end
end

# During timestepping
function timestep_with_async_coupling!(mesh, state, dt)
    ctx = Jexpresso.coupling_ctx

    if !isnothing(ctx) && is_coupling_root(ctx)
        # Start async send/recv
        velocity_data = extract_velocity(state)
        start_buffered_exchange!(ctx, :velocity, velocity_data)

        # Do local computation while communication happens
        compute_interior_nodes!(mesh, state, dt)

        # Wait for exchange to complete
        received_pressure = finish_buffered_exchange!(ctx, :velocity)
    else
        received_pressure = zeros(n_interface)
    end

    # Broadcast to all ranks
    if !isnothing(ctx)
        broadcast_to_local!(ctx, received_pressure)
    end

    # Apply boundary conditions
    apply_boundary_conditions!(mesh, state, received_pressure)
end
```

### Pattern 3: Synchronization Between Codes

```julia
# Global barrier across all codes
synchronize_coupling(ctx)

# Barrier only within Jexpresso (not with external codes)
synchronize_local(ctx)
```

### Pattern 4: Collective Operations Across Codes

```julia
# Broadcast from one code to all others (coupling roots only)
if is_coupling_root(ctx)
    data = [1.0, 2.0, 3.0]
    broadcast_from_code!(ctx, data, source_code_id=1)  # Code 1 broadcasts to all
end

# Allreduce across all codes
if is_coupling_root(ctx)
    local_sum = [compute_local_sum()]
    allreduce_across_codes!(ctx, local_sum, MPI.SUM)  # Sum across all codes
end
```

## Architecture

After coupling initialization, three communicators are created:

```
Before: MPI.COMM_WORLD = [Rank0, Rank1, Rank2, Rank3, Rank4, Rank5]

After:
  Code 1 (Jexpresso):        Code 2 (External):
    comm_local: [R0, R1, R2]   comm_local: [R3, R4, R5]
    root: R0                    root: R3

    comm_inter: [R0, R3]  ← Only roots communicate between codes
```

- **`comm_world`**: Original communicator with all ranks
- **`comm_local`**: Communicator for ranks within this code only
- **`comm_inter`**: Communicator connecting coupling roots (rank 0 of each code)

Jexpresso automatically uses `comm_local` for all internal operations when coupling is enabled.

## API Reference

See `JexpressoCoupling.jl` for complete API documentation. Key functions:

**Initialization:**
- `initialize_coupling(comm, code_id, n_codes; code_name)` - Set up coupling
- `finalize_coupling(ctx)` - Clean up (called automatically)
- `is_coupling_root(ctx)` - Check if this rank handles inter-code communication

**Synchronous Exchange:**
- `exchange_field_data!(ctx, send, recv, partner_id; tag)` - Blocking send/recv
- `send_field_to_partner!(ctx, data, partner_id; tag)` - Blocking send
- `recv_field_from_partner!(ctx, data, partner_id; tag)` - Blocking receive

**Asynchronous Exchange:**
- `register_exchange_buffer!(ctx, name, size, partner_id; tag)` - Pre-allocate buffer
- `start_buffered_exchange!(ctx, name, send_data)` - Begin async exchange
- `finish_buffered_exchange!(ctx, name)` - Complete and retrieve data

**Local Distribution:**
- `broadcast_to_local!(ctx, data)` - Share data with all local ranks
- `scatter_to_local!(ctx, global_data, local_data)` - Distribute to local ranks
- `gather_from_local!(ctx, local_data, global_data)` - Collect from local ranks

**Synchronization:**
- `synchronize_coupling(ctx)` - Barrier across all codes
- `synchronize_local(ctx)` - Barrier within this code only

## Troubleshooting

**Problem: "Root ranks not available on this rank"**
- Solution: This function is only available on coupling roots or after broadcasting

**Problem: Deadlock during exchange**
- Solution: Ensure both codes call exchange functions in the same order
- Use `exchange_field_data!` instead of separate send/recv for safety

**Problem: Data size mismatch**
- Solution: Ensure send and receive buffers are the same size on both sides
- Check with explicit size synchronization if needed

**Problem: Wrong number of codes detected**
- Solution: Verify all codes are launched with the same `--n-codes` value
- Check MPI configuration allows inter-job communication

## Example

See `src/mpi/example_coupling.jl` for a complete working example.
