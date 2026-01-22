# Known Issues and Limitations with Coupling

## Current Status

✅ **Working:**
- Coupling initialization and communicator splitting
- Data exchange between codes (send/recv)
- Mesh partitioning with `with_mpi(comm)`
- Rank management and identification

⚠️ **Potential Issues:**

### 1. I/O Functions Hardcode MPI.COMM_WORLD

**Location:** `src/io/write_output.jl` (multiple functions)

**Issue:**
Several output functions hardcode `comm = MPI.COMM_WORLD` instead of using the coupling-aware communicator:

```julia
function write_output(...)
    comm = MPI.COMM_WORLD  # ← Should use local comm when coupling
    rank = MPI.Comm_rank(comm)
    ...
end
```

**Impact:**
- Output files may have conflicts if both coupled codes write simultaneously
- Rank numbering in output files may be incorrect (world ranks vs local ranks)
- Parallel I/O operations may fail or produce unexpected results

**Workaround:**
When coupling, each code should write to separate output directories:
```julia
# In user_inputs.jl or initialization
if !isnothing(Jexpresso.coupling_ctx)
    ctx = Jexpresso.coupling_ctx
    OUTPUT_DIR = "output_$(ctx.code_name)"  # e.g., "output_Jexpresso-A"
else
    OUTPUT_DIR = "output"
end
```

**Permanent Fix (TODO):**
Modify I/O functions to accept `comm` parameter and use the coupling-aware communicator.

### 2. Timer Functions Use MPI.COMM_WORLD

**Location:** `src/auxiliary/timing.jl`, `src/macros/je_macros.jl`

**Issue:**
Timing functions default to `MPI.COMM_WORLD`:
```julia
function MPIFunctionTimer(comm::MPI.Comm=MPI.COMM_WORLD; skip_first_n::Int=1)
```

**Impact:**
- Timing statistics may be incorrect (averaged over all codes instead of just Jexpresso)
- Performance profiles may be confusing

**Workaround:**
Timing should still work but will show aggregated statistics. For accurate timing of just Jexpresso:
```julia
# Explicitly pass local comm if available
if !isnothing(Jexpresso.coupling_ctx)
    timer = MPIFunctionTimer(Jexpresso.coupling_ctx.comm_local)
end
```

### 3. Auxiliary Functions May Use Wrong Rank

**Location:** `src/auxiliary/auxiliary_functions.jl:2`

**Issue:**
Some helper functions directly query `MPI.COMM_WORLD`:
```julia
rank = MPI.Comm_rank(MPI.COMM_WORLD)
```

**Impact:**
- Conditional logic based on rank may behave incorrectly
- Functions that should run on "rank 0" may not execute as expected

**Workaround:**
Most code goes through the `rank` variable which is correctly set. Only direct queries are affected.

## Testing Checklist

When testing coupling, verify:

- [ ] Mesh loads and partitions correctly
- [ ] Each code writes to its own output directory
- [ ] No file conflicts or corruption
- [ ] Rank-conditional code executes as expected
- [ ] Timing output makes sense (may be aggregated)
- [ ] No MPI deadlocks or communication errors

## Recommended Improvements

### Priority 1: I/O Functions
Refactor `write_output.jl` functions to accept `comm` parameter:

```julia
# Before
function write_output(SD, sol, mesh, OUTPUT_DIR, ...)
    comm = MPI.COMM_WORLD
    ...
end

# After
function write_output(SD, sol, mesh, OUTPUT_DIR, ...; comm=MPI.COMM_WORLD)
    rank = MPI.Comm_rank(comm)
    ...
end
```

Then call with:
```julia
write_output(...; comm=Jexpresso.comm)  # Uses coupling-aware comm
```

### Priority 2: Global Communicator Access

Create a global function to get the correct communicator:

```julia
# In run.jl or Jexpresso.jl
function get_mpi_comm()
    return Jexpresso.comm  # Returns comm_local when coupling, COMM_WORLD otherwise
end

function get_mpi_rank()
    return Jexpresso.rank  # Returns local_rank when coupling, world_rank otherwise
end

function get_mpi_size()
    return Jexpresso.nparts  # Returns local_size when coupling, world_size otherwise
end
```

Then replace hardcoded MPI calls:
```julia
# Instead of
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Use
comm = Jexpresso.get_mpi_comm()
rank = Jexpresso.get_mpi_rank()
```

### Priority 3: Separate Output Directories

Automatically set output directories based on coupling status:

```julia
# In run.jl after coupling initialization
if !isnothing(coupling_ctx)
    OUTPUT_DIR = joinpath(OUTPUT_DIR, coupling_ctx.code_name)
    mkpath(OUTPUT_DIR)
end
```

## Contributing Fixes

If you fix any of these issues:

1. Add tests in `test_jexpresso_coupling.jl`
2. Update this document
3. Document the API changes
4. Submit a pull request

## Questions?

See:
- `COUPLING_WITH_EXTERNAL_CODE.md` - Integration guide
- `COUPLING_TESTS.md` - Testing procedures
- `src/mpi/COUPLING_USAGE.md` - API documentation
