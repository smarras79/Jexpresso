# Coupling Tests for Jexpresso

This directory contains test scripts to verify the MPI coupling functionality.

## Test Scripts

### 1. `test_simple_coupling.jl` - Basic Coupling Test (Recommended First)

**Purpose**: Tests the coupling infrastructure without running full Jexpresso simulations.

**What it does**:
- Splits MPI ranks into two mock "codes" (Jexpresso-A and Jexpresso-B)
- Initializes coupling contexts
- Displays rank assignments and communicator information
- Performs a simple data exchange between the two codes
- Verifies data was correctly exchanged and broadcast

**Usage**:
```bash
# Run with 4 ranks (2 per code)
mpirun -np 4 julia --project=. test_simple_coupling.jl

# Run with 6 ranks (3 per code)
mpirun -np 6 julia --project=. test_simple_coupling.jl

# Run with just 2 ranks (minimal - 1 per code)
mpirun -np 2 julia --project=. test_simple_coupling.jl
```

**Expected Output**:
```
======================================================================
COUPLING TEST STARTED
======================================================================
[Rank  0] I am Jexpresso-A (code_id=1)
[Rank  1] I am Jexpresso-B (code_id=2)
...
[ Info: [Jexpresso-A] Coupling initialized: local_size=2, world_rank=0, ...
[ Info: [Jexpresso-B] Coupling initialized: local_size=2, world_rank=2, ...
...
[Jexpresso-A Root] AFTER exchange:
  Send: [10.0, 20.0, 30.0, 40.0, 50.0]
  Recv: [100.0, 200.0, 300.0, 400.0, 500.0]
  ✓ Data exchange SUCCESSFUL!
...
COUPLING TEST COMPLETED SUCCESSFULLY!
```

**What to verify**:
- [ ] Both codes initialize successfully
- [ ] Each code shows correct local_size (half of world_size)
- [ ] Coupling roots are identified (world_rank 0 and world_size/2)
- [ ] Data exchange completes without deadlock
- [ ] Received data matches expected values
- [ ] All ranks receive data after broadcast

---

### 2. `test_jexpresso_coupling.jl` - Full Jexpresso Coupling Test

**Purpose**: Tests coupling with two actual Jexpresso instances running simultaneously.

**What it does**:
- Launches two complete Jexpresso simulations with coupling enabled
- Each instance has a unique name tag (Jexpresso-A, Jexpresso-B)
- Both run the same test case but in separate communicators
- Can be extended to exchange actual simulation data

**Usage**:
```bash
# Run with 4 ranks (2 per Jexpresso instance)
mpirun -np 4 julia --project=. test_jexpresso_coupling.jl

# Run with more ranks for realistic simulation
mpirun -np 8 julia --project=. test_jexpresso_coupling.jl
```

**Before running**: Make sure you have a valid test case set up. By default, it runs `CompEuler/theta`. You can edit the script to use a different case.

**Expected behavior**:
- Two separate Jexpresso instances start
- Each shows its own rank assignments
- Both run their simulations independently
- No conflicts in communicators or file I/O

**What to verify**:
- [ ] Both Jexpresso instances start without errors
- [ ] Each instance shows correct local MPI configuration
- [ ] No MPI communicator conflicts
- [ ] Both simulations complete successfully
- [ ] Output files are properly separated (if applicable)

---

## Troubleshooting

### Problem: "Need at least 2 MPI ranks"
**Solution**: Coupling requires at least 2 ranks. Use `mpirun -np 2` or higher.

### Problem: "UndefVarError: JexpressoCoupling not defined"
**Solution**: Make sure you're running from the Jexpresso root directory and the module is properly installed.

### Problem: Deadlock (program hangs)
**Solution**:
- Ensure both codes are calling exchange functions
- Check that rank splitting is correct
- Verify MPI configuration allows inter-communicator operations

### Problem: "Data exchange FAILED"
**Solution**:
- Check MPI installation supports inter-communicator communication
- Verify both codes have the same buffer sizes
- Try with fewer ranks first (2 or 4)

### Problem: Different data on non-root ranks
**Solution**: This likely means `broadcast_to_local!` isn't working. Check MPI configuration.

---

## Next Steps

Once these tests pass:

1. **Test with your external code**: Use `launch_coupled.jl` template
2. **Add data exchange to Jexpresso**: Modify boundary conditions or user files
3. **Performance testing**: Try with realistic processor counts
4. **Production runs**: Launch your coupled simulations

---

## Test Checklist

Use this checklist when validating your coupling setup:

- [ ] `test_simple_coupling.jl` with 2 ranks: PASS
- [ ] `test_simple_coupling.jl` with 4 ranks: PASS
- [ ] `test_simple_coupling.jl` with 6+ ranks: PASS
- [ ] `test_jexpresso_coupling.jl` with 4 ranks: PASS
- [ ] `test_jexpresso_coupling.jl` with 8+ ranks: PASS
- [ ] Both Jexpresso instances show distinct names: PASS
- [ ] Data exchange completes without deadlock: PASS
- [ ] All ranks receive correct data: PASS

---

## For Developers

### Adding Custom Tests

To add your own coupling tests:

1. Copy `test_simple_coupling.jl` as a template
2. Modify the data exchange logic for your needs
3. Add verification checks for your specific data
4. Document expected behavior

### Debugging Tips

```julia
# Add these to your test for detailed MPI info:
println("My communicator: $(MPI.Comm_rank(ctx.comm_local))")
println("Communicator size: $(MPI.Comm_size(ctx.comm_local))")
println("World rank: $(MPI.Comm_rank(MPI.COMM_WORLD))")

# Force synchronization to check for deadlocks:
MPI.Barrier(MPI.COMM_WORLD)
println("Checkpoint reached on rank $world_rank")
```

---

## Support

If tests fail, check:
1. MPI installation: `mpirun --version`
2. Julia MPI.jl configuration: `julia -e 'using MPI; MPI.versioninfo()'`
3. Jexpresso installation: `julia --project=. -e 'using Jexpresso'`

For coupling issues, see `src/mpi/COUPLING_USAGE.md`.
