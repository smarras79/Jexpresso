# ReducedGravityDSGS

This case uses the wind-driven reduced-gravity DoubleGyre configuration from
`problems/ShallowWater/ReducedGravity`, replacing its constant artificial
viscosity with the `DSGS_SW()` stabilization configured in
`problems/ShallowWater/SoliWaveIslandDSGS`.

The physics, initial condition, boundary conditions, primitive variables, and
mesh are loaded directly from the ReducedGravity sibling. Only
`user_inputs.jl` is case-specific.

The DynSGS wave-speed cap uses the reduced gravity `g' = 0.015 m/s²`, and the
wet/dry regularization uses `Hmin = 10⁻³ m`, matching the ReducedGravity flux.

Run from the repository root with:

```bash
julia --project=. src/Jexpresso.jl ShallowWater ReducedGravityDSGS
```

The production deck integrates six 360-day model years with a 300-second
timestep, so use a shortened copy of the deck for a quick smoke test.

To write restart checkpoints, add a positive interval (in seconds) to the deck:

```julia
:restart_time => 30 * day,
```

The latest checkpoint is saved in `restart/` inside the run's output directory,
independently of the VTK diagnostic output. It contains `t.h5`, the conservative
fields and reference state in `var_1_<rank>.h5` through `var_3_<rank>.h5`.

To resume, set these inputs and keep `:tend` later than the checkpoint time:

```julia
:lrestart => true,
:restart_input_file_path => "/absolute/path/to/previous/run/restart",
```

The checkpoint sets `:tinit` automatically. Checkpoint writing is disabled by default.

To restart from VTK output instead, point at the previous run's output directory
(the directory containing `simulation.pvd`):

```julia
:lrestart_vtk => true,
:restart_vtk_input_dir => "/absolute/path/to/previous/run",
# Optional: select a particular snapshot; otherwise use the latest PVD entry.
# :restart_vtk_iout => 2,
```

Set only one of `:lrestart` and `:lrestart_vtk`. VTK restart reads the physical
time from the selected snapshot, independently of the new diagnostic schedule.
New VTK files include conservative and reference fields; existing ReducedGravity
VTK files can use their `H`, `Hu`, and `Hv` fields with the case's lake-at-rest
reference. Use the same mesh, polynomial order, MPI rank count, timestep, and
model settings. AMR restart remains unsupported for this case.

Both formats deliberately reinitialize DynSGS history from the restored
solution, including the startup hold. No DynSGS history is written or read, so
continuation can differ from an uninterrupted run. Existing history sidecars
are ignored.

The regression test checks both formats, restored states and times, fresh
DynSGS history, continuation, checkpoint cadence, and warm-up behavior:

```bash
julia --project=. test/test_shallow_water_restart.jl
```
