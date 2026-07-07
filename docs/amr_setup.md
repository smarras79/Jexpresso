# Setting up and running AMR cases (2D and 3D)

This is a practical guide to enabling adaptive mesh refinement (AMR) in a
Jexpresso case, based on the two minimal reference cases
(`problems/CompEuler/theta_amr` for 2D, `problems/CompEuler/3d_amr` for 3D)
and the most feature-complete AMR case in the tree
(`problems/CompEuler/giga_les_MOST_amr`, which additionally uses pre-adaptation
and VTK/p4est-forest restart).

AMR is built on `GridapP4est`/p4est: the mesh is a forest of octrees, and
refining/coarsening an element replaces it with a p4est octree operation
rather than a full remesh.

## 1. Minimal setup: `:lamr`

The smallest AMR-enabled case needs three keys in `user_inputs.jl`:

```julia
:lamr          => true,   # turns AMR on; this also flips the internal :ladapt flag
:amr_freq      => 200,    # re-evaluate refine/coarsen flags every N base timesteps
:amr_max_level => 1,      # maximum octree refinement level (each level halves Δx)
```

`amr_freq` and `Δt` together define the AMR cadence: internally
`Δt_amr = amr_freq * Δt`, and the ODE integrator is re-armed with a new
`tspan` of length `Δt_amr` after every AMR pass (`build_tspan` in
`src/io/mod_inputs.jl`, `amr_strategy!` in
`src/kernel/Adaptivity/Projection.jl`). A smaller `amr_freq` means more
frequent (and more expensive) remeshing.

`:Δt` in `user_inputs.jl` is the *base* (unrefined, level-0) timestep, not the
actual one used at runtime. The actual `dt` is rescaled every step by the
mesh's current maximum refinement level (`ad_lvl_max`, MPI-reduced across all
ranks):

```julia
dt = Δt / 2^ad_lvl_max
```

i.e. `ad_lvl_max = 0` → `dt = Δt`, `ad_lvl_max = 1` → `dt = Δt/2`,
`ad_lvl_max = 2` → `dt = Δt/4`, and so on (`src/kernel/solvers/TimeIntegrators.jl`,
computed right before each CFL check). This keeps the CFL number stable as
elements shrink under refinement — set `:Δt` for your coarsest (level-0)
resolution and let AMR shrink it automatically as `amr_max_level`/
`preadapt_max_level` push the mesh finer.

You also need a genuinely `N`-dimensional coarse mesh — set `:gmsh_filename`
to a mesh whose vertices carry exactly `Dc` coordinates for a `Dc`-dimensional
problem. If `GmshDiscreteModel` reports `Dp != Dc` (check with the snippet in
the FAQ), AMR will fail with an `AssertionError` deep in `GridapP4est`, not a
useful message about the mesh.

## 2. Required user code: `user_get_adapt_flags!`

Every AMR case must define, in its `initialize.jl`, the function that decides
which elements refine, coarsen, or stay put. Its signature is fixed by the
caller (`do_adapt!` in `src/kernel/Adaptivity/Projection.jl:3079-3093`) and is
the **same for 2D and 3D** — it always receives the full set of microphysics
fields even if your case doesn't use moisture (`Tabs, qn, qc, qi, qr, qs, qg,
Pr, Ps, Pg, S_micro, qsatt` are just unused parameters in a dry case):

```julia
function user_get_adapt_flags!(adapt_flags, inputs, old_ad_lvl, q, qe,
                               Tabs, qn, qc, qi, qr,
                               qs, qg, Pr, Ps, Pg,
                               S_micro, qsatt,
                               connijk, nelem, ngl,
                               coords,
                               max_level)
    # ... set adapt_flags[iel] = refine_flag / coarsen_flag / nothing_flag
end
```

Inside, loop over elements and set `adapt_flags[iel]` per element using your
own refinement criterion, gated by `old_ad_lvl[iel] < max_level` (don't refine
past the ceiling) and `old_ad_lvl[iel] > <floor>` (don't coarsen below it).
`refine_flag`, `coarsen_flag`, and `nothing_flag` are Jexpresso constants; any
element left unset defaults to `nothing_flag`.

**2D vs 3D differences** (see `theta_amr/initialize.jl` vs
`giga_les_MOST_amr/initialize.jl`):
- 2D: `connijk[iel, i, j]` indexed with two loop indices `i,j = 1:ngl`
  (`ngl*ngl` nodes per element).
- 3D: `connijk[iel, i, j, k]` indexed with three loop indices
  (`ngl*ngl*ngl` nodes per element).
- `coords` is `mesh.coords`, an `(npoin, Dc)` array — `coords[:,1]`/`[:,2]` are
  x/y in both cases, `coords[:,3]` is z in 3D only.

`theta_amr`'s version is the simplest reference: it refines where a
potential-temperature perturbation `θ - θ_ref` exceeds a tolerance, and
coarsens everywhere else. `giga_les_MOST_amr`'s version is more elaborate
(refines on cloud/precipitation fields `qn`, with height-dependent coarsening
floors) but follows the exact same contract.

## 3. Two different kinds of "pre-refinement"

These are easy to confuse — they run at different times and solve different
problems.

### `:linitial_refine` — static uniform refinement, no user function needed

```julia
:linitial_refine => true,
:init_refine_lvl => 1,   # uniform refinement levels applied once, before t=0
```

Uniformly refines the whole coarse mesh `init_refine_lvl` times before the
run starts. No refine/coarsen criterion involved — every element gets
refined the same number of times. Cheap way to start at a finer base
resolution without writing a finer mesh file.

### `:lpreadapt` — criterion-driven pre-refinement, needs `user_get_preadapt_flags!`

```julia
:lpreadapt         => true,
:preadapt_max_level => 2,
```

Unlike `:linitial_refine`, pre-adaptation is *spatially selective*: before
t=0, Jexpresso repeatedly calls a user-supplied flag function and refines
only the flagged elements, iterating until no new refinement occurs or
`preadapt_max_level` is reached (`mod_mesh_mesh_driver` in
`src/kernel/mesh/mesh.jl:5163-5199`). This is how `giga_les_MOST_amr` puts
extra resolution in, say, the boundary layer and lower troposphere from the
very first timestep, without wasting resolution near the domain top.

Requires a second user function, also in `initialize.jl`:

```julia
function user_get_preadapt_flags!(adapt_flags, inputs, mesh, old_ad_lvl, connijk, nelem, ngl, max_level)
    for iel = 1:nelem
        for i = 1:ngl, j = 1:ngl, k = 1:ngl   # k loop only for 3D
            ips = connijk[iel, i, j, k]
            x, y, z = mesh.x[ips], mesh.y[ips], mesh.z[ips]
            # ... set adapt_flags[iel] based on geometry alone (no solution
            #     fields exist yet — this runs before initialize())
        end
    end
end
```

Note the criterion here can only depend on **geometry** (`mesh.x/y/z`), since
there is no solution field yet at pre-adaptation time — that's the key
difference from `user_get_adapt_flags!`, which sees `q`/`qe`.

`:lamr`, `:lpreadapt`, and `:linitial_refine` are independent switches and can
combine (e.g. pre-adapt to bias the initial mesh, then let runtime AMR take
over), but a case using `:lpreadapt` must define
`user_get_preadapt_flags!`, and any case using `:lamr` must define
`user_get_adapt_flags!` — both are simple `UndefVarError`/`MethodError`s at
runtime if missing, not caught at setup time.

## 4. Restart with AMR: `:lrestart_amr`

```julia
:lrestart_amr => true,
# :restart_vtk_iout => 12,   # optional: explicit iter_N to resume from;
                             # if omitted, restarts from the LAST entry in
                             # simulation.pvd (read_pvd_last_entry)
```

This is the mechanism that actually resumes AMR runs correctly. On write
(every diagnostic-output step, gated by `:lamr`), Jexpresso saves the p4est
forest alongside the VTK output (`write_p4est_checkpoint`, called from the
`affect!` callback in `src/kernel/solvers/TimeIntegrators.jl:413-418`,
producing `iter_N/iter_N.p4est`). On restart:

1. `sem_setup.jl` (`src/kernel/infrastructure/sem_setup.jl:138-169`) loads the
   forest via `load_p4est_checkpoint_model` and rebuilds the exact adapted
   mesh — not just the coarse one — using `mod_mesh_mesh_driver` with
   all-`nothing_flag` (no further adaptation, just reconstruction), then
   restores `mesh.ad_lvl` from the loaded forest.
2. Your case's `initialize.jl` must then read the matching solution data back
   in. Add this at the end of `initialize()`, after building/initializing
   `q` as normal:

   ```julia
   if get(inputs, :lrestart_amr, false)
       read_vtk_amr_restart!(q, mesh, inputs; output_dir=OUTPUT_DIR,
                        varnames=["ρ", "u", "w", "θ", "p"])  # your case's qoutvars
   end
   ```

   `read_vtk_amr_restart!` (in `src/io/write_output.jl`) is a case-agnostic
   utility — it takes the point-data variable names to read as `varnames`
   (defaulting to giga_les's 7-variable moist set if omitted), so pass your
   own case's `qoutvars` list here. It reads the VTK output at the resolved
   `iter_N` and interpolates the requested fields onto the just-rebuilt mesh.

3. You must also define `user_read_vtu_point_data!` in the case's
   `user_primitives.jl` — it's how `read_vtk_amr_restart!` turns the primitive
   VTK fields back into your conserved state `q.qn`:

   ```julia
   function user_read_vtu_point_data!(q, vars, ip_map, mesh)
       # vars is the Dict returned by read_vtu_point_data — keyed by the
       # varnames you passed to read_vtk_amr_restart!.
       # ip_map maps mesh point index -> VTK array index (handles any
       # reordering between the write and this restart's rank/partition).
       for ip = 1:mesh.npoin
           j = ip_map[ip]
           ρ = vars["ρ"][j]
           q.qn[ip, 1] = ρ
           q.qn[ip, 2] = ρ * vars["u"][j]
           # ... reconstruct the rest of your conserved variables ...
           q.qn[ip, end] = vars["p"][j]
       end
   end
   ```

   `read_vtk_amr_restart!` raises a clear error naming this function if it's
   missing, rather than failing deep inside VTK parsing.

`theta_amr`/`3d_amr` don't use `:lrestart_amr` at all — they're the
"AMR runs, no restart" reference. `theta_amr_hang` (2D) and
`giga_les_MOST_amr` (3D) are the references for the full VTK/p4est restart
path.

## 5. Quick-reference: AMR-related `user_inputs.jl` keys

| Key | Default | Meaning |
|---|---|---|
| `:lamr` | `false` | Master AMR switch; also sets `:ladapt = true` internally |
| `:amr_freq` | `0` | Base timesteps between AMR passes |
| `:amr_max_level` | `0` | Max refinement level `user_get_adapt_flags!` may request |
| `:amr_start_time` | `0.0` | Delay before the first AMR pass (added into the first `tspan`) |
| `:linitial_refine` | `false` | Uniform static refinement before t=0, no user function |
| `:init_refine_lvl` | `0` | How many uniform refinement passes for `:linitial_refine` |
| `:lpreadapt` | `false` | Criterion-driven static refinement before t=0; also sets `:ladapt = true` |
| `:preadapt_max_level` | `0` | Max level for `:lpreadapt`; also used as the floor `user_get_adapt_flags!` coarsens back down to |
| `:lrestart_amr` | `false` | VTK + p4est forest restart (resumes the real adapted mesh) |
| `:restart_vtk_iout` | auto | Explicit `iter_N` to resume from; omit to restart from the last entry in `simulation.pvd` |

## 6. Testing locally

Invoke as a script, not `-e '...include(...)'` — the latter leaves
`PROGRAM_FILE` unset and silently skips the whole driver (see the
`_LOADED_CASE_DIR`/lazy-loader ordering notes in `src/Jexpresso.jl` if this
ever regresses):

```bash
mpirun -np 1 julia --project=. src/Jexpresso.jl CompEuler theta_amr
mpirun -np 1 julia --project=. src/Jexpresso.jl CompEuler 3d_amr
```

Make sure the `mpirun`/`mpiexec` on your `PATH` matches the MPI
`MPIPreferences`/`MPI.jl` is actually bound to — see
[FAQ.md](../FAQ.md#a-parallel-run-hangs-forever-at-mpiinit--first-collective).
A mismatch here can surface as a confusing, unrelated-looking MPI error deep
inside a `p4est_*` call (e.g. `internal_Comm_size` / "before initializing or
after finalizing") rather than an obvious launcher error.
