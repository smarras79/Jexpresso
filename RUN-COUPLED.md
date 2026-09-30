# Running Jexpresso coupled with Alya (Fortran)

This guide explains how to run **Jexpresso (Julia)** coupled with the
**Alya proxy (Fortran)**. The single most important requirement for a coupled
run is covered first, because almost every coupled-launch failure comes from
getting it wrong:

> ## ⚠️ Both codes MUST use the *same* MPI library
>
> Jexpresso and Alya are launched as **one** MPI job (MPMD) that share a single
> `MPI_COMM_WORLD` — see [COUPLING-ALGORITHM.md](COUPLING-ALGORITHM.md), §1.
> The two executables therefore have to be built and run against the **same MPI
> implementation, the same version, and the same ABI**:
>
> - **Fortran side (Alya):** linked at compile time by `mpif90`.
> - **Julia side (Jexpresso):** bound by `MPI.jl` via `MPIPreferences`.
>
> If the Julia `MPI.jl` is bound to a *different* MPI than the one `mpif90` used
> to build `Alya.x` (for example OpenMPI on one side and MPICH on the other, or
> even two different OpenMPI builds), the job will fail to start, **hang forever
> inside `MPI_Init`**, or silently corrupt the exchanged data. There is no
> partial credit here — they must match exactly.

The rest of this document walks you through (1) identifying which MPI each side
uses, (2) making them match, (3) compiling Alya, and (4) launching the coupled
run. For general (non-coupled) MPI setup of Jexpresso, see
[INSTALL.md, Section 5 — "Running in parallel with MPI"](INSTALL.md#5-running-in-parallel-with-mpi).

> **Shortcut: MPI.jl on its default binding, the MPICH bundled with Julia
> (`MPICH_jll`).** Then there is nothing to reconcile and no system MPI to
> install — only `gfortran`. `compilef90.sh` builds Alya against that very MPICH,
> and `run_coupled.sh` launches with the `mpiexec` that ships with it:
>
> ```bash
> cd AlyaProxy && bash compilef90.sh && cd ..
> ./run_coupled.sh 2 2
> ```
>
> Section 2 shows how to check the binding (`binary = MPICH_jll`), and Section 4
> what the two scripts do. On macOS the
> [hostname fix](INSTALL.md#55-macos-hostname-fix-mpich-and-mpich_jll-only) for
> MPICH applies to it as well.

---

## 1. Identify the MPI used by the Fortran side (Alya)

`Alya.x` is compiled with the `mpif90` wrapper (see
[`AlyaProxy/compilef90.sh`](AlyaProxy/compilef90.sh)). Whichever MPI owns the
`mpif90` on your `PATH` is the MPI that Alya — and therefore the **whole coupled
job** — must use. Find out exactly what it is:

```bash
which mpif90              # which wrapper will be used
mpif90 --version         # implementation + version (e.g. "gfortran ..."/OpenMPI/MPICH banner)
mpif90 -show             # the underlying compile/link line — shows -I/-L paths and -lmpi
mpiexec --version        # the matching launcher's implementation + version
```

From this, write down three things — you will match all three on the Julia side:

| What | How to read it |
|------|----------------|
| **Implementation** | OpenMPI or MPICH (printed by `mpif90 --version` / `mpiexec --version`) |
| **Version** | e.g. `5.0.6` (OpenMPI) or `4.2.x` (MPICH) |
| **Library path** | the `-L/.../lib` shown by `mpif90 -show`, e.g. `/opt/homebrew/lib` |

> **Tip — make sure the right `mpif90` is first on `PATH`.** If you have several
> MPIs installed, `which mpif90`, `which mpiexec`, and `which mpirun` should all
> resolve to the **same** installation. If they don't, fix your `PATH` (or use
> absolute paths everywhere) before going any further.

---

## 2. Identify the MPI used by the Julia side (MPI.jl)

Ask `MPI.jl` what it is currently bound to:

```bash
julia --project=. -e '
  using MPIPreferences; println("binary = ", MPIPreferences.binary)
  using MPI;            println(MPI.identify_implementation())'
```

- `binary = "system"` means MPI.jl uses a system MPI; `identify_implementation()`
  then prints which one (e.g. `(OpenMPI, v"5.0.6")`).
- `binary = "MPICH_jll"` — MPI.jl's default — means MPI.jl uses the MPICH that
  ships with Julia. That works for coupling as it is: `compilef90.sh` builds
  Alya against that same MPICH and `run_coupled.sh` launches with its own
  `mpiexec` (Section 4). Skip to Section 5.
- `binary = "MPItrampoline_jll"` or `"OpenMPI_jll"` means another bundled MPI,
  which the coupling scripts do not support. Bind MPI.jl to a system MPI
  (Section 3, Option A) or back to `MPICH_jll`.

**For a system MPI, this must report the same implementation and version you
wrote down in Section 1.** If it doesn't, fix it in the next step.

---

## 3. Compare the two sides — and reconcile them if they differ

Put the output of Section 1 (Fortran) next to Section 2 (Julia). Both must agree
on all three rows:

| Must match | Fortran (Alya) | Julia (Jexpresso) |
|------------|----------------|-------------------|
| Implementation | `mpif90 --version` | `MPI.identify_implementation()` |
| Version | `mpiexec --version` | `MPI.identify_implementation()` |
| Library / launcher | `mpif90 -show` → `-L.../lib` | `MPIPreferences.binary == "system"` + same `mpiexec` |

**If all three rows agree, skip to Section 5.** For example, the verified setup

```
Julia:   binary = system   ("OpenMPI", v"5.0.8")
Fortran: mpif90 -show  →  -L/opt/homebrew/Cellar/open-mpi/5.0.8/lib
```

already matches (OpenMPI 5.0.8 on both), so no reconciliation is needed.

### If they do NOT match

Pick **one** of the two directions below — whichever lets both sides land on the
same implementation **and** version. Changing one side means leaving the other
alone.

#### Option A *(usual)* — rebind the Julia side to the Fortran side's MPI

Use this when `mpif90` already points at the MPI you want for the whole job
(e.g. the cluster/Homebrew default). Point `MPI.jl` at the **same** library that
`mpif90 -show` reported, then rebuild and precompile:

```bash
# Use the SAME lib dir that `mpif90 -show` printed.
# Homebrew OpenMPI on Apple Silicon:
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=["/opt/homebrew/lib"])'
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

> The path is the `-L.../lib` from `mpif90 -show` — e.g.
> `/opt/homebrew/Cellar/open-mpi/5.0.8/lib`, of which `/opt/homebrew/lib` is the
> Homebrew symlink. If the MPI is in a standard location (`/usr/bin`,
> `/usr/local/bin`), `MPIPreferences.use_system_binary()` with no `extra_paths`
> suffices. Full options:
> [INSTALL.md, Section 5.2](INSTALL.md#52-point-mpipreferences-at-your-mpi).

#### Option B — rebuild the Fortran side against the Julia side's MPI

Use this when you want to keep the MPI that `MPI.jl` is already bound to (e.g.
you don't want to re-precompile Julia, or MPI.jl points at the MPI you prefer).
Recompile `Alya.x` with the `mpif90` that belongs to **that** MPI — pin it
explicitly via the `MPIF90` variable so the right wrapper is used regardless of
`PATH` (see Section 5):

```bash
cd AlyaProxy
# Point MPIF90 at the wrapper of the MPI that MPI.jl reports. Example:
MPIF90=/opt/homebrew/Cellar/open-mpi/5.0.8/bin/mpif90 bash compilef90.sh
cd ..
```

If MPI.jl is on the MPICH bundled with Julia (`MPICH_jll`), there is no system
wrapper to find and none is needed: `compilef90.sh` builds against that MPICH
directly (Section 4). For the other bundled MPIs (`MPItrampoline_jll`,
`OpenMPI_jll`) install a matching **system** MPI and use Option A instead.

### Confirm the match

Re-run the Section 2 check and re-read `mpif90 -show`; the three rows above must
now line up. **Only proceed once they do.**

---

## 4. Choose a consistent MPI for coupling

For a coupled run, pick **one** of the following — the same choice applies to
*both* codes:

- **OpenMPI for both** *(recommended on Linux/HPC).* Build Alya with OpenMPI's
  `mpif90`; bind MPI.jl with `use_system_binary` to the same OpenMPI `lib`.
- **MPICH for both** *(recommended on macOS, where OpenMPI 5 can deadlock in
  `MPI_Init`).* Build Alya with MPICH's `mpif90`; bind MPI.jl to the same MPICH
  `lib`. On macOS, also apply the
  [hostname `/etc/hosts` fix](INSTALL.md#55-macos-hostname-fix-mpich-and-mpich_jll-only)
  required for any MPICH-based MPI.
- **The MPICH bundled with Julia (`MPICH_jll`) for both** *(MPI.jl's default;
  no system MPI needed, only `gfortran`).* Leave MPI.jl as it is.
  `compilef90.sh` sees the binding and builds `Alya.x` against that very MPICH;
  `run_coupled.sh` launches with the `mpiexec` that ships with it, started with
  the library path MPI.jl gives it. The macOS hostname fix above applies here
  too.

> **Same installation, not just the same family.** Two builds of MPICH in one
> job — a Homebrew MPICH for Alya and Julia's `MPICH_jll`, say — are two
> installations, and they are not guaranteed to form a working world. So with
> MPI.jl on `MPICH_jll`, Alya is built against `MPICH_jll` itself, and
> `run_coupled.sh` refuses to launch an `Alya.x` that loads its `libmpi` from
> anywhere else. The JLL's own `bin/mpifort` cannot do this build: it still
> points at the directory the JLL was built in, `/workspace/destdir`. So
> `compilef90.sh` gives `gfortran` the JLL's include and library directories
> directly and sets rpaths, so that `Alya.x` runs with no environment set.
> The JLL ships no `mpi_f08` module; the build uses `mpif.h`, which the proxy
> supports.

---

## 5. Compile the Alya proxy with the chosen MPI

`Alya.x` is built from `myAlya.f90` (a symlink to `alya_all2all_time_loop.f90`)
by [`AlyaProxy/compilef90.sh`](AlyaProxy/compilef90.sh), which calls `mpif90`.
The only thing that matters is that this `mpif90` belongs to the **same MPI**
you bound MPI.jl to in Section 3.

**In the normal case, this is the whole step:**

```bash
cd AlyaProxy
bash compilef90.sh
cd ..
```

With no arguments the script asks Julia what MPI.jl is bound to and builds
against that same implementation — which is exactly the rule from Section 1, so
there is nothing for you to keep in sync by hand. It then finds that
implementation's `mpif90` **by prefix, not by `PATH`** (Homebrew can only *link*
one MPI at a time, so `which mpif90` may well be the other one's), verifies the
wrapper really is what it claims, picks the Fortran binding the MPI actually
provides, and checks the finished binary against MPI.jl — exiting non-zero if
they disagree. When MPI.jl is on the MPICH bundled with Julia (`MPICH_jll`),
there is no wrapper to find: it builds against that MPICH directly and checks
that `Alya.x` loads that very `libmpi` (Section 4). If Julia cannot be asked,
the script prints Julia's error rather than guessing.

To build against a specific implementation instead, name it:

```bash
bash compilef90.sh mpich
bash compilef90.sh openmpi
```

`bash compilef90.sh --help` lists everything, including `MPIF90=` to pin an
exact wrapper, `USE_MPIF08=0|1` to force the Fortran binding, `FFLAGS_EXTRA=`
for extra compiler flags, and `SKIP_JULIA_CHECK=1` for machines with no
Jexpresso environment.

§5a and §5b below spell out the per-implementation details — mainly what to
check afterwards, and the `PATH` setup the *launcher* still needs.

> **`ld: library not found for -lSystem` (macOS).** A toolchain failure, not an
> MPI one — compilation succeeded and only the link failed. gfortran passes the
> linker `-syslibroot $SDKROOT`, and the SDK it was built against has moved or
> gone, which is what an Xcode / Command Line Tools update does. The script
> exports `SDKROOT=$(xcrun --show-sdk-path)` for you; if it still fails, check
> that `xcrun --show-sdk-path` returns a directory that exists, then
> `sudo xcode-select --reset && xcode-select --install`, then
> `brew reinstall gcc` to rebuild gfortran against the current SDK.

### 5a. Building against OpenMPI

Use this when Section 3 left MPI.jl bound to OpenMPI. Build:

```bash
cd AlyaProxy
bash compilef90.sh openmpi
cd ..
```

The script locates OpenMPI's own `mpif90` itself (via `brew --prefix open-mpi`
on macOS, `/usr/lib/*/openmpi/bin` and friends on Linux). If yours is somewhere
it cannot find, point at it directly:

```bash
cd AlyaProxy
MPIF90=/path/to/openmpi/bin/mpif90 bash compilef90.sh
cd ..
```

Confirm the result links OpenMPI:

```bash
otool -L AlyaProxy/Alya.x | grep -i mpi
```

You want `libmpi.40.dylib` under the OpenMPI prefix. If you see
`libmpi.12.dylib` or a path containing `mpich`, the wrong wrapper was used —
rebuild with `MPIF90` set as above.

Put OpenMPI's `bin` first on `PATH` for the shell you launch from, so the
`mpirun` that starts the job belongs to the same OpenMPI:

```bash
OMPI_PREFIX=$(brew --prefix open-mpi)
export PATH="$OMPI_PREFIX/bin:$PATH"
hash -r
which mpiexec mpirun mpif90
```

All three must resolve inside `$OMPI_PREFIX`. On Linux use that
installation's `bin` directly, e.g.
`export PATH=/usr/lib/x86_64-linux-gnu/openmpi/bin:$PATH`.

> **OpenMPI 5 on macOS.** This combination is prone to hanging inside
> `MPI_Init` before anything runs — see Section 7. If that happens, switch to
> MPICH (§5b) rather than fighting it.

### 5b. Building against MPICH

Use this when Section 3 left MPI.jl bound to MPICH. MPICH is the recommended
implementation on macOS (Section 4), and the one to move to if OpenMPI hangs in
`MPI_Init`. Bind the Julia side first
([INSTALL.md §5.2, case 2](INSTALL.md#52-point-mpipreferences-at-your-mpi)),
then build:

```bash
cd AlyaProxy
bash compilef90.sh mpich
cd ..
```

The script locates MPICH's own `mpif90` itself (via `brew --prefix mpich` on
macOS, `/usr/lib/*/mpich/bin` and friends on Linux). If yours is somewhere it
cannot find, point at it directly:

```bash
cd AlyaProxy
MPIF90=/path/to/mpich/bin/mpif90 bash compilef90.sh
cd ..
```

Confirm the result links MPICH and **not** OpenMPI:

```bash
otool -L AlyaProxy/Alya.x | grep -i mpi
```

You want `libmpi.12.dylib` / `libmpifort.12.dylib` under the MPICH prefix. If
you instead see `libmpi.40.dylib` or any path containing `open-mpi`, the wrong
wrapper was used — rebuild with `MPIF90` set as above.

Put MPICH's `bin` first on `PATH` for the shell you launch from:

```bash
MPICH_PREFIX=$(brew --prefix mpich)
export PATH="$MPICH_PREFIX/bin:$PATH"
hash -r
which mpiexec mpirun mpif90
```

All three must resolve inside `$MPICH_PREFIX`. On Linux use that
installation's `bin` directly, e.g.
`export PATH=/usr/lib/x86_64-linux-gnu/mpich/bin:$PATH`.

> **`Cannot open module file 'mpi_f08.mod'`.** Some MPICH builds ship no
> Fortran 2008 module (and a `gfortran` newer than the one that built MPI
> cannot read the module even when it is present). This is not fatal: the proxy
> also compiles against the legacy `include 'mpif.h'` bindings, and
> `compilef90.sh` detects the situation and switches automatically. The MPI
> calls are identical — only the handle types differ. Force either binding with
> `USE_MPIF08=0` or `USE_MPIF08=1` if you need to.
>
> On that legacy path the script also adds `-fallow-argument-mismatch`. With
> `mpif.h` there are no interfaces, so every MPI routine is an external
> procedure, and gfortran ≥ 10 refuses a file that calls one external with two
> different argument types:
> ```
> Error: Type mismatch between actual argument at (1) and actual argument
> at (2) (REAL(8)/INTEGER(4)).
> ```
> That is `MPI_Bcast` being called on an `INTEGER` buffer in one place and a
> `REAL(8)` buffer in another — ordinary, correct MPI usage that the F08 module
> would have typed properly. The script probes for the flag (falling back to
> `-Wno-argument-mismatch` on gfortran 9). For a compiler needing something
> else, pass it through: `FFLAGS_EXTRA=-some-flag bash compilef90.sh`.

> **MPICH on macOS needs two more things.** Apply the
> [`/etc/hosts` hostname fix](INSTALL.md#55-macos-hostname-fix-mpich-and-mpich_jll-only),
> and if `MPI_Init`/`MPI_Finalize` fails with an OFI/libfabric error naming a
> VPN interface (`nic=utun7`) or no interface at all (`nic=(n/a)`), pin the
> fabric — `export FI_PROVIDER=tcp; export FI_TCP_IFACE=lo0`. Both are covered
> in [INSTALL.md §5.7](INSTALL.md#57-troubleshooting).

> **MPICH does not accept OpenMPI's `-x` flag.** Its launcher is Hydra, which
> uses `-env VAR VAL` / `-genv VAR VAL` and errors with
> `match_arg: unrecognized argument x` if given `-x`. This is why Section 6
> tells you to `export JEXPRESSO_COUPLED=1` instead — that form works with both
> launchers.

### 5c. Verify both sides agree

§5a and §5b each checked `Alya.x` in isolation. This step puts it next to the
Julia side, which is the comparison that actually matters — it catches a wrong
wrapper even when `-show` looked fine.

Print what the Fortran side linked. On macOS:

```bash
otool -L AlyaProxy/Alya.x | grep -i mpi
```

On Linux:

```bash
ldd AlyaProxy/Alya.x | grep -i mpi
```

Print what the Julia side loads:

```bash
julia --project=. -e 'using MPI; println(MPI.API.libmpi); println(MPI.identify_implementation())'
```

The two must name the **same installation**, not merely the same family:

| You see | That is | Must be paired with |
|---|---|---|
| `libmpi.40.dylib`, path contains `open-mpi` | OpenMPI | MPI.jl bound to the same OpenMPI |
| `libmpi.12.dylib`, `libmpifort.12.dylib`, path contains `mpich` | MPICH | MPI.jl bound to the same MPICH |

`bash tools/check_mpi_setup.sh` compares these two sides for you (Stage 5) and
refuses to draw conclusions from an Alya run when they disagree.

This produces `./AlyaProxy/Alya.x`, linked against the MPI you selected in
Section 4. If you later switch the Julia side to a different MPI, you **must
recompile `Alya.x`** with the matching `mpif90` — the two always travel together.

---

## 6. Launch the coupled run (MPMD)

A coupled run is a single MPMD launch: one launcher starts Alya ranks **and**
Julia ranks, separated by `:`, so they share `MPI_COMM_WORLD`. Use the launcher
that belongs to the MPI both codes were built against.

The `JEXPRESSO_COUPLED=1` environment variable tells Jexpresso it is running
coupled. The example case is `CompEuler 3dAlya` (3D), which matches the proxy as
shipped: same domain (`hexa_TFI_10x1x10.msh` spans exactly the proxy's
`rem_min`/`rem_max`) and the same 2000 steps.

`CompEuler thetaAlya` (2D) needs a 2D proxy first. In
`AlyaProxy/alya_all2all_time_loop.f90` set `ndime = 2`,
`rem_min = [-5000.0, 0.0, 0.0]`, `rem_max = [5000.0, 10000.0, 0.0]` and
`rem_nx = [10, 10, 1]` (the domain of its `hexa_TFI_10x10.msh`), then rebuild
`Alya.x`. Against the 3D proxy Jexpresso stops at setup with `nfields mismatch`:
it sends two velocity components and the proxy expects `ndime = 3`.

> **Export it, don't pass it with `-x` after the `:`.** Not every launcher
> honours `-x VAR=value` inside a *secondary* MPMD app context. When it is
> silently dropped, Jexpresso takes its **standalone** path, runs its
> collectives on the full `MPI_COMM_WORLD` (Alya's ranks included) and
> deadlocks with no error message at all — the job just sits there. Exporting
> the variable for the whole job cannot fail this way, and Alya ignores it.
> Verify your launcher's behaviour with `bash tools/check_mpi_setup.sh`
> (Stage 6).

### OpenMPI or MPICH (system MPI) — `mpirun`

```bash
export JEXPRESSO_COUPLED=1

# 3D: 2 Alya ranks + 2 Jexpresso ranks
mpirun -np 2 ./AlyaProxy/Alya.x : -np 2 \
    julia --project=. ./src/Jexpresso.jl CompEuler 3dAlya
```

```bash
# 2D (only with the 2D proxy described above)
export JEXPRESSO_COUPLED=1
mpirun -np 2 ./AlyaProxy/Alya.x : -np 2 \
    julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya
```

> **`prterun` error (OpenMPI 5)?** OpenMPI's launcher sometimes cannot resolve
> the `juliaup` shim. Pass the real `julia` binary explicitly:
>
> ```bash
> JULIA_BIN=$(julia -e 'print(joinpath(Sys.BINDIR, "julia"))')
> echo "$JULIA_BIN"   # sanity check, e.g. ~/.julia/juliaup/julia-1.11.x+.../bin/julia
>
> export JEXPRESSO_COUPLED=1
> mpirun -np 2 ./AlyaProxy/Alya.x \
>      : -np 2 "$JULIA_BIN" --project=. ./src/Jexpresso.jl CompEuler 3dAlya
> ```

### The MPICH bundled with Julia (`MPICH_jll`)

With MPI.jl on `MPICH_jll` and `Alya.x` built against it (Section 4),
`run_coupled.sh` needs nothing more: it launches with MPI.jl's own `mpiexec`.
By hand, the launcher comes from MPI.jl, and the env var is passed with MPICH's
`-env`:

```bash
julia --project=. -e '
  using MPI
  run(`$(mpiexec()) -n 2 ./AlyaProxy/Alya.x : -n 2 -env JEXPRESSO_COUPLED 1 $(Base.julia_cmd()) --project=. src/Jexpresso.jl CompEuler 3dAlya`)'
```

### Convenience script

[`run_coupled.sh`](run_coupled.sh) wraps the launch above. Its Julia ranks
load Jexpresso as a package (`using Jexpresso; Jexpresso.run_case(...)`), and
it brings the package image up to date before it launches, so that the ranks
never start compiling it all at once inside the job:

```bash
./run_coupled.sh                        # 2 Alya ranks + 2 Julia ranks (default)
./run_coupled.sh 2 4                    # 2 Alya ranks + 4 Julia ranks
PRECOMPILE_COUPLED=1 ./run_coupled.sh   # bake this case into the package image first
SCRIPT=1 ./run_coupled.sh               # run src/Jexpresso.jl as a script instead
```

- `PRECOMPILE_COUPLED=1` is explained in
  [Precompiling a coupled case](#precompiling-a-coupled-case).
- `SCRIPT=1` compiles everything at every launch. It suits an edit-and-run
  session, since nothing has to be precompiled after each edit.
- `REBUILD_SYSIMAGE=1` builds a sysimage of the dependencies only (see the same
  section).

It also passes `--startup-file=no` and per-rank output tagging (`-prepend-rank`
on MPICH, `--output tag` on OpenMPI 5, `--tag-output` on OpenMPI 4; disable
with `TAG_OUTPUT=0`). Without
tagging, a parallel run's stdout is block-buffered and interleaved, which is
the main reason a working run can look like a hang.

For the full performance recipe on macOS — sysimage, rank counts against
performance cores, thread settings, and what to avoid — see
[INSTALL.md §6.1](INSTALL.md#61-fastest-parallel-execution-on-macos-with-mpich).

### What a coupled run actually costs

Do not calibrate your expectations on `Alya.x` run alone. With no Jexpresso
ranks in the world it receives nothing and posts no receives, so its time loop
is empty and finishes in about a second. Coupled, the wall-clock time is
Jexpresso's: a full SEM right-hand side plus one coupling exchange **per
timestep**, for `(tend - tinit) / Δt` steps — 2000 of them with the `3dAlya`
defaults.

Alya receives exactly its own number of steps, and nothing in the protocol
checks that the two agree
([COUPLING-ALGORITHM.md §4.3](COUPLING-ALGORITHM.md#43-the-step-contract)).
Check it in the log. Alya prints `Steps: 2000` at startup, and Jexpresso prints
this at the end of its solve:

```
 # Coupling: 2000 exchanges sent to Alya
```

The two numbers must be equal.

**Telling "slow" from "deadlocked".** The `t=` lines only appear at
`:diagnostics_at_times`, which for `3dAlya` is every 100 time units — 200
timesteps apart. A run that is merely slow looks exactly like one that is stuck.
Turn on the per-step heartbeat, which prints the first five steps and then every
hundredth:

```bash
export JEXPRESSO_STEP_HEARTBEAT=1
./run_coupled.sh 2 2
```

If the `#   step N   t = ...` lines keep coming, it is running and you are
measuring throughput. If they stop, it is genuinely blocked — and the step
number tells you where.

Each exchange interpolates the solution from the Jexpresso SEM mesh to every
Alya grid point. Locating those points — bin lookup, bounding-box test, Newton
solve for the reference coordinates — depends only on the two geometries, so on
a static mesh Jexpresso now does it **once** and reuses the result, leaving one
dot product per point per equation each step. You will see this line on the
first exchange, with totals over all Jexpresso ranks:

```
[coupling] interpolation cache built: 1000/1000 points located in elements ...
```

The two numbers should be equal. A point that is not located takes the value of
its nearest mesh node instead of the interpolant.

The cache is disabled automatically when the mesh can change under it
(`:ladapt` or `:lamr`), and can be turned off with
`:lcouple_cache_interp => false` in `user_inputs.jl`.

If it is still slower than you want, the exchange frequency is the next lever —
but note that it currently happens on **every** step regardless of the
`:Δt_couple` entry in the example `user_inputs.jl`, which nothing reads yet.

For a shorter run while experimenting, change the number of steps on **both**
sides: `:Δt`/`:tend` in the case's `user_inputs.jl` **and** `dt`/`tend` in
`AlyaProxy/alya_all2all_time_loop.f90`, then rebuild `Alya.x`. Changing only
Jexpresso's side breaks the one-exchange-per-Alya-step contract: with fewer
steps Alya waits forever for the rest, and with more the extra messages have no
receive.

### Where the time goes, and running from the REPL

These timings are for `3dAlya` (2000 steps) with 2 Jexpresso ranks, run back to
back on one 4-core Linux machine. Your absolute times will differ; the
proportions are what carry over.

**The coupling itself costs about 5 %.** The same launch without Alya took
421 s:

```bash
mpiexec -n 2 julia --project=. src/Jexpresso.jl CompEuler 3dAlya   # JEXPRESSO_COUPLED unset
```

The coupled launch, in the same form, took 442 s. The coupling adds one
interpolation to Alya's
1000 points and one small message per step. Alya's ranks also busy-wait
inside MPI, so each keeps a core at 100 % for the whole run. That is harmless
with a core per rank, and expensive if the ranks outnumber your cores.

**What makes every launch slow is compilation.** Each rank is a new Julia
process, and it compiles Jexpresso before its first step. A coupled launch
(ranks started as in the REPL snippet below, caches in place) spent 60 % of its
time before the first step:

| phase | ends at |
|---|---|
| Julia started, packages loaded, `MPI_Init` | 26 s |
| driver compiled (`Entering with_mpi block … DONE`) | 112 s |
| mesh and SEM set up from `.jexpresso_cache` | 151 s |
| parameters, `Precompile warm-up`, integrator warm-up | 239 s |
| 2000 time steps, about 82 ms each | 403 s |

In a REPL, a second `Jexpresso.run_case` in the same session reuses all of that
compiled code, which is why a standalone case feels fast there. A coupled job
cannot do the same. `mpiexec` starts the ranks, so the REPL process is never one
of them, and `Alya.x` runs once per job.

**Launching from the REPL** still works, but it starts new processes:

```julia
julia> using MPI    # the REPL started as `julia --project=.` in the repository root
julia> ranks = `$(Base.julia_cmd()) --project=. --startup-file=no -e 'using Jexpresso; Jexpresso.run_case("CompEuler", "3dAlya")'`
julia> withenv("JEXPRESSO_COUPLED" => "1") do
           run(`$(mpiexec()) -prepend-rank -n 2 ./AlyaProxy/Alya.x : -n 2 $ranks`)
       end
```

- `-prepend-rank` is the MPICH flag. On Open MPI, use the tagging flag listed
  under the convenience script above.
- `withenv` keeps `JEXPRESSO_COUPLED` out of the REPL's own environment, so a
  `run_case` typed there afterwards still runs standalone.

These ranks load the precompiled package (`using Jexpresso`) instead of
compiling the script `src/Jexpresso.jl`. On its own, that only gets each rank
going about 25 s sooner (26 s instead of 51 s). The compilation itself goes
away once the case is precompiled for coupled runs, as described next.

### Precompiling a coupled case

Jexpresso can compile a coupled case ahead of time, into its package image, so
that a launch starts computing almost at once. Precompilation runs in a single
process with no Alya, so it runs the case *dry*. Every coupling function runs as
in a real coupled run, which gives the image the code and types it has to hold.
Alya is simulated in the same process: a virtual world with Alya at ranks 0 and
1, Alya's grid made of a few points inside Jexpresso's mesh, and no messages.
The case runs twice, for three steps each time. The first run builds the
case's mesh and SEM caches if they do not exist yet; the second loads them,
which is the path every later launch takes.

Turn it on once per checkout:

```bash
PRECOMPILE_COUPLED=1 ./run_coupled.sh 2 2
```

From Julia, `Jexpresso.set_coupled_precompile!("CompEuler", "3dAlya")` does the
same, and `Jexpresso.set_coupled_precompile!(nothing)` turns it off. The choice
is kept in `LocalPreferences.toml`, so every later precompile bakes in the same
case, wherever it is triggered: `Pkg.precompile()`, a `using Jexpresso` in the
REPL, or `run_coupled.sh`. It only helps launches whose ranks load the package:
`run_coupled.sh` and the REPL recipe above do, while `SCRIPT=1` and
`julia src/Jexpresso.jl` do not. A coupled launch that runs without it prints
`Not precompiled for coupled runs`.

On `3dAlya` with 2 + 2 ranks, on the same machine as above:

| | case not precompiled | case precompiled |
|---|---|---|
| first time step starts at | 239 s | 73 s |
| whole launch | 402 s | 232 s |

The results are unchanged: Alya's output files were byte-identical.

What it costs:

- Precompiling Jexpresso takes several minutes instead of about one (6.5
  minutes on the machine above). It happens again whenever Jexpresso's source
  changes, before the next launch or at the next `using Jexpresso`.
- The package image grows from 17 MB to 150–180 MB, and every
  `using Jexpresso`, including those of standalone runs, loads it about 8 s
  more slowly.
- Editing the case's own files does not make the image stale. The edited file
  is reloaded at launch instead. For `user_inputs.jl` that costs nothing, but
  an edited hook that the right-hand side calls (`user_flux.jl`,
  `user_source.jl`, …) is recompiled at every launch until you rebuild with
  `PRECOMPILE_COUPLED=1 ./run_coupled.sh`. That flag always rebuilds.

Jexpresso's structs used to carry their array sizes in their types, for
example `St_metrics{Float64, (50, 5, 5, 5), …}`, so every partition size was a
different type. A rank of a 2-rank launch holds 50 of the 100 elements, so it
could not use the code the single-process precompile had compiled for all 100,
and its warm-ups compiled for about 30 s. The sizes are now ordinary fields (see
the top of `src/kernel/globalStructs.jl`), and a launch compiles for about a
second. Measured back to back on another, slower machine, with the case
precompiled:

- The first time step starts at 25 s instead of 52 s, about 3 s after Julia has
  loaded.
- The whole launch takes 123 s instead of 154 s.
- The time per step is unchanged (48–50 ms), and so is Alya's output.

What remains is loading the package image, about 22 s.

The sysimage (`REBUILD_SYSIMAGE=1`, `create_Jexpresso_sysimage.jl`) holds only
the dependencies. Its trace script, `precompile_jexpresso.jl`, includes
`src/Jexpresso.jl`, which starts a case only when it is the program being run.
The sysimage makes packages load faster, but it contains none of Jexpresso's
compiled code.

To compare with a standalone run in the REPL:

```julia
julia> using Jexpresso
julia> @time Jexpresso.run_case("CompEuler", "3dAlya")   # compiles, like every coupled launch
julia> @time Jexpresso.run_case("CompEuler", "3dAlya")   # reuses the compiled code
```

Both calls run on a single process, so each time step takes longer than it does
with two ranks.

---

## 7. Troubleshooting

### "It hangs forever and never starts" — find out where, first

Do **not** start by re-reading the coupling code. Run the automated ladder,
which times out each stage instead of freezing:

```bash
bash tools/check_mpi_setup.sh
```

It answers, in order: what is on `PATH`; what MPI.jl is bound to; whether the
launcher and MPI.jl are the same MPI; whether `MPI_Init` completes with 1 rank;
whether it completes with 2 ranks; whether `Alya.x` links the same MPI; and
whether `JEXPRESSO_COUPLED` actually reaches the Julia ranks.

If you prefer to do it by hand, the decisive test is this six-line program —
it contains no Jexpresso code at all:

```bash
JULIA_BIN=$(julia -e 'print(joinpath(Sys.BINDIR, "julia"))')
mpirun -np 2 "$JULIA_BIN" --project=. -e '
  using MPI; MPI.Init()
  println("rank ", MPI.Comm_rank(MPI.COMM_WORLD), " of ", MPI.Comm_size(MPI.COMM_WORLD))
  flush(stdout); MPI.Finalize()'
```

Read the result as follows:

| What you see | What it means |
|---|---|
| `rank 0 of 2` and `rank 1 of 2` | MPI is healthy. The hang is above this layer — see the next two bullets. |
| nothing, forever | `MPI_Init` never completes. Launcher and MPI.jl are different MPIs (Sections 1–3), or on macOS + MPICH the hostname does not resolve ([INSTALL.md §5.5](INSTALL.md#55-macos-hostname-fix-mpich-and-mpich_jll-only)). |
| `rank 0 of 1` printed twice | Each rank came up as its own world — same mismatch, non-blocking flavour. Jexpresso will then hang on its first collective. |

**Narrow it further by removing the coupling from the picture.** If

```bash
mpirun -np 2 "$JULIA_BIN" --project=. ./src/Jexpresso.jl CompEuler 3dAlya
```

also hangs (with `JEXPRESSO_COUPLED` unset), the problem is *not* in the
coupling — it is in the MPI setup or in Julia startup, and the two sections
below apply.

### Two things that look like a hang but are not

- **First-time precompilation.** A cold project takes many minutes to
  precompile, and launching *N* ranks at once makes them serialise on the
  precompile lock. Do it once, serially, before any `mpirun`:
  ```bash
  julia --project=. -e 'using Pkg; Pkg.precompile()'
  ```
- **Buffered output.** Under `mpirun`, stdout is a pipe rather than a
  terminal, so Julia block-buffers it: the run can be making progress with
  nothing on screen. Tag and flush per rank —
  OpenMPI: `mpirun --output tag …` (OpenMPI 4: `--tag-output`);
  MPICH: `mpiexec -prepend-rank …`.

### The run starts, steps, then stops at the first diagnostic output

Symptom: the `t=` line for the first `:diagnostics_at_times` entry prints, and
nothing follows — no CFL lines, no `.pvtu` write, no further steps. Alya keeps
whatever it had and waits.

Cause: some Jexpresso code on the diagnostic path ran an **MPI collective on
`MPI.COMM_WORLD`**. In a coupled run `COMM_WORLD` also contains Alya's ranks,
which never execute Jexpresso's diagnostics, so every Jexpresso rank blocks in
that collective forever.

> **Rule for anyone adding MPI calls to Jexpresso:** never use
> `MPI.COMM_WORLD` for a Jexpresso-internal collective. Use `get_mpi_comm()`,
> which is the Jexpresso-only communicator under coupling and is `COMM_WORLD`
> in standalone runs, so it is always correct. `MPI.COMM_WORLD` is right only
> for the coupling handshake itself, where Alya participates too — those call
> sites live in `src/kernel/coupling/couplingStructs.jl` and take `world`
> explicitly.

This is easy to introduce accidentally: an `Allreduce` added to a routine that
used to compute a local maximum will pass every standalone test and deadlock
only when coupled.

### Everything else

- **Job hangs in `MPI_Init`, or one side never reaches the handshake.** The
  classic symptom of mismatched MPIs. Redo Sections 1–3 and confirm the table
  matches; recompile `Alya.x` (Section 5) if you changed the Julia binding.
- **The job starts, Alya prints its coupling labels, then everything stops.**
  `JEXPRESSO_COUPLED` did not reach the Julia ranks, so Jexpresso is running
  standalone on a `COMM_WORLD` that contains Alya. `export` it (Section 6)
  rather than passing `-x` after the `:`.
- **Jexpresso finishes its solve but the job never exits**, and Alya never
  prints `Alya: time loop complete`. The two time loops disagree on the number
  of steps: compare `# Coupling: N exchanges sent to Alya` with Alya's
  `Steps:` line, and make `:Δt`/`:tend` match the proxy's `dt`/`tend`
  (Section 6, "What a coupled run actually costs").
- **`prterun` / launcher cannot find Julia.** Use the explicit `$JULIA_BIN`
  form shown in Section 6.
- **Stale or conflicting MPI binding on the Julia side.** Remove the recorded
  preference and reconfigure from Section 3:
  ```bash
  rm -f LocalPreferences.toml
  ```
- **macOS + MPICH `gethostbyname` failure at `MPI_Init`.** Apply the
  [hostname `/etc/hosts` fix](INSTALL.md#55-macos-hostname-fix-mpich-and-mpich_jll-only).
- **Wrong `mpif90`/`mpiexec` picked up.** Confirm `which mpif90`, `which mpiexec`,
  and `which mpirun` all point at the same installation, or use absolute paths.

See also the [FAQ.md](FAQ.md) for common run and installation errors.

For the full description of what the two codes exchange once they are running,
see [COUPLING-ALGORITHM.md](COUPLING-ALGORITHM.md).
