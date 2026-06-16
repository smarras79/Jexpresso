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
- `binary = "MPItrampoline_jll"` (or another `*_jll`) means MPI.jl uses the
  **bundled native MPI** — this is **not** what you want for coupling unless you
  also build Alya against that exact JLL (see the note in Section 4).

**This must report the same implementation and version you wrote down in
Section 1.** If it doesn't, fix it in the next step.

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

If that MPI provides no `mpif90` at all (for instance MPI.jl is on the bundled
`MPItrampoline_jll`), you cannot easily build Alya against it — install a matching
**system** MPI and use Option A instead. See the
[note on the JLL route in Section 4](#4-choose-a-consistent-mpi-for-coupling).

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

> **About the `MPICH_jll` *native* route.** The bundled-JLL route in
> [INSTALL.md Section 5, Route C](INSTALL.md#5-running-in-parallel-with-mpi) is
> for **standalone** Jexpresso. It is **not** appropriate for coupling in the
> general case, because Alya is a separately compiled Fortran binary linked
> against a *system* MPI — and that system MPI is almost never ABI-identical to
> the JLL. Use a **system** OpenMPI or MPICH for coupled runs so both codes link
> the very same library. (The only exception: if you deliberately build `Alya.x`
> against the *same* MPICH that is ABI-compatible with `MPItrampoline_jll`, you
> may launch with MPI.jl's bundled `mpiexec` as shown in Section 6 — advanced,
> macOS-only.)

---

## 5. Compile the Alya proxy with the chosen MPI

`Alya.x` is built from `myAlya.f90` (a symlink to `alya_all2all_time_loop.f90`)
by [`AlyaProxy/compilef90.sh`](AlyaProxy/compilef90.sh), which calls `mpif90`.
The only thing that matters is that this `mpif90` belongs to the **same MPI**
you bound MPI.jl to in Section 3.

**First make sure the right wrapper is selected.** `compilef90.sh` uses the
`mpif90` first on your `PATH` by default, but you can pin a specific one with the
`MPIF90` variable (recommended when several MPIs are installed):

```bash
cd AlyaProxy

# Option A — the mpif90 on PATH is already the right one (verify, then build):
which mpif90 && mpif90 -show
bash compilef90.sh

# Option B — pin the exact wrapper that matches MPI.jl:
MPIF90=/opt/homebrew/Cellar/open-mpi/5.0.8/bin/mpif90 bash compilef90.sh

cd ..
```

The script prints the wrapper's `-show` line so you can eyeball that the
include/lib paths point at the same MPI as the Julia side.

**Worked example (the verified Homebrew OpenMPI 5.0.8 setup).** If
`MPI.identify_implementation()` reports `("OpenMPI", v"5.0.8")` and `mpif90
-show` reports `-L/opt/homebrew/Cellar/open-mpi/5.0.8/lib`, the two already
match — just build:

```bash
cd AlyaProxy && bash compilef90.sh && cd ..
```

**Verify what `Alya.x` actually linked against** (catches a wrong wrapper even
when `-show` looked fine):

```bash
# macOS:
otool -L AlyaProxy/Alya.x | grep -i mpi
#   → .../open-mpi/5.0.8/lib/libmpi.40.dylib   (must be the SAME MPI as MPI.jl)

# Linux:
ldd AlyaProxy/Alya.x | grep -i mpi
```

This produces `./AlyaProxy/Alya.x`, linked against the MPI you selected in
Section 4. If you later switch the Julia side to a different MPI, you **must
recompile `Alya.x`** with the matching `mpif90` — the two always travel together.

---

## 6. Launch the coupled run (MPMD)

A coupled run is a single MPMD launch: one launcher starts Alya ranks **and**
Julia ranks, separated by `:`, so they share `MPI_COMM_WORLD`. Use the launcher
that belongs to the MPI both codes were built against.

The `JEXPRESSO_COUPLED=1` environment variable tells Jexpresso it is running
coupled. The example case is `CompEuler thetaAlya` (2D) or `CompEuler 3dAlya` (3D).

### OpenMPI or MPICH (system MPI) — `mpirun`

```bash
# 2D: 2 Alya ranks + 2 Jexpresso ranks
mpirun -np 2 ./AlyaProxy/Alya.x : -np 2 -x JEXPRESSO_COUPLED=1 \
    julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya
```

```bash
# 3D
mpirun -np 2 ./AlyaProxy/Alya.x : -np 2 -x JEXPRESSO_COUPLED=1 \
    julia --project=. ./src/Jexpresso.jl CompEuler 3dAlya
```

> **`prterun` error (OpenMPI 5)?** OpenMPI's launcher sometimes cannot resolve
> the `juliaup` shim. Pass the real `julia` binary explicitly:
>
> ```bash
> JULIA_BIN=$(julia -e 'print(joinpath(Sys.BINDIR, "julia"))')
> echo "$JULIA_BIN"   # sanity check, e.g. ~/.julia/juliaup/julia-1.11.x+.../bin/julia
>
> mpirun -np 2 ./AlyaProxy/Alya.x \
>      : -np 2 -x JEXPRESSO_COUPLED=1 "$JULIA_BIN" --project=. ./src/Jexpresso.jl CompEuler 3dAlya
> ```

### MPICH_jll bundled launcher (advanced, macOS)

Only valid if `Alya.x` was built against an MPICH that is ABI-compatible with
`MPItrampoline_jll` (see the note in Section 4). Here the launcher comes from
MPI.jl, and the env var is passed with MPICH's `-env`:

```bash
julia --project=. -e '
  using MPI
  run(`$(mpiexec()) -n 2 ./AlyaProxy/Alya.x : -n 2 -env JEXPRESSO_COUPLED 1 $(Base.julia_cmd()) --project=. src/Jexpresso.jl CompEuler 3dAlya`)'
```

### Convenience script

[`run_coupled.sh`](run_coupled.sh) wraps the system-MPI launch above, and also
builds a PackageCompiler sysimage to remove Julia's cold-start JIT cost:

```bash
./run_coupled.sh               # 1 Alya rank + 2 Julia ranks (default)
./run_coupled.sh 1 4           # 1 Alya rank + 4 Julia ranks
REBUILD_SYSIMAGE=1 ./run_coupled.sh   # force sysimage rebuild first
```

---

## 7. Troubleshooting

- **Job hangs in `MPI_Init`, or one side never reaches the handshake.** The
  classic symptom of mismatched MPIs. Redo Sections 1–3 and confirm the table
  matches; recompile `Alya.x` (Section 5) if you changed the Julia binding.
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
