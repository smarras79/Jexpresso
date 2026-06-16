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

## 3. Make the Julia side match the Fortran side

Point `MPI.jl` at the **same** system MPI that `mpif90` uses, using the library
path from Section 1, then rebuild and precompile:

```bash
# Use the SAME lib dir that `mpif90 -show` reported (here: Homebrew on macOS).
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=["/opt/homebrew/lib"])'
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

> If the system MPI is already in a standard location (`/usr/bin`,
> `/usr/local/bin`), `MPIPreferences.use_system_binary()` with no `extra_paths`
> is enough. The detailed configuration options are in
> [INSTALL.md, Section 5.2](INSTALL.md#52-point-mpipreferences-at-your-mpi).

**Confirm the match** by running the Section 2 check again and comparing it
side-by-side with Section 1. Both must agree on all three rows:

| Must match | Fortran (Alya) | Julia (Jexpresso) |
|------------|----------------|-------------------|
| Implementation | `mpif90 --version` | `MPI.identify_implementation()` |
| Version | `mpiexec --version` | `MPI.identify_implementation()` |
| Library / launcher | `mpif90 -show` → `-L.../lib` | `MPIPreferences.binary == "system"` + same `mpiexec` |

Only proceed once these line up.

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

With the correct `mpif90` first on your `PATH` (verify with `which mpif90`),
build the proxy:

```bash
cd AlyaProxy
bash compilef90.sh        # runs: mpif90 -cpp -DUSEMPIF08 myAlya.f90 -o Alya.x
cd ..
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

For the full description of what the two codes exchange once they are running,
see [COUPLING-ALGORITHM.md](COUPLING-ALGORITHM.md).
