# Radiative Transfer — Planned Features

Work items for PRs following the `yt/Radiation_merged` branch merge.
Items marked **high-priority** should be tackled first.

---

## 1. BOMEX and Lasher-Trapp 2001 test cases  ⭐ high-priority

Add two dynamics test cases that drive the RT solver:

- **BOMEX** — shallow cumulus benchmark (Siebesma et al. 2003).
- **Lasher-Trapp 2001** — mixed-phase deep convection case.

Each case should live under `problems/` as a self-contained problem directory and call `build_radiative_transfer_problem` through the standard driver interface.

---

## 2. Shortwave and longwave reflection verification  ⭐ high-priority

Confirm that boundary reflection mechanics for both SW and LW are physically correct:

- Specular / Lambertian surface reflection.
- Top-of-atmosphere incoming SW.
- Comparison against analytic or reference solutions.

---

## 3. Warm-start persistence across adaptive mesh steps

Currently warm starts are only valid when the mesh does not change. When the angular or spatial mesh adapts, the solution vector changes size and layout, making the stored warm start stale. Implement interpolation of the previous solution onto the new mesh so that warm starts remain usable across refinement steps.

---

## 4. Preconditioning memory footprint reduction (MUMPS and global_ilu)

The global ILU and MUMPS preconditioners currently assemble the full global matrix on every rank, which is memory-intensive at scale. Investigate and implement a reduced-footprint strategy, e.g. distributed assembly with local solve or overlap-based partitioning.

---

## 5. Angular refinement interior hanging node performance

The current handling of interior hanging nodes for angular refinement is correct but slow. Profile and optimize the bottleneck — likely the constraint matrix assembly or the ghost communication phase for interior nodes.

---

## 6. p-refinement in angle

Add capability to increase the polynomial order per angular element (p-refinement) as an alternative or complement to h-refinement. This requires:

- Per-element nop tracking in the angular mesh (already partially in place via `extra_meshes_extra_nops`).
- Adaptation criterion selecting between h- and p-refinement.
- Updated constraint matrix construction for p-nonconforming interfaces.

---

## 7. MPI support for the 2D+1D radiation case (single extra dimension)

The 2D spatial + 1D angular (`build_rad_2d.jl`) solver currently treats MPI as an afterthought: the angular connectivity and constraint assembly in `adaptive_spatial_angular_numbering_2D_1D!` runs serially on each rank without ghost communication. Port the full MPI-aware pattern from the 3D solver:

- Ghost layer exchange for hanging angular nodes across partition boundaries.
- Parallel restriction/prolongation matrices for the non-conforming angular mesh.
- `Allreduce`-based global numbering analogous to `setup_global_numbering_extra_dim`.

---

## 8. Spatial AMR for the 2D radiation case + MPI

Extend `build_rad_2d.jl` to support spatially non-conforming (hanging) meshes in 2D, mirroring the `SpatialAMRCache` / `spatial_constraint_matrices.jl` / `spatial_ghost_comms.jl` infrastructure already in place for the 3D solver. MPI-correctness (ghost exchange, rank-aware constraint application) must be included from the start.

---

## 10. Multi-level angular refinement (2+ levels)

Currently the angular adaptivity framework supports a single level of h-refinement: a coarse parent element can be subdivided once into children. Support two or more successive refinement levels, where a child element can itself become a parent and be further subdivided.

Key work items:

- **Recursive refinement data structures**: `ref_level` and `neighbors` arrays currently encode a single parent–child relationship. Extend to a tree (or depth-tagged flat array) so that grandparent → parent → child chains are tracked.
- **Multi-level constraint matrices**: At a 2-level interface the fine-mesh node is constrained through two interpolation steps (grandparent basis → parent basis → child basis). Build the composed constraint matrix `R_2 = R_child * R_parent` and integrate it into the existing constraint assembly in `extra_amr_matrices.jl`.
- **Hanging node numbering**: `setup_global_numbering_adaptive_angular_scalable` (Geom.jl) currently assumes at most one level of refinement. Generalise the interior-node exclusion and ghost communication logic for multi-level chains.
- **Ghost communication**: The 10-phase ghost exchange in `angular_comms.jl` must propagate constraint data across potentially two hops (child → parent → grandparent), which may span rank boundaries. Extend the Isend/Irecv pattern accordingly.
- **Adaptivity criterion**: The existing criterion in `compute_adaptivity_criterion` / `compute_adaptivity_criterion3D_2D` fires per element independently. Add a max-depth cap (`max_ref_level` input key) and prevent refinement beyond it.
- **Coarsening (derefinement)**: Multi-level refinement without coarsening causes unbounded growth of the angular DOF count. Implement derefinement when the criterion falls below a lower threshold, collapsing a set of children back to their parent.

---

## 9. Adaptivity for 3D extra dimensions (Vlasov–Poisson / Maxwell / spectral radiation)

Generalize the extra-dimension adaptivity framework from 2D angular meshes (`extra_mesh` with θ, ϕ) to fully 3D extra-dimension meshes. Target use cases:

- **Vlasov–Poisson / Maxwell**: 3D velocity-space adaptivity (vx, vy, vz).
- **Spectral integration of radiation**: frequency as a third extra dimension, enabling spectrally-resolved RT without a gray-atmosphere assumption.

Key work items:
- Extend `extra_mesh` structures and `connijk_spa` to carry a third extra coordinate.
- Generalize hanging-node constraint construction (currently 2D Lagrange tensor products) to 3D.
- Update ghost communication phases and global numbering for the higher-dimensional case.
