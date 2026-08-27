#=============================================================================
 precond_api.jl -- the extension point for the IMEX3D stage-solve
                   preconditioner.

 WHAT A PRECONDITIONER HAS TO BE HERE
 ------------------------------------
 `gmres_solve!` (krylov.jl) asks for exactly one thing: a callable `precon!(Z)`
 that overwrites the `npoin x nimp` field `Z` with `M^-1 Z`. It never inspects
 the preconditioner, never asks for a matrix, and never needs `M` to be linear
 in any way it can check. Everything else in this scheme -- the operator, the
 column topology, the distributed inner product, the ARK driver -- is already
 blind to which `M` that is: `IMEX3DCache` is parametric in the preconditioner
 type (`PCT`), so a new one changes no field type and no method signature
 anywhere upstream.

 So the machinery was always open. What was NOT open was the DECK: `:imex_precond`
 accepted `:column` or `:none` and raised on anything else, which meant a new
 preconditioner could not be reached without editing `build_imex3d` and
 `build_imex3d_schur`. This file is the third value, `:custom`, and the three
 generic functions a custom object answers.

 THE CONTRACT, IN FULL
 ---------------------
 A user supplies a BUILDER through `:imex_precond_build`, called once at setup
 with the context below, and adds methods to `imex_precond_apply!` (required)
 and `imex_precond_refresh!` / `imex_precond_describe` (both optional) for
 whatever type the builder returns. Nothing else. `user_inputs.jl` is
 `include`d into the `Jexpresso` module, so those methods can be written
 straight into the case directory with no prefix and no `import`.

 IN PLACE MEANS IN PLACE, AND THAT IS THE ONE MISTAKE WORTH CATCHING
 -------------------------------------------------------------------
 `precon!(Z)` is called for its EFFECT ON `Z`; the return value is discarded by
 the Arnoldi loop. An `imex_precond_apply!` that computes `M^-1 Z` correctly
 and returns it as a fresh array therefore acts as the IDENTITY -- the run
 still converges, to the right answer, at the iteration count of no
 preconditioner at all, and the only symptom is that the scheme is
 inexplicably slow. `imex_precond_selfcheck` below turns that into an error at
 setup instead, which is why `:imex_verify` is worth leaving on while a new
 preconditioner is being written.

 WHAT MAKES A PRECONDITIONER LEGAL IN PARALLEL
 ---------------------------------------------
 Two rules, both inherited from the iteration rather than from this file:

   * it must leave the field SINGLE-VALUED across ranks -- every holder of a
     shared node agreeing -- because the multiplicity-weighted inner product in
     krylov.jl is an inner product only on such fields. The built-in column
     preconditioner ends in `column_scatter!` for this reason; a node-local
     preconditioner (Jacobi, say) gets it for free by acting identically on
     every copy of a node.

   * every collective inside it must be reached by every rank on every
     application. The iteration itself is bit-identical across ranks (see
     `imex3d_solve!`), so this holds automatically unless the preconditioner
     branches on a rank-local quantity.

 A preconditioner that CHANGES between applications (an inner Krylov solve with
 its own tolerance, say) is flexible-GMRES territory and this iteration is not
 flexible: it applies `M^-1` again at the end of each cycle to form the update
 (`krylov.jl`, `precon!(C)`), which assumes the same `M^-1` it used inside the
 cycle. Keep the application deterministic, or set `:imex_restart` to
 `:imex_maxiter` so there is only ever one cycle.
=============================================================================#

"""
    imex_precond_apply!(pc, V, params, gdt) -> V

Apply `M^-1` **in place** to `V`, an `npoin x nimp` field: `nimp = 5` on the
five-field stage solve (`ρ, ρu, ρv, ρw, ρθ`, in `op.vars` order) and `nimp = 1`
on the scalar Schur system (`P = β·Θ`).

`gdt` is `γΔt` for the stage being solved. It is not constant over a run --
Jexpresso passes every diagnostic time as a `tstop` and the integrator shortens
the step that lands on one -- so a preconditioner that depends on it must
notice when it moves. Both built-ins compare against the value they were built
at and refactorise; see `imex3d_precond!`.

REQUIRED for a custom preconditioner. The default method exists only to name
the omission.

IT IS DECLARED WITH NO ARGUMENT TYPES AT ALL, deliberately. A fallback that
constrained `V::AbstractMatrix` or `gdt::Real` would be AMBIGUOUS with the
obvious user method `imex_precond_apply!(pc::MyPrecond, V, params, gdt)` --
more specific in its first argument, less specific in the others, so Julia can
order neither and raises a method ambiguity from inside the Arnoldi loop.
Untyped, ANY method a user writes is strictly more specific and simply wins.
"""
function imex_precond_apply! end

imex_precond_apply!(pc, V, params, gdt) =
    error("IMEX3D: no `imex_precond_apply!` method for a preconditioner of type ",
          typeof(pc), ". `:imex_precond_build` returned that object, and the ",
          "stage solve has no way to apply it. Define\n\n",
          "    Jexpresso.imex_precond_apply!(pc::", typeof(pc),
          ", V::AbstractMatrix, params, gdt) = ...\n\n",
          "which must overwrite V with M^-1 V IN PLACE and return it. See ",
          "src/kernel/solvers/hevi/precond_api.jl.")

"""
    imex_precond_refresh!(pc, params, gdt) -> pc

Rebuild whatever the preconditioner holds against the operator's current
coefficients. Called from `ark_relinearize!` under `:imex_linearization => :PS`,
after `β`, `θ̄`, `grad(θ̄)` and the eddy viscosity have been refreshed on the
3D operator and before the next step.

The default is a NO-OP, which is right for a preconditioner that reads the
operator's coefficient arrays live and wrong for one that copied them at setup.
Under `:RS` (the default linearisation) this is never called at all.

`gdt` here is the run's NOMINAL `γΔt`, not necessarily the value the
preconditioner was last applied at -- the integrator shortens the step that
lands on a `tstop`, and this hook fires between steps rather than inside a
stage. A preconditioner that tracks its own `gdt` should rebuild at the value
it is holding and ignore this one; the built-ins do exactly that
(`refactorize!(..., fac.gdt[])`). It is passed because a preconditioner that
tracks nothing has no other way to know what to build against.

Refreshing the operator and not the preconditioner is a SILENT cost: the band
still converges, just to a worse iteration count. See `ark_relinearize!`.

Untyped for the reason on `imex_precond_apply!`: a constrained fallback would
be ambiguous with the natural user method rather than beaten by it.
"""
imex_precond_refresh!(pc, params, gdt) = pc

"""
    imex_precond_describe(pc) -> String

One line for the setup report, after `preconditioner:`. Say what `M` is and
what it cost to build; the report is where a user finds out that the thing they
selected is the thing that ran.
"""
imex_precond_describe(pc) = string(nameof(typeof(pc)), " (:imex_precond => :custom)")

"""
    imex_precond_context(; params, inputs, op, topo, comm, gdt, nimp, schur,
                           lwall_flux) -> NamedTuple

The single argument handed to a `:imex_precond_build` builder.

| field | what it is |
|---|---|
| `params` | the full parameter bundle: `mesh`, `metrics`, `Minv`, `qp.qe`, … |
| `inputs` | the deck `Dict`, so a builder can read its own keys |
| `op` | the 3D acoustic operator `A`. `hevi_apply_A!(W, V, params, op)` applies it; `op.beta`, `op.thetabar`, `op.vars`, `op.nimp` are readable |
| `topo` | the `ColumnTopology`: who owns which vertical column and how many levels it has |
| `comm` | the MPI communicator the run is on |
| `gdt` | `γΔt` at setup — the value to build against |
| `nimp` | width of the field this will be applied to: 5, or 1 on the Schur path |
| `schur` | `true` when this preconditions the scalar Schur operator `H` rather than `I - γΔt·A` |
| `lwall_flux` | whether the implicit vertical mass flux is zeroed at floor and lid |

FIELDS MAY BE ADDED HERE. Destructure by name (`ctx.topo`), never by position.
"""
imex_precond_context(; params, inputs, op, topo, comm, gdt, nimp, schur,
                       lwall_flux) =
    (params = params, inputs = inputs, op = op, topo = topo, comm = comm,
     gdt = Float64(gdt), nimp = Int(nimp), schur = Bool(schur),
     lwall_flux = Bool(lwall_flux))

"""
    imex_precond_selfcheck(pc, params, npoin, nimp, gdt, comm)

Apply the preconditioner once to a deterministic field and refuse the two
failures that are otherwise silent for a whole run.

`NOT IN PLACE` -- the field came back bit-identical, so the Arnoldi loop would
be preconditioning with the identity. See the header.

`NOT FINITE` -- a NaN or Inf out of `M^-1` propagates into the Krylov basis and
the run dies several stages later inside the operator, where nothing points
back here.

The probe is built from the GLOBAL point id rather than from `rand()`, so every
holder of a shared node feeds in the same value and a parallel preconditioner
is exercised on a field it is allowed to see.
"""
function imex_precond_selfcheck(pc, params, npoin::Int, nimp::Int, gdt::Real,
                                comm::MPI.Comm)

    gip = params.mesh.ip2gip
    Z = zeros(Float64, npoin, nimp)
    @inbounds for q = 1:nimp, ip = 1:npoin
        Z[ip, q] = sinpi(1.0e-3 * (Float64(gip[ip]) + 29.0 * q))
    end
    Z0 = copy(Z)

    imex_precond_apply!(pc, Z, params, Float64(gdt))

    nbad = count(!isfinite, Z)
    nbad = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(nbad, MPI.SUM, comm) : nbad
    nbad == 0 ||
        error("IMEX3D: the custom preconditioner (", typeof(pc), ") produced ",
              nbad, " non-finite values. A NaN here enters the Krylov basis and ",
              "the run fails several stages later inside the operator, a long ",
              "way from the cause.")

    same = Z == Z0
    same = MPI.Comm_size(comm) > 1 ?
           (MPI.Allreduce(Int(same), MPI.MIN, comm) == 1) : same
    same &&
        error("IMEX3D: the custom preconditioner (", typeof(pc), ") left the ",
              "field bit-identical. `imex_precond_apply!` is called for its ",
              "EFFECT ON V and its return value is discarded, so a method that ",
              "computes M^-1 V into a NEW array acts as the identity: the run ",
              "converges to the right answer at the iteration count of no ",
              "preconditioner at all, and nothing else says so. Write into V. ",
              "If the identity really is what you meant, ask for it by name: ",
              ":imex_precond => :none.")

    return nothing
end

"""
    build_custom_precond(ctx) -> pc

Resolve `:imex_precond_build`, call it once, and prove the object it returned
can be applied.

The builder is anything callable with the context: a named function defined in
`user_inputs.jl`, a closure, a callable struct. It runs on EVERY rank -- it is
free to use collectives, and must not branch on a rank-local quantity in a way
that changes how many it reaches.

`imex_precond_selfcheck` runs here UNCONDITIONALLY, not under `:imex_verify`:
it costs exactly one preconditioner application, which is by definition one
Krylov iteration's worth of work, and the two failures it catches are the two
that otherwise cost a whole run before anyone notices.
"""
function build_custom_precond(ctx)

    builder = get(ctx.inputs, :imex_precond_build, nothing)
    builder === nothing &&
        error("IMEX3D: :imex_precond => :custom needs :imex_precond_build => f, ",
              "a callable that takes the build context and returns the ",
              "preconditioner object. See src/kernel/solvers/hevi/precond_api.jl ",
              "and the worked example in README_IMEX3D.md.")
    applicable(builder, ctx) ||
        error("IMEX3D: :imex_precond_build is a ", typeof(builder),
              ", which is not callable with the single build-context argument. ",
              "The signature is `f(ctx) -> pc`; `ctx` is the NamedTuple ",
              "documented on `imex_precond_context`.")

    pc = builder(ctx)
    pc === nothing &&
        error("IMEX3D: :imex_precond_build returned `nothing`. Return the ",
              "preconditioner object, or select :imex_precond => :none if the ",
              "identity is what is wanted.")

    imex_precond_selfcheck(pc, ctx.params, Int(ctx.params.mesh.npoin),
                           ctx.nimp, ctx.gdt, ctx.comm)
    return pc
end
