#=============================================================================
 test_sgs_alloc.jl -- every allocate_SGS method must tolerate the keywords the
 call site passes, whether or not that model uses them.

 THE BUG THIS PINS. params_setup.jl has ONE call site for allocate_SGS and it
 passes the union of what any model might want. When the dynamic-SGS work added
 nelem/neqs/SD/C1/C2 there unconditionally, every deck using SMAG or VREM began
 failing at setup with "no method matching allocate_SGS", because:

   * `allocate_SGS(::Any, ::Any, ::Any, ::Any, ::Any; kwargs...)` looks like it
     absorbs anything, but it does NOT get the chance -- `::SMAG` is MORE
     SPECIFIC, so dispatch commits to the SMAG method and the keyword sorter
     throws there.
   * DSGS names those keywords, so it worked; AV has no method of its own, so
     it reached the fallback and worked. Only the two models with a specific
     method AND a narrow keyword list broke -- which is why it looked like a
     case-specific failure rather than a signature problem.

 So this does not test one signature: it calls EVERY tag with the full keyword
 set the call site can produce.

     julia --project=<env> test/physics/test_sgs_alloc.jl
=============================================================================#

using Test
using Jexpresso: allocate_SGS, PhysicalConst, AV, SMAG, VREM, DSGS, NSD_2D, NSD_3D,
                 AbstractSGSModel
using KernelAbstractions

const PC = PhysicalConst{Float64}()
const NP = 64

# Exactly what params_setup.jl passes on the branch that added them. Keep this
# in step with that call site: it is the point of the test.
callsite_kwargs(SD) = (C_s = 0.23, nelem = 8, neqs = 5, SD = SD,
                       C1 = 1.0, C2 = 0.5)

@testset "allocate_SGS tolerates the call site's keywords" begin
    for SD in (NSD_2D(), NSD_3D())
        for tag in (AV(), SMAG(), VREM(), DSGS())
            nm = string(typeof(tag))
            got = try
                allocate_SGS(NP, Float64, CPU(), PC, tag; callsite_kwargs(SD)...)
            catch e
                @error "allocate_SGS($nm) rejected the call site's keywords" exception=e
                :THREW
            end
            # THE INVARIANT: no tag may throw on the call site's keywords.
            # That is what broke, and it is branch-independent.
            @test got !== :THREW

            # The RETURN: AV has no method of its own by design and must come
            # back as nothing. SMAG, VREM and DSGS all build real models --
            # DSGS's allocator arrived with the dynamic-SGS merge, which is
            # also what put the extra keywords at the call site in the first
            # place, so the two belong in the same assertion.
            # AV has no method of its own by design. DSGS's struct is 3D-only
            # (`SD === NSD_3D() || return nothing` in its allocator), so in 2D
            # it legitimately returns nothing and the older in-rhs path runs
            # instead -- asserting a model there fails for a correct reason.
            if tag isa AV || (tag isa DSGS && SD == NSD_2D())
                @test got === nothing
            else
                @test got isa AbstractSGSModel
            end
        end
    end

    # And the minimal call -- only C_s -- must still work, because that is what
    # OTHER branches' call sites pass.
    for tag in (SMAG(), VREM())
        @test allocate_SGS(NP, Float64, CPU(), PC, tag; C_s = 0.23) isa AbstractSGSModel
    end
end
