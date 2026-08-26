#!/usr/bin/env julia
#
# diff_user_inputs.jl -- prove that a refactor of a case's user_inputs.jl did
#                        not change a single value it produces.
#
#     git show HEAD:problems/<EQS>/<CASE>/user_inputs.jl > /tmp/ui_old.jl
#     julia --startup-file=no tools/diff_user_inputs.jl \
#           /tmp/ui_old.jl problems/<EQS>/<CASE>/user_inputs.jl
#
# TWO THINGS HERE ARE LOAD-BEARING, and both are ways this check can pass while
# testing nothing:
#
#   * the two versions are evaluated in SEPARATE MODULES. Both define
#     `user_inputs`, so including them into one namespace would call the second
#     one twice and report agreement no matter what the first said.
#   * the difference counter is `global`. Written as a soft-scope local it is
#     never incremented, the final verdict reads `bad == 0` and is true by
#     construction -- which is exactly what the first version of this script
#     did. Verified by planting a one-value change and confirming it fails.
#
# Compare the input Dict produced by the pre- and post-refactor user_inputs().
# No Jexpresso load: the deck only NAMES a few constructors, so stubbing them
# exercises the real function body at a fraction of the cost. The stubs record
# their arguments, so a change of tableau or viscosity model still shows up.
module Stub
    struct Tagged; name::Symbol; args::Any; end
    Base.show(io::IO, t::Tagged) = print(io, t.name, t.args)
    for f in (:IMEX_ARK, :HEVI_ARK, :CarpenterKennedy2N54, :AV, :SMAG, :VREMAN,
              :NCL, :CL, :TOTAL, :PERT, :Inexact, :Exact, :ContGal, :LagGal)
        @eval $f(a...) = Tagged($(QuoteNode(f)), a)
        @eval export $f
    end
end
function load(path, name)
    m = Module(name)
    Core.eval(m, :(using Main.Stub))
    Base.include(m, path)
    Base.invokelatest(getfield(m, :user_inputs))
end
length(ARGS) == 2 ||
    error("usage: julia tools/diff_user_inputs.jl <old user_inputs.jl> <new user_inputs.jl>")
old = load(abspath(ARGS[1]), :Old)
new = load(abspath(ARGS[2]), :New)
ko, kn = Set(keys(old)), Set(keys(new))
isempty(setdiff(ko,kn)) || println("REMOVED: ", setdiff(ko,kn))
isempty(setdiff(kn,ko)) || println("ADDED:   ", setdiff(kn,ko))
bad = 0
for k in sort(collect(ko ∩ kn), by=string)
    if string(old[k]) != string(new[k])
        println("DIFFERS: ", k, "\n   old = ", old[k], "\n   new = ", new[k]); global bad += 1
    end
end
println(bad == 0 && ko == kn ? "IDENTICAL: $(length(kn)) keys, same values" : "*** $bad DIFFERENCE(S) ***")
