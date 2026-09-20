# Precompile execution script for PackageCompiler.
# PackageCompiler runs this once (single MPI process) to trace which
# code paths and type specializations to bake into the sysimage.
#
# Requirements before running create_sysimage():
#   - Comment out `using Revise` and `using SnoopCompile` in src/Jexpresso.jl
#     (they are incompatible with sysimage builds)
#

push!(empty!(ARGS), "CompEuler", "3d")
include("src/Jexpresso.jl")
