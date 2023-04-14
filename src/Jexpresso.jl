module Jexpresso

using Dates
using Revise

const TInt   = Int64
const TFloat = Float64

#using DocStringExtensions

include(joinpath("problems", "AbstractProblems.jl"))

include(joinpath("kernel", "abstractTypes.jl"))

include(joinpath("kernel", "globalStructs.jl"))

include(joinpath("kernel", "physics", "globalConstantsPhysics.jl"))

include(joinpath("kernel", "physics", "constitutiveLaw.jl"))

include(joinpath("kernel", "infrastructure", "sem_setup.jl"))

include(joinpath("kernel", "boundaryconditions", "BCs.jl"))

include(joinpath("kernel", "operators", "rhs.jl"))

include(joinpath("kernel", "solvers", "TimeIntegrators.jl"))

include(joinpath("kernel", "solvers", "Axb.jl"))

include(joinpath("io", "mod_inputs.jl"))

include(joinpath("io", "write_output.jl"))

include("./run.jl")

end
