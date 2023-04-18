module Jexpresso

using Dates
using Revise

const TInt   = Int64
const TFloat = Float64

#using DocStringExtensions

include("./problems/AbstractProblems.jl")

include("./kernel/abstractTypes.jl")

include("./kernel/globalStructs.jl")

include("./kernel/infrastructure/sem_setup.jl")

include("./kernel/boundaryconditions/BCs.jl")

include("./kernel/operators/rhs.jl")

include("./kernel/solvers/TimeIntegrators.jl")

include("./kernel/solvers/Axb.jl")

include("./io/mod_inputs.jl")

include("./io/write_output.jl")

include("./run.jl")

end
