abstract type AbstractProblem end

struct Wave1D <: AbstractProblem end
struct AdvDiff <: AbstractProblem end
struct LinearCLaw <: AbstractProblem end
struct NS <: AbstractProblem end
struct SW <: AbstractProblem end
