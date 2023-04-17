abstract type AbstractProblem end

struct AdvDiff <: AbstractProblem end
struct LinearCLaw <:AbstractProblem end
struct Burgers <: AbstractProblem end
struct Elliptic <: AbstractProblem end
struct ShallowWater <: AbstractProblem end
struct SoilTopo <: AbstractProblem end
struct Helmholtz <: AbstractProblem end
struct Euler <: AbstractProblem end
