abstract type AbstractEquations end

struct AdvDiff <: AbstractEquations end
struct LinearCLaw <:AbstractEquations end
struct Burgers <: AbstractEquations end
struct Elliptic <: AbstractEquations end
struct ShallowWater <: AbstractEquations end
struct SoilTopo <: AbstractEquations end
struct Helmholtz <: AbstractEquations end
struct CompEuler <: AbstractEquations end
