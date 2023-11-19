#
# General abstract not tied to any specific problem
#

#
# Quadrature rules 
#
abstract type AbstractIntegrationType end
struct Exact <: AbstractIntegrationType end
struct Inexact <: AbstractIntegrationType end
#
# Monolithic/tensor-product
#
abstract type AbstractMatrixType end
struct Monolithic <: AbstractMatrixType end
struct TensorProduct <: AbstractMatrixType end

#
# Space dimensions
#
abstract type AbstractSpaceDimensions end
struct NSD_1D <: AbstractSpaceDimensions end
struct NSD_2D <: AbstractSpaceDimensions end
struct NSD_3D <: AbstractSpaceDimensions end

#
# Space discretization
#
abstract type AbstractDiscretization end
struct ContGal <: AbstractDiscretization end
struct DiscGal <: AbstractDiscretization end

abstract type AbstractPointsType end
struct LG <: AbstractPointsType end
struct LGL <: AbstractPointsType end
struct CG <: AbstractPointsType end
struct CGL <: AbstractPointsType end
struct LGR <: AbstractPointsType end
#
# System of reference
#
abstract type AbstractMetricForm end
struct COVAR <: AbstractMetricForm end
struct CNVAR <: AbstractMetricForm end


#
# Coservation vs non-conservation formulation
#
abstract type AbstractLaw end
struct CL <: AbstractLaw end
struct NCL <: AbstractLaw end

#
# Solve for perturbation vs not perturbation variables
# ex. solve for either α or for α' = α-αref:
#
abstract type AbstractPert end
struct PERT  <: AbstractPert end
struct TOTAL <: AbstractPert end

abstract type AbstractOutFormat end
struct PNG <: AbstractOutFormat end
struct ASCII <: AbstractOutFormat end
struct VTK <: AbstractOutFormat end
struct HDF5 <: AbstractOutFormat end

#
# Boundary flags/conditions
#
abstract type AbstractBC end
struct PERIODIC1D_CG <: AbstractBC end
struct DefaultBC <: AbstractBC end
struct LinearClaw_1 <: AbstractBC end
struct LinearClaw_KopRefxmax <: AbstractBC end
struct DirichletExample <: AbstractBC end
struct bc_space_function <: AbstractBC end
