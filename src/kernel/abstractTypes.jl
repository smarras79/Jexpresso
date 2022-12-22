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
struct CG <:  AbstractDiscretization end


#
# Boundary flags/conditions
#
abstract type AbstractBC end
struct PERIODIC1D_CG <: AbstractBC end
