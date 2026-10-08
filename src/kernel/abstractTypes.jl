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
struct FD <: AbstractDiscretization end
# Finite volumes as order-zero DG (operators/fv.jl); mod_inputs maps it to DiscGal() + :lfv.
struct FV <: AbstractDiscretization end

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
struct THETA <: AbstractPert end

#
# viscosity type
#
abstract type AbstractVT end
struct AV    <: AbstractVT end
struct SMAG  <: AbstractVT end
struct DSMAG <: AbstractVT end
struct VREM  <: AbstractVT end
struct WALE  <: AbstractVT end
struct DSGS  <: AbstractVT end
# Marras-Nazarov residual-based Dynamic SGS for the 2D ideal GLM-MHD
# system (9 fields). Kept as its own tag rather than folded into DSGS()
# because the residual set, the equation-of-state and the wave speed all
# differ from the Euler-theta system DSGS() is written for.
struct DSGS_MHD <: AbstractVT end
# Marras-Nazarov residual-based Dynamic SGS for the 2D non-linear
# shallow-water system q = (H, Hu, Hv): residual of the three equations,
# wave speed |v| + sqrt(g H), one kinematic coefficient on every slot
# (problems/ShallowWater/SoliWaveIslandDSGS).
struct DSGS_SW <: AbstractVT end


abstract type AbstractVolumeFlux end
struct ranocha <: AbstractVolumeFlux end
struct artiano_ec <: AbstractVolumeFlux end
struct artiano_etec <: AbstractVolumeFlux end
struct artiano_tec <: AbstractVolumeFlux end
struct kennedy_gruber <: AbstractVolumeFlux end
struct central_euler <: AbstractVolumeFlux end
struct central_theta <: AbstractVolumeFlux end

#
# Numerical (interface) fluxes for DG
#
abstract type AbstractNumericalFlux end
struct upwind_flux <: AbstractNumericalFlux end
struct rusanov_flux <: AbstractNumericalFlux end

# Approximate Riemann solvers (dg_fluxes.jl) for DG and FV.
# ScalarLaw: scalar law, signal speeds from user_wave_speed_bounds (default ±λ).
# EulerIdealGas(γ): 2D Euler, energy form (ρ, ρu, ρv, ρE).
abstract type AbstractFluxPhysics end
struct ScalarLaw <: AbstractFluxPhysics end
struct EulerIdealGas <: AbstractFluxPhysics
    γ :: Float64
end
EulerIdealGas() = EulerIdealGas(1.4)

# HLL (Harten, Lax & van Leer 1983)
struct hll_flux{P <: AbstractFluxPhysics} <: AbstractNumericalFlux
    phys :: P
end
hll_flux() = hll_flux(ScalarLaw())

# HLLC (Toro, Spruce & Speares 1994); HLL for a scalar law
struct hllc_flux{P <: AbstractFluxPhysics} <: AbstractNumericalFlux
    phys :: P
end
hllc_flux() = hllc_flux(ScalarLaw())

# Roe (1981) with Harten's entropy fix; Murman–Roe for a scalar law
struct roe_flux{P <: AbstractFluxPhysics} <: AbstractNumericalFlux
    phys :: P
end
roe_flux() = roe_flux(ScalarLaw())

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

abstract type AbstractOutFormat end
struct PNG <: AbstractOutFormat end
struct ASCII <: AbstractOutFormat end
struct VTK <: AbstractOutFormat end
struct HDF5 <: AbstractOutFormat end
struct NETCDF <: AbstractOutFormat end
struct NONE <: AbstractOutFormat end
