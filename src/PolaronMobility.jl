"""
PolaronMobility.jl - https://github.com/jarvist/PolaronMobility.jl
  Codes by Jarvist Moore Frost, 2017-2023; Bradley A.A. Martin 2020-2023
  PolaronMobility.jl is a Julia package which calculates the temperature-dependent polaron
  mobility for a material. It implements the Feynman variational approach with
  Osaka/Hellwarth free-energy for the polaron system.  
"""
module PolaronMobility

export Polaron, polaron                                   # Type to hold the polaron data.
export Material, material                                 # Type to hold material specific data.
export frohlichalpha, ϵ_ionic_mode, F, A, B, C, feynmanvw # Frohlich alpha, energy and minimisation functions.              
export HellwarthBScheme, HellwarthAScheme                 # Hellwarth effective mode schemes.
export polaron_memory_function, D, S                      # Polaron memory functions
export polaron_complex_impedence, polaron_complex_conductivity, optical_absorption        # Response functions.
export polaron_mobility, Hellwarth_mobility, Kadanoff_mobility_lowT, FHIP_mobility_lowT   # Mobility functions.
export save_polaron, load_polaron # Functions to save and load Polaron types to/from .jld format.
export reduce_array
export polaron_effective_mass
export holstein_B, holstein_energy, holsteinvw
export addunits!

export puconvert, punit, pustrip, m0_pu, e_pu, ħ_pu, k_pu, ω0_pu, a0_pu, E0_pu, β0_pu, T0_pu, μ0_pu, t0_pu

##### load in library routines... #####
# stdlib
using LinearAlgebra
using Printf
using JLD

#Key numerics used:
# One-dimensional numerical integration in Julia using adaptive Gauss-Kronrod quadrature
import QuadGK.quadgk
# Special function arising from cutoff k-space integrals (c.f. Holstein model).
import SpecialFunctions.erf
# Using the powerful Julia Optim package to optimise the variational parameters
using Optim

# Import required Unitful units, for mutability (generating new polaron units)
using Unitful
using Unitful: @unit, Dimension, Dimensions, NoDims, NoUnits, Units, dimension, uconvert, ustrip

include("FeynmanTheory.jl")     # Frohlich actions + variational functions.
include("HolsteinTheory.jl")    # Holstein actions + variational functions.
include("HellwarthTheory.jl")   # multimode -> equivalent mode.
include("MemoryFunction.jl")    # Memory function X calculation.
include("EffectiveMass.jl")     # Effective mass (Feynman's ansatz).
include("ResponseFunctions.jl") # Linear reponse functions for polaron.
include("Material.jl")          # Material type and constructors.
include("Polaron.jl")           # Polaron type and constructors.
include("PolaronUnits.jl")      # Implements internal Feynman 'polaron units', and SI conversions

# Register newly defined units with Unitful
Unitful.register(PolaronMobility)
__init__() = Unitful.register(PolaronMobility)

# QOL function for removing singleton dimensions.
reduce_array(a) = length(a) == 1 ? only(a) : dropdims(a, dims=tuple(findall(size(a) .== 1)...))

end # module

