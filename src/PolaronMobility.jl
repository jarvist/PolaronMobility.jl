"""
PolaronMobility.jl - https://github.com/jarvist/PolaronMobility.jl
  Codes by Jarvist Moore Frost, 2017-2023; Bradley A.A. Martin 2020-2023
  PolaronMobility.jl is a Julia package which calculates the temperature-dependent polaron
  mobility for a material. It implements the Feynman variational approach with
  Osaka/Hellwarth free-energy for the polaron system.  
"""
module PolaronMobility

 # Type to hold the polaron data.
export Polaron, polaron                                  

 # Type to hold material specific data.
export Material, material                                

# Frohlich alpha, energy and minimisation functions.   
export frohlichalpha, ϵ_ionic_mode, F, A, B, C, feynmanvw

# Hellwarth effective mode schemes.
export HellwarthBScheme, HellwarthAScheme                 

# Polaron memory functions
export polaron_memory_function, D, S, polaron_memory_function_athermal, polaron_memory_function_thermal                      

# Response functions
export polaron_complex_impedence, polaron_complex_conductivity, optical_absorption      

# Mobility functions
export polaron_mobility, Hellwarth_mobility, Kadanoff_mobility_lowT, FHIP_mobility_lowT 

# Functions to save and load Polaron types to/from .jld format
export save_polaron, load_polaron

# Polaron effective mass
export polaron_effective_mass

# Generalised polaron functions
export general_memory_function, vw_variation, polaron_propagator, phonon_propagator, electron_energy

# Holstein functions
export holstein_interaction_energy_integrand, holstein_interaction_energy, holstein_energy, holstein_structure_factor, holstein_memory_function, holstein_mobility

# K-space functions
export cartesian_k_integrand, spherical_k_integral, holstein_coupling, frohlich_coupling, holstein_interaction_energy_integrand_k_space, frohlich_interaction_energy_integrand_k_space, holstein_interaction_energy_k_space, frohlich_interaction_energy_k_space, holstein_energy_k_space, frohlich_energy_k_space, holstein_structure_factor_k_space, frohlich_structure_factor_k_space, holstein_memory_function_k_space, holstein_mobility_k_space, frohlich_memory_function_k_space, frohlich_mobility_k_space

##### load in library routines... #####
# stdlib
using LinearAlgebra, Printf, JLD

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

# Unitful functions
export puconvert, punit, pustrip, m0_pu, e_pu, ħ_pu, k_pu, ω0_pu, a0_pu, E0_pu, β0_pu, T0_pu, μ0_pu, t0_pu, addunits!

include("FeynmanTheory.jl")     # Frohlich actions + variational functions.
include("HolsteinTheory.jl")    # Holstein actions + variational functions.
include("KSpaceTheory.jl")      # K-space dependent integral code.
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
reduce_array(a) = length(a) == 1 ? only(a) : Array(dropdims(a, dims=tuple(findall(size(a) .== 1)...)))

# Utility functions
export reduce_array

end # module
