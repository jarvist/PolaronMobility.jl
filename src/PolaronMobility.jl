"""
PolaronMobility.jl - https://github.com/jarvist/PolaronMobility.jl
  Codes by Jarvist Moore Frost, 2017-2022; Bradley A.A. Martin 2020-2022
  PolaronMobility.jl is a Julia package which calculates the temperature-dependent polaron
  mobility for a material. It implements the Feynman variational approach with
  Osaka/Hellwarth free-energy for the polaron system.  
"""
module PolaronMobility

export Polaron, polaron                                   # Type to hold the polaron data.
export Material, material                                 # Type to hold material specific data.
export frohlichalpha, ϵ_ionic_mode, F, A, B, C, feynmanvw # Frohlich alpha, energy and minimisation functions.              
export HellwarthBScheme, HellwarthAScheme                 # Hellwarth effective mode schemes.
export polaron_memory_function                            # Polaron memory functions
export polaron_complex_impedence, polaron_complex_conductivity, optical_absorption        # Response functions.
export polaron_mobility, Hellwarth_mobility, Kadanoff_mobility_lowT, FHIP_mobility_lowT   # Mobility functions.
export save_polaron, load_polaron # Functions to save and load Polaron types to/from .jld format.
export reduce_array


##### load in library routines... #####
# stdlib
using LinearAlgebra
using Printf
using JLD

#Key numerics used:
# one-dimensional numerical integration in Julia using adaptive Gauss-Kronrod quadrature
import QuadGK.quadgk
# Using the powerful Julia Optim package to optimise the variational parameters
using Optim

include("FeynmanTheory.jl")     # Actions + variational functions.
include("HellwarthTheory.jl")   # multimode -> equivalent mode.
include("MemoryFunction.jl")    # Memory function X calculation.
include("ResponseFunctions.jl") # Linear reponse functions for polaron.
include("Material.jl")          # Material type and constructors.
include("Polaron.jl")           # Polaron type and constructors.

# Physical Constants
"Planck's constant, (kgm²s⁻¹)."
const hbar = const ħ = 1.054571817e-34
"Electron charge, (kgm²s⁻²)."
const eV = const q = const ElectronVolt = 1.602176634e-19
"Electron mass, (kg)."
const me = MassElectron = 9.1093837015e-31
"Boltzmann's constant, (kgm²K⁻¹)."
const Boltzmann = const kB = 1.380649e-23
"Permittivity of free space, (C²N⁻¹m⁻²)."
const ε_0 = ϵ_0 = 8.85418682e-12
"Speed of light, (ms⁻¹)."
const c = 299792458
"Atomic mass unit, (kg)"
const amu = 1.660_539_066_60e-27

reduce_array(a) = length(a) == 1 ? only(a) : dropdims(a, dims=tuple(findall(size(a) .== 1)...))

end # module

