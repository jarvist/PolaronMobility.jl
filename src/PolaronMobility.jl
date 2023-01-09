"""
PolaronMobility.jl - https://github.com/jarvist/PolaronMobility.jl
  Codes by Jarvist Moore Frost, 2017-2022; Bradley A.A. Martin 2020-2022
  PolaronMobility.jl is a Julia package which calculates the temperature-dependent polaron
  mobility for a material. It implements the Feynman variational approach with
  Osaka/Hellwarth free-energy for the polaron system.  
"""
module PolaronMobility

export Polaron, NewPolaron # Type to hold the data
export frohlichalpha, feynmanvw, F, polaronmobility, savepolaron
export HellwarthBScheme, HellwarthAScheme
export polaron_memory_function  # Polaron memory functions
export optical_absorption       # Polaron optical absorption
export ϵ_ionic_mode, multi_frohlichalpha, extended_feynmanvw, multi_F, polaron_mobility, polaron_complex_impedence, polaron_complex_conductivity
export frohlichPartial, IRtoDielectric, IRtoalpha, DielectricFromIRmode
export Hellwarth1999mobilityRHS
export make_polaron, save_polaron, load_polaron, combine_polarons


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

include("Types.jl")            # Polaron types
include("FeynmanTheory.jl")    # Actions + variational functions
include("HellwarthTheory.jl")  # multimode -> equivalent mode.
include("MobilityTheories.jl") # Main polaronmobility function
include("MemoryFunction.jl")   # Memory function X calculation.
include("Susceptibility.jl")   # ImX calculation
include("OedipusRex.jl")       # Optical Absorption

end # module

