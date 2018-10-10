# PolaronMobility.jl - https://github.com/jarvist/PolaronMobility.jl
# Codes by Jarvist Moore Frost, 2017-2018
# Calculate Polaron Mobility - by a Osaka/Hellwarth variational solution to the Feynman model
module PolaronMobility

# These codes were developed with Julia 0.5.0 - Julia 0.6.2, and require the Optim and Plots packages.

export Polaron # Type to hold the data
export frohlichalpha, feynmanvw, F, polaronmobility, savepolaron, plotpolaron
export HellwarthBScheme, HellwarthAScheme
export ImX

##### load in library routines... #####
# stdlib
if VERSION > v"0.7.0-DEV.2005"
    using LinearAlgebra
    using Printf
end
# one-dimensional numerical integration in Julia using adaptive Gauss-Kronrod quadrature
import QuadGK.quadgk

# Using the powerful Julia Optim package to optimise the variational parameters
using Optim

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space

include("types.jl")            # Polaron types
include("FeynmanTheory.jl")    # Actions + variational functions  
include("HellwarthTheory.jl")  # multimode -> equivalent mode.
include("MobilityTheories.jl") # Main polaronmobility function 
include("Susceptibility.jl")   # ImX calculation
include("OedipusRex.jl")       # Optical Absorption

end # module

