push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility
using Base.Test

include("FeynmanAlpha.jl") # Simple test of Frohlich alpha vs. literature values
include("FeynmanAthermal.jl") # Athermal Feynman tests
include("HellwarthMultipleBranches.jl") # Test Hellwarth et al. 1999 PRB 'B' multiple branch reduction scheme
include("FrostPolaronMobility2017.jl") # Reproduce values published in Frost 2017 PRB 

println("\nThat's me! If I finished without interupting, all tests have passed.")

