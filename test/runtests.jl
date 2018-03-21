push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility
using Base.Test

include("FeynmanAlpha.jl") # Simple test of Frohlich alpha vs. literature values
include("HellwarthMultipleBranches.jl") # Test Hellwarth et al. 1999 PRB 'B' multiple branch reduction scheme
include("FrostPolaronMobility2017.jl") # Reproduce values published in Frost 2017 PRB 

println("That's me! If I finished without interupting, let's assume the tests were A.OK.")

