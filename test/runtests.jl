push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility

if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "PolaronMobility" begin

include("FrohlichAlpha.jl") # Simple test of Frohlich alpha vs. literature values
include("FeynmanAthermal.jl") # Athermal Feynman tests
include("HellwarthMultipleBranches.jl") # Test Hellwarth et al. 1999 PRB 'B' multiple branch reduction scheme
include("FrostPolaronMobility2017.jl") # Reproduce values published in Frost 2017 PRB 

end

println("\nThat's me! If I finished without interupting, all tests have passed.")

