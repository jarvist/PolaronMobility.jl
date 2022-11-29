push!(LOAD_PATH, "../src/") # load module from local directory

using PolaronMobility
using Printf # used for some stdout

using Test

@testset "PolaronMobility" begin

    include("FrohlichAlpha.jl") # Simple test of Frohlich alpha vs. literature values
    include("FeynmanAthermal.jl") # Athermal Feynman tests
    include("HellwarthEffectiveFrequency.jl") # Test Hellwarth et al. 1999 PRB 'B' multiple branch reduction scheme
    include("FrostPolaronMobility2017.jl") # Reproduce values published in Frost 2017 PRB 
    include("MultipleBranches.jl") # Test explicit Oct 2019:-> explicit phonon branches
    include("MultiplePhonons.jl") # Test 2021 work on multiple phonon modes
end

println("\nThat's me! If I finished without interupting, all tests have passed.")

