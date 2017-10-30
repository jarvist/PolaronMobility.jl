push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility
using Base.Test

# write your own tests here
#@test 1 == 2

include("FeynmanAlpha.jl")
include("HellwarthMultipleBranches.jl")
include("FrostPolaronMobility2017.jl")
