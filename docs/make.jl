push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility, Documenter

makedocs()

deploydocs(repo="github.com/jarvist/PolaronMobility.jl.git")

