module SpeedyPolaronMobility

using PolaronMobility

export speedymakepolaron

include("speedymakepolaron.jl")

# Physical constants
const ħ = 1.05457162825e-34         # kg m2 / s
const eV = 1.602176487e-19          # kg m2 / s2
const me = 9.10938188e-31           # kg
const kB = 1.3806504e-23            # kg m2 / K s2
const ϵ_0 = 8.854E-12               # permittivity of free space C2N−1m−2, 
const amu = 1.660_539_066_60e-27    # kg

end # module
