# Organic materials - polaron mobility 

push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility 

##### load in library routines... #####
# Plot figures with Plots, which defaults to Pyplot backend
#using Plots
#default(grid=false) # No silly dotted grid lines
#default(size=(400,300)) # A good small size for two-column EPS output
#default(size=(800,600)) # Nice size for small-ish PNGs for slides

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space

#####
# Call simulation

# PCBM: 4.0, 5.0, ??? , effective-mass=1.0
PCBM=polaronmobility(300, 4.0, 5.0, 2.0E12, 1.00)
savepolaron("PCBM",PCBM)

println("That's me!")

