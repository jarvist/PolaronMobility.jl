# HalidePerovskites
# Uses the newly forked off FeynmanKadanoffOsakaHellwarth.jl package
# Codes by Jarvist Moore Frost, 2017
# These codes were developed with Julia 0.5.0, and requires the Optim and Plots packages.

# This file, when run under Julia, should regenerate polaron data associated with Arxiv paper:
# https://arxiv.org/abs/1708.04158
# Slow cooling of hot polarons in halide perovskite solar cells
# Jarvist Moore Frost, Lucy D. Whalley and Aron Walsh
# (Submitted on 14 Aug 2017 [v1])

push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility 
using PlotPolaron # Plots dependency

##### load in library routines... #####
# Plot figures with Plots, which defaults to Pyplot backend
using Plots
pyplot() # PyPlot (matplotlib) backend, to Plots
#gr() # GR backend, to Plots
default(grid=false) # No silly dotted grid lines
default(size=(400,300)) # A good small size for two-column EPS output
#default(size=(800,600)) # Nice size for small-ish PNGs for slides

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space
# Units
Å=1E-10

# --> 0.12 for electrons and 0.15 for holes, in MAPI. See 2014 PRB.
# MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=
MAPIe=polaronmobility(10:20:1000, 4.5, 24.1, 2.25E12, 0.12)
savepolaron("MAPI-electron",MAPIe)
plotpolaron("MAPI-electron", MAPIe)
MAPIh=polaronmobility(10:20:1000, 4.5, 24.1, 2.25E12, 0.15)
plotpolaron("MAPI-hole", MAPIh)
savepolaron("MAPI-hole",MAPIh)

## Polaron radius vs. Temperature
# More complex figure for the Thermal Pathways paper;
# https://github.com/WMD-group/2017-03-ThermalPathways/
p=MAPIe 

plot(p.T,p.rfsi./Å, markersize=2,marker=:rect,
    label="Polaron radius",xlab="Temperature (K)",ylab="Polaron Radius (Angstrom)",ylims=(0,Inf))
plot!(p.T,p.rfsmallalpha./Å,label="T=0 Schultz small alpha polaron radius")
savefig("ThermalPathways-radius.png")
savefig("ThermalPathways-radius.pdf")

plot(p.T,p.v./p.w,label="v/w",markersize=2,marker=:circle,
    xlab="Temperature (K)", ylab="\hbar\omega")
savefig("ThermalPathways-vwratio.png")
savefig("ThermalPathways-vwratio.pdf")

    plot!(p.T,p.v,label="v",markersize=2,marker=:uptriangle)
plot!(p.T,p.w,label="w",markersize=2,marker=:downtriangle)
savefig("ThermalPathways-vw.png")
savefig("ThermalPathways-vw.pdf")

println("That's me!")

