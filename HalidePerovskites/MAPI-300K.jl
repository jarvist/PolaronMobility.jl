# HalidePerovskites
# Uses the newly forked off FeynmanKadanoffOsakaHellwarth.jl package
# Codes by Jarvist Moore Frost, 2017
# These codes were developed with Julia 0.5.0, and requires the Optim and Plots packages.

# This file, when run under Julia, should regenerate all data (and plots) associated with Arxiv paper:
# https://arxiv.org/abs/1704.05404
# Polaron mobility in halide perovskites
# Jarvist Moore Frost
# (Submitted on 18 Apr 2017 [v1])

push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility

##### load in library routines... #####
# Plot figures with Plots, which defaults to Pyplot backend
using Plots
#pyplot()
gr()
#default(grid=false) # No silly dotted grid lines
#default(size=(400,300)) # A good small size for two-column EPS output
#gr()
#default(size=(800,600)) # Nice size for small-ish PNGs for slides

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space

#####
# Call simulation

# CsSnX3 X={Cl,Br,I}
# L. Huang, W.Lambrecht - PRB 88, 165203 (2013)
# Dielectric consts, from TABLE VII
# Effective masses from TABLE VI, mh*
const cm1=2.997e10 # cm-1 to Herz

#Ts,Kμs, Hμs, FHIPμs, ks, Ms, As, Bs, Cs, Fs, Taus
#effectivemass=0.12 # the bare-electron band effective-mass. 
# --> 0.12 for electrons and 0.15 for holes, in MAPI. See 2014 PRB.
# MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=
println("OK, solving Polaron problem...")
MAPIe = polaronmobility("MAPI-electron", [150,300], 4.5, 24.1, 2.25E12, 0.12; verbose=true, figures=false)
MAPIh = polaronmobility("MAPI-hole",     [150,300], 4.5, 24.1, 2.25E12, 0.15; verbose=true, figures=false)

MAPIe = polaronmobility("MAPI-electron", [150,300], 4.5, 24.1, 1.25E12, 0.12; verbose=true, figures=false)


s=ImX(0.1:0.1:6,MAPIe.v[2],MAPIe.w[2],MAPIe.βred[2], MAPIe.α[1],MAPIe.ω[1],MAPIe.mb[1])
plot( s.nu,s.ImX,label="ImX",
    markersize=3,marker=:downtriangle, xlab="nu (units Omega)",ylab="ImX")
savefig("MAPIe-ImX.png")
plot( s.nu,s.μ,label="mu",
    markersize=3,marker=:uptriangle, xlab="nu (units Omega)",ylab="Mob")
savefig("MAPIe-mu.png")


println("That's me!")

