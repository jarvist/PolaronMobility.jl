# HalidePerovskites
# Uses the newly forked off FeynmanKadanoffOsakaHellwarth.jl package
# Codes by Jarvist Moore Frost, 2017
# These codes were developed with Julia 0.5.0, and requires the Optim and Plots packages.

# This file, when run under Julia, should regenerate all data (and plots) associated with Arxiv paper:
# https://arxiv.org/abs/1704.05404
# Polaron mobility in halide perovskites
# Jarvist Moore Frost
# (Submitted on 18 Apr 2017 [v1])

push!(LOAD_PATH,"./") # load module from local directory

using FeynmanKadanoffOsakaHellwarth

##### load in library routines... #####
# Plot figures with Plots, which defaults to Pyplot backend
using Plots
#default(size=(800,600)) # For the .png file output

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space

" Copy and pasted out of a Jupyter notebook; this calculates 'alpha' parameters
for various materials, as a comparison to the literature used when figuring out
the oft-quoted units. "
function checkalpha()
	println(" Alpha-parameter, Cross check 'feynmanalpha()' fn vs. literature values.\n")
	print("\t NaCl Frohlich paper α=",feynmanalpha(2.3, 5.6, (4.9E13/(2*pi)), 1.0)) 
    println(" should be ~about 5 (Feynman1955)")
	print("\t CdTe  α=",feynmanalpha(7.1,   10.4,  5.08E12, 0.095)) 
    println(" Stone 0.39 / Devreese 0.29 ")
	print("\t GaAs  α=",feynmanalpha(10.89, 12.9,  8.46E12, 0.063)) 
    println(" Devreese 0.068 ")

    println()
    println("Guess at PCBM: 4.0, 6.0 ; α=",feynmanalpha(4.0,6.0, 1E12, 50))
    println("MAPI:")
    println("MAPI  4.5, 24.1, 9THz ; α=",feynmanalpha(4.5,   24.1,  9.0E12,    0.12))
    println("MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=",feynmanalpha(4.5,   24.1,  2.25E12,    0.12))
    println("MAPI  6.0, 25.7, 9THz ; α=",feynmanalpha(6.0,   25.7,  9.0E12,    0.12))
    println("MAPI  6.0, 36.0, 9THz ; α=",feynmanalpha(6.0,   36,  9.0E12,    0.12))
    println("MAPI  6.0, 36.0, 1THz ; α=",feynmanalpha(6.0,   36,  1.0E12,    0.12))
end
checkalpha()

#####
# Call simulation

#Ts,Kμs, Hμs, FHIPμs, ks, Ms, As, Bs, Cs, Fs, Taus
#effectivemass=0.12 # the bare-electron band effective-mass. 
# --> 0.12 for electrons and 0.15 for holes, in MAPI. See 2014 PRB.
# MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=
Ts,eKμs, eHμs, eFHIPμs, eks, eMs, eAs, eBs, eCs, eFs, eTaus = polaronmobility("MAPI-electron", 4.5, 24.1, 2.25E12, 0.12)
Ts,hKμs, hHμs, hFHIPμs, hks, hMs, hAs, hBs, hCs, hFs, hTaus = polaronmobility("MAPI-hole",     4.5, 24.1, 2.25E12, 0.15)

# PCBM: 4.0, 5.0, 2.25Thz, effective-mass=1.0
#polaronmobility("PCBM",     4.0, 5.0, 2.25E12, 1.00)

# CsPbI3: 3.82, 8.27, 2.83 THz (highest freq mode in cubic) 
#    ^-- v. quick Kmesh=3x3x3 PBESol Calc, JMF 2017-04-12. Files: jarvist@titanium:~/phonopy-work/2017-03-Scott-Distortions/1000-Dielectric/3x3x3-kmesh/
# Horribly non-converged! 
#                         Kmesh=6x6x6; 7.2/12.1
# Kmesh=9x9x9x, Ediff=10^-9;           6.1/12.0, 2.57 THz
polaronmobility("CsPbI3-electron",     6.1,6.1+12.0, 2.57E12, 0.12)


function SendnerCrosscheck()
    cm1=2.997e10 # cm-1 to Herz

    # Rob's / Sendner's paper - values extracted form IR measures
    # https://doi.org/10.1039%2Fc6mh00275g
    #Ts,a,MAPI=polaronmobility("Rob-MAPI", 5.0, 33.5, 40*cm1, 0.104)
    #Ts,a,MAPBr=polaronmobility("Rob-MAPBr", 4.7, 32.3, 51*cm1, 0.117)
    #Ts,a,MAPCl=polaronmobility("Rob-MAPCl", 4.0, 29.8, 70*cm1, 0.2)

    # Private communication. It is not well described in the paper but they
    # used one the Hellwarth effective-mode method to reduce the observed IR
    # oscillators down to a single mode. Rob provided these values by email (as
    # below). The Effetive Masses and dielectric constants are from the paper.
    # Combined, this reproduces their mobilities, and the internal w and
    # v parameters (again, email from Rob).  
    # For omega we used: MAPbI/Br/Cl = 112.9/149.4/214.0
    Ts,a,MAPI=polaronmobility("Rob-MAPI", 5.0, 33.5, 112.9*cm1, 0.104)
    Ts,a,MAPBr=polaronmobility("Rob-MAPBr", 4.7, 32.3,149.4*cm1, 0.117)
    Ts,a,MAPCl=polaronmobility("Rob-MAPCl", 4.0, 29.8, 214.0*cm1, 0.2)

    plot(Ts,MAPI,label="(Rob's values) MAPI",markersize=2,marker=:uptriangle,ylim=(0,400))
    plot!(Ts,MAPBr,label="(Rob's values) MAPBr",markersize=2,marker=:diamond)
    plot!(Ts,MAPCl,label="(Rob's values) MAPCl",markersize=2,marker=:diamond)
    savefig("Rob-comparison.png")
    savefig("Rob-comparison.eps")
end




#####
## Expt. data to compare against
# Milot/Herz 2015 Time-Resolved-Microwave-Conductivity mobilities
# Data from table in SI of: DOI: 10.1002/adfm.201502340
# Absolute values possibly dodge due to unknown yield of charge carriers; but hopefully trend A.OK!
Milot= [
8 184
40 321
80 143
120 62
140 40
160 52
180 44
205 41
230 39
265 26
295 35
310 24
320 24
330 19
340 16
355 15 
]

# IV estimated mobilities (?!) from large single crystals, assumed ambient T
# Nature Communications 6, Article number: 7586 (2015)
# doi:10.1038/ncomms8586
Saidaminov = 
[ 300 67.2 ]

#Semonin2016,
#  doi = {10.1021/acs.jpclett.6b01308},
Semonin = 
[ 300 115 ] # +- 15 cm^2/Vs, holes+electrons

#####
## Calculated mobilities vs. expt
plot(Milot[:,1],Milot[:,2],label="Milot T-dep TRMC Polycrystal",marker=:hexagon,markersize=3,
xlab="Temperature (K)",ylab="Mobility (cm\$^2\$/Vs)", ylims=(0,400) )
plot!(Saidaminov[:,1],Saidaminov[:,2],label="Saidaminov JV Single Crystal", markersize=6,marker=:utriangle)
plot!(Semonin[:,1],Semonin[:,2],label="Semonin Single Crystal TRMC", markersize=6,marker=:hexagon)

#plot!(Ts,eKμs,label="(electrons) Kadanoff Polaron mobility",marker=2)
plot!(Ts,eHμs,label="(electrons) Hellwarth1999 Polaron mobility",markersize=3,marker=:diamond)

#plot!(Ts,hKμs,label="(holes) Kadanoff Polaron mobility",marker=2)
plot!(Ts,hHμs,label="(holes) Hellwarth1999 Polaron mobility",markersize=3,marker=:dtriangle)

savefig("MAPI-eh-mobility-calculated-experimental.png")
savefig("MAPI-eh-mobility-calculated-experimental.eps")

plot!(Ts,MAPI,label="(Rob's values) MAPI",markersize=2,marker=:rect)
savefig("MAPI-eh-mobility-calculated-experimental-Rob.png")
savefig("MAPI-eh-mobility-calculated-experimental-Rob.eps")


println("That's me!")



