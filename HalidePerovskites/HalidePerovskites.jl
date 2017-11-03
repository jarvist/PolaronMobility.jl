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
using PlotPolaron # Plots dependency

##### load in library routines... #####
#"Rob-MAPCl",  Plot figures with Plots, which defaults to Pyplot backend
using Plots
pyplot()
default(grid=false) # No silly dotted grid lines
default(size=(400,250)) # A good small size for two-column EPS output

#default(size=(800,600)) # Nice size for small-ish PNGs for slides

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

# CsSnX3 X={Cl,Br,I}
# L. Huang, W.Lambrecht - PRB 88, 165203 (2013)
# Dielectric consts, from TABLE VII
# Effective masses from TABLE VI, mh*
const cm1=2.997e10 # cm-1 to Herz

function CsSn()
    CsSnCl3=polaronmobility(10:10:1000, 4.80, 29.4,  243cm1, 0.140, figures=true) # alpha= 1.386311
    savepolaron("CsSnCl3",CsSnCl3)

    CsSnBr3=polaronmobility(10:10:1000, 5.35, 32.4,  183cm1, 0.082, figures=true) # alpha= 1.094468 
    savepolaron("CsSnBr3",CsSnBr3)
 
    CsSnI3=polaronmobility(10:10:1000, 6.05, 48.2,  152cm1, 0.069, figures=true) # alpha= 1.020355 
    savepolaron("CsSnI3",CsSnI3)
 
end

CsSn() # generate all that lovely tin data

# Ts, βreds, Kμs, Hμs, FHIPμs, vs, ws, ks, Ms, As, Bs, Cs, Fs, Taus
# CsSnBr3.dat:300 0.877382 511.513 356.358 874.499 8.13686 7.31511 12.6976 0.23729 -3.78992 2.15682 0.84436 0.788735 0.185396
# CsSnCl3.dat:300 1.165048 212.034 147.36 272.993 6.50772 5.58522 11.1557 0.357614 -3.5556 2.42611 0.947833 0.181657 0.143969
# CsSnI3.dat:300 0.728755 703.502 487.085 1448.02 9.64022 8.80508 15.4044 0.198691 -4.01438 2.18761 0.859401 0.967371 0.207864

# So for CsSnSi3, a hole-mobility of 487(Hellwarth)-703(Kadanoff). Compares to expt. 400 (transport) -  585 (Hall effect). Ace!

#Ts,Kμs, Hμs, FHIPμs, ks, Ms, As, Bs, Cs, Fs, Taus
#effectivemass=0.12 # the bare-electron band effective-mass. 
# --> 0.12 for electrons and 0.15 for holes, in MAPI. See 2014 PRB.
# MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=
# MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=
MAPIe=polaronmobility(10:10:1000, 4.5, 24.1, 2.25E12, 0.12)
plotpolaron("MAPI-electron", MAPIe)
savepolaron("MAPI-electron",MAPIe)

MAPIh=polaronmobility(10:10:1000, 4.5, 24.1, 2.25E12, 0.15)
plotpolaron("MAPI-hole", MAPIh)
savepolaron("MAPI-hole",MAPIh)

# PCBM: 4.0, 5.0, 2.25Thz, effective-mass=1.0
#polaronmobility(4.0, 5.0, 2.25E12, 1.00)

# CsPbI3: 3.82, 8.27, 2.83 THz (highest freq mode in cubic) 
#    ^-- v. quick Kmesh=3x3x3 PBESol Calc, JMF 2017-04-12. Files: jarvist@titanium:~/phonopy-work/2017-03-Scott-Distortions/1000-Dielectric/3x3x3-kmesh/
# Horribly non-converged! 
#                         Kmesh=6x6x6; 7.2/12.1
# Kmesh=9x9x9x, Ediff=10^-9;           6.1/12.0, 2.57 THz
function CsPbI()
    CsPbI=polaronmobility(10:10:1000, 6.1,6.1+12.0, 2.57E12, 0.12)
    plotpolaron("CsPbI3-electron",CsPbI)
    savepolaron("CsPbI3-electron",CsPbI)
end
CsPbI() # generate CsPbI3 data for later plotting

function SendnerCrosscheck()
    const cm1=2.997e10 # cm-1 to Herz

    # Rob's / Sendner's paper - values extracted form IR measures
    # https://doi.org/10.1039%2Fc6mh00275g

    # These data mainly from Private communication. 
    # It is not well described in the paper but they
    # used one the Hellwarth effective-mode method to reduce the observed IR
    # oscillators down to a single mode. Rob provided these values by email (as
    # below). The Effetive Masses and dielectric constants are from the paper.
    # Combined, this reproduces their mobilities, and the internal w and
    # v parameters (again, data from email from Rob):
    #  Paper: "For omega we used: MAPbI/Br/Cl = 112.9/149.4/214.0"
    RobMAPI =polaronmobility(10:10:400, 5.0, 33.5, 112.9*cm1, 0.104)
    RobMAPBr=polaronmobility(10:10:400, 4.7, 32.3,149.4*cm1, 0.117)
    RobMAPCl=polaronmobility(10:10:400, 4.0, 29.8, 214.0*cm1, 0.2)

    savepolaron("Rob-MAPI",RobMAPI)
    savepolaron("Rob-MAPBr",RobMAPBr)
    savepolaron("Rob-MAPCl",RobMAPCl)

    plot(RobMAPI.T,RobMAPI.Hμ,label="(Rob's values) MAPI",markersize=2,marker=:uptriangle,ylim=(0,400))
    plot!(RobMAPBr.T,RobMAPBr.Hμ,label="(Rob's values) MAPBr",markersize=2,marker=:diamond)
    plot!(RobMAPCl.T,RobMAPCl.Hμ,label="(Rob's values) MAPCl",markersize=2,marker=:diamond)
    savefig("Rob-comparison.png")
#    savefig("Rob-comparison.eps")
end

SendnerCrosscheck() # cross-check data for Rob

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
plot!(Saidaminov[:,1],Saidaminov[:,2],label="",markersize=6,marker=:utriangle)
plot!(Semonin[:,1],Semonin[:,2],label="",markersize=6,marker=:square)
#plot!(Saidaminov[:,1],Saidaminov[:,2],label="Saidaminov JV Single Crystal", markersize=6,marker=:utriangle)
#plot!(Semonin[:,1],Semonin[:,2],label="Semonin Single Crystal TRMC", markersize=6,marker=:hexagon)


#plot!(MAPIe.T,MAPIe.Kμs,label="(electrons) Kadanoff Polaron mobility",marker=2)
#plot!(Ts,hHμs,label="Calculated (hole) mobility",markersize=3,marker=:dtriangle)

plot!(MAPIe.T,MAPIe.Hμ,label="Calculated (electron) mobility",markersize=3,marker=:dtriangle)
plot!(MAPIh.T,MAPIh.Hμ,label="Calculated (hole) mobility",markersize=3,marker=:utriangle)

savefig("MAPI-eh-mobility-calculated-experimental.png")
savefig("MAPI-eh-mobility-calculated-experimental.eps")

plot!(ylims=(),yscale=:log10)
plot!(ylims=(),xscale=:log10) # on log-log axes; powerlaw becomes straight
plot!(x->1e5*x^(-3/2),20:1000,yscale=:log10) # Trend line for the T^-3/2 powerlaw

savefig("MAPI-eh-mobility-calculated-experimental-log10.png")
savefig("MAPI-eh-mobility-calculated-experimental-log10.eps")

# OK; fit some stuff.
using LsqFit # Least squares fit, orig from Optim package.

function fitPowerLaw()
    model(x,p)=p[1]*x.^(p[2])
    p0=[1.0,-3.0/2]

    # Data we're fitting against
    plot(MAPIe.T,MAPIe.Hμ,label="Calculated (electron) mobility",markersize=3,marker=:dtriangle,ylims=(1,400))

    # They're nasty 'Any' data types; force to Float64
    T=convert(Array{Float64},MAPIe.T)
    μ=convert(Array{Float64},MAPIe.Hμ)

    # Ta da! Least squares fit.
    fit=curve_fit(model,T,μ,p0)
    p=fit.param
    #plot!(T,model(T,p))

    # Fit just from 100K onwards
    fit=curve_fit(model,T[10:end],μ[10:end],p0)
    p=fit.param
    plot!(T[10:end],model(T[10:end],p), label="Calc [10:end] Powerlaw: $p")

    # Fit all from 50K onwards - very dubious of quality of fit
    #fit=curve_fit(model,T[5:end],μ[5:end],p0)
    #p=fit.param
    #plot!(T[5:end],model(T[5:end],p), label="Calc [5:end] Powerlaw: $p")

    # TRMC Polycrystal data for comparison
    plot!(Milot[:,1],Milot[:,2],label="Milot T-dep TRMC Polycrystal",marker=:hexagon,markersize=3,
    xlab="Temperature (K)",ylab="Mobility (cm\$^2\$/Vs)", ylims=(0,400) )

    # Fit all data points except for spurious first one
    MT=convert(Array{Float64},Milot[2:end,1])
    Mu=convert(Array{Float64},Milot[2:end,2])
    fit=curve_fit(model,MT,Mu,p)
    p=fit.param
    plot!(T[5:end], model(T[5:end],p), label="Milot Powerlaw Fit 2:end: $p")

    # Fit only sensible data points > 120 K and up
    MT=convert(Array{Float64},Milot[4:end,1])
    Mu=convert(Array{Float64},Milot[4:end,2])
    fit=curve_fit(model,MT,Mu,p)
    p=fit.param
    plot!(T[5:end], model(T[5:end],p), label="Milot Powerlaw Fit 4:end: $p")

    savefig("MAPI-fit-lin.png")
    savefig("MAPI-fit-lin.pdf")

    plot!(ylims=(1,), yscale=:log10)# log-log; power-law --> straight line
    plot!(xlims=(1,), xscale=:log10)
    savefig("MAPI-fit-log.png")
    plot!(legend=false) # no key
    savefig("MAPI-fit-log.pdf")
end

fitPowerLaw()

#plot!(RobMAPI.T,RobMAPI.Hμ,label="(Rob's values) MAPI",markersize=2,marker=:rect)
#savefig("MAPI-eh-mobility-calculated-experimental-Rob.png")
#savefig("MAPI-eh-mobility-calculated-experimental-Rob.eps")


println("That's me!")



