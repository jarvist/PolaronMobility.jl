# ((freq THz)) ((IR Activity / e^2 amu^-1))
# These data from MAPbI3-Cubic_PeakTable.csv
# https://github.com/WMD-group/Phonons/tree/master/2015_MAPbI3/SimulatedSpectra
# Data published in Brivio2015 (PRB)
# https://doi.org/10.1103/PhysRevB.92.144308
MAPI= [
96.20813558773261 0.4996300522819191
93.13630357703363 1.7139631746083817
92.87834578121567 0.60108592692181
92.4847918585963 0.0058228799414729
92.26701437594754 0.100590086574602
89.43972834606603 0.006278895133832249
46.89209141511332 0.2460894564364346
46.420949316788 0.14174282581124137
44.0380222871706 0.1987196948553428
42.89702947649343 0.011159939465770681
42.67180170168193 0.02557751102757614
41.46971205834201 0.012555230726601503
37.08982543385215 0.00107488277468418
36.53555265689563 0.02126940080871224
30.20608114002676 0.009019481779712388
27.374810898415028 0.03994453721421388
26.363055017011728 0.05011922682554448
9.522966890022039 0.00075631870522737
4.016471586720514 0.08168931020200264
3.887605410774121 0.006311654262282101
3.5313112232401513 0.05353548710183397
2.755392921480459 0.021303020776321225
2.4380741812443247 0.23162784335484837
2.2490917637719408 0.2622203718355982
2.079632190634424 0.23382298607799906
2.0336707697261187 0.0623239656843172
1.5673011873879714 0.0367465760261409
1.0188379384951798 0.0126328938653956
1.0022960504442775 0.006817361620021601
0.9970130778462072 0.0103757951973341
0.9201781906386209 0.01095811116040592
0.800604081794174 0.0016830270365341532
0.5738689505255512 0.00646428491253749
#0.022939578929507105 8.355742795827834e-05   # Acoustic modes!
#0.04882611767873102 8.309858592685e-06
#0.07575149723846182 2.778248540373041e-05
]

# Change to SI, but not actually needed as units cancel everywhere

#MAPI_SI = [ MAPI_orig[:,1].*10^12*2*π MAPI_orig[:,2].*1 ]

# OK, black magic here - perhaps our units of oscillator strength are not what we need? maybe already effectively 'squared'?
#MAPI = [ MAPI[:,1] MAPI[:,2].^0.5]

MAPI_low=MAPI[19:33,:] # Just inorganic components, everything below 10THz; modes 3-18


# Hellwarth Table II - BiSiO frequencies
HellwarthII = [
    106.23 8.86
    160.51 9.50
    180.33 20.85
    206.69 10.05
    252.76 27.00
    369.64 61.78
    501.71 52.87
    553.60 86.18
    585.36 75.41
    607.29 98.15
    834.53 89.36
]

# Most simple scheme
# Hellwarth (58), assuming further typo on LHS, actually should be W_e
function HellwarthBscheme(LO)
    println("Hellwarth B Scheme... (athermal)") 
    H58 = sum( (LO[:,2].^2)./ LO[:,1].^2 )
    println("Hellwarth (58) summation: ",H58)

    H59 = sum( LO[:,2].^2 ) # sum of total ir activity squarred
    println("Hellwarth (59) summation (total ir activity ^2): ", H59)
    println("Hellwarth (59) W_e (total ir activity ): ", sqrt(H59))


    omega = sqrt(H59 / H58)
    println("Hellwarth (61) Omega (freq): ",omega)
end

HellwarthBscheme(HellwarthII)
println(" should agree with values given in Hellwarth(60) W_e=196.9 cm^-1 and Hellwarth(61) Ω_e=500 cm^-1")
println("\t MAPI: (all values)")
HellwarthBscheme(MAPI)
println("\t MAPI: (low-frequency, non molecular IR)")
HellwarthBscheme(MAPI_low)


# More complex scheme, involving thermodynamic Beta
# Hellwarth(50), RHS
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 

function HellwarthAscheme(LO,T=295)
    println("Hellwarth A scheme...T=$T K") 
    β=LO[:,1].*2*pi*1E12*ħ/(kB*T) #assuming units of THz
    H50 = sum( ((LO[:,2].^2).*coth.(β))./LO[:,1] )
    println("Hellwarth (50) summation: ",H50)

    H51= sum( LO[:,2].^2 ) # sum of total ir activity squarred
    println("Hellwarth (51) summation (total ir activity ^2): ", H51)
    println("Hellwarth (51) W_e (total ir activity ): ", sqrt(H51))

    # OK; so this is deriving Omega / coth(Beta/2)
    omegacoth=H51/H50
    println("omegacoth: ",omegacoth)

    # NOT FINISHED - need to somehow decouple Omega from both sides of the eqn. 
    for freq in 0.1:0.1:20
        pseudo_omega=omegacoth*coth(freq * 2*pi*1E12*ħ/(2*kB*T))
        if freq>pseudo_omega
            println("freq: $freq pseudo-omega: $pseudo_omega")
            break
        end
    end
end

HellwarthAscheme(HellwarthII)
HellwarthAscheme( [HellwarthII[:,1].*0.02998 HellwarthII[:,2]] ) # convert data to Thz
println(" should agree with values given in Hellwarth\n TableII: H50sum= 91.34 cm^-1, \n W_e=196.9 cm^-1 and Hellwarth(53) Ω_e=504 cm^-1")

println("\t MAPI: (all values)")
HellwarthAscheme(MAPI)
println("\t MAPI: (low-frequency, non molecular IR)")
HellwarthAscheme(MAPI_low)


# Integrate through Lorentz oscillators to get dielectric fn
# Should give 'extra' contribution from these modes, extrapolated to zero omega
function integrate_dielectric(LO,V0)
    summate=sum( (LO[:,2])./(LO[:,1].^2) )
    summate*4*π/V0
end

Å=1E-10 # angstrom in metres
r=6.29Å # Sensible cubic cell size
V0=(r)^3
println("volume: $V0")
const amu=1.66054e-27
const ε0=8.854187817E-12

MAPI_SI = [ MAPI[:,1].*10^12*2*π MAPI[:,2]./(q^2/amu) ]

println(" MAPI: ",integrate_dielectric(MAPI,1.0))
println(" MAPI_low: ",integrate_dielectric(MAPI_low,1.0))
println(" MAPI_SI: ",integrate_dielectric(MAPI_SI,V0))
println(" MAPI_SI: fudged epislon0 ",integrate_dielectric(MAPI_SI,V0)*ε0/(4*π))
println(" MAPI_SI_low: fudged epislon0 ",integrate_dielectric(MAPI_SI[19:33,:],V0)*ε0/(4*π))


println()
#println("From ε_S-ε_Inf, expect this to be: ",ε_S-ε_Inf)

