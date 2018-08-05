# types.jl

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space


# Structure to store data of polaron solution + other parameters, for each temperature
struct Polaron
    T
    # Mobilities
    Kμ; Hμ; FHIPμ
    # Spring constant and renormalised (phonon-drag) mass
    k; M
    # Osaka free energy components (A,B,C) and total (F). See Hellwarth et al. 1999 PRB Part IV
    A; B; C; F
    # Relaxation time from Kadanoff Boltzmann transport equation
    Tau
    # Raw variational parameters
    v; w
    # Reduced thermodynamic beta
    βred
    # Feynman polaron radius (Schultz), in SI units. Then also the small-alpha asymptotic approx
    rfsi; rfsmallalpha
    # Setup of simulation. These parameters are sent to the function.
    # Alpha = Frohlich alpha
    α
    # Band effective mass
    mb
    # Effective dielectric frequency
    ω
end
Polaron()=Polaron([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]) # structure initialisation

