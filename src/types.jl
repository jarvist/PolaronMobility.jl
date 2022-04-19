# types.jl

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space
const c = 3e8


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

struct NewPolaron
    α      # Frohlich alpha (unitless)
    T      # Temperature (K)
    β      # Reduced Thermodynamic beta (unitless)
    ω      # Phonon frequency (rad THz)
    v      # Variational parameter (s^-1)
    w      # Variational parameter (s^-1)
    κ      # Fictitious spring constant (multiples of m_e) (kg / s^2)
    M      # Fictitious particle (multiples of m_e) (kg)
    F      # Free energy (meV)
    Ω      # Electric field frequencies (multiples of phonon frequency ω) (s^-1)
    Z      # Complex impedence (V/A)
    σ      # Complex conductivity (A/V)
    μ      # Mobility (cm^2/Vs)
end

# Broadcast Polaron data.
function Base.show(io::IO, x::NewPolaron)
    flush(stdout)
    print(io, "-------------------------------------------------\n Polaron Information: \n-------------------------------------------------\n", "Fröhlich coupling | α = ", round.(x.α, digits = 3), " | sum(α) = ", round(sum(x.α), digits = 3),"\nTemperatures | T = ", round.(x.T, digits = 3), " K \nReduced thermodynamic | β = ", round.(x.β, digits = 3), "\nPhonon frequencies | ω = ", round.(x.ω, digits = 3), " 2π THz\nVariational parameters | v = ", round.(x.v, digits = 3), " ω | w = ", round.(x.w, digits = 3), " ω\nFictitious spring constant | κ = ", round.(x.κ, digits = 3), " kg/s²\nFictitious mass | M = ", round.(x.M, digits = 3), " kg\nFree energy | F = ", round.(x.F, digits = 3), " meV\nElectric field frequency | Ω = ", round.(Float64.(x.Ω), digits = 3),  " 2π THz\nComplex impedance | Z = ", x.Z .|> y -> round.(ComplexF64.(y), digits = 3), " V/A\nComplex conductivity | σ = ", x.σ .|> y -> round.(ComplexF64.(y), digits = 3), " A/V\nMobility | μ = ", round.(Float64.(x.μ), digits = 3), " cm²/Vs")
end
