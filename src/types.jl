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
    v      # Variational parameter (s^-1)
    w      # Variational parameter (s^-1)
    κ      # Fictitious spring constant (multiples of m_e) (kg / s^2)
    M      # Fictitious particle (multiples of m_e) (kg)
    F      # Free energy (meV)
    Ω      # Electric field frequencies (multiples of phonon frequency ω) (s^-1)
    Z      # Complex impedence (cm^2/Vs)
    σ      # Complex conductivity (cm^-1)
    μ      # Mobility (cm^2/Vs)
end

# Broadcast Polaron data.
function Base.show(io::IO, x::NewPolaron)
    flush(stdout)
    print(io, "---------------------- \n Polaron Information: \n----------------------\n", "α = ", round(x.α, digits = 3), "\nT = ", round.(x.T, digits = 3), " K \nβ = ", round.(x.β, digits = 3), "\nv = ", round.(x.v, digits = 3), " s^-1\nw = ", round.(x.w, digits = 3), " s^-1\nκ = ", round.(x.κ, digits = 3), " kg/s^2\nM = ", round.(x.M, digits = 3), " kg\nF = ", round.(x.F, digits = 3), " meV\nΩ = ", round.(Float64.(x.Ω), digits = 3),  " s^-1\nZ = ", x.Z .|> y -> round.(ComplexF64.(y), digits = 3), " cm^2/Vs\nσ = ", x.σ .|> y -> round.(ComplexF64.(y), digits = 3), " cm^-1\nμ = ", round.(Float64.(x.μ), digits = 3))
end