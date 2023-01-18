# types.jl

# Physical Constants

"Planck's constant, (kgm²s⁻¹)."
const hbar = const ħ = 1.054571817e-34
"Electron charge, (kgm²s⁻²)."
const eV = const q = const ElectronVolt = 1.602176634e-19
"Electron mass, (kg)."
const me = MassElectron = 9.1093837015e-31
"Boltzmann's constant, (kgm²K⁻¹)."
const Boltzmann = const kB = 1.380649e-23
"Permittivity of free space, (C²N⁻¹m⁻²)."
const ε_0 = ϵ_0 = 8.85418682e-12
"Speed of light, (ms⁻¹)."
const c = 299792458
"Atomic mass unit, (kg)"
const amu = 1.660_539_066_60e-27

# Structures

"Polaron. Structure to store data of polaron solution and other parameters, for each temperature or frequency."
# struct Polaron
#     "T, temperature (K)."
#     T
#     "Kμ, Kadanoff mobility (cm²V⁻¹s⁻¹)."
#     Kμ
#     "Hμ, Hellwarth mobility (cm²V⁻¹s⁻¹)."
#     Hμ
#     "FHIPμ, FHIP mobility (cm²V⁻¹s⁻¹)."
#     FHIPμ
#     "k, spring constant."
#     k
#     "M, renormalised (phonon-drag) mass (mₑ)."
#     M
#     "Osaka free energy components (A,B,C) and total (F) (unitless). See Hellwarth et al. 1999 PRB Part IV."
#     A
#     B
#     C
#     F
#     "Tau, relaxation time from Kadanoff Boltzmann transport equation (s)."
#     Tau
#     "v and w, raw variational parameters (unitless)."
#     v
#     w
#     "βred, reduced thermodynamic beta (unitless)."
#     βred
#     "rfsi, Feynman polaron radius (Schultz) (m)."
#     rfsi
#     "rfsmallalpha, small-alpha asymptotic approximation (unitless)."
#     rfsmallalpha

#     # Setup of simulation. These parameters are sent to the function.

#     "α, Fröhlich alpha coupling parameter (unitless)."
#     α
#     "mb, Band effective mass (mₑ)."
#     mb
#     "ω, effective dielectric frequency (2π THz)."
#     ω
# end

struct Polaron
    α
    αsum
    T
    ω
    β
    Ω
    v
    w
    F
    κ
    M
    z
    σ
    μ
end

# structure initialisation
# Polaron() = Polaron([], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [])

"NewPolaron. Structure to store data of polaron solution and other parameters, for each temperature or frequency."
struct NewPolaron
    "α, Fröhlich alpha coupling parameter (unitless)"
    α
    "T, temperature (K)."
    T
    "β, reduced Thermodynamic beta (unitless)."
    β
    "ω, phonon frequency (2π THz)."
    ω
    "v, variational parameters (unitless)."
    v
    "w, variational parameters (unitless)."
    w
    "κ, fictitious spring constant (multiples of mₑ) (kgs⁻²)."
    κ
    "M, fictitious particle (multiples of mₑ) (kg)."
    M
    "F, free energy (meV)."
    F
    "Ω, electric field frequencies (multiples of phonon frequency ω) (2π THz)."
    Ω
    "Z, complex impedence (V/A)."
    Z
    "σ, complex conductivity (A/V)."
    σ
    "μ, FHIP mobility (cm²V⁻¹s⁻¹)."
    μ
end

# Broadcast Polaron data.
function Base.show(io::IO, x::Polaron)
    flush(stdout)
    println("-------------------------------------------------")
    println(" Polaron Information: ")
    println("-------------------------------------------------") 
    println(IOContext(stdout, :limit => true), "Fröhlich coupling | α = ", round.(x.α, digits=3), " | sum(α) = ", round.(x.αsum, digits=3)) 
    println(IOContext(stdout, :limit => true), "Temperatures | T = ", round.(x.T, digits=3), " K")
    println(IOContext(stdout, :limit => true), "Phonon frequencies | ω = ", round.(x.ω, digits=3), " 2π THz")
    println(IOContext(stdout, :limit => true), "Reduced thermodynamic | β = ", round.(x.β, digits=3))
    println(IOContext(stdout, :limit => true), "Variational parameters | v = ", round.(x.v, digits=3), " ω | w = ", round.(x.w, digits=3), " ω")
    println(IOContext(stdout, :limit => true), "Fictitious spring constant | κ = ", round.(x.κ, digits=3), " kg/s²")
    println(IOContext(stdout, :limit => true), "Fictitious mass | M = ", round.(x.M, digits=3), " kg")
    println(IOContext(stdout, :limit => true), "Free energy | F = ", round.(x.F, digits=3), " meV")
    println(IOContext(stdout, :limit => true), "Electric field frequency | Ω = ", round.(Float64.(x.Ω), digits=3), " 2π THz")
    println(IOContext(stdout, :limit => true), "Complex impedance | z = ", x.z .|> y -> round.(ComplexF64.(y), digits=3), " V/A")
    println(IOContext(stdout, :limit => true), "Complex conductivity | σ = ", x.σ .|> y -> round.(ComplexF64.(y), digits=3), " A/V")
    println(IOContext(stdout, :limit => true), "Mobility | μ = ", round.(Float64.(x.μ), digits=3), " cm²/Vs")
    println("-------------------------------------------------")
end

