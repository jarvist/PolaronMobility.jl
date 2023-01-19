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
struct OldPolaron
    "T, temperature (K)."
    T
    "Kμ, Kadanoff mobility (cm²V⁻¹s⁻¹)."
    Kμ
    "Hμ, Hellwarth mobility (cm²V⁻¹s⁻¹)."
    Hμ
    "FHIPμ, FHIP mobility (cm²V⁻¹s⁻¹)."
    FHIPμ
    "k, spring constant."
    k
    "M, renormalised (phonon-drag) mass (mₑ)."
    M
    "Osaka free energy components (A,B,C) and total (F) (unitless). See Hellwarth et al. 1999 PRB Part IV."
    A
    B
    C
    F
    "Tau, relaxation time from Kadanoff Boltzmann transport equation (s)."
    Tau
    "v and w, raw variational parameters (unitless)."
    v
    w
    "βred, reduced thermodynamic beta (unitless)."
    βred
    "rfsi, Feynman polaron radius (Schultz) (m)."
    rfsi
    "rfsmallalpha, small-alpha asymptotic approximation (unitless)."
    rfsmallalpha

    # Setup of simulation. These parameters are sent to the function.

    "α, Fröhlich alpha coupling parameter (unitless)."
    α
    "mb, Band effective mass (mₑ)."
    mb
    "ω, effective dielectric frequency (2π THz)."
    ω
end

struct Polaron
    α       # Fröhlich coupling, unitless
    αeff    # Effective Fröhlich coupling summed for multiple modes, unitless
    T       # Temperature, K
    ω       # Phonon mode frequencies, 2π⋅THz
    # mb      # Particle mass
    β       # Reduced thermodynamic beta ħω₀/kBT, unitless
    Ω       # Electric field frequency, 2π⋅THz
    v0       # Variational parameter v, unitless
    w0       # Variational parameter w, unitless
    F0       # Polaron free energy, meV
    A0       # Bare electron free energy, meV
    B0       # ⟨S⟩ₜ interaction energy, meV
    C0       # ⟨Sₜ⟩ₜ free energy of trial system, meV
    v       # Variational parameter v, unitless
    w       # Variational parameter w, unitless
    F       # Polaron free energy, meV
    A       # Bare electron free energy, meV
    B       # ⟨S⟩ₜ interaction energy, meV
    C       # ⟨Sₜ⟩ₜ free energy of trial system, meV
    # Fw      # Weak coupling energy approximation, meV
    # Fs      # Strong coupling energy approximation, meV
    κ       # Fictitious spring constant, unitless (multiples of mₑω₀², electron band-mass and (2π⋅THz)²)
    M       # Fictitious mass, unitless (multiples of mₑ electron band-mass)
    R       # Schultz polaron radius, unitless
    z       # Complex impedence, V/A
    σ       # Complex conductivity, A/V
    μ       # DC mobility, cm²/Vs
    # μK      # Kadanoff DC mobility (cm²/Vs)
    # μH      # Hellwarth DC mobility (cm²/Vs)
    # τ       # relaxation time from Kadanoff Boltzmann transport equation

    function Polaron(x...)
        reduce_array(a) = length(a) == 1 ? only(a) : dropdims(a, dims = tuple(findall(size(a) .== 1)...))
        new(reduce_array.(x)...)
    end
end 

# structure initialisation
OldPolaron() = OldPolaron([], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [])

# Broadcast Polaron data.
function Base.show(io::IO, x::Polaron)
    flush(stdout)
    IOContext(stdout, :limit=>true, :displaysize=>(100,80))
    println("\033[K-----------------------------------------------------------------------")
    println("\033[K                         Polaron Information:                          ")
    println("\033[K-----------------------------------------------------------------------") 

    println(IOContext(stdout, :limit => true, :compact => true), "\033[KPhonon frequencies         | ω = ", x.ω, " 2π THz")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KFröhlich coupling          | α = ", x.α, " | sum(α) = ", x.αeff) 

    println("\033[K-----------------------------------------------------------------------") 
    println("\033[K                       Ground State Information:                       ")
    println("\033[K-----------------------------------------------------------------------") 

    println(IOContext(stdout, :limit => true, :compact => true), "\033[KGS variational parameter   | v₀ = ", x.v0, " ω")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KGS variational parameter   | w₀ = ", x.w0, " ω")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KGS Energy                  | E₀ = ", x.F0, " meV")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KGS Electron energy         | A₀ = ", x.A0, " meV")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KGS Interaction energy      | B₀ = ", x.B0, " meV")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KGS Trial energy            | C₀ = ", x.C0, " meV")

    println("\033[K-----------------------------------------------------------------------") 
    println("\033[K                    Finite Temperature Information:                    ")
    println("\033[K-----------------------------------------------------------------------") 

    println(IOContext(stdout, :limit => true, :compact => true), "\033[KTemperatures               | T = ", x.T, " K")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KReduced thermodynamic      | β = ", x.β)
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KVariational parameter      | v = ", x.v, " ω")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KVariational parameter      | w = ", x.w, " ω")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KFree energy                | F = ", x.F, " meV")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KElectron energy            | A = ", x.A, " meV")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KInteraction energy         | B = ", x.B, " meV")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KTrial energy               | C = ", x.C, " meV")

    println("\033[K-----------------------------------------------------------------------") 
    println("\033[K                       Trial System Information:                       ")
    println("\033[K-----------------------------------------------------------------------") 

    println(IOContext(stdout, :limit => true, :compact => true), "\033[KFictitious spring constant | κ = ", x.κ, " kg/s²")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KFictitious mass            | M = ", x.M, " mₑ")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KPolaron radius             | R = ", x.R, " rₚ")

    println("\033[K-----------------------------------------------------------------------") 
    println("\033[K                     Linear Reponse Information:                       ")
    println("\033[K-----------------------------------------------------------------------") 

    println(IOContext(stdout, :limit => true, :compact => true), "\033[KElectric field frequency   | Ω = ", x.Ω, " 2π THz")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KComplex impedance          | z = ", x.z, " V/A")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KComplex conductivity       | σ = ", x.σ, " A/V")
    println(IOContext(stdout, :limit => true, :compact => true), "\033[KMobility                   | μ = ", x.μ, " cm²/Vs")

    println("\033[K-----------------------------------------------------------------------") 
end

