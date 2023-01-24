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

reduce_array(a) = length(a) == 1 ? only(a) : dropdims(a, dims=tuple(findall(size(a) .== 1)...))

struct Material
    optical # Optical dielectric constant
    static  # Static dielectric constant
    ionic   # Ionic dielectric contributions
    mb      # Effective band mass
    α       # Fröhlich coupling, unitless
    freqs   # Phonon frequencies, THz
    ir      # Infrared activities
    volume  # Unit cell volumes
    function Material(x...)
        new(reduce_array.(x)...)
    end
end

function material(ϵ_optical, ϵ_static, m_eff, phonon_freq)
    ϵ_ionic = ϵ_static - ϵ_optical
    α = frohlichalpha(ϵ_optical, ϵ_static, phonon_freq, m_eff)
    return Material(ϵ_optical, ϵ_static, ϵ_ionic, m_eff, α, phonon_freq, 1, 1)
end

function material(ϵ_optic, ϵ_static, m_eff, phonon_freqs, ir_activity, volume)
    ϵ_ionic = ϵ_ionic_mode.(phonon_freqs, ir_activity, volume)
    α = frohlichalpha.(ϵ_optic, ϵ_ionic, sum(ϵ_ionic), phonon_freqs, m_eff)
    return Material(ϵ_optic, ϵ_static, ϵ_ionic, m_eff, α, phonon_freqs, ir_activity, volume)
end

function Base.show(io::IO, ::MIME"text/plain", x::Material)
    flush(io)
    println("\e[K------------------------------------------")
    println("\e[K           Material Information           ")
    println("\e[K------------------------------------------")
    println(IOContext(io, :compact => true, :limit => true), "\e[KOptic dielectric   | ϵ∞ = ", x.optical, " ")
    println(IOContext(io, :compact => true, :limit => true), "\e[KStatic dielectric  | ϵ0 = ", x.static, " ")
    println(IOContext(io, :compact => true, :limit => true), "\e[KIonic dielectric   | ϵᵢ = ", x.ionic, " ")
    println(IOContext(io, :compact => true, :limit => true), "\e[KBand mass          | mb = ", x.mb, " mₑ")
    println(IOContext(io, :compact => true, :limit => true), "\e[KFröhlich coupling  | α = ", x.α)
    println(IOContext(io, :compact => true, :limit => true), "\e[KPhonon frequencies | f = ", x.freqs, " THz")
    println(IOContext(io, :compact => true, :limit => true), "\e[KIR activities      | IR = ", x.ir, " ")
    println(IOContext(io, :compact => true, :limit => true), "\e[KUnit cell volume   | V₀ = ", x.volume, " m³")
    println("\e[K-------------------------------------------")
end

struct Polaron
    α       # Fröhlich coupling, unitless
    αeff    # Effective Fröhlich coupling summed for multiple modes, unitless
    T       # Temperature, K
    ω       # Phonon mode frequencies, 2π⋅THz
    β       # Reduced thermodynamic beta ħω₀/kBT, unitless
    Ω       # Electric field frequency, 2π⋅THz
    v0      # Variational parameter v, unitless
    w0      # Variational parameter w, unitless
    F0      # Polaron free energy, meV
    A0      # Bare electron free energy, meV
    B0      # ⟨S⟩ₜ interaction energy, meV
    C0      # ⟨Sₜ⟩ₜ free energy of trial system, meV
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
        new(reduce_array.(x)...)
    end
end

# structure initialisation
OldPolaron() = OldPolaron([], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [])

# Broadcast Polaron data.
function Base.show(io::IO, ::MIME"text/plain", x::Polaron)
    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                         Polaron Information:                          ")
    println("\e[K-----------------------------------------------------------------------")

    println(IOContext(io, :limit => true, :compact => true), "\e[KPhonon frequencies         | ω = ", x.ω, " 2π THz")
    println(IOContext(io, :limit => true, :compact => true), "\e[KFröhlich coupling          | α = ", x.α, " | sum(α) = ", x.αeff)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                       Ground State Information:                       ")
    println("\e[K-----------------------------------------------------------------------")

    println(IOContext(io, :limit => true, :compact => true), "\e[KGS variational parameter   | v₀ = ", x.v0, " ω")
    println(IOContext(io, :limit => true, :compact => true), "\e[KGS variational parameter   | w₀ = ", x.w0, " ω")
    println(IOContext(io, :limit => true, :compact => true), "\e[KGS Energy                  | E₀ = ", x.F0, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KGS Electron energy         | A₀ = ", x.A0, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KGS Interaction energy      | B₀ = ", x.B0, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KGS Trial energy            | C₀ = ", x.C0, " meV")

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                    Finite Temperature Information:                    ")
    println("\e[K-----------------------------------------------------------------------")

    println(IOContext(io, :limit => true, :compact => true), "\e[KTemperatures               | T = ", x.T, " K")
    println(IOContext(io, :limit => true, :compact => true), "\e[KReduced thermodynamic      | β = ", x.β)
    println(IOContext(io, :limit => true, :compact => true), "\e[KVariational parameter      | v = ", x.v, " ω")
    println(IOContext(io, :limit => true, :compact => true), "\e[KVariational parameter      | w = ", x.w, " ω")
    println(IOContext(io, :limit => true, :compact => true), "\e[KFree energy                | F = ", x.F, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KElectron energy            | A = ", x.A, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KInteraction energy         | B = ", x.B, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KTrial energy               | C = ", x.C, " meV")

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                       Trial System Information:                       ")
    println("\e[K-----------------------------------------------------------------------")

    println(IOContext(io, :limit => true, :compact => true), "\e[KFictitious spring constant | κ = ", x.κ, " kg/s²")
    println(IOContext(io, :limit => true, :compact => true), "\e[KFictitious mass            | M = ", x.M, " mₑ")
    println(IOContext(io, :limit => true, :compact => true), "\e[KPolaron radius             | R = ", x.R, " rₚ")

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                     Linear Reponse Information:                       ")
    println("\e[K-----------------------------------------------------------------------")

    println(IOContext(io, :limit => true, :compact => true), "\e[KElectric field frequency   | Ω = ", x.Ω, " 2π THz")
    println(IOContext(io, :limit => true, :compact => true), "\e[KComplex impedance          | z = ", x.z, " V/A")
    println(IOContext(io, :limit => true, :compact => true), "\e[KComplex conductivity       | σ = ", x.σ, " A/V")
    println(IOContext(io, :limit => true, :compact => true), "\e[KMobility                   | μ = ", x.μ, " cm²/Vs")

    println("\e[K-----------------------------------------------------------------------")
end

