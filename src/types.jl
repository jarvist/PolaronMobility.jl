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
const ε_0 = 8.85418682e-12
"Speed of light, (ms⁻¹)."
const c = 299792458

# Structures

"Polaron. Structure to store data of polaron solution and other parameters, for each temperature or frequency."
struct Polaron
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

# structure initialisation
Polaron() = Polaron([], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [])

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

function combine_polarons(x::Matrix{NewPolaron})

    Ωlen, Tlen = size(x)
    N = length(x[1,1].v)

    α = x[1,1].α
    T = vcat([x[1, j].T for j in 1:Tlen]...)
    β = hcat([x[1, j].β for j in 1:Tlen]...)
    ω = x[1,1].ω
    v = [x[1, j].v[k] for j in 1:Tlen, k in 1:N]
    w = [x[1, j].w[k] for j in 1:Tlen, k in 1:N]
    κ = [x[1, j].κ[k] for j in 1:Tlen, k in 1:N]
    M = [x[1, j].M[k] for j in 1:Tlen, k in 1:N]
    F = [x[1, j].F for j in 1:Tlen]
    Ω = [x[i, 1].Ω for i in 1:Ωlen]
    Z = [x[i, j].Z for i in 1:Ωlen, j in 1:Tlen]
    σ = [x[i, j].σ for i in 1:Ωlen, j in 1:Tlen]
    μ = [x[1, j].μ for j in 1:Tlen]

    return NewPolaron(α, T, β, ω, v, w, κ, M, F, Ω, Z, σ, μ)
end

function combine_polarons(x::Array{NewPolaron, 3})

    Ωlen, Tlen, αlen = size(x)
    N = length(x[1,1,1].v)

    α = vcat([x[1, 1, k].α for k in 1:αlen]...)
    T = vcat([x[1, j, 1].T for j in 1:Tlen]...)
    β = vcat([x[1, j, 1].β for j in 1:Tlen]...)
    ω = x[1,1,1].ω
    v = [x[1, j, k].v[n] for j in 1:Tlen, k in 1:αlen, n in 1:N]
    w = [x[1, j, k].w[n] for j in 1:Tlen, k in 1:αlen, n in 1:N]
    κ = [x[1, j, k].κ[n] for j in 1:Tlen, k in 1:αlen, n in 1:N]
    M = [x[1, j, k].M[n] for j in 1:Tlen, k in 1:αlen, n in 1:N]
    F = [x[1, j, k].F for j in 1:Tlen, k in 1:αlen]
    Ω = [x[i, 1, 1].Ω for i in 1:Ωlen]
    Z = [x[i, j, k].Z for i in 1:Ωlen, j in 1:Tlen, k in 1:αlen]
    σ = [x[i, j, k].σ for i in 1:Ωlen, j in 1:Tlen, k in 1:αlen]
    μ = [x[i, j, 1].μ for i in 1:Ωlen, j in 1:Tlen]

    return NewPolaron(α, T, β, ω, v, w, κ, M, F, Ω, Z, σ, μ)
end

# Broadcast Polaron data.
function Base.show(io::IO, x::NewPolaron)
    flush(stdout)
    print(io, "-------------------------------------------------\n Polaron Information: \n-------------------------------------------------\n", "Fröhlich coupling | α = ", round.(x.α, digits=3), " | sum(α) = ", round.(sum(x.α), digits=3), "\nTemperatures | T = ", round.(x.T, digits=3), " K \nReduced thermodynamic | β = ", round.(x.β, digits=3), "\nPhonon frequencies | ω = ", round.(x.ω, digits=3), " 2π THz\nVariational parameters | v = ", round.(x.v, digits=3), " ω | w = ", round.(x.w, digits=3), " ω\nFictitious spring constant | κ = ", round.(x.κ, digits=3), " kg/s²\nFictitious mass | M = ", round.(x.M, digits=3), " kg\nFree energy | F = ", round.(x.F, digits=3), " meV\nElectric field frequency | Ω = ", round.(Float64.(x.Ω), digits=3), " 2π THz\nComplex impedance | Z = ", x.Z .|> y -> round.(ComplexF64.(y), digits=3), " V/A\nComplex conductivity | σ = ", x.σ .|> y -> round.(ComplexF64.(y), digits=3), " A/V\nMobility | μ = ", round.(Float64.(x.μ), digits=3), " cm²/Vs")
end

