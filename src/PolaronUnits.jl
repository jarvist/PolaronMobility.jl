# The following five constants are used as the “base” atomic units

# Physical constants

"Planck's constant, (kgm²s⁻¹)."
const hbar = ħ = 1.054571817e-34
const h = ħ * 2π
"Electron charge, (kgm²s⁻²)."
const eV = q = ElectronVolt = 1.602176634e-19
"Electron mass, (kg)."
const me = MassElectron = 9.1093837015e-31
"Boltzmann's constant, (kgm²K⁻¹)."
const Boltzmann = kB = 1.380649e-23
"Permittivity of free space, (C²N⁻¹m⁻²)."
const ε_0 = ϵ_0 = 8.85418682e-12
"Speed of light, (ms⁻¹)."
const c = 299792458
"Atomic mass unit, (kg)"
const amu = 1.660_539_066_60e-27

"""
    PolaronUnits.m0_pu
A unit equal to the electron rest mass mₑ ≈ 9.109,383,7015 × 10^-31 kg.
Printed as "mₑ".
`Unitful.me` is a quantity (with units `kg`) whereas `PolaronUnits.m0_pu` is a unit equal to
`Unitful.me`.
Dimension: [`Unitful.𝐌`](@ref).
See also: `Unitful.me`, `Unitful.kg`.
"""
@unit m0_pu "m₀"  PolaronMassScale      Unitful.me                     false

"""
    PolaronUnits.e_pu
A unit equal to the elementary charge e = 1.602,176,634 × 10^-19 C.
Printed as "e".
`Unitful.q` is a quantity (with units `C`) whereas `PolaronUnits.e_pu` is a unit equal to
`Unitful.q`.
Dimension: 𝐈 𝐓.
See also: `Unitful.q`, `Unitful.C`.
"""
@unit e_pu  "e"   ElementaryCharge      Unitful.q                      false

"""
    PolaronUnits.ħ_pu
A unit equal to the reduced Planck constant ħ = h / 2π ≈ 1.054,571,8176 × 10^-34 J × s.
Printed as "ħ".
`Unitful.ħ` is a quantity (with units `J × s`) whereas `PolaronUnits.ħ_pu` is a unit equal to
`Unitful.ħ`.
Dimension: 𝐋^2 𝐌 𝐓^-1.
See also: `Unitful.ħ`, `Unitful.J`, `Unitful.s`.
"""
@unit ħ_pu "ħ" ReducedPlanckConstant 1Unitful.ħ false

"""
    PolaronUnits.k_pu
A unit equal to the Boltzmann constant k = 1.380,649 × 10^-23 J / K.
Printed as "k".
`Unitful.k` is a quantity (with units `J / K`) whereas `PolaronUnits.k_pu` is a unit equal to
`Unitful.k`.
Dimension: 𝐋^2 𝐌 𝚯^-1 𝐓^-2.
See also: `Unitful.k`, `Unitful.J`, `Unitful.K`.
"""
@unit k_pu "k" BoltzmannConstant 1Unitful.k false

"""
    PolaronUnits.ω0_pu
A unit equal to the reduced Planck constant ħ = h / 2π ≈ 1.054,571,8176 × 10^-34 J × s.
Printed as "ħ".
`Unitful.ħ` is a quantity (with units `J × s`) whereas `PolaronUnits.ħ_pu` is a unit equal to
`Unitful.ħ`.
Dimension: 𝐋^2 𝐌 𝐓^-1.
See also: `Unitful.ħ`, `Unitful.J`, `Unitful.s`.
"""
@unit ω0_pu "ω₀" PolaronAngularFrequency 1Unitful.THz2π false

# Polaron radius is derived from the base polaron units
"""
    PolaronUnits.a0_pu
A unit equal to the characteristic polaron radius
```
a₀ = sqrt(ħ / m₀ / ω₀)
```
Printed as "a₀".
Dimension: 𝐋.
See also: `Unitful.ε0`, `Unitful.ħ`, `Unitful.me`, `Unitful.q`, `Unitful.m`.
"""
@unit a0_pu "a₀" PolaronRadius √(1ħ_pu/1m0_pu/1ω0_pu) false

# Polaron energy is derived from the base polaron units
"""
    PolaronUnits.E0_pu
A unit equal to the phonon energy
```
E₀ = ħω₀
```
Printed as "E₀".
Dimension: 𝐋^2 𝐌 𝐓^-2.
See also: `Unitful.me`, `Unitful.q`, `Unitful.ε0`, `Unitful.ħ`, `Unitful.J`, `Unitful.eV`, [`UnitfulAtomic.Ry`](@ref).
"""
@unit E0_pu "E₀" PolaronEnergy 1ħ_pu*1ω0_pu false

"""
    punit(x::Unitful.Quantity)
    punit(x::Unitful.Units)
    punit(x::Unitful.Dimensions)
Returns the appropriate polaron unit (a `Unitful.Units` object) for the dimension of `x`.
"""
punit(x) = punit(dimension(x))

# `punit` for `Dimension` types
punit(x::Dimension{:Length})      = (a0_pu)^x.power
punit(x::Dimension{:Mass})        = (m0_pu)^x.power
punit(x::Dimension{:Time})        = (ħ_pu/E0_pu)^x.power
punit(x::Dimension{:Current})     = (e_pu*E0_pu/ħ_pu)^x.power
punit(x::Dimension{:Temperature}) = (E0_pu/k_pu)^x.power

# For dimensions not specified above, there is no polaron unit.
punit(::Dimension{D}) where D = throw(ArgumentError("No polaron unit defined for dimension $D."))

# `punit` for `Dimensions` types
@generated punit(::Dimensions{N}) where N = prod(punit, N)
punit(::typeof(NoDims)) = NoUnits

# Simplifications for some derived dimensions, so that e.g. `punit(u"J")` returns `E₀`
# instead of `a₀^2 mₑ E₀^2 ħ^-2`. The following units/dimensions are considered:
#   * Energy: E₀
#   * Momentum: ħ/a₀
#   * Action/angular momentum: ħ
#   * Force: E₀/a₀
#   * E-field: E₀/(e*a₀)
#   * B-field: ħ/(e*a₀^2)
#   * Voltage/electric potential: E₀/e
for unit in (:(E0_pu), :(e_pu), :(ω0_pu), :(m0_pu), :(ħ_pu/a0_pu), :(ħ_pu), :(E0_pu/a0_pu),
    :(E0_pu/(e_pu*a0_pu)), :(ħ_pu/(e_pu*a0_pu^2)), :(E0_pu/e_pu), :(e_pu^2/(a0_pu*E0_pu)))
    @eval punit(::typeof(dimension($unit))) = $unit
end

"""
    puconvert(x::Unitful.Quantity)
Convert a quantity to the appropriate polaron unit.
"""
puconvert(x) = uconvert(punit(x), x)

"""
    puconvert(u::Unitful.Units, x::Number)
Interpret `x` as a quantity given in polaron units and convert it to the unit `u`.
"""
puconvert(u::Units, x::Number) = uconvert(u, x*punit(u))

"""
    pustrip(x::Unitful.Quantity)
Returns the value of the quantity converted to polaron units as a number type (i.e., with the
units removed). This is equivalent to `Unitful.ustrip(puconvert(x))`.
"""
pustrip(x) = ustrip(puconvert(x))

"""
    addunits!(polaron::Polaron)
"""
function addunits!(polaron::FrohlichPolaron; unit="pu")
    polaron.Fs = pustrip.(polaron.Fs .* polaron.ωeff) .* punit(1Unitful.meV)
    polaron.Fl = pustrip.(polaron.Fl .* polaron.ωeff) .* punit(1Unitful.meV)
    polaron.Ms = pustrip.(polaron.Ms .* polaron.mb) .* punit(1Unitful.me)
    polaron.Ml = pustrip.(polaron.Ml .* polaron.mb) .* punit(1Unitful.me)
    polaron.Rs = pustrip.(polaron.Rs ./ sqrt(2 * polaron.ωeff * polaron.mb)) .* punit(1Unitful.Å)
    polaron.Rl = pustrip.(polaron.Rl ./ sqrt(2 * polaron.ωeff * polaron.mb)) .* punit(1Unitful.Å)
    polaron.ΩFC = pustrip.(polaron.ΩFC) .* punit(1Unitful.THz2π)
    polaron.F0 = pustrip.(polaron.F0) .* punit(1Unitful.meV)
    polaron.A0 = pustrip.(polaron.A0) .* punit(1Unitful.meV)
    polaron.B0 = pustrip.(polaron.B0) .* punit(1Unitful.meV)
    polaron.C0 = pustrip.(polaron.C0) .* punit(1Unitful.meV)
    polaron.v0 = pustrip.(polaron.v0) .* punit(1Unitful.THz2π)
    polaron.w0 = pustrip.(polaron.w0) .* punit(1Unitful.THz2π)
    polaron.κ0 = pustrip.(polaron.κ0 .* polaron.mb) .* punit(1u"μN/m")
    polaron.M0 = pustrip.(polaron.M0 .* polaron.mb) .* punit(1Unitful.me)
    polaron.M0a = pustrip.(polaron.M0a .* polaron.mb) .* punit(1Unitful.me)
    polaron.M0r = pustrip.(polaron.M0r .* polaron.mb) .* punit(1Unitful.me)
    polaron.R0 = pustrip.(polaron.R0 ./ sqrt(2 * polaron.ωeff * polaron.mb)) .* punit(1Unitful.Å)
    polaron.F = pustrip.(polaron.F) .* punit(1Unitful.meV)
    polaron.A = pustrip.(polaron.A) .* punit(1Unitful.meV)
    polaron.B = pustrip.(polaron.B) .* punit(1Unitful.meV)
    polaron.C = pustrip.(polaron.C) .* punit(1Unitful.meV)
    polaron.v = pustrip.(polaron.v) .* punit(1Unitful.THz2π)
    polaron.w = pustrip.(polaron.w) .* punit(1Unitful.THz2π)
    polaron.κ = pustrip.(polaron.κ .* polaron.mb) .* punit(1u"N/m")
    polaron.M = pustrip.(polaron.M .* polaron.mb) .* punit(1Unitful.me)
    polaron.Ma = pustrip.(polaron.Ma .* polaron.mb) .* punit(1Unitful.me)
    polaron.Mr = pustrip.(polaron.Mr .* polaron.mb) .* punit(1Unitful.me)
    polaron.R = pustrip.(polaron.R ./ sqrt(2 * polaron.ωeff * polaron.mb)) .* punit(1Unitful.Å)
    polaron.μ = pustrip.(polaron.μ ./ polaron.mb) .* punit(1u"cm^2/V/s")
    polaron.μFHIP = pustrip.(polaron.μFHIP ./ polaron.mb) .* punit(1u"cm^2/V/s")
    polaron.μD = pustrip.(polaron.μD ./ polaron.mb) .* punit(1u"cm^2/V/s")
    polaron.μK = pustrip.(polaron.μK ./ polaron.mb) .* punit(1u"cm^2/V/s")
    polaron.μH = pustrip.(polaron.μH ./ polaron.mb) .* punit(1u"cm^2/V/s")
    polaron.μH0 = pustrip.(polaron.μH0 ./ polaron.mb) .* punit(1u"cm^2/V/s")
    polaron.τ = pustrip.(polaron.τ) .* punit(1u"ns")
    polaron.χ = pustrip.(polaron.χ .* polaron.mb) .* punit(u"Ω") .* punit(1Unitful.THz2π)
    polaron.z = pustrip.(polaron.z .* polaron.mb) .* punit(u"Ω")
    polaron.σ = pustrip.(polaron.σ ./ polaron.mb) .* punit(u"S")
    polaron.T = pustrip.(polaron.T) .* punit(1Unitful.K)
    polaron.β = pustrip.(polaron.β ./ Unitful.ħ) .* punit(1 / 1Unitful.meV)
    polaron.Ω = pustrip.(polaron.Ω) .* punit(1Unitful.THz2π)
    polaron.ω = pustrip.(polaron.ω) .* punit(1Unitful.THz2π)
    polaron.ωeff = pustrip.(polaron.ωeff) .* punit(1Unitful.THz2π)
    if unit == "su"
        suconvert!(polaron)
    end
end

function addunits!(polaron::HolsteinPolaron; unit="pu")
    polaron.F0 = pustrip.(polaron.F0) .* punit(1Unitful.meV)
    polaron.A0 = pustrip.(polaron.A0) .* punit(1Unitful.meV)
    polaron.B0 = pustrip.(polaron.B0) .* punit(1Unitful.meV)
    polaron.C0 = pustrip.(polaron.C0) .* punit(1Unitful.meV)
    polaron.v0 = pustrip.(polaron.v0) .* punit(1Unitful.THz2π)
    polaron.w0 = pustrip.(polaron.w0) .* punit(1Unitful.THz2π)
    polaron.κ0 = pustrip.(polaron.κ0 .* polaron.mb) .* punit(1u"N/m")
    polaron.M0 = pustrip.(polaron.M0 .* polaron.mb) .* punit(1Unitful.me)
    polaron.M0a = pustrip.(polaron.M0a .* polaron.mb) .* punit(1Unitful.me)
    polaron.M0r = pustrip.(polaron.M0r .* polaron.mb) .* punit(1Unitful.me)
    polaron.R0 = pustrip.(polaron.R0 ./ sqrt(2 * polaron.ωeff * polaron.mb)) .* punit(1Unitful.Å)
    polaron.F = pustrip.(polaron.F) .* punit(1Unitful.meV)
    polaron.A = pustrip.(polaron.A) .* punit(1Unitful.meV)
    polaron.B = pustrip.(polaron.B) .* punit(1Unitful.meV)
    polaron.C = pustrip.(polaron.C) .* punit(1Unitful.meV)
    polaron.v = pustrip.(polaron.v) .* punit(1Unitful.THz2π)
    polaron.w = pustrip.(polaron.w) .* punit(1Unitful.THz2π)
    polaron.κ = pustrip.(polaron.κ .* polaron.mb) .* punit(1u"N/m")
    polaron.M = pustrip.(polaron.M .* polaron.mb) .* punit(1Unitful.me)
    polaron.Ma = pustrip.(polaron.Ma .* polaron.mb) .* punit(1Unitful.me)
    polaron.Mr = pustrip.(polaron.Mr .* polaron.mb) .* punit(1Unitful.me)
    polaron.R = pustrip.(polaron.R ./ sqrt(2 * polaron.ωeff * polaron.mb)) .* punit(1Unitful.Å)
    polaron.μ = pustrip.(polaron.μ ./ polaron.mb) .* punit(1u"cm^2/V/s")
    polaron.χ = pustrip.(polaron.χ .* polaron.mb) .* punit(u"Ω") .* punit(1Unitful.THz2π)
    polaron.z = pustrip.(polaron.z .* polaron.mb) .* punit(u"Ω")
    polaron.σ = pustrip.(polaron.σ ./ polaron.mb) .* punit(u"S")
    polaron.T = pustrip.(polaron.T) .* punit(1Unitful.K)
    polaron.β = pustrip.(polaron.β ./ Unitful.ħ) .* punit(1 / 1Unitful.meV)
    polaron.Ω = pustrip.(polaron.Ω) .* punit(1Unitful.THz2π)
    polaron.ω = pustrip.(polaron.ω) .* punit(1Unitful.THz2π)
    polaron.ωeff = pustrip.(polaron.ωeff) .* punit(1Unitful.THz2π)
    polaron.J = pustrip.(polaron.J) .* punit(1Unitful.meV) 
    polaron.a = pustrip.(polaron.a) .* punit(1Unitful.Å) 
    if unit == "su"
        suconvert!(polaron)
    end
end

"""
    addunits!(polaron::Material)
"""
function addunits!(material::Material)
    material.feff = material.feff .* ω0_pu
    material.f = material.f .* ω0_pu
    material.ϵi = material.ϵi .* punit(u"ϵ0")
    material.mb = material.mb .* m0_pu
    material.ϵo = material.ϵo .* punit(u"ϵ0")
    material.ϵs = material.ϵs .* punit(u"ϵ0")
end


function suconvert!(polaron::FrohlichPolaron)
    polaron.ω = polaron.ω .|> Unitful.THz2π
    polaron.ωeff = polaron.ωeff .|> Unitful.THz2π
    polaron.Fs = polaron.Fs .|> Unitful.meV
    polaron.Fl = polaron.Fl .|> Unitful.meV
    polaron.Ms = polaron.Ms .|> Unitful.kg
    polaron.Ml = polaron.Ml .|> Unitful.kg
    polaron.Rs = polaron.Rs .|> Unitful.Å
    polaron.Rl = polaron.Rl .|> Unitful.Å
    polaron.ΩFC = polaron.ΩFC .|> Unitful.THz2π
    polaron.F0 = polaron.F0 .|> Unitful.meV
    polaron.A0 = polaron.A0 .|> Unitful.meV
    polaron.B0 = polaron.B0 .|> Unitful.meV
    polaron.C0 = polaron.C0 .|> Unitful.meV
    polaron.v0 = polaron.v0 .|> Unitful.THz2π
    polaron.w0 = polaron.w0 .|> Unitful.THz2π
    polaron.κ0 = polaron.κ0 .|> u"N/m"
    polaron.M0 = polaron.M0 .|> Unitful.kg
    polaron.M0a = polaron.M0a .|> Unitful.kg
    polaron.M0r = polaron.M0r .|> Unitful.kg
    polaron.R0 = polaron.R0 .|> Unitful.Å
    polaron.T = polaron.T .|> Unitful.K
    polaron.β = polaron.β .|> u"meV^-1"
    polaron.F = polaron.F .|> Unitful.meV
    polaron.A = polaron.A .|> Unitful.meV
    polaron.B = polaron.B .|> Unitful.meV
    polaron.C = polaron.C .|> Unitful.meV
    polaron.v = polaron.v .|> Unitful.THz2π
    polaron.w = polaron.w .|> Unitful.THz2π
    polaron.κ = polaron.κ .|> u"N/m"
    polaron.M = polaron.M .|> Unitful.kg
    polaron.Ma = polaron.Ma .|> Unitful.kg
    polaron.Mr = polaron.Mr .|> Unitful.kg
    polaron.R = polaron.R .|> Unitful.Å
    polaron.μ = polaron.μ .|> u"cm^2/V/s"
    polaron.μFHIP = polaron.μFHIP .|> u"cm^2/V/s"
    polaron.μD = polaron.μD .|> u"cm^2/V/s"
    polaron.μK = polaron.μK .|> u"cm^2/V/s"
    polaron.μH = polaron.μH .|> u"cm^2/V/s"
    polaron.μH0 = polaron.μH0 .|> u"cm^2/V/s"
    polaron.τ = polaron.τ .|> Unitful.ns
    polaron.Ω = polaron.Ω .|> Unitful.THz2π
    polaron.χ = polaron.χ .|> u"kΩ" .* Unitful.THz2π
    polaron.z = polaron.z .|> Unitful.kΩ
    polaron.σ = polaron.σ .|> Unitful.mS
end

function suconvert!(polaron::HolsteinPolaron)
    polaron.ω = polaron.ω  .|> Unitful.THz2π
    polaron.ωeff = polaron.ωeff .|> Unitful.THz2π
    polaron.J = polaron.J .|> Unitful.meV
    polaron.a = polaron.a .|> Unitful.Å
    polaron.F0 = polaron.F0 .|> Unitful.meV
    polaron.A0 = polaron.A0 .|> Unitful.meV
    polaron.B0 = polaron.B0 .|> Unitful.meV
    polaron.C0 = polaron.C0 .|> Unitful.meV
    polaron.v0 = polaron.v0 .|> Unitful.THz2π
    polaron.w0 = polaron.w0 .|> Unitful.THz2π
    polaron.κ0 = polaron.κ0 .|> u"N/m"
    polaron.M0 = polaron.M0 .|> Unitful.kg
    polaron.M0a = polaron.M0a .|> Unitful.kg
    polaron.M0r = polaron.M0r .|> Unitful.kg
    polaron.R0 = polaron.R0 .|> Unitful.Å
    polaron.T = polaron.T .|> Unitful.K
    polaron.β = polaron.β .|> u"meV^-1"
    polaron.F = polaron.F .|> Unitful.meV
    polaron.A = polaron.A .|> Unitful.meV
    polaron.B = polaron.B .|> Unitful.meV
    polaron.C = polaron.C .|> Unitful.meV
    polaron.v = polaron.v .|> Unitful.THz2π
    polaron.w = polaron.w .|> Unitful.THz2π
    polaron.κ = polaron.κ .|> u"N/m"
    polaron.M = polaron.M .|> Unitful.kg
    polaron.Ma = polaron.Ma .|> Unitful.kg
    polaron.Mr = polaron.Mr .|> Unitful.kg
    polaron.R = polaron.R .|> Unitful.Å
    polaron.μ = polaron.μ .|> u"cm^2/V/s"
    polaron.Ω = polaron.Ω .|> Unitful.THz2π
    polaron.χ = polaron.χ .|> u"kΩ" .* Unitful.THz2π
    polaron.z = polaron.z .|> Unitful.kΩ
    polaron.σ = polaron.σ .|> Unitful.mS
end