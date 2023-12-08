# The following five constants are used as the “base” atomic units

# Physical constants

"Planck's constant, (kgm²s⁻¹)."
const hbar = ħ = 1.054571817e-34
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
@unit ħ_pu "ħ" ReducedPlanckConstant Unitful.ħ false

"""
    PolaronUnits.k_pu
A unit equal to the Boltzmann constant k = 1.380,649 × 10^-23 J / K.
Printed as "k".
`Unitful.k` is a quantity (with units `J / K`) whereas `PolaronUnits.k_pu` is a unit equal to
`Unitful.k`.
Dimension: 𝐋^2 𝐌 𝚯^-1 𝐓^-2.
See also: `Unitful.k`, `Unitful.J`, `Unitful.K`.
"""
@unit k_pu "k" BoltzmannConstant Unitful.k false

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

"""
    PolaronUnits.J0_pu
Transfer integral energy unit
`Unitful.ħ` is a quantity (with units `J × s`) whereas `PolaronUnits.ħ_pu` is a unit equal to
`Unitful.ħ`.
Dimension: 𝐋^2 𝐌 𝐓^-1.
See also: `Unitful.ħ`, `Unitful.J`, `Unitful.s`.
"""
@unit J0_pu "J₀" PolaronTransferIntegral 1Unitful.meV false
@unit ωh0_pu "ωₕ₀" PolaronAdiabticity 1J0_pu/1ħ_pu false
@unit mh0_pu "mₕ₀" PolaronHolsteinMass 1ħ_pu^2/1J0_pu/1Unitful.Å^2 false

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
@unit ah0_pu "aₕ₀" PolaronHolsteinRadius √(1ħ_pu/1mh0_pu/1ωh0_pu) false

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
@unit Eh0_pu "Eₕ₀" PolaronHolsteinEnergy 1J0_pu false

# Polaron thermodynamic temperature is derived from the base polaron units
"""
    PolaronUnits.β_pu
A unit equal to the phonon energy
```
β₀ = ħ/k
```
Printed as "β₀".
Dimension: 𝐋^2 𝐌 𝐓^-2.
See also: `Unitful.me`, `Unitful.q`, `Unitful.ε0`, `Unitful.ħ`, `Unitful.J`, `Unitful.eV`, [`UnitfulAtomic.Ry`](@ref).
"""
@unit T0_pu "T₀" PolaronTemperature 1Unitful.K false
@unit β0_pu "β₀" PolaronBeta 1ħ_pu/1k_pu false
@unit βh0_pu "βₕ₀" PolaronHolsteinBeta 1J0_pu/1k_pu false

@unit μ0_pu "μ₀" PolaronMobility 1e_pu/1m0_pu/1ω0_pu false
@unit μh0_pu "μₕ₀" PolaronHolsteinMobility 1e_pu/1mh0_pu/1ωh0_pu false
@unit t0_pu "t₀" PolaronTime 1Unitful.ns false

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
function addunits!(polaron::FrohlichPolaron)
    polaron.ω = polaron.ω .* ω0_pu
    polaron.Fs = polaron.Fs .* E0_pu
    polaron.Fl = polaron.Fl .* E0_pu
    polaron.Ms = polaron.Ms .* m0_pu
    polaron.Ml = polaron.Ml .* m0_pu
    polaron.Rs = polaron.Rs .* a0_pu
    polaron.Rl = polaron.Rl .* a0_pu
    polaron.ΩFC = polaron.ΩFC .* ω0_pu
    polaron.F0 = polaron.F0 .* E0_pu
    polaron.A0 = polaron.A0 .* E0_pu
    polaron.B0 = polaron.B0 .* E0_pu
    polaron.C0 = polaron.C0 .* E0_pu
    polaron.κ0 = polaron.κ0 .* m0_pu * ω0_pu^2
    polaron.M0 = polaron.M0 .* m0_pu
    polaron.M0a = polaron.M0a .* m0_pu
    polaron.M0r = polaron.M0r .* m0_pu
    polaron.R0 = polaron.R0 .* a0_pu
    polaron.T = polaron.T .* T0_pu
    polaron.β = polaron.β .* β0_pu
    polaron.F = polaron.F .* E0_pu
    polaron.A = polaron.A .* E0_pu
    polaron.B = polaron.B .* E0_pu
    polaron.C = polaron.C .* E0_pu
    polaron.κ = polaron.κ .* m0_pu * ω0_pu^2
    polaron.M = polaron.M .* m0_pu
    polaron.Ma = polaron.Ma .* m0_pu
    polaron.Mr = polaron.Mr .* m0_pu
    polaron.R = polaron.R .* a0_pu
    polaron.μ = polaron.μ .* μ0_pu
    polaron.μFHIP = polaron.μFHIP .* μ0_pu
    polaron.μD = polaron.μD .* μ0_pu
    polaron.μK = polaron.μK .* μ0_pu
    polaron.μH = polaron.μH .* μ0_pu
    polaron.μH0 = polaron.μH0 .* μ0_pu
    polaron.τ = polaron.τ .* t0_pu
    polaron.Ω = polaron.Ω .* ω0_pu
    polaron.χ = polaron.χ .* ω0_pu
    polaron.z = polaron.z .* punit(u"Ω")
    polaron.σ = polaron.σ .* punit(u"S")
end

function addunits!(polaron::HolsteinPolaron)
    polaron.ω = polaron.ω .* ωh0_pu
    polaron.F0 = polaron.F0 .* Eh0_pu
    polaron.A0 = polaron.A0 .* Eh0_pu
    polaron.B0 = polaron.B0 .* Eh0_pu
    polaron.C0 = polaron.C0 .* Eh0_pu
    polaron.κ0 = polaron.κ0 .* mh0_pu * ωh0_pu^2
    polaron.M0 = polaron.M0 .* mh0_pu
    polaron.M0a = polaron.M0a .* mh0_pu
    polaron.M0r = polaron.M0r .* mh0_pu
    polaron.R0 = polaron.R0 .* ah0_pu
    polaron.T = polaron.T .* T0_pu
    polaron.β = polaron.β .* βh0_pu
    polaron.F = polaron.F .* Eh0_pu
    polaron.A = polaron.A .* Eh0_pu
    polaron.B = polaron.B .* Eh0_pu
    polaron.C = polaron.C .* Eh0_pu
    polaron.κ = polaron.κ .* mh0_pu * ωh0_pu^2
    polaron.M = polaron.M .* mh0_pu
    polaron.Ma = polaron.Ma .* mh0_pu
    polaron.Mr = polaron.Mr .* mh0_pu
    polaron.R = polaron.R .* ah0_pu
    polaron.μ = polaron.μ .* μh0_pu
    polaron.χ = polaron.χ .* ωh0_pu
    polaron.z = polaron.z .* punit(u"Ω")
    polaron.σ = polaron.σ .* punit(u"S")
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

