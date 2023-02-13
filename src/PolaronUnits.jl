module PolaronUnits

    import Unitful
    using Unitful: @unit, Dimension, Dimensions, NoDims, NoUnits, Units, dimension, uconvert, ustrip
    
    export puconvert, punit, pustrip, m0_pu, e_pu, ħ_pu, k_pu, ω0_pu, a0_pu, E0_pu
    
    # The following five constants are used as the “base” atomic units
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
    @unit ħ_pu  "ħ"   ReducedPlanckConstant Unitful.ħ                      false

    """
        PolaronUnits.k_pu
    A unit equal to the Boltzmann constant k = 1.380,649 × 10^-23 J / K.
    Printed as "k".
    `Unitful.k` is a quantity (with units `J / K`) whereas `PolaronUnits.k_pu` is a unit equal to
    `Unitful.k`.
    Dimension: 𝐋^2 𝐌 𝚯^-1 𝐓^-2.
    See also: `Unitful.k`, `Unitful.J`, `Unitful.K`.
    """
    @unit k_pu  "k"   BoltzmannConstant     Unitful.k                      false

    """
        PolaronUnits.ω_pu
    A unit equal to the reduced Planck constant ħ = h / 2π ≈ 1.054,571,8176 × 10^-34 J × s.
    Printed as "ħ".
    `Unitful.ħ` is a quantity (with units `J × s`) whereas `PolaronUnits.ħ_pu` is a unit equal to
    `Unitful.ħ`.
    Dimension: 𝐋^2 𝐌 𝐓^-1.
    See also: `Unitful.ħ`, `Unitful.J`, `Unitful.s`.
    """
    @unit ω0_pu  "ω₀"   PolaronAngularFrequency 1Unitful.THz2π               false

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
    @unit a0_pu "a₀" PolaronRadius √(1ħ_pu/1m0_pu/1ω0_pu)                       false
    
    # Polaron energy is derived from the base polaron units
    """
        PolaronUnits.E0_pu
    A unit equal to the phonon energy
    ```
    E₀ = ħω₀
    ```
    Printed as "Eₕ".
    Dimension: 𝐋^2 𝐌 𝐓^-2.
    See also: `Unitful.me`, `Unitful.q`, `Unitful.ε0`, `Unitful.ħ`, `Unitful.J`, `Unitful.eV`, [`UnitfulAtomic.Ry`](@ref).
    """
    @unit E0_pu "E₀"  PolaronEnergy         1ħ_pu*1ω0_pu                       false

    """
        punit(x::Unitful.Quantity)
        punit(x::Unitful.Units)
        punit(x::Unitful.Dimensions)
    Returns the appropriate polaron unit (a `Unitful.Units` object) for the dimension of `x`.
    # Examples
    ```jldoctest
    julia> punit(2.3u"cm")
    a₀
    julia> punit(u"T")
    a₀^-2 e^-1 ħ
    ```
    """
    punit(x) = punit(dimension(x))

    # `aunit` for `Dimension` types
    punit(x::Dimension{:Length})      = (a0_pu)^x.power
    punit(x::Dimension{:Mass})        = (m0_pu)^x.power
    punit(x::Dimension{:Time})        = (ħ_pu/E0_pu)^x.power
    punit(x::Dimension{:Current})     = (e_pu*E0_pu/ħ_pu)^x.power
    punit(x::Dimension{:Temperature}) = (E0_pu/k_pu)^x.power

    # For dimensions not specified above, there is no atomic unit.
    punit(::Dimension{D}) where D = throw(ArgumentError("No polaron unit defined for dimension $D."))

    # `aunit` for `Dimensions` types
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
    Convert a quantity to the appropriate atomic unit.
    # Examples
    ```jldoctest
    julia> puconvert(13.6u"eV")
    0.4997907858599377 E₀
    julia> puconvert(20u"nm")
    377.94522509156565 a₀
    ```
    """
    puconvert(x) = uconvert(punit(x), x)

    """
        puconvert(u::Unitful.Units, x::Number)
    Interpret `x` as a quantity given in atomic units and convert it to the unit `u`.
    # Examples
    ```jldoctest
    julia> puconvert(u"eV", 1)  # convert 1 Eₕ to eV
    27.211386013449417 eV
    julia> puconvert(u"m", 1)   # convert 1 a₀ to m
    5.2917721067e-11 m
    ```
    """
    puconvert(u::Units, x::Number) = uconvert(u, x*punit(u))

    """
        pustrip(x::Unitful.Quantity)
    Returns the value of the quantity converted to polaron units as a number type (i.e., with the
    units removed). This is equivalent to `Unitful.ustrip(puconvert(x))`.
    # Examples
    ```jldoctest
    julia> pustrip(13.6u"eV")
    0.4997907858599377
    julia> pustrip(20u"nm")
    377.94522509156565
    ```
    """
    pustrip(x) = ustrip(puconvert(x))      

    Unitful.register(PolaronUnits)

    __init__() = Unitful.register(PolaronUnits)

end