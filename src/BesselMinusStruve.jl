# bessel_minus_struve.jl

module bessel_minus_struve
export BesselI_minus_StruveL

using ArbNumerics

"""
Implementation of I[n, z] - L[-n, z] where I[n, z] is the modified Bessel function of the first kind, and L[-n, z] is the modified Struve function. n is an integer ≥ 0 and the order of the functions, and z∈ℜ>0 a positive real argument. Notice the order are opposite signs.

The expansion of the BesselI and StruveL functions are:

    I[n, z] = ∑_{k=0}^∞ (z/2)^{2k+n} / (Γ(k+1)⋅Γ(k+n+1))
    L[-n, z] = ∑_{k=0}^∞ (z/2)^{2k-n+1} / (Γ(k+3/2)⋅Γ(k-n+3/2))

with Γ(x) the Gamma function. The values of I[n, z] and L[-n, z] can be incredibly close, so ArbReal types are used as they allow for accurate arbitrary precision floats so that the difference, I[n, z] - L[-n, z], can be distinguished. ArbReal types are faster than BigFloat types, and the error of a float can be tracked using ball(::ArbReal)[2].

Unfortunately, as ArbReal(bits = 64) * ArbReal(bits = 128) = ArbReal(bits = 64), if the starting precision is too small to accurately produce I[n, z] - L[-n, z], then the entire summation has to be restarted at a higher precision otherwise rounding errors will be too large, which slows the algorithm sometimes. Therefore, if the precision is found to be too low, the precision is just doubled (rather than increased by some set amount) to reach a required precision quickly.

The general structure of this arbitrary precision summation algorithm is used elsewhere too since you just change the arguments and term appropriately.
"""

function BesselI_minus_StruveL(n, z, a; prec = 64) # z > 0 & n >= 0

    # Initialise precision of ArbReal to prec.
    p = prec
    setextrabits(0)
    setprecision(ArbReal, p) # ArbReal(bit = 64 + 8) has same precision as Float64 (which is why we add 8 bits).

    n = ArbReal("$n")
    z = ArbReal("$z")
    a = ArbReal("$a")

    k = ArbReal("0.0")
    result = ArbReal("0.0")
    term = ArbReal("1.0")
    err = eps(result)  # Machine accuracy of specified precision = prec.

    while abs(midpoint(term)) > err * abs(midpoint(result))
        term = ArbReal((abs(z) * a / 2)^(2 * k + n + 1) / (ArbNumerics.gamma(k + 1) * ArbNumerics.gamma(k + n + 2)) - (abs(z) * a / 2)^(2 * k - n) / (ArbNumerics.gamma(k + 3//2) * ArbNumerics.gamma(k - n + 1//2)))
        result += term
        # println("term: k = ", k, "\nterm value: ", ball(term), "\ncumulant result: ", ball(result), "\n")
        k += 1
        # Double precision if rounding error in result exceeds accuracy specified by prec.
        while radius(result) > err * abs(midpoint(result))
            p *= 2
            setprecision(ArbReal, p)
            # variables have to be re-parsed into the higher precision ArbReal type.
            # println("Not precise enough. Error = ", radius(result), " > ", ArbReal("$err"), ". Increasing precision to ", p, " bits.\n")
            n = ArbReal("$n")
            z = ArbReal("$z")
            a = ArbReal("$a")
            k = ArbReal("0")
            result = ArbReal("0.0")
            term = ArbReal("1.0")
        end
    end
    # println("z: ", ArbReal(z, bits = prec), ". Final result: ", ArbReal(result, bits = prec))
    ArbReal(result * √π * ArbNumerics.gamma(-n - 1/2) * (abs(z) / 2)^n * z / 2, bits = prec) # return to specified precision.
end
# @time bsi = BesselI_minus_StruveL(3, 55, 1)
# @show(bsi)

# β = ArbReal("2.0")
# α = ArbReal("7.0")
# v = ArbReal("5.8")
# w = ArbReal("1.6")
# R = ArbReal((v^2 - w^2) / (w^2 * v))
# a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
# z = ArbReal("90.61")
# n = ArbReal("12.0")
# c = BesselI_minus_StruveL(n, z, a; prec = 24)

end # end of module
