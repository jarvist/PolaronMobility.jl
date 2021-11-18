# hyperbolic_expansion.jl

module hypergeometric_expansion
export hypergeom_exp

include("arb_hypgeom.jl")
import .arb_hypgeom
using ArbNumerics

function arb_binomial(x, y; prec = 64)
    setextrabits(0)
    setprecision(ArbReal, prec)
    x = ArbReal("$x")
    y = ArbReal("$y")
    one = ArbReal("1")
    ArbNumerics.gamma(x + one) / (ArbNumerics.gamma(y + one) * ArbNumerics.gamma(x - y + one))
end

"""
This function implements the binomial expansion of the denominator of the integrals:

    ∫_0^1 cosh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt
    ∫_0^1 sinh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt

which is:

    ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 cosh(zt) t^{2m} dt
    ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 sinh(zt) t^{2m} dt.

and also substitutes the integrals for:

    ∫_0^1 cosh(zt) t^{2m} dt = ∑_{t=0}^∞ z^{2t} / ((2t+2m+1)⋅(2t)!)
    ∫_0^1 sinh(zt) t^{2m} dt = ∑_{t=0}^∞ z^{2t+1} / ((2t+2m+2)⋅(2t+1)!)

by replacing those integrals with one_F_two_fast() from arb_hypgeom.jl. h=0 gives cosh version, h=1 gives sinh version.
"""
function hypergeom_exp(z, n, β, a, h; prec = 64)
# h = 0 gives cosh, h = 1 gives sinh.

    # Initialise precision of ArbReal to prec.
    p = prec
    setextrabits(0)
    setprecision(ArbReal, p)

    z = ArbReal("$z")
    n = ArbReal("$n")
    β = ArbReal("$β")
    a = ArbReal("$a")

    m = ArbReal("0")
    result = ArbReal("0.0")
    term = ArbReal("1.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    while abs(midpoint(term)) > err * abs(midpoint(result))

        term = ArbReal(β * arb_binomial(-n - 3/2, m; prec = p) * (-1)^m * (β / (2 * a))^(2 * m) * arb_hypgeom.one_f_two_fast(z, m, h; prec = p) / a^(n + 2))
        result += term
        # println("term: m = ", m, "\nterm value: ", ball(ArbReal(term, bits = prec)), "\ncumulant result: ", ball(ArbReal(result, bits = prec)), "\n")
        m += ArbReal("1")

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if radius(result) > err * abs(midpoint(result))
            p *= 2
            setprecision(ArbReal, p)
            # println("Not precise enough. Error = ", abs(radius(result)/midpoint(result)), " > ", err, ". Increasing precision to ", p, " bits.\n")

            n = ArbReal("$n")
            z = ArbReal("$z")
            β = ArbReal("$β")
            a = ArbReal("$a")

            m = ArbReal("0")
            result = ArbReal("0.0")
            term = ArbReal("1.0")
        end
    end
    # println("z: ", ArbReal(z, bits = prec), ". Final result: ", ArbReal(result, bits = prec))
    ArbReal(result, bits = prec)
end

# β = ArbReal("2.0")
# α = ArbReal("7.0")
# v = ArbReal("5.8")
# w = ArbReal("1.6")
# R = ArbReal((v^2 - w^2) / (w^2 * v))
# a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
# z = ArbReal("90.61")
# n = ArbReal("12.0")
# @time c = hypergeom_exp(z, n, β, a, 0; prec = 64)
# @show(c)
end # end module
