# hyperbolic_expansion.jl

module hypergeometric_expansion
export hypergeom_exp

include("arb_hypgeom.jl")
import .arb_hypgeom
using ArbNumerics

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
    setextrabits(0)
    setprecision(ArbReal, prec + 8)

    z = ArbReal("$z")
    n = ArbReal("$n")
    β = ArbReal("$β")
    a = ArbReal("$a")

    m_binomial_coeff(m) = ArbReal(ArbNumerics.gamma(-n - 1//2) / (ArbNumerics.gamma(m + 1) * ArbNumerics.gamma(-n - m - 1//2)))

    S(m) = ArbReal(arb_hypgeom.one_f_two_fast(z, m, h; prec = prec))

    m = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    while true
        term = ArbReal(m_binomial_coeff(m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3) * S(m))

        # Break loop if term smaller than accuracy of result.
        if abs(term) < err
            break
        end
        result += term
        m += ArbReal("1")

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if ball(result)[2] > err
            setprecision(ArbReal, precision(result) * 2)
            n = ArbReal("$n")
            z = ArbReal("$z")
            β = ArbReal("$β")
            a = ArbReal("$a")

            m = ArbReal("0")
            result = ArbReal("0.0")

            m_binomial_coeff(m) = ArbReal(ArbNumerics.gamma(-n - 1//2) / (ArbNumerics.gamma(m + 1) * ArbNumerics.gamma(-n - m - 1//2)))

            S(m) = ArbReal(arb_hypgeom.one_f_two_fast(z, m, h; prec = precision(result)))
        end
    end
    ArbReal(result, bits = prec + 8)
end

end # end module
