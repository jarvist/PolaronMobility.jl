# hyperbolic_integral_expansion.jl

"""
Note: This is meant to be the full expansion for the hyperbolic integral using the correct arbitrary precision methods, but it currently does not converge to the correct values, so made a typo/error somewhere. It's not massively efficient either, so optimisation will be important!
"""

include("arb_hypgeom.jl")
import .one_F_two_fast

include("hypergeometric_expansion.jl")
import .hypergeom_exp

using ArbNumerics
using QuadGK

"""
Solving the hyperbolic integral with expansions and summations.
"""
function arb_binomial(x, y)
    x = ArbReal("$x")
    y = ArbReal("$y")
    one = ArbReal("1")
    ArbNumerics.gamma(x + one) / (ArbNumerics.gamma(y + one) * ArbNumerics.gamma(x - y + one))
end

function hyperbolic_expansion(Ω, β, α, v, w; prec = 64)

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec + 8)

    Ω = ArbReal("$Ω")
    β = ArbReal("$β")
    α = ArbReal("$α")
    v = ArbReal("$v")
    w = ArbReal("$w")

    R = ArbReal((v^2 - w^2) / (w^2 * v))
    a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
    b = ArbReal(R * β / sinh(β * v / 2))

    coefficient = ArbReal(α * β^(3/2) * v^3 / (3 * a * √π * w^3 * sinh(β / 2)))

    M_c(z, n) = ArbReal(hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 0; prec = prec))
    M_s(z, n) = ArbReal(hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 1; prec = prec))

    c = ArbReal(cosh(Ω * β / 2))
    s = ArbReal(sinh(Ω * β / 2))

    n = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(ArbReal(0, bits = prec))  # Machine accuracy of specified precision prec.

    while true

        previous_result = ArbReal("$result")

        if mod(n, 2) == 0
            even_term = ArbReal(arb_binomial(n, n/2) * M_c(1, n))
            @fastmath @inbounds @simd for z in [ArbReal(Ω + 1), ArbReal(Ω - 1)]
                even_term += ArbReal(arb_binomial(n, n/2) * (s * M_s(z, n) - c * M_c(z, n)) / 2)
            end
        else
            even_term = ArbReal("0.0")
        end

        k_term = ArbReal("0.0")
        @fastmath @inbounds @simd for k in ArbReal("0"):ArbReal(floor(n / 2 - 1 / 2))
            @fastmath @inbounds @simd for z in [ArbReal(1 + v * (n - 2 * k)), ArbReal(1 - v * (n - 2 * k))]
                k_term += ArbReal(arb_binomial(n, k) * M_c(z, n))
            end
            @fastmath @inbounds @simd for z in [ArbReal(Ω + 1 + v * (n - 2 * k)), ArbReal(Ω - 1 + v * (n - 2 * k)), ArbReal(Ω + 1 - v * (n - 2 * k)), ArbReal(Ω - 1 - v * (n - 2 * k))]
                k_term += ArbReal(arb_binomial(n, k) * (s * M_s(z, n) - c * M_c(z, n)) / 2)
            end
        end

        term = ArbReal(arb_binomial(-3/2, n) * (-b / (2 * a))^n * (even_term
         + k_term))
        # Break loop if term smaller than accuracy of result.
        if abs(term) < err
            break
        end
        result += ArbReal(term)
        println("term: n = ", n, "\nterm value: ", coefficient * term, "\ncurrent result: ", coefficient * result, "\n")
        n += ArbReal("1")

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if ball(result)[2] > err
            setprecision(ArbReal, precision(result) * 2)

            println("Not precise enough. Required error < ", err, ". Increasing precision to ", precision(result) * 2, " bits.\n")

            Ω = ArbReal("$Ω")
            β = ArbReal("$β")
            α = ArbReal("$α")
            v = ArbReal("$v")
            w = ArbReal("$w")

            R = ArbReal((v^2 - w^2) / (w^2 * v))
            a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
            b = ArbReal(R * β / sinh(β * v / 2))

            coefficient = ArbReal(α * β^(3/2) * v^3 / (3 * a * √π * w^3 * sinh(β / 2)))

            M_c(z, n) = ArbReal(hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 0; prec = precision(result)))
            M_s(z, n) = ArbReal(hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 1; prec = precision(result)))

            c = ArbReal(cosh(Ω * β / 2))
            s = ArbReal(sinh(Ω * β / 2))

            n = ArbReal("$(n - 1)")
            result = ArbReal("$previous_result")
        end
    end
    result *= ArbReal(coefficient)
    println("Frequency: ", ArbReal(Ω, bits = prec + 8), ". Final result: ", ArbReal(result, bits = prec + 8))
    ArbReal(result, bits = prec + 8)
end

"""
Solving the hyperbolic integral with quadgk
"""
function hyperbolic_integral(Ω, β, α, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    integrand(x) = (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2)

    integral = quadgk(x -> integrand(x), 0, β / 2)
    println(Ω, " ", coefficient * integral[1])
    return coefficient * integral[1]
end

"""
Plot results for comparison.
"""

Ω_range = 7.01
β = 2.0
α = 7.0
v = 5.8
w = 1.6
R_int = [hyperbolic_integral(Ω, β, α, v, w) for Ω in Ω_range]
# R_osc = [oscillatory_integral_expansion(Ω, β, α, v, w; prec = 32) for Ω in Ω_range]
R_exp = [hyperbolic_expansion(Ω, β, α, v, w; prec = 32) for Ω in Ω_range]
# @show(R_exp)
# p = plot(Ω_range, abs.(R_exp), yaxis = :log)
# plot!(Ω_range, abs.(R_int))
# display(p)
