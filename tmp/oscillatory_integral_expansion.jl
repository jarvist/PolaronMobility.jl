"""
So, ReX = oscillatory integral + hyperbolic integral. I expanded the first integral ages ago, but it was not an arbitrary precision algorithm. I tested solving ReX with this old expansion for the oscillatory integral added to the hyperbolic integral solved by arbitrary precision quadgk, and it still diverged when digits beyond Float64 became important. Thus, the oscillatory integral needs to be rewritten in the arbitrary precision summation algorithm (the algorithm I used for BesselI-StruveL). I will do this here and also test it against brute forcing the integral with arb. prec. quadgk provided quadgk manages to solve this!
"""

include("bessel_minus_struve.jl")
import .bessel_minus_struve

using QuadGK
using Plots
using ArbNumerics
plotly()
Plots.PlotlyBackend()

function oscillatory_integral(Ω, β, α, v, w)

    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 * sinh(Ω * β / 2) / (3 * √π * w^3 * sinh(β / 2))

    integrand(x) = (sin(Ω * x) * cos(x) / (x^2 + a^2 - b * cos(v * x))^(3 / 2))

    integral = QuadGK.quadgk(x -> integrand(x), 0, Inf)[1]
    @show(coefficient * integral)
    BigFloat(coefficient * integral)
end

function arb_binomial(x, y)
    x = ArbReal("$x")
    y = ArbReal("$y")
    one = ArbReal("1")
    ArbNumerics.gamma(x + one) / (ArbNumerics.gamma(y + one) * ArbNumerics.gamma(x - y + one))
end

function oscillatory_integral_expansion(Ω, β, α, v, w; prec = 64)

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

    θ(z, n) = ArbReal(√π * ArbNumerics.gamma(-n - 1/2) * abs(z)^n * z * bessel_minus_struve.BesselI_minus_StruveL(n + 1, a * abs(z); prec = prec) / 2)

    s = ArbReal(sinh(Ω * β / 2))

    n = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(ArbReal(0, bits = prec + 8))  # Machine accuracy of specified precision prec.

    while true

        previous_result = ArbReal("$result")

        if mod(n, 2) == 0
            even_term = ArbReal("0.0")
            @fastmath @inbounds @simd for z in [ArbReal(Ω + 1), ArbReal(Ω - 1)]
                even_term += ArbReal(arb_binomial(n, n/2) * s * θ(z, n) / 2)
            end
        else
            even_term = ArbReal("0.0")
        end

        k_term = ArbReal("0.0")
        @fastmath @inbounds @simd for k in ArbReal("0"):ArbReal(floor(n / 2 - 1 / 2))
            k_coeff = ArbReal(arb_binomial(n, k))
            @fastmath @inbounds @simd for z in [ArbReal(Ω + 1 + v * (n - 2 * k)), ArbReal(Ω - 1 + v * (n - 2 * k)), ArbReal(Ω + 1 - v * (n - 2 * k)), ArbReal(Ω - 1 - v * (n - 2 * k))]
                k_term += ArbReal(k_coeff * s * θ(z, n) / 2)
            end
        end

        term = ArbReal(arb_binomial(-3/2, n) * (-b / (4 * a))^n * (even_term
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
            setprecision(ArbReal, precision(result) + 8)

            println("Not precise enough. Required error < ", err, ". Increasing precision to ", precision(result) + 8, " bits.\n")

            Ω = ArbReal("$Ω")
            β = ArbReal("$β")
            α = ArbReal("$α")
            v = ArbReal("$v")
            w = ArbReal("$w")

            R = ArbReal((v^2 - w^2) / (w^2 * v))
            a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
            b = ArbReal(R * β / sinh(β * v / 2))

            coefficient = ArbReal(α * β^(3/2) * v^3 / (3 * a * √π * w^3 * sinh(β / 2)))

            θ(z, n) = ArbReal(√π * ArbNumerics.gamma(-n - 1/2) * abs(z)^n * z * bessel_minus_struve.BesselI_minus_StruveL(n + 1, a * abs(z); prec = precision(result)) / 2)

            s = ArbReal(sinh(Ω * β / 2))

            n = ArbReal("$(n - 1)")
            result = ArbReal("$previous_result")
        end
    end

    result *= ArbReal(coefficient)
    println("Frequency: ", ArbReal(Ω, bits = prec + 8), ". Final result: ", ArbReal(result, bits = prec + 8))
    ArbReal(result, bits = prec + 8)
end

Ω_range = 0.01:0.5:20.01
β = 2.0
α = 7.0
v = 5.8
w = 1.6
# For Ω_range = 3.1:
# 9.828299 seconds (46.17 M allocations: 1.859 GiB, 12.89% gc time)
@time oscill_expansion = [oscillatory_integral_expansion(Ω, β, α, v, w; prec = 64) for Ω in Ω_range]
# oscill_expansion = fill(36.009791933858445279331705890513148915488)
# @show(oscill_expansion)
# 273.411207 seconds (931.57 M allocations: 32.002 GiB, 5.19% gc time)
@time oscill_integral = [oscillatory_integral(Ω, β, α, v, w) for Ω in Ω_range]

diff = oscill_expansion .- oscill_integral
@show(diff)
# oscill_integral = fill(-5.81034244862553995092)
# @show(oscill_integral)
# p = plot(Ω_range, abs.(oscill_expansion), yaxis = :log, label = "Oscillatory", xlabel = "Ω", ylabel = "ReX")
# plot!(Ω_range, abs.(oscill_integral), label = "Integral")
# display(p)
