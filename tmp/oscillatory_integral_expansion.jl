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

function oscillatory_integral(Ω, β, α, v, w; prec = 64)

    setprecision(BigFloat, prec)

    Ω = BigFloat("$Ω")
    β = BigFloat("$β")
    α = BigFloat("$α")
    v = BigFloat("$v")
    w = BigFloat("$w")

    R = BigFloat("$((v^2 - w^2) / (w^2 * v))")
    a = BigFloat("$(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))")
    b = BigFloat("$(R * β / sinh(β * v / 2))")

    coefficient = BigFloat("$(2 * α * β^(3 / 2) * v^3 / (3 * √π * w^3 * sinh(β / 2)))")

    integrand(x) = BigFloat("$((1 - cos(Ω * x) * cosh(Ω * β / 2)) * cos(x) / (x^2 + a^2 - b * cos(v * x))^(3 / 2))")

    integral = QuadGK.quadgk(x -> integrand(x), BigFloat("0"), BigFloat("$Inf"))[1]

    BigFloat("$(coefficient * integral)")
end

function oscillatory_integral_expansion(Ω, β, α, v, w; prec = 64) # z > 0 & n >= 0

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec + 8) # ArbReal(bit = 64 + 8) has same precision as Float64 (which is why we add 8 bits).

    Ω = ArbReal("$Ω")
    β = ArbReal("$β")
    α = ArbReal("$α")
    v = ArbReal("$v")
    w = ArbReal("$w")

    R = ArbReal((v^2 - w^2) / (w^2 * v))
    a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
    b = ArbReal(R * β / sinh(β * v / 2))

    coefficient = ArbReal(2 * α * v^3 * β^(3 / 2) * sinh(Ω * β / 2) / (3 * sqrt(π) * sinh(β / 2) * w^3))

    n = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision = prec.

    while true

        n_coefficient = ArbReal(-π * (-b / (4 * a))^n / (4 * a))

        if mod(n, 2) == 0
            even_term = ArbReal("0.0")
            for z in [ArbReal(Ω + 1), ArbReal(Ω - 1)]
                even_term += ArbReal(abs(z)^(n + 1) * sign(z) * bessel_minus_struve.BesselI_minus_StruveL(n + 1, a * abs(z); prec = precision(result)))
            end
            even_term = ArbReal(even_term / (ArbNumerics.gamma(n / 2 + 1)^2))
        else
            even_term = ArbReal("0.0")
        end

        k_term = ArbReal("0.0")
        k = ArbReal("0")
        for k in ArbReal("0"):ArbReal(floor((n - 1) / 2))

            k_coeff = ArbReal(ArbNumerics.gamma(k + 1) * ArbNumerics.gamma(n - k + 1))

            if isnan(k_coeff) || k_coeff == 0.0
                k_term = ArbReal("0.0")
            else
                for z in [ArbReal(Ω + 1 + v * (n - 2 * k)), ArbReal(Ω - 1 + v * (n - 2 * k)), ArbReal(Ω + 1 - v * (n - 2 * k)), ArbReal(Ω - 1 - v * (n - 2 * k))]
                    k_term += ArbReal(abs(z)^(n + 1) * sign(z) * bessel_minus_struve.BesselI_minus_StruveL(n + 1, a * abs(z); prec = precision(result)))
                end
                k_term = k_term / k_coeff
            end
        end

        term = n_coefficient * (even_term + k_term)

        # Break loop if term smaller than accuracy of result. (I.e. indistinguishable at set precison).
        if abs(term) < err
            break
        end
        # @show(n, term, result)
        result += term
        n += ArbReal("1")

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if ball(result)[2] > err
            setprecision(ArbReal, precision(result) * 2)
            Ω = ArbReal("$Ω")
            β = ArbReal("$β")
            α = ArbReal("$α")
            v = ArbReal("$v")
            w = ArbReal("$w")

            R = ArbReal((v^2 - w^2) / (w^2 * v))
            a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
            b = ArbReal(R * β / sinh(β * v / 2))

            n = ArbReal("0")
            result = ArbReal("0.0")
        end
    end
    result *=  coefficient
    ArbReal(result, bits = prec + 8) # return to specified precision.
end

Ω_range = 0.001:0.01:3.01
# For Ω_range = 3.1:
# 9.828299 seconds (46.17 M allocations: 1.859 GiB, 12.89% gc time)
@time oscill_expansion = [oscillatory_integral_expansion(Ω, 100, 5, 4.0, 2.8; prec = 64) for Ω in Ω_range]
# oscill_expansion = fill(36.009791933858445279331705890513148915488)
@show(oscill_expansion)
# 273.411207 seconds (931.57 M allocations: 32.002 GiB, 5.19% gc time)
# @time oscill_integral = [oscillatory_integral(Ω, 2.3, 5, 4.0, 2.8; prec = 128) for Ω in Ω_range]
# oscill_integral = fill(-5.81034244862553995092)
# @show(oscill_integral)
p = plot(Ω_range, abs.(oscill_expansion), yaxis = :log, label = "Oscillatory", xlabel = "Ω", ylabel = "ReX")
# plot!(Ω_range, oscill_integral, label = "Integral")
display(p)
