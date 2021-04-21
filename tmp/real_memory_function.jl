
include("hypergeometric_expansion.jl")
import .hypergeom_exp

include("bessel_minus_struve.jl")
import .BesselI_minus_StruveL

using ArbNumerics
using Plots
plotly()
Plots.PlotlyBackend()

function arb_binomial(x, y)
    x = ArbReal("$x")
    y = ArbReal("$y")
    one = ArbReal("1")
    ArbNumerics.gamma(x + one) / (ArbNumerics.gamma(y + one) * ArbNumerics.gamma(x - y + one))
end

function ℜχ(Ω, β, α, v, w; prec = 24)

    # Initialise precision of ArbReal to prec.
    p = prec
    setextrabits(0)
    setprecision(ArbReal, prec)

    Ω = ArbReal("$Ω")
    β = ArbReal("$β")
    α = ArbReal("$α")
    v = ArbReal("$v")
    w = ArbReal("$w")

    R = ArbReal((v^2 - w^2) / (w^2 * v))
    a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
    b = ArbReal(R * β / sinh(β * v / 2))

    coefficient = ArbReal(α * β^(3/2) * v^3 / (3 * a * √π * w^3 * sinh(β / 2)))

    M_c(z, n) = hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 0; prec = p)
    M_s(z, n) = hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 1; prec = p)
    θ(z, n) = bessel_minus_struve.BesselI_minus_StruveL(n, z, a; prec = p)

    c = ArbReal(cosh(Ω * β / 2))
    s = ArbReal(sinh(Ω * β / 2))

    n = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    while true

        previous_result = result

        if mod(n, 2) == 0
            even_term = ArbReal(arb_binomial(n, n/2) * M_c(1, n))
            z = [ArbReal(Ω + 1), ArbReal(Ω - 1)]
            even_term += ArbReal(sum(arb_binomial(n, n/2) .* (s .* (M_s.(z, n) .+ θ.(z, n)) .- c .* M_c.(z, n)) ./ 2))
        else
            even_term = ArbReal("0.0")
        end

        k_term = ArbReal("0.0")
        for k in ArbReal("0"):ArbReal(floor(n / 2 - 1 / 2))
            z = [ArbReal(1 + v * (n - 2 * k)), ArbReal(1 - v * (n - 2 * k))]
            k_term += ArbReal(sum(arb_binomial.(n, k) .* M_c.(z, n)))

            z = [ArbReal(Ω + 1 + v * (n - 2 * k)), ArbReal(Ω - 1 + v * (n - 2 * k)), ArbReal(Ω + 1 - v * (n - 2 * k)), ArbReal(Ω - 1 - v * (n - 2 * k))]
            k_term += ArbReal(sum(arb_binomial.(n, k) .* (s .* (M_s.(z, n) + θ.(z, n)) .- c .* M_c.(z, n)) ./ 2))
        end

        term = ArbReal(coefficient * arb_binomial(-3/2, n) * (-b / (2 * a))^n * (even_term + k_term))

        result += term
        # Break loop if term smaller than accuracy of result.
        if abs(term/result) < err
            break
        end
        println("term: n = ", n, "\nterm value: ", ball(ArbReal(term, bits = prec)), "\ncumulant result: ", ball(ArbReal(result, bits = prec)), "\n")
        n += ArbReal("1")

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if radius(result) > err
            p += precision(result)
            setprecision(ArbReal, p)

            println("Not precise enough. Error = ", radius(result), " > ", err, ". Increasing precision to ", p, " bits.\n")

            Ω = ArbReal("$Ω")
            β = ArbReal("$β")
            α = ArbReal("$α")
            v = ArbReal("$v")
            w = ArbReal("$w")

            R = ArbReal((v^2 - w^2) / (w^2 * v))
            a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
            b = ArbReal(R * β / sinh(β * v / 2))

            coefficient = ArbReal(α * β^(3/2) * v^3 / (3 * a * √π * w^3 * sinh(β / 2)))

            n = ArbReal("$(n - 1)")

            M_c(z, n) = hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 0; prec = p)
            M_s(z, n) = hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 1; prec = p)
            θ(z, n) = bessel_minus_struve.BesselI_minus_StruveL(n, z, a; prec = p)

            result = ArbReal("$previous_result")

            c = ArbReal(cosh(Ω * β / 2))
            s = ArbReal(sinh(Ω * β / 2))
        end
    end
    println("Frequency: ", ArbReal(Ω, bits = prec), ". Final result: ", ArbReal(result, bits = prec))
    ArbReal(result, bits = prec)
end

Ω_range = 20.01
β = 10.0
α = 7.0
v = 5.8
w = 1.6
@time ReX = [ℜχ(Ω, β, α, v, w; prec = 24) for Ω in Ω_range]
@show(ReX)
p = plot(Ω_range, ReX, xlabel = "Ω", ylabel = "ReX")
display(p)
