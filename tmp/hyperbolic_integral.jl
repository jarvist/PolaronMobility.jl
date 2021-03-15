# hyperbolic_integral.jl

include("arb_hypgeom.jl")
import .arb_hypgeom

using ArbNumerics
using Plots
using QuadGK
plotly()
Plots.PlotlyBackend()

# Original integral
function hyperbolic_integral(z, n, β, a, h)
    if h == 0
        integrand = function (t)
            cosh(z * t) / (4 * a^2 / β^2 - t^2)^(n + 3 // 2)
        end
    elseif h == 1
        integrand = function (t)
            sinh(z * t) / (4 * a^2 / β^2 - t^2)^(n + 3 // 2)
        end
    end
    integral = quadgk(t -> integrand(t), 0, 1)[1]
end

function hyperbolic_integral_mexpansion(z, n, β, a, h; prec = 64)
# h = 0 gives cosh, h = 1 gives sinh.

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec + 8)

    z = ArbReal("$z")
    n = ArbReal("$n")
    β = ArbReal("$β")
    a = ArbReal("$a")

    m = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    m_binomial_coeff(m) = ArbReal("$(ArbNumerics.gamma(-n - 1//2) / (ArbNumerics.gamma(m + 1) * ArbNumerics.gamma(-n - m - 1//2)))")

    if h == 0
        integrand = function (m, t)
            cosh(z * t) * t^(2 * m)
        end
    elseif h == 1
        integrand = function (m, t)
            sinh(z * t) * t^(2 * m)
        end
    end

    integral(m) = ArbReal("$(quadgk(t -> integrand(m, t), 0, 1)[1])")

    while true
        term = ArbReal("$(m_binomial_coeff(m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3) * integral(m))")

        # Break loop if term smaller than accuracy of result.
        if abs(term) < eps(result)
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
        end
    end
    ArbReal(result, bits = prec + 8)
end

function hyperbolic_hypergeom_mexpansion(z, n, β, a, h; prec = 64)
# h = 0 gives cosh, h = 1 gives sinh.

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec + 8)

    z = ArbReal("$z")
    n = ArbReal("$n")
    β = ArbReal("$β")
    a = ArbReal("$a")

    m_binomial_coeff(m) = ArbReal("$(ArbNumerics.gamma(-n - 1//2) / (ArbNumerics.gamma(m + 1) * ArbNumerics.gamma(-n - m - 1//2)))")

    S(m) = ArbReal("$(arb_hypgeom.one_f_two_fast(z, m, h; prec = prec))")

    m = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    while true
        term = ArbReal("$(m_binomial_coeff(m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3) * S(m))")

        # Break loop if term smaller than accuracy of result.
        if abs(term) < eps(result)
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

            m_binomial_coeff(m) = ArbReal("$(ArbNumerics.gamma(-n - 1//2) / (ArbNumerics.gamma(m + 1) * ArbNumerics.gamma(-n - m - 1//2)))")

            S(m) = ArbReal("$(arb_hypgeom.one_f_two_fast(z, m, h; prec = (precision(result)-8) * 2))")

            m = ArbReal("0")
            result = ArbReal("0.0")
        end
    end
    ArbReal(result, bits = prec + 8)
end

function hyperbolic_hypergeom_mexpansion2(z, n, β, a, h; prec = 64)
# h = 0 gives cosh, h = 1 gives sinh.

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec + 8)

    z = ArbReal("$z")
    n = ArbReal("$n")
    β = ArbReal("$β")
    a = ArbReal("$a")

    m = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    m_binomial_coeff(m) = ArbReal("$(ArbNumerics.gamma(-n - 1//2) / (ArbNumerics.gamma(m + 1) * ArbNumerics.gamma(-n - m - 1//2)))")

    if h == 0
        one_f_two = function (m)
            ArbReal("$(arb_hypgeom.one_f_two(m + 1/2, (1/2, m + 3/2), (z / 2)^2; prec = prec) / (2 * m + 1))")
        end
    elseif h == 1
        one_f_two = function (m)
            ArbReal("$(z * arb_hypgeom.one_f_two(m + 1, (3/2, m + 2), (z / 2)^2; prec = prec) / (2 * m + 2))")
        end
    end

    while true
        term = ArbReal("$(m_binomial_coeff(m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3) * one_f_two(m))")

        # Break loop if term smaller than accuracy of result.
        if abs(term) < eps(result)
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
        end
    end
    ArbReal(result, bits = prec + 8)
end

n = 3
a = 6.23
β = 3.2
z_range = 0.01:1:200
h = 1

@time original_integral = [hyperbolic_integral(z, n, β, a, h) for z in z_range]
@time m_expansion = [hyperbolic_integral_mexpansion(z, n, β, a, h) for z in z_range]
@time hypgeom_expansion = [hyperbolic_hypergeom_mexpansion(z, n, β, a, h; prec = 64) for z in z_range]
@time hypgeom_expansion2 = [hyperbolic_hypergeom_mexpansion2(z, n, β, a, h; prec = 64) for z in z_range]
p = plot(z_range, hypgeom_expansion, yaxis=:log, label="1f2 expansion")
plot!(z_range, m_expansion, label="m expansion")
plot!(z_range, original_integral, label="original")
display(p)
@show(original_integral, m_expansion, hypgeom_expansion)
