
include("arb_hypgeom.jl")
import .arb_hypgeom

using ArbNumerics
using QuadGK

"""
This file checks the binomial expansion of the denominators of the integrals:

    ∫_0^1 cosh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt
    ∫_0^1 sinh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt

which should give:

    ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 cosh(zt) t^{2m} dt
    ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 sinh(zt) t^{2m} dt

These are terms in hyperbolic_integral_five() from check_hyperbolic_integral.jl that depend on the iterator m of the second binomial expansion in check_hyperbolic_integral.jl (the positive integer n is the iterator from the first binomial expansion). This is reproduced here as we then want to check if substituting:

    ∫_0^1 cosh(zt) t^{2m} dt = 1F2(m+1/2; {1/2, m+3/2}; z^2/4)/(2m+1)
                             = ∑_{t=0}^∞ z^{2t} / ((2t+2m+1)⋅(2t)!)

    ∫_0^1 sinh(zt) t^{2m} dt = z ⋅ 1F2(m+1; {3/2, m+2}; z^2/4)/(2m+2)
                             = ∑_{t=0}^∞ z^{2t+1} / ((2t+2m+2)⋅(2t+1)!)

i.e. one_f_two() or one_f_two_fast() from arb_hypgeom.jl, accurately numerically reproduces the top two integrals. The reason for doing this over the expansions already done in check_hyperbolic_integral.jl (which is only as accurate as Float64) is that we want to ensure these integrals are evaluated accurately to precisions greater than Float64 (beyond the first ~16 digits) using the ArbReal summation algorithm first used in bessel_minus_struve.jl, and to hopefully produce cleaner, more digestible code by separating out all these nested sums rather than writing them out explicitly as one massive confusing functions.

Also, one_f_two() versus one_f_two_fast() are compared for efficiency (@time) and accuracy to the top two integrals.
"""

"""
This function uses quadgk to evalulate the integrals:

    ∫_0^1 cosh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt
    ∫_0^1 sinh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt.

In this function, h=0 gives the cosh version and h=1 gives the sinh version.
"""
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

"""
This function then implements the binomial expansion of the denominator of the integrals:

    ∫_0^1 cosh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt
    ∫_0^1 sinh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt

i.e.:

    ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 cosh(zt) t^{2m} dt
    ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 sinh(zt) t^{2m} dt.

Here h=0 is the cosh version, h=1 is the sinh version.
"""
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

    m_binomial_coeff(m) = ArbReal(ArbNumerics.gamma(-n - 1//2) / (ArbNumerics.gamma(m + 1) * ArbNumerics.gamma(-n - m - 1//2)))

    if h == 0
        integrand = function (m, t)
            cosh(z * t) * t^(2 * m)
        end
    elseif h == 1
        integrand = function (m, t)
            sinh(z * t) * t^(2 * m)
        end
    end

    integral(m) = ArbReal(quadgk(t -> integrand(m, t), 0, 1)[1])

    while true
        term = ArbReal(m_binomial_coeff(m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3) * integral(m))

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

            if h == 0
                integrand = function (m, t)
                    cosh(z * t) * t^(2 * m)
                end
            elseif h == 1
                integrand = function (m, t)
                    sinh(z * t) * t^(2 * m)
                end
            end
        end
    end
    ArbReal(result, bits = prec + 8)
end

"""
This function then implements the binomial expansion of the denominator of the integrals:

    ∫_0^1 cosh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt
    ∫_0^1 sinh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt

i.e.:

    ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 cosh(zt) t^{2m} dt
    ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 sinh(zt) t^{2m} dt.

but also substitutes the integrals:

    ∫_0^1 cosh(zt) t^{2m} dt = ∑_{t=0}^∞ z^{2t} / ((2t+2m+1)⋅(2t)!)
    ∫_0^1 sinh(zt) t^{2m} dt = ∑_{t=0}^∞ z^{2t+1} / ((2t+2m+2)⋅(2t+1)!)

by replacing those integrals in the previous function with one_F_two_fast() from arb_hypgeom.jl. h=0 gives cosh version, h=1 gives sinh version.
"""
function hyperbolic_hypergeom_mexpansion(z, n, β, a, h; prec = 64)
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

"""
This function then implements the binomial expansion of the denominator of the integrals:

    ∫_0^1 cosh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt
    ∫_0^1 sinh(zt) / (4a^2/β^2 - t^2)^(n+3/2) dt

i.e.:

    ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 cosh(zt) t^{2m} dt
    ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 sinh(zt) t^{2m} dt.

but also substitutes the integrals:

    ∫_0^1 cosh(zt) t^{2m} dt = 1F2(m+1/2; {1/2, m+3/2}; z^2/4)/(2m+1)
    ∫_0^1 sinh(zt) t^{2m} dt = z ⋅ 1F2(m+1; {3/2, m+2}; z^2/4)/(2m+2)

by replacing those integrals with one_F_two() from arb_hypgeom.jl. h=0 gives cosh version, h=1 gives sinh version.
"""
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

    m_binomial_coeff(m) = ArbReal(ArbNumerics.gamma(-n - 1//2) / (ArbNumerics.gamma(m + 1) * ArbNumerics.gamma(-n - m - 1//2)))

    if h == 0
        one_f_two = function (m)
            ArbReal(arb_hypgeom.one_f_two(m + 1/2, (1/2, m + 3/2), (z / 2)^2; prec = prec) / (2 * m + 1))
        end
    elseif h == 1
        one_f_two = function (m)
            ArbReal(z * arb_hypgeom.one_f_two(m + 1, (3/2, m + 2), (z / 2)^2; prec = prec) / (2 * m + 2))
        end
    end

    while true
        term = ArbReal(m_binomial_coeff(m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3) * one_f_two(m))

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

            if h == 0
                one_f_two = function (m)
                    ArbReal(arb_hypgeom.one_f_two(m + 1/2, (1/2, m + 3/2), (z / 2)^2; prec = precision(result)) / (2 * m + 1))
                end
            elseif h == 1
                one_f_two = function (m)
                    ArbReal(z * arb_hypgeom.one_f_two(m + 1, (3/2, m + 2), (z / 2)^2; prec = precision(result)) / (2 * m + 2))
                end
            end
        end
    end
    ArbReal(result, bits = prec + 8)
end

"""
Plot and time each of these stages/implementations of the integrals and compare.
"""
using Plots
plotly()
Plots.PlotlyBackend()

n = 3
a = 6.23
β = 3.2
z_range = 0.01:1:200 # z=0 diverges
h = 1

# For n=3, a=6.23, β=3.2, z_range=0.01:1:200, h=1 I get:
# 0.161909 seconds (564.27 k allocations: 23.958 MiB, 5.98% gc time) for me.
@time original_integral = [hyperbolic_integral(z, n, β, a, h) for z in z_range]
# 15.588500 seconds (69.64 M allocations: 3.085 GiB, 13.72% gc time) for me.
@time m_expansion = [hyperbolic_integral_mexpansion(z, n, β, a, h) for z in z_range]
# Below uses one_F_two_fast()
# 150.126938 seconds (652.59 M allocations: 28.175 GiB, 15.40% gc time) for me.
@time hypgeom_expansion = [hyperbolic_hypergeom_mexpansion(z, n, β, a, h; prec = 64) for z in z_range]
# Below uses one_F_two()
# 322.089672 seconds (988.47 M allocations: 66.587 GiB, 7.72% gc time) for me.
@time hypgeom_expansion2 = [hyperbolic_hypergeom_mexpansion2(z, n, β, a, h; prec = 64) for z in z_range]
p = plot(z_range, hypgeom_expansion, yaxis=:log, label="1f2 expansion")
plot!(z_range, m_expansion, label="m expansion")
plot!(z_range, original_integral, label="original")
display(p) # plots all look identical.
@show(original_integral, m_expansion, hypgeom_expansion, hypgeom_expansion2)
# original_integral: [2.839772599386279e-8 ... 4.443248342101801e78]
# m_expansion: ArbReal{72}[2.8397725993862784615e-8 ... 4.44324834210182985933e+78]
# hypgeom_expansion = ArbReal{72}[2.8397725993862784742e-8 ... 4.44324834210182176148e+78]
# hypgeom_expansion2 = ArbReal{72}[2.8397725993862784741e-8 ... 4.4432483421018217615e+78]
