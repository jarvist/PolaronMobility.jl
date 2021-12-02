# check_hyperbolic_integrals.jl

# module check_hyperbolic_integral
# export hyperbolic_integral_one, hyperbolic_integral_two, hyperbolic_integral_three, hyperbolic_integral_four, hyperbolic_integralfive, hyperbolic_integral_six

using QuadGK
using Plots
using SpecialFunctions
plotly()
Plots.PlotlyBackend()

"""
This file tests a mathematical expansion of the integral:

    ∫_0^{β/2} (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2) dx

using some change of variables, two binomial expansions, a product of coshines expansion, and an identity that relates simpler resultant integrals to 1F2(a, {b1, b2}; z) hypergeometric functions.

Code is far from tidy or efficient. It is more a proof of implementation.
"""

"""
Original integral:

    ∫_0^{β/2} (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2) dx
"""
function hyperbolic_integral_one(Ω, β, α, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    integrand(x) = (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2)

    integral = quadgk(x -> integrand(x), 0, β / 2)

    return coefficient * integral[1]
end

"""
Transform integral to limits (0, 1):

    β/2 ∫_0^1 (1 - cosh(Ωβ(1-x)/2)) cosh(xβ/2) / (a^2 - β^2 x^2/4 - bcosh(vβx/2))^{3/2} dx
"""
function hyperbolic_integral_two(Ω, β, α, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    integrand(x) = (1 - cosh(Ω * β * (1 - x) / 2)) * cosh(x * β / 2) / (a^2 - β^2 * x^2 / 4 - b * cosh(v * β * x / 2))^(3 / 2)

    integral = quadgk(x -> integrand(x), 0, 1)

    return coefficient * β * integral[1] / 2
end

"""
First binomial expansion of denominator in integrand:

    ∑_{n=0}^∞ C(-3/2, n) (2/β)^{2n + 2} (-b)^n ∫_0^1 (1 - cosh(Ωβ(1-x)/2)) cosh(βx/2) cosh^n(βvx/2) / (4a^2/β^2 - x^2)^{n + 3/2} dx.

C(n, k) are binomial coefficients.
"""
function hyperbolic_integral_three(Ω, β, α, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    integrand(x, n) = (1 - cosh(Ω * β * (1 - x) / 2)) * cosh(x * β / 2) * cosh(β * v * x / 2)^n / (1 - x^2 + 4 * R * coth(β * v / 2) / β)^(n + 3 / 2)

    n_binomial_coeff(n) = -2 * √π / (gamma(-n - 1 / 2) * gamma(n + 1))

    n_coefficient(n) = n_binomial_coeff(n) * (2 / β)^(2 * n + 2) * (-b)^n

    total_sum = 0.0

    for n in 0:8 # Seems to converge very quickly.
        integral = quadgk(x -> integrand(x, n), 0, 1)[1]
        total_sum += n_coefficient(n) * integral
    end

    return coefficient * total_sum
end

"""
Second binomial expansion of the denominator in the integrand:

    ∑_{n=0}^∞ C(-3/2, n) (2/β)^{2n + 2} (-b)^n ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} ∫_0^1 (1 - cosh(Ωβ(1-x)/2)) cosh(βx/2) cosh^n(βvx/2) x^{2m} dx
"""
function hyperbolic_integral_four(Ω, β, α, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    integrand(x, n, m) = (1 - cosh(Ω * β * (1 - x) / 2)) * cosh(x * β / 2) * cosh(β * v * x / 2)^n * (x^2)^m

    n_binomial_coeff(n) = -2 * √π / (gamma(-n - 1 / 2) * gamma(n + 1))

    m_binomial_coeff(n, m) = gamma(-n - 1 / 2) / (gamma(m + 1) * gamma(-m - n - 1 / 2))

    n_coefficient(n) = n_binomial_coeff(n) * (2 / β)^(2 * n + 2) * (-b)^n

    m_coefficient(n, m) = m_binomial_coeff(n, m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3)

    total_sum_n = 0.0

    for n in 0:8

        total_sum_m = 0.0

        for m in 0:10
            integral = quadgk(x -> integrand(x, n, m), 0, 1)[1]
            total_sum_m += m_coefficient(n, m) * integral
        end
        total_sum_n += n_coefficient(n) * total_sum_m
    end

    return coefficient * total_sum_n
end

"""
Expand product of coshines in the integrand:

    ∑_{n=0}^∞ C(-3/2, n) (2/β)^{2n + 2} (-b/2)^n ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} {

    C(n, n/2)(1 - nmod2)[∫_0^1 cosh(βx/2) x^{2m} dx - (1/2) ∑_{z_2} (cosh(Ωβ/2) ∫_0^1 cosh(βx(z_2)/2) x^{2m} dx - sinh(Ωβ/2) ∫_0^1 sinh(βx(z_2)/2) x^{2m} dx)]

    + ∑_{k=0}^{floor(n/2-1/2)} C(n, k) [∑_{z_3} ∫_0^1 cosh(βx(z_3)/2) x^{2m} dx - (1/2) ∑_{z_4} (cosh(Ωβ/2) ∫_0^1 cosh(βx(z_4)/2) x^{2m} dx - sinh(Ωβ/2) ∫_0^1 sinh(βx(z_4)/2) x^{2m} dx)]
    }

where:

    z_2 ∈ {Ω + 1, Ω - 1}
    z_3 ∈ {1 + v(n - 2k), 1 - v(n - 2k)}
    z_4 ∈ {Ω + 1 + v(n - 2k), Ω - 1 + v(n - 2k), Ω + 1 - v(n - 2k), Ω - 1 - v(n - 2k)}.
"""
function hyperbolic_integral_five(Ω, β, α, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    n_binomial_coeff_one(n) = -2 * √π / (gamma(-n - 1 / 2) * gamma(n + 1))

    m_binomial_coeff(n, m) = gamma(-n - 1 / 2) / (gamma(m + 1) * gamma(-m - n - 1 / 2))

    n_coefficient(n) = n_binomial_coeff_one(n) * (2 / β)^(2 * n + 2) * (-b)^n

    m_coefficient(n, m) = m_binomial_coeff(n, m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3)

    ############################################################################

    cosh_integral(z, m) = quadgk(x -> cosh(β * x * z / 2) * (x^2)^m, 0, 1)[1]

    cosh_sinh_integral(z, m) = cosh(Ω * β / 2) * quadgk(x -> cosh(β * x * z / 2) * (x^2)^m, 0, 1)[1] - sinh(Ω * β / 2) * quadgk(x -> sinh(β * x * z / 2) * (x^2)^m, 0, 1)[1]

    k_binomial_coeff(n, k) = gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1))

    n_binomial_coeff_two(n) = 2^n * gamma((n + 1) / 2) / (√π * gamma(n / 2 + 1))

    z_1_args = 1
    z_2_args = [Ω + 1, Ω - 1]
    z_3_args(n, k) = [1 + v * (n - 2 * k), 1 - v * (n - 2 * k)]
    z_4_args(n, k) = [Ω + 1 + v * (n - 2 * k), Ω - 1 + v * (n - 2 * k), Ω + 1 - v * (n - 2 * k), Ω - 1 - v * (n - 2 * k)]

    ############################################################################

    total_sum = 0.0
    for n in 0:8

        n_coeff = n_coefficient(n)

        for k in 0:Int64(floor((n - 1) / 2))

            k_coeff = k_binomial_coeff(n, k)
            z_3, z_4 = z_3_args(n, k), z_4_args(n, k)

            for m in 0:10

                integral_3 = sum([cosh_integral(z, m) for z in z_3])
                integral_4 = sum([cosh_sinh_integral(z, m) for z in z_4])

                total_sum += n_coeff * m_coefficient(n, m) * k_coeff * (integral_3 - integral_4 / 2) / 2^n
            end
        end

        n_bin_coeff_2 = n_binomial_coeff_two(n)

        if mod(n, 2) == 0

            for m in 0:10

                integral_1 = cosh_integral(z_1_args, m)
                integral_2 = sum([cosh_sinh_integral(z, m) for z in z_2_args])

                total_sum += n_coeff * m_coefficient(n, m) * n_bin_coeff_2 * (integral_1 - integral_2 / 2) / 2^n
            end
        end
    end

    return coefficient * total_sum
end

"""
Substitute remaining integrals for equivalent summations and extract exponential term:

    e^{Ωβ/2}/2 ∑_{n=0}^∞ C(-3/2, n) (2/β)^{2n + 2} (-b/2)^n ∑_{m=0}^∞ C(-n-3/2, m) (-1)^m (β/2a)^{2n+2m+3} {

    C(n, n/2)(1 - nmod2)[2e^{-Ωβ/2} J(1, 2t, 0) - (1/2) ∑_{z_2} (J(z_2, 2t, 1) - J(z_2, 2t+1, -1))]

    + ∑_{k=0}^{floor(n/2-1/2)} C(n, k) [2e^{-Ωβ/2} ∑_{z_3} J(z_3, 2t, 0) - (1/2) ∑_{z_4} (J(z_4, 2t, 1) - J(z_4, 2t+1, -1))]
    }

where:

    z_2 ∈ {Ω + 1, Ω - 1}
    z_3 ∈ {1 + v(n - 2k), 1 - v(n - 2k)}
    z_4 ∈ {Ω + 1 + v(n - 2k), Ω - 1 + v(n - 2k), Ω + 1 - v(n - 2k), Ω - 1 - v(n - 2k)}

and:

    J(z, t, c) ≡ ∑_{t=0}^∞ (βz/2)^t (1-c⋅e^{-Ωβ}) / ((2m + t + 1)⋅(t!)).

The integral identities are:

    ∫_0^1 cosh(zx) x^{2m} dx = 1F2(m+1/2; {1/2, m+3/2}; z^2/4)/(2m+1)
                             = ∑_{t=0}^∞ z^{2t} / ((2t+2m+1)⋅(2t)!)

    ∫_0^1 sinh(zx) x^{2m} dx = z ⋅ 1F2(m+1; {3/2, m+2}; z^2/4)/(2m+2)
                             = ∑_{t=0}^∞ z^{2t+1} / ((2t+2m+2)⋅(2t+1)!)

where 1F2 is the hypergeometric function pFq with p = 1 and q = 2.
"""
function hyperbolic_integral_six(Ω, β, α, v, w; prec = 64)

    setprecision(BigFloat, prec)

    Ω = BigFloat("$Ω")
    β = BigFloat("$β")
    α = BigFloat("$α")
    v = BigFloat("$v")
    w = BigFloat("$w")

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    n_binomial_coeff_one(n) = -2 * √π / (gamma(-n - 1 / 2) * gamma(n + 1))

    m_binomial_coeff(n, m) = gamma(-n - 1 / 2) / (gamma(m + 1) * gamma(-m - n - 1 / 2))

    n_coefficient(n) = n_binomial_coeff_one(n) * (2 / β)^(2 * n + 2) * (-b)^n

    m_coefficient(n, m) = m_binomial_coeff(n, m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3)

    k_binomial_coeff(n, k) = gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1))

    n_binomial_coeff_two(n) = 2^n * gamma((n + 1) / 2) / (√π * gamma(n / 2 + 1))

    z_1_args = 1
    z_2_args = [Ω + 1, Ω - 1]
    z_3_args(n, k) = [1 + v * (n - 2 * k), 1 - v * (n - 2 * k)]
    z_4_args(n, k) = [Ω + 1 + v * (n - 2 * k), Ω - 1 + v * (n - 2 * k), Ω + 1 - v * (n - 2 * k), Ω - 1 - v * (n - 2 * k)]

    J(z, t, c, m) = (β * z / 2)^t * (1 + c * exp(-Ω * β)) / ((2 * m + t + 1) * gamma(t + 1))

    total_sum = 0.0
    for n in 0:3, m in 0:3, t in 0:100

        n_coeff = n_coefficient(n)
        m_coeff = m_coefficient(n, m)

        for k in 0:Int64(floor((n - 1) / 2))

            k_coeff = k_binomial_coeff(n, k)
            z_3, z_4 = z_3_args(n, k), z_4_args(n, k)

            integral_3 = sum([2 * J(z, 2 * t, 0, m) * exp(-Ω * β / 2) for z in z_3])
            integral_4 = sum([J(z, 2 * t, 1, m) - J(z, 2 * t + 1, -1, m) for z in z_4])
            total_sum += n_coeff * m_coeff * k_coeff * (integral_3 - integral_4 / 2) / 2^n
        end

        n_bin_coeff_2 = n_binomial_coeff_two(n)

        if mod(n, 2) == 0

            integral_1 = 2 * J(z_1_args, 2 * t, 0, m) * exp(-Ω * β / 2)
            integral_2 = sum([J(z, 2 * t, 1, m) - J(z, 2 * t + 1, -1, m) for z in z_2_args])
            total_sum += n_coeff * m_coeff * n_bin_coeff_2 * (integral_1 - integral_2 / 2) / 2^n
        end
    end

    return coefficient * total_sum * exp(Ω * β / 2) / 2
end

# end # end of module

"""
Plotting of the above six functions to compare them. They should all produce the same graph (which they do).
"""
Ω_range = 0.01:0.1:20
@time hyper_int_1 = [hyperbolic_integral_one(Ω, 10, 7, 5.8, 1.6) for Ω in Ω_range]
@time hyper_int_2 = [hyperbolic_integral_two(Ω, 10, 7, 5.8, 1.6) for Ω in Ω_range]
@time hyper_int_3 = [hyperbolic_integral_three(Ω, 10, 7, 5.8, 1.6) for Ω in Ω_range]
@time hyper_int_4 = [hyperbolic_integral_four(Ω, 10, 7, 5.8, 1.6) for Ω in Ω_range]
@time hyper_int_5 = [hyperbolic_integral_five(Ω, 10, 7, 5.8, 1.6) for Ω in Ω_range]
@time hyper_int_6 = [hyperbolic_integral_six(Ω, 10, 7, 5.8, 1.6) for Ω in Ω_range]

p = plot(Ω_range, abs.(hyper_int_1), yaxis=:log, label = "hyper_int_1")
plot!(Ω_range,  abs.(hyper_int_2), label = "hyper_int_2")
plot!(Ω_range,  abs.(hyper_int_3), label = "hyper_int_3")
plot!(Ω_range,  abs.(hyper_int_4), label = "hyper_int_4")
plot!(Ω_range,  abs.(hyper_int_5), label = "hyper_int_5")
plot!(Ω_range,  abs.(hyper_int_6), label = "hyper_int_6")
display(p)

"""
Mini tests of some individual parts of the above six functions. Such as the coshine product expansion, the 1F2 hypergeometric function, etc.
"""

# function inner_hyp_integral(x, n, β, v, w)
#
#     # Initialise constants.
#     R = (v^2 - w^2) / (w^2 * v)
#     a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
#     b = R * β / sinh(β * v / 2)
#
#     m_binomial_coeff(n, m) = gamma(-n - 1 / 2) / (gamma(m + 1) * gamma(-m - n - 1 / 2))
#
#     normal_term = (4 * a^2 / β^2 - x^2)^(-n - 3 / 2)
#
#     total_sum_m = 0.0
#     for m in 0:120
#         m_coefficient = m_binomial_coeff(n, m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3)
#         total_sum_m += x^(2 * m) * m_coefficient
#     end
#
#     return normal_term, total_sum_m
# end
#
# # m = inner_hyp_integral(0.99, 0, 200, 5.8, 1.6)
# # @show(m)
#
# function cosh_expansion(x, n)
#
#     normal = cosh(x)^n
#
#     k_binomial_coeff(n, k) = gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1))
#
#     n_binomial_coeff_two(n) = 2^n * gamma((n + 1) / 2) / (√π * gamma(n / 2 + 1))
#
#     total_sum_k = 0.0
#     for k in 0:Int64(floor((n - 1) / 2))
#
#         total_sum_k += n_binomial_coeff_two(n) * (1 - mod(n, 2)) / 2^n + k_binomial_coeff(n, k) * cosh(x * (n - 2 * k)) / 2^(n - 1)
#     end
#
#     return normal, total_sum_k
# end
#
# # x_range = 0.01:0.01:1.0
# # cosh_exp = [cosh_expansion(x, 2) for x in x_range]
# # cosh_exp_norm = [x[1] for x in cosh_exp]
# # cosh_exp_sum = [x[2] for x in cosh_exp]
# # p = plot(x_range, cosh_exp_norm, yaxis = :log, label = "cosh_exp_norm")
# # plot!(x_range, cosh_exp_sum, label = "cosh_exp_sum")
# # display(p)
#
# function hyperbolic_expansion_small(β, v, n, m)
#
#     integrand(x, n, m) = cosh(x * β / 2) * cosh(β * v * x / 2)^n * (x^2)^m
#     integral = quadgk(x -> integrand(x, n, m), 0, 1)[1]
#
#     ##########################################################################
#
#     z_3(n, k) = [1 + v * (n - 2 * k), 1 - v * (n - 2 * k)]
#
#     cosh_integral(z, m) = quadgk(x -> cosh(β * x * z / 2) * (x^2)^m, 0, 1)[1]
#
#     k_binomial_coeff(n, k) = gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1))
#
#     n_binomial_coeff_two(n) = 2^n * gamma((n + 1) / 2) / (√π * gamma(n / 2 + 1))
#
#     total_sum_k = 0.0
#     for k in 0:Int64(floor((n - 1) / 2))
#
#         integral_3 = sum([cosh_integral(z, m) for z in z_3(n, k)])
#
#         total_sum_k += k_binomial_coeff(n, k) * integral_3
#     end
#
#     if mod(n, 2) == 0
#
#         integral_1 = cosh_integral(1, m)
#
#         total_sum_k += n_binomial_coeff_two(n) * integral_1
#     end
#
#     return integral, total_sum_k / 2^n
# end
#
# # β_range = 0.01:0.1:1
# # hyper_exp_s = [hyperbolic_expansion_small(β, 5.8, 1, 1) for β in β_range]
# # hyper_exp_s_int = [x[1] for x in hyper_exp_s]
# # hyper_exp_s_exp = [x[2] for x in hyper_exp_s]
# # p = plot(β_range, hyper_exp_s_int, yaxis = :log, label = "hyper_exp_s_int")
# # plot!(β_range, hyper_exp_s_exp, label = "hyper_exp_s_exp")
# # display(p)
#
# function hyperbolic_expansion_large(Ω, β, v, n, m)
#
#     integrand(x, n, m) = cosh(Ω * β * (1 - x) / 2) * cosh(x * β / 2) * cosh(β * v * x / 2)^n * (x^2)^m
#     integral = quadgk(x -> integrand(x, n, m), 0, 1)[1]
#
#     ##########################################################################
#
#     z_2 = [1 + Ω, -1 + Ω]
#     z_4(n, k) = [1 + Ω + v * (n - 2 * k), -1 + Ω - v * (n - 2 * k), 1 + Ω - v * (n - 2 * k), -1 + Ω + v * (n - 2 * k)]
#
#     cosh_sinh_integral(z, m) = cosh(Ω * β / 2) * quadgk(x -> cosh(β * x * z / 2) * (x^2)^m, 0, 1)[1] - sinh(Ω * β / 2) * quadgk(x -> sinh(β * x * z / 2) * (x^2)^m, 0, 1)[1]
#
#     k_binomial_coeff(n, k) = gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1))
#
#     n_binomial_coeff_two(n) = 2^n * gamma((n + 1) / 2) / (√π * gamma(n / 2 + 1))
#
#     total_sum_k = 0.0
#     for k in 0:Int64(floor((n - 1) / 2))
#
#         integral_4 = sum([cosh_sinh_integral(z, m) for z in z_4(n, k)])
#
#         total_sum_k += k_binomial_coeff(n, k) * integral_4
#     end
#
#     if mod(n, 2) == 0
#
#         integral_2 = sum([cosh_sinh_integral(z, m) for z in z_2])
#
#         total_sum_k += n_binomial_coeff_two(n) * integral_2
#     end
#
#     return integral, total_sum_k / 2^(n + 1)
# end
#
# # Ω_range = 0.01:0.1:1
# # hyper_exp_l = [hyperbolic_expansion_large(Ω, 2, 5.8, 1, 1) for Ω in Ω_range]
# # hyper_exp_l_int = [x[1] for x in hyper_exp_l]
# # hyper_exp_l_exp = [x[2] for x in hyper_exp_l]
# # p = plot(Ω_range, hyper_exp_l_int, label = "hyper_exp_l_int")
# # plot!(Ω_range, hyper_exp_l_exp, label = "hyper_exp_l_exp")
# # display(p)
#
# function final_cosh_expansion(z, β, m)
#
#     cosh_integral = quadgk(x -> cosh(β * x * z / 2) * (x^2)^m, 0, 1)[1]
#
#     ############################################################################
#
#     J(z, t, c, m) = (β * z / 2)^t / ((2 * m + t + 1) * gamma(t + 1))
#
#     total_sum = 0.0
#     for t in 0:100
#         total_sum += J(z, 2 * t, 0, m)
#     end
#
#     return cosh_integral, total_sum
# end
#
# # z_range = 0.01:0.1:10
# # f_cosh_exp = [final_cosh_expansion(z, 6, 0) for z in z_range]
# # f_cosh_exp_int = [x[1] for x in f_cosh_exp]
# # f_cosh_exp_exp = [x[2] for x in f_cosh_exp]
# # p = plot(z_range, f_cosh_exp_int, yaxis=:log, label = "f_int_exp_int")
# # plot!(z_range, f_cosh_exp_exp, label = "f_int_exp_exp")
# # display(p)
#
# function final_hyperbolic_expansion(Ω, z, β, m)
#
#     cosh_sinh_integral = cosh(Ω * β / 2) * quadgk(x -> cosh(β * x * z / 2) * (x^2)^m, 0, 1)[1] - sinh(Ω * β / 2) * quadgk(x -> sinh(β * x * z / 2) * (x^2)^m, 0, 1)[1]
#
#     ############################################################################
#
#     J(z, t, c, m) = (β * z / 2)^t * (1 + c * exp(-Ω * β)) / ((2 * m + t + 1) * gamma(t + 1))
#
#     total_sum = 0.0
#     for t in 0:100
#         total_sum += J(z, 2 * t, 1, m) - J(z, 2 * t + 1, -1, m)
#     end
#
#     return cosh_sinh_integral, total_sum * exp(Ω * β / 2) / 2
# end
#
# # Ω_range = 0.01:0.1:10
# # f_hyp_exp = [final_hyperbolic_expansion(Ω, 6, 6, 10) for Ω in Ω_range]
# # f_hyp_exp_int = [x[1] for x in f_hyp_exp]
# # f_hyp_exp_exp = [x[2] for x in f_hyp_exp]
# # p = plot(Ω_range, f_hyp_exp_int, yaxis=:log, label = "f_int_exp_int")
# # plot!(Ω_range, f_hyp_exp_exp, label = "f_int_exp_exp")
# # display(p)
