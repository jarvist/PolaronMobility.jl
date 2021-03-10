using QuadGK
using Plots
using SpecialFunctions
plotly()
Plots.PlotlyBackend()

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

    for n in 0:8
        integral = quadgk(x -> integrand(x, n), 0, 1)[1]
        total_sum += n_coefficient(n) * integral
    end

    return coefficient * total_sum
end

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

    for n in 0:30

        total_sum_m = 0.0

        for m in 0:100
            integral = quadgk(x -> integrand(x, n, m), 0, 1)[1]
            total_sum_m += m_coefficient(n, m) * integral
        end
        total_sum_n += n_coefficient(n) * total_sum_m
    end

    return coefficient * total_sum_n
end

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

    ############################################################################

    z_1 = 1
    z_2 = [1 + Ω, 1 - Ω]
    z_3(n, k) = [1 + v * (n - 2 * k), 1 - v * (n - 2 * k)]
    z_4(n, k) = [1 + Ω + v * (n - 2 * k), 1 - Ω + v * (n - 2 * k), 1 + Ω - v * (n - 2 * k), 1 - Ω - v * (n - 2 * k)]

    total_sum_n = 0.0
    for n in 0:30

        total_sum_m = 0.0
        for m in 0:100

            total_sum_k = 0.0
            for k in 0:Int64(floor((n - 1) / 2))

                integral_3 = sum([cosh_integral(z, m) for z in z_3(n, k)])
                integral_4 = sum([cosh_sinh_integral(z, m) for z in z_4(n, k)])

                total_sum_k += k_binomial_coeff(n, k) * (integral_3 - integral_4 / 2)
            end

            if mod(n, 2) == 0

                integral_1 = cosh_integral(z_1, m)
                integral_2 = sum([cosh_sinh_integral(z, m) for z in z_2])

                total_sum_k += n_binomial_coeff_two(n) * (integral_1 - integral_2 / 2)
            end

            total_sum_m += m_coefficient(n, m) * total_sum_k
        end

        total_sum_n += n_coefficient(n) * total_sum_m / 2^n
    end

    return coefficient * total_sum_n
end

Ω_range = 0.01:0.1:1
hyper_int_1 = [(hyperbolic_integral_one(Ω, 2, 7, 5.8, 1.6)) for Ω in Ω_range]
hyper_int_2 = [(hyperbolic_integral_two(Ω, 2, 7, 5.8, 1.6)) for Ω in Ω_range]
hyper_int_3 = [(hyperbolic_integral_three(Ω, 2, 7, 5.8, 1.6)) for Ω in Ω_range]
hyper_int_4 = [(hyperbolic_integral_four(Ω, 2, 7, 5.8, 1.6)) for Ω in Ω_range]
hyper_int_5 = [(hyperbolic_integral_five(Ω, 2, 7, 5.8, 1.6)) for Ω in Ω_range]

p = plot(Ω_range, hyper_int_1, label = "hyper_int_1")
plot!(Ω_range, hyper_int_2, label = "hyper_int_2")
plot!(Ω_range, hyper_int_3, label = "hyper_int_3")
plot!(Ω_range, hyper_int_4, label = "hyper_int_4")
plot!(Ω_range, hyper_int_5, label = "hyper_int_5")
display(p)

function inner_hyp_integral(x, n, β, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    m_binomial_coeff(n, m) = gamma(-n - 1 / 2) / (gamma(m + 1) * gamma(-m - n - 1 / 2))

    normal_term = (4 * a^2 / β^2 - x^2)^(-n - 3 / 2)

    total_sum_m = 0.0
    for m in 0:120
        m_coefficient = m_binomial_coeff(n, m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3)
        total_sum_m += x^(2 * m) * m_coefficient
    end

    return normal_term, total_sum_m
end

# m = inner_hyp_integral(0.99, 0, 200, 5.8, 1.6)
# @show(m)

function cosh_expansion(x, n)

    normal = cosh(x)^n

    k_binomial_coeff(n, k) = gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1))

    n_binomial_coeff_two(n) = 2^n * gamma((n + 1) / 2) / (√π * gamma(n / 2 + 1))

    total_sum_k = 0.0
    for k in 0:Int64(floor((n - 1) / 2))

        total_sum_k += n_binomial_coeff_two(n) * (1 - mod(n, 2)) / 2^n + k_binomial_coeff(n, k) * cosh(x * (n - 2 * k)) / 2^(n - 1)
    end

    return normal, total_sum_k
end

# x_range = 0.01:0.01:1.0
# cosh_exp = [cosh_expansion(x, 2) for x in x_range]
# cosh_exp_norm = [x[1] for x in cosh_exp]
# cosh_exp_sum = [x[2] for x in cosh_exp]
# p = plot(x_range, cosh_exp_norm, yaxis = :log, label = "cosh_exp_norm")
# plot!(x_range, cosh_exp_sum, label = "cosh_exp_sum")
# display(p)
