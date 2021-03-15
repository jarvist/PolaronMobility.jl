using QuadGK
using Plots
using SpecialFunctions
using ArbNumerics
plotly()
Plots.PlotlyBackend()

"""
Original hyperbolic integral to be solved via an expansion.
"""
function hyperbolic_integral(Ω, β, α, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    integrand(x) = (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2)

    integral = quadgk(x -> integrand(x), 0, β / 2)
    # @show(integral[1])
    return coefficient * integral[1]
end

"""
This now implements infinite loops that are broken upon convergence criteria (i.e. next term being less than a fraction of the total result so far). Issue is to high β or orders of loops, some values overflow and return NaN => need ArbFloats next.
"""
function hyperbolic_integral_loop(Ω, β, α, v, w)

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

    result = 0.0
    n = 0

    eps = 1e-3

    while true
        next_term = 0.0
        n_coeff = n_coefficient(n) / 2^n

        k_result = 0.0
        for k in 0:Int64(floor((n - 1) / 2))

            k_coeff = k_binomial_coeff(n, k)
            z_3, z_4 = z_3_args(n, k), z_4_args(n, k)

            m = 0
            m_result = 0.0
            while true

                m_coeff = m_coefficient(n, m)

                t = 0
                t_result = 0.0
                while true
                    integral_3 = sum([2 * J(z, 2 * t, 0, m) * exp(-Ω * β / 2) for z in z_3])
                    integral_4 = sum([J(z, 2 * t, 1, m) - J(z, 2 * t + 1, -1, m) for z in z_4])
                    t_result += (integral_3 - integral_4 / 2)
                    t += 1
                    if abs(integral_3 - integral_4 / 2) < eps^2 * abs(t_result)
                        break
                    end
                end

                m_result += t_result * m_coeff
                m += 1

                if abs(t_result * m_coeff) < eps * abs(m_result)
                    break
                end
            end
            k_result += m_result * k_coeff
        end

        n_bin_coeff_2 = n_binomial_coeff_two(n)

        even_term = 0.0
        if mod(n, 2) == 0

            m = 0
            m_result = 0.0
            while true

                m_coeff = m_coefficient(n, m)

                t = 0
                t_result = 0.0
                while true
                    integral_1 = 2 * J(z_1_args, 2 * t, 0, m) * exp(-Ω * β / 2)
                    integral_2 = sum([J(z, 2 * t, 1, m) - J(z, 2 * t + 1, -1, m) for z in z_2_args])
                    t_result +=  (integral_1 - integral_2 / 2)
                    t += 1

                    if abs(integral_1 - integral_2 / 2) < eps^2 * abs(t_result)
                        break
                    end
                end

                m_result += t_result * m_coeff
                m += 1

                if abs(t_result * m_coeff) < eps * abs(m_result)
                    break
                end
            end
            even_term += m_result * n_bin_coeff_2
        end

        result += (k_result + even_term) * n_coeff
        n += 1

        if abs((k_result + even_term) * n_coeff) < eps * abs(result)
            break
        end
    end

    return coefficient * result * exp(Ω * β / 2) / 2
end

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

    @show(coefficient * total_sum * exp(Ω * β / 2) / 2)

    return coefficient * total_sum * exp(Ω * β / 2) / 2
end

function hyperbolic_integral_ArbFloat(Ω, β, α, v, w)

    # Initialise constants.
    R = ArbFloat("$((v^2 - w^2) / (w^2 * v))")
    a = ArbFloat("$(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))")
    b = ArbFloat("$(R * β / sinh(β * v / 2))")

    coefficient = ArbFloat("$(2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3))")

    n_binomial_coeff_one(n) = ArbFloat("$(-2 * √π / (gamma(-n - 1 / 2) * gamma(n + 1)))")

    m_binomial_coeff(n, m) = ArbFloat("$(gamma(-n - 1 / 2) / (gamma(m + 1) * gamma(-m - n - 1 / 2)))")

    n_coefficient(n) = ArbFloat("$(n_binomial_coeff_one(n) * (2 / β)^(2 * n + 2) * (-b)^n)")

    m_coefficient(n, m) = ArbFloat("$(m_binomial_coeff(n, m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3))")

    k_binomial_coeff(n, k) = ArbFloat("$(gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1)))")

    n_binomial_coeff_two(n) = ArbFloat("$(2^n * gamma((n + 1) / 2) / (√π * gamma(n / 2 + 1)))")

    z_1_args = ArbFloat("1")
    z_2_args = [ArbFloat("$(Ω + 1)"), ArbFloat("$(Ω - 1)")]
    z_3_args(n, k) = [ArbFloat("$(1 + v * (n - 2 * k))"), ArbFloat("$(1 - v * (n - 2 * k))")]
    z_4_args(n, k) = [ArbFloat("$(Ω + 1 + v * (n - 2 * k))"), ArbFloat("$(Ω - 1 + v * (n - 2 * k))"), ArbFloat("$( Ω + 1 - v * (n - 2 * k))"), ArbFloat("$(Ω - 1 - v * (n - 2 * k))")]

    J(z, t, c, m) = ArbFloat("$((β * z / 2)^t * (1 + c * exp(-Ω * β)) / ((2 * m + t + 1) * gamma(t + 1)))")

    result = ArbFloat("0.0")
    n = ArbFloat("0")

    eps = 1e-3

    while n < ArbFloat("3")
        next_term = ArbFloat("0.0")
        n_coeff = ArbFloat("$(n_coefficient(n) / 2^n)")

        m = ArbFloat("0")
        m_result = ArbFloat("0.0")
        while m < ArbFloat("3")

            m_coeff = ArbFloat("$(m_coefficient(n, m))")

            k_result = ArbFloat("0.0")
            for k in ArbFloat("0"):floor(ArbFloat("$((n - 1) / 2)"))

                k_coeff = ArbFloat("$(k_binomial_coeff(n, k))")
                z_3, z_4 = z_3_args(n, k), z_4_args(n, k)

                t = ArbFloat("0")
                t_result = ArbFloat("0.0")
                while t < ArbFloat("500")
                    integral_3 = ArbFloat("$(sum([ArbFloat("$(2 * J(z, 2 * t, 0, m) * exp(-Ω * β / 2))") for z in z_3]))")
                    integral_4 = ArbFloat("$(sum([ArbFloat("$(J(z, 2 * t, 1, m) - J(z, 2 * t + 1, -1, m))") for z in z_4]))")
                    t_result = ArbFloat("$(t_result + integral_3 - integral_4 / 2)")
                    t = ArbFloat("$(t + 1)")
                end

                k_result = ArbFloat("$(k_result + t_result * k_coeff * n_coeff)")
            end

            m_result = ArbFloat("$(m_result + k_result * m_coeff)")
            m = ArbFloat("$(m + 1)")
        end

        other_term = ArbFloat("$(m_result)")

        n_bin_coeff_2 = ArbFloat("$(n_binomial_coeff_two(n))")

        even_term = ArbFloat("0.0")
        if mod(n, 2) == 0

            m = ArbFloat("0")
            m_result = ArbFloat("0.0")
            while m < ArbFloat("3")

                m_coeff = ArbFloat("$(m_coefficient(n, m))")

                t = ArbFloat("0")
                t_result = ArbFloat("0.0")
                while t < ArbFloat("500")
                    integral_1 = ArbFloat("$(2 * J(z_1_args, 2 * t, 0, m) * exp(-Ω * β / 2))")
                    integral_2 = ArbFloat("$(sum([ArbFloat("$(J(z, 2 * t, 1, m) - J(z, 2 * t + 1, -1, m))") for z in z_2_args]))")
                    t_result = ArbFloat("$(t_result + (integral_1 - integral_2 / 2))")
                    t = ArbFloat("$(t + 1)")
                end

                m_result = ArbFloat("$(m_result + t_result * m_coeff)")
                m = ArbFloat("$(m + 1)")
            end
            even_term = ArbFloat("$(even_term + m_result * n_bin_coeff_2 * n_coeff)")
        end

        result = ArbFloat("$(result + (other_term + even_term))")
        # @show(result)

        n = ArbFloat("$(n + 1)")
    end
    @show(result)
    return ArbFloat("$(coefficient * result * exp(Ω * β / 2) / 2)")
end

"""
Original hyperbolic integral with ArbReal
"""
function hyperbolic_integral_arb(Ω, β, α, v, w)

    # Initialise constants.
    R = ArbFloat("$((v^2 - w^2) / (w^2 * v))")
    a = ArbFloat("$(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))")
    b = ArbFloat("$(R * β / sinh(β * v / 2))")

    coefficient = ArbFloat("$(2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3))")

    @show(coefficient)

    integrand(x) = ArbFloat("$((1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2))")

    integral = quadgk(x -> integrand(x), ArbFloat("0"), ArbFloat("$(β / 2)"))
    @show(integral[1])
    return ArbFloat("$(coefficient * integral[1])")
end

# Ω_range = 0.01:1:20
# hyper_int = [(abs(hyperbolic_integral_arb(Ω, 10, 5, 4.0, 2.8))) for Ω in Ω_range]
# # hyper_int_loop = [(hyperbolic_integral_loop(Ω, 3, 5, 4.0, 2.8)) for Ω in Ω_range]
# # hyper_int_big = [(-hyperbolic_integral_ArbFloat(Ω, 20, 5, 4.0, 2.8)) for Ω in Ω_range]
#
# p = plot(Ω_range, hyper_int, yaxis=:log, label = "hyper_int")
# # plot!(Ω_range, hyper_int_loop, label = "hyper_int_loop")
# # plot!(Ω_range, hyper_int_big, label = "hyper_int_2")
# display(p)

function oscillatory_integral_arb(Ω, β, α, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    integrand(x) = sin(Ω * x) * cos(x) / (a^2 + x^2 - b * cos(v * x))^(3 / 2)

    integral = quadgk(x -> integrand(x), 0, Inf)
    @show(integral[1] * coefficient * sinh(Ω * β / 2))
    return coefficient * integral[1] * sinh(Ω * β / 2)
end

function M(n, a, z, k_c)
    x = ArbReal("$(a * abs(z))")
    bessel = besseli(ArbReal("$(n + 1)"), x)
    digits = Int64(ceil(log10(abs(bessel)))) + 64
    setextrabits(0)
    setprecision(ArbReal, digits = digits)
    setworkingprecision(ArbReal, digits = digits)
    z = ArbReal("$(z)")
    x = ArbReal("$(a * abs(z))")
    total_sum = besseli(ArbReal("$(n + 1)"), x)
    current = ArbReal("$(total_sum + 1.0)")
    k = ArbReal("0")
    while abs(current) > ArbReal("$(abs(total_sum) * 1e-17)")
        # println("Current ", abs(current), "\nTotal sum ", ArbReal("$(total_sum * 1e-3)"))
        current = ArbReal("$((x / 2)^(2 * k - ArbReal("$(n)")) / (gamma(k + 3 / 2) * gamma(k - ArbReal("$(n)") + 1 / 2)))")
        total_sum = ArbReal("$(total_sum - current)")
        k = ArbReal("$(k + 1)")
    end
    println("Arb ", total_sum)
    return ArbReal("$(total_sum * sign(z) * abs(z)^(n + 1) * k_c)")
end

function BesselI_minus_StruveL(n, a, z, k_c)

    x = a * abs(z)

    # Get input precision as we'll return answer to this precision.
    p = Int64(precision(x))
    setprecision(ArbFloat, bits = p)

    # Calculate number of bits needed to fully represent the first term of sum then add input precision to ensure final answer is at least as precise as precision(x).
    new_p = digits2bits(Int64(ceil(log10(abs(besseli(ArbFloat("$(n + 1)"), ArbReal("$x"))))) + 1)) + 2 * p

    # Set working precision to bits. NB that precision(ArbReal) + 8 = precison(BigFloat).
    setextrabits(0)
    setprecision(ArbFloat, new_p)
    setprecision(BigFloat, new_p)

    # Turn inputs and first term (BesselI(n + 1, x)) into ArbReal types because only ArbNumerics allows arbitrary precision bessel calculations.
    x = ArbFloat("$x")
    n = ArbFloat("$n")
    z = ArbFloat("$z")
    k_c = ArbFloat("$k_c")
    total_sum = ArbFloat("$(ArbNumerics.besseli(n + 1, x))")

    # Setup summation terms as ArbReal types. NB I set current (as in the current term in the sum) to 1.0 here instead of 0.0 so the while loop starts.
    current = ArbFloat("1.0")
    k = ArbFloat("0")

    # Loop over each kth term in the gamma expansion of the Struve L function at negative order -(n + 1), and subtract from besselI(n + 1) term (which was set as the first term of total_sum).
    while ArbFloat("$(abs(current))") > ArbFloat("$(abs(total_sum) * 1e-300)")
        current = (x / ArbFloat("2"))^(ArbFloat("2") * k - n) / (ArbNumerics.gamma(k + ArbFloat("3") / ArbFloat("2")) * ArbNumerics.gamma(k - n + ArbFloat("1") / ArbFloat("2")))
        total_sum = ArbFloat("$(total_sum - current)")
        k = ArbFloat("$(k + 1)")
    end

    # Change precision back to input precision.
    setprecision(ArbFloat, p)
    setprecision(BigFloat, p)
    # @show(ArbFloat("$(total_sum * sign(z) * abs(z)^(n + 1) * k_c)"))
    # Return as BigFloat type with precision p.
    return ArbFloat("$(total_sum * sign(z) * abs(z)^(n + 1) * k_c)")
end

function oscillitory_integral_expansion(Ω, β, α, v, w)

    # Initialise constants.
    R = ArbFloat("$((v^2 - w^2) / (w^2 * v))")
    a = ArbFloat("$(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))")
    b = ArbFloat("$(R * β / sinh(β * v / 2))")

    # Coefficient independent of summation variables n & k.
    coefficient = ArbFloat("$(2 * α * v^3 * β^(3 / 2) * exp(Ω * β / 2) / (3 * sqrt(π) * sinh(β / 2) * w^3))")

    # Initialise total value of double summation as a BigFloat.
    total_sum = ArbFloat("0.0")
    next_current = ArbFloat("1.0")
    n =  ArbFloat("0")

    while ArbFloat("$(abs(next_current))") > ArbFloat("$(abs(total_sum) * 1e-3)") # Infinite summation from the Binomial expansion of (x^2 + a^2 - b * cos(v * x))^(-3/2). |b * cos(v * x) / (x^2 + a^2)| < 1 for v > 0 and β > 0. Here just limit summation to when next term adds negiglible amount to total sum.
        current = ArbFloat("0.0")

        # Coefficient that depends on n summation variable.
        n_coefficient = ArbFloat("$(-π * (-b)^n * (1 - exp(-Ω * β))  / (4 * a)^(n + 1) / 2)")

        # Finite sum over k for even n from (cos(v * x))^n expansion.
        for k = ArbFloat("-1"):ArbFloat("$(floor((n - 1) / 2))")

            # Coefficients dependent upon k.
            k_coeff = vcat(
                repeat(
                    [ArbFloat("$(1 / (SpecialFunctions.gamma(k + 1) * SpecialFunctions.gamma(n - k + 1)))")],
                    4,
                ),
                repeat(
                    [ArbFloat("$((1 + mod(n, 2)) / SpecialFunctions.gamma(n / 2 + 1)^2)")],
                    2,
                ),
            )

            # Arguments of the BesselI - StruveL functions.
            z_args = [
                ArbFloat("$(Ω + 1 + v * (n - 2 * k))"),
                ArbFloat("$(Ω - 1 + v * (n - 2 * k))"),
                ArbFloat("$(Ω + 1 - v * (n - 2 * k))"),
                ArbFloat("$(Ω - 1 - v * (n - 2 * k))"),
                ArbFloat("$(Ω + 1)"),
                ArbFloat("$(Ω - 1)"),
            ]

            # Sum over all arguments with their respective k-dependent coefficients. Function M(n, x) gives 'BesselI[n, x] - StruveL[-n, x]'.
            for (z, k_c) in zip(z_args, k_coeff)
                if isnan(k_c)
                    k_c = ArbFloat("0.0")
                end
                current = ArbFloat("$(current + BesselI_minus_StruveL(n, a, z, k_c))")
            end
        end


        # Add nth contribution to the total sum.
        total_sum += ArbFloat("$(n_coefficient * current)")

        # Remember nth contribution for convergence comparison.
        next_current = ArbFloat("$(n_coefficient * current)")

        @show(n, next_current)

        # Move to next nth term.
        n = ArbFloat("$(n + 1)")
    end

    # Return final value obtained from double summation and hyperbolic integral.
    println("Frequency = $Ω \n Total sum = $(ArbFloat("$(coefficient * total_sum)")) \n\n")
    return ArbFloat("$(coefficient * total_sum)")
end

Ω_range = 0.01:3:60
# osc_int = [(abs(oscillatory_integral_arb(Ω, 10, 5, 4.0, 2.8))) for Ω in Ω_range]
# osc_int_exp = [abs(oscillitory_integral_expansion(Ω, 10, 5, 4.0, 2.8)) for Ω in Ω_range]
hyper_int_big = [(-hyperbolic_integral_ArbFloat(Ω, 20, 5, 4.0, 2.8)) for Ω in Ω_range]

# p = plot(Ω_range, osc_int, yaxis=:log, label = "osc_int")
# p = plot(Ω_range, osc_int_exp, label = "osc_int_exp")
plot!(Ω_range, hyper_int_big, label = "hyper_int_2")
display(p)
