# OedipusRex.jl
# - Polaron optical absorption
#     A WORK IN PROGRESS!

"""
    ReX(v=20,w=20)

Codes to implement 'Optical Absorption of Polarons in the
Feynman-Hellwarth-Iddings-Platzman Approximation',
https://doi.org/10.1103/PhysRevB.5.2367 ;

In particular we are trying to calculate Re[χ] , the real part of the polaron optical absorption.
# (13) - defn of Re[χ]

"""
# function ReX(v=20,w=20)
#     R=(v^2-w^2)/(w^2*v)
#
#     Reχintegrand(u,Ω) =  (1-cos(Ω*u)*exp(im*u))/(R*(1-exp(im*v*u))-im*u)^(3/2)
#
#     Ω=1
#     #Reχ = quadgk(u->Reχintegrand(u,Ω), 0.0, Inf)
# # OK, problematic as this is a complex (multi valued integration!)
#
#     [  Reχintegrand(u,Ω) for u=0:20 ]
#     # Oh yikes, and it explodes as u->0 takes it to ->0 in the denominator
# end

"""
--------------------------------------------------------------------------------
Finite temperature implementation for ℜχ using BesselI and StruveL functions.
--------------------------------------------------------------------------------

Calculated in a similar nature to ℑχ, however compared to the resultant contour integral for ℑχ, the ℜχ version interchanges a cosine for a sine (which results in using BesselI and StruveL functions instead of just BesselK functions) and also introduces an extra integral of hyperbolic functions. This makes it trickier to evaluate than ℑχ.
"""

"""
ℜχ(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the real part of χ(Ω) in a zero temperature approximation (equation (16) in Devreese's et al.) for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function ℜχ(Ω, α, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)

    # Integrand from equation (16).
    integrand(x, n) =
        (
            (n + 1 / 2) * x^(n - 1 / 2) * exp(-R * x) -
            R * x^(n + 1 / 2) * exp(-R * x)
        ) * log(abs((1 + n * v + x)^2 / (Ω^2 - (1 + n * v + x)^2))^(1 / 2))

    # Initialise total sum, current sum (for convergence comparison) and zeroth term counter.
    total_sum = 0.0
    next_current = 0.0
    n = 0

    # Infinite summation. Here just limit summation to when next term adds negiglible amount to total sum.
    while abs(next_current) >= abs(total_sum) * 1e-3
        coef = R^n
        inte = QuadGK.quadgk(x -> integrand(x, n), 0.0, Inf)[1]
        # Calculate nth term of expansion from equation (16).
        current = (-1)^n * (1 / (SpecialFunctions.gamma(n + 3 / 2) * SpecialFunctions.gamma(-n-1/2) * SpecialFunctions.gamma(n + 1))) * inte

        # Remember nth term and add it to the total sum.
        next_current = current
        total_sum += next_current * coef

        # Go to next term.
        n += 1
    end

    # Return the total sum.
    return 2 * α * (v / w)^3 / 3 * total_sum
end

"""
BesselI_minus_StruveL(n::Int, z::Float64)

    Calculates BesselI[n, z] - StruveL[-n, z] using a series expansion of Gamma functions. n is the order and x the argument to the functions. Adaptive precision is used to produce the correct value, starting at the precision of the input argument x. Used for evaluating ℜχ at finite temperatures.
"""
function BesselI_minus_StruveL(n, z; prec = 64) # z > 0 & n >= 0

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec + 8)

    n = ArbReal("$n")
    z = ArbReal("$z")

    k = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    while true
        bessel_term = ArbReal((z / 2)^(2 * k + n) / (ArbNumerics.gamma(k + 1) * ArbNumerics.gamma(k + n + 1)))
        struve_term = ArbReal((z / 2)^(2 * k - n + 1) / (ArbNumerics.gamma(k + 3//2) * ArbNumerics.gamma(k - n + 3//2)))
        term = bessel_term - struve_term

        # Break loop if term smaller than accuracy of result.
        if abs(term) < eps(result)
            break
        end

        result += term
        k += ArbReal("1")

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if ball(result)[2] > err
            setprecision(ArbReal, precision(result) * 2)
            n = ArbReal("$n")
            z = ArbReal("$z")
            k = ArbReal("0")
            result = ArbReal("0.0")
        end
    end
    ArbReal(result, bits = prec + 8)
end

"""
Quicknote: this is the expansion of the second diverging integral in ReX. Still WIP.
"""
function hyperbolic_integral(Ω, β, α, v, w; prec = 64)

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
    for n in 0:1, m in 0:1, t in 0:700

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

    return total_sum * exp(Ω * β / 2) / 2
end

"""
ℜχ(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64, N::Int)

    Calculate the real part of χ(Ω) at finite temperatures for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function ℜχ(Ω, β, α, v, w)

    # Set arguments to BigFloat precision. Without this the calculations break down due to large values of hyperbolic functions.
    Ω = BigFloat("$Ω")
    β = BigFloat("$β")
    α = BigFloat("$α")
    v = BigFloat("$v")
    w = BigFloat("$w")

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    # Coefficient independent of summation variables n & k.
    coefficient = 2 * α * v^3 * β^(3 / 2) * exp(Ω * β / 2) / (3 * sqrt(π) * sinh(β / 2) * w^3)

    # Initialise total value of double summation as a BigFloat.
    total_sum = BigFloat("0.0")
    next_current = BigFloat("1.0")
    n = BigFloat("0")

    while n < 11 # Infinite summation from the Binomial expansion of (x^2 + a^2 - b * cos(v * x))^(-3/2). |b * cos(v * x) / (x^2 + a^2)| < 1 for v > 0 and β > 0. Here just limit summation to when next term adds negiglible amount to total sum.
        current = BigFloat("0.0")

        # Coefficient that depends on n summation variable.
        n_coefficient = -π * (-b)^n * (1 - exp(-Ω * β))  / (4 * a)^(n + 1) / 2

        # Finite sum over k for even n from (cos(v * x))^n expansion.
        for k = -1:floor((n - 1) / 2)

            # Coefficients dependent upon k.
            k_coeff = vcat(
                repeat(
                    [1 / (SpecialFunctions.gamma(k + 1) * SpecialFunctions.gamma(n - k + 1))],
                    4,
                ),
                repeat(
                    [(1 + mod(n, 2)) / SpecialFunctions.gamma(n / 2 + 1)^2],
                    2,
                ),
            )

            # Arguments of the BesselI - StruveL functions.
            z_args = [
                Ω + 1 + v * (n - 2 * k),
                Ω - 1 + v * (n - 2 * k),
                Ω + 1 - v * (n - 2 * k),
                Ω - 1 - v * (n - 2 * k),
                Ω + 1,
                Ω - 1,
            ]

            # Sum over all arguments with their respective k-dependent coefficients. Function M(n, x) gives 'BesselI[n, x] - StruveL[-n, x]'.
            for (z, k_c) in zip(z_args, k_coeff)
                current += BigFloat("$(BesselI_minus_StruveL(n + 1, a * abs(z); prec = 200))") * sign(z) * abs(z)^(n + 1) * k_c
            end
        end


        # Add nth contribution to the total sum.
        total_sum += n_coefficient * current

        # Remember nth contribution for convergence comparison.
        next_current = n_coefficient * current

        # Move to next nth term.
        n += 1
    end

    # Hyperbolic integral that makes up the second term of ℜχ.
    integrand(x) = (1 - cosh(Ω * β * (1 - x) / 2)) * cosh(x * β / 2) / (a^2 - β^2 * x^2 / 4 - b * cosh(v * β * x / 2))^(3 // 2)

    I1 = β * quadgk(x -> integrand(x), BigFloat("0.0"), BigFloat("1.0"))[1] / (2 * exp(Ω * β / 2))

    # I1 = hyperbolic_integral_six(Ω, β, α, v, w; prec = 200) / exp(Ω * β / 2)

    # Return final value obtained from double summation and hyperbolic integral.
    println("Frequency = $Ω \n Total sum = $(total_sum) \n Integral = $(I1) \n Difference = $(total_sum + I1) \n\n")
    return total_sum, I1, (total_sum + I1)
end

"""
In-file tests and plotting for ReX
"""

using ArbNumerics
using SpecialFunctions
using QuadGK
using Plots
plotly()
Plots.PlotlyBackend()
setprecision(BigFloat, 128)
Ω_range = 0.01:1:20
RealX = [ℜχ(Ω, 4, 7, 5.8, 1.6) for Ω in Ω_range]
ReX1 = [abs(x[1]) for x in RealX]
ReX2 = [abs(x[2]) for x in RealX]
ReX3 = [abs(x[3]) for x in RealX]
# p = plot(Ω_range, abs.(ReX1), yaxis=:log, label = "Struve expansion", xlabel = "Ω")
# plot!(p, Ω_range, abs.(ReX2), label = "Hyperbolic int")
# plot!(p, Ω_range, abs.(ReX3), label = "ReX")
# display(p)
q = plot(Ω_range, ReX1, yaxis = :log, label = "Expansion", xlabel = "Ω", ylabel = "ReX")
plot!(Ω_range, ReX2, label="Integral")
plot!(Ω_range, ReX3, label="Sum")
# plot!(Ω_range, [abs(ℜχ(Ω, 7, 5.8, 1.6)) for Ω in Ω_range], label="Struve expansion")

display(q)

"""
----------------------------------------------------------------------
Polaron absorption coefficient Γ(Ω).
----------------------------------------------------------------------
"""

"""
optical_absorption(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the absorption coefficient Γ(Ω) for the polaron at zero-temperature (equation (11a) in Devreese's et al.) for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function optical_absorption(Ω, α, v, w)
    Reχ = ℜχ(Ω, α, v, w)
    Imχ = ℑχ(Ω, α, v, w)
    Ω * Imχ / (Ω^4 - 2 * Ω^2 * Reχ + Reχ^2 + Imχ^2)
end

"""
optical_absorption(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the absorption coefficient Γ(Ω) for the polaron at finite temperatures (equation (11b) in Devreese's et al.) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function optical_absorption(Ω, β, α, v, w)
    ((Ω^2 - w^2)^2 / (Ω^5 - Ω^3 * v^2)^2) * ℑχ(Ω, β, α, v, w)
end
