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
        current = (-1)^n * (1 / (gamma(n + 3 / 2) * gamma(-n-1/2) * gamma(n + 1))) * inte

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
M(n::Int, x::Float64)

    Calculates BesselI[n, x] - StruveL[-n, x] using a series expansion of Gamma functions. n is the order and x the argument to the functions. Adaptive precision is used to produce the correct value, starting at the precision of the input argument x. Used for evaluating ℜχ at finite temperatures.
"""
# function M(n, a, z, k_c)
#     x = ArbReal(a * abs(z))
#     bessel = besseli(ArbReal(n + 1), x)
#     digits = Int64(ceil(log10(abs(bessel)))) + 64
#     setextrabits(0)
#     setprecision(ArbReal, digits = digits)
#     z = ArbReal(z)
#     x = ArbReal(a * abs(z))
#     total_sum = besseli(ArbReal(n + 1), x)
#     current = ArbReal(1.0)
#     k = ArbReal(0)
#     while abs(current) > 1e-128
#         current = (ArbReal(x) / ArbReal(2))^(ArbReal(2) * k - ArbReal(n)) / (gamma(k + ArbReal(3) / ArbReal(2)) * gamma(k - ArbReal(n) + ArbReal(1) / ArbReal(2)))
#         total_sum -= current
#         k += ArbReal(1)
#     end
#     return total_sum * sign(z) * abs(z)^ArbReal(n + 1) * ArbReal(k_c)
# end

"""
ℜχ(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64, N::Int)

    Calculate the real part of χ(Ω) at finite temperatures for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
# function ℜχ(Ω, β, α, v, w)
#
#     # Set arguments to BigFloat precision. Without this the calculations break down due to large values of hyperbolic functions.
#     Ω = BigFloat(Ω)
#     β = BigFloat(β)
#     α = BigFloat(α)
#     v = BigFloat(v)
#     w = BigFloat(w)
#
#     # Initialise constants.
#     R = (v^2 - w^2) / (w^2 * v)
#     a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
#     b = R * β / sinh(β * v / 2)
#
#     # Coefficient independent of summation variables n & k.
#     coefficient = 2 * α * v^3 * β^(3 / 2) * exp(Ω * β / 2) / (3 * sqrt(π) * sinh(β / 2) * w^3)
#
#     # Initialise total value of double summation as a BigFloat.
#     total_sum = BigFloat(0.0)
#     next_current = BigFloat(1.0)
#     n = 0
#
#     while abs(next_current) > abs(total_sum) * 1e-20 # Infinite summation from the Binomial expansion of (x^2 + a^2 - b * cos(v * x))^(-3/2). |b * cos(v * x) / (x^2 + a^2)| < 1 for v > 0 and β > 0. Here just limit summation to when next term adds negiglible amount to total sum.
#         current = BigFloat(0.0)
#
#         # Coefficient that depends on n summation variable.
#         n_coefficient = -π * (-b)^n * (1 - exp(-Ω * β))  / (4 * a)^(n + 1) / 2
#
#         # Finite sum over k for even n from (cos(v * x))^n expansion.
#         for k = -1:floor((n - 1) / 2)
#
#             # Coefficients dependent upon k.
#             k_coeff = vcat(
#                 repeat(
#                     [1 / (SpecialFunctions.gamma(k + 1) * SpecialFunctions.gamma(n - k + 1)),],
#                     4,
#                 ),
#                 repeat(
#                     [(1 + mod(n, 2)) / SpecialFunctions.gamma(n / 2 + 1)^2,],
#                     2,
#                 ),
#             )
#
#             # Arguments of the BesselI - StruveL functions.
#             z_args = [
#                 Ω + 1 + v * (n - 2 * k),
#                 Ω - 1 + v * (n - 2 * k),
#                 Ω + 1 - v * (n - 2 * k),
#                 Ω - 1 - v * (n - 2 * k),
#                 Ω + 1,
#                 Ω - 1,
#             ]
#
#             # Sum over all arguments with their respective k-dependent coefficients. Function M(n, x) gives 'BesselI[n, x] - StruveL[-n, x]'.
#             for (z, k_c) in zip(z_args, k_coeff)
#                 current += M(n, a, z, k_c)
#             end
#         end
#
#         digits = Int64(ceil(log10(abs(current)))) + 64
#         setextrabits(0)
#         setworkingprecision(ArbReal, digits = digits)
#         current = ArbReal(current)
#
#         # Add nth contribution to the total sum.
#         total_sum += ArbReal(n_coefficient * current)
#
#         # Remember nth contribution for convergence comparison.
#         next_current = ArbReal(n_coefficient * current)
#
#         # Move to next nth term.
#         n += 1
#     end
#
#     digits = Int64(ceil(log10(abs(total_sum)))) + 64
#     setextrabits(0)
#     setworkingprecision(ArbReal, digits = digits)
#     total_sum = ArbReal(total_sum)
#
#     # Hyperbolic integral that makes up the second term of ℜχ.
#     I1 = quadgk(
#         x -> (1 - cosh(Ω * x)) * cosh(x - β / 2)/ exp(Ω * β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2),
#         0, β / 2, atol = 1e-20)
#
#     # x, w = gauss(1000, 0, β / 2)
#     # f(x) = (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2)
#
#     # I1 = ArbReal(sum(f.(x) .* w))
#     I1 = ArbReal(I1[1])
#     coefficient = ArbReal(coefficient)
#
#     # Return final value obtained from double summation and hyperbolic integral.
#     println("Frequency = $Ω \n Total sum = $(coefficient * total_sum) \n Integral = $(coefficient * I1) \n Difference = $(coefficient * total_sum + coefficient * I1) \n\n")
#     return coefficient * total_sum, coefficient * I1, coefficient * (total_sum + I1)
# end

function BesselI_minus_StruveL(n, x)

    # Get input precision as we'll return answer to this precision.
    p = Int64(precision(x))
    setprecision(ArbReal, p + 8)

    # Calculate number of bits needed to fully represent the first term of sum then add input precision to ensure final answer is at least as precise as precision(x).
    new_p = digits2bits(Int64(ceil(log10(abs(besseli(ArbReal("$n") + 1, ArbReal("$x"))))) + 1)) + 2 * p

    # Set working precision to bits. NB that precision(ArbReal) + 8 = precison(BigFloat).
    setextrabits(0)
    setprecision(BigFloat, new_p)
    setprecision(ArbReal, new_p)

    # Turn inputs and first term (BesselI(n + 1, x)) into ArbReal types because only ArbNumerics allows arbitrary precision bessel calculations.
    x = ArbReal("$x")
    n = ArbReal("$n")
    total_sum = ArbReal("$(ArbNumerics.besseli(n + 1, x))")

    # Setup summation terms as ArbReal types. NB I set current (as in the current term in the sum) to 1.0 here instead of 0.0 so the while loop starts.
    current = ArbReal("1.0")
    k = ArbReal("0")

    # Loop over each kth term in the gamma expansion of the Struve L function at negative order -(n + 1), and subtract from besselI(n + 1) term (which was set as the first term of total_sum).
    while ArbReal("$(abs(current))") > ArbReal("$(10.0^(-p - 8))")
        current = (ArbReal("$x") / ArbReal("2"))^(ArbReal("2") * ArbReal("$k") - ArbReal("$n")) / (ArbNumerics.gamma(ArbReal("$k") + ArbReal("3") / ArbReal("2")) * ArbNumerics.gamma(ArbReal("$k") - ArbReal("$n") + ArbReal("1") / ArbReal("2")))
        total_sum = ArbReal("$(total_sum - current)")
        k = ArbReal("$(k + 1)")
    end

    # Change precision back to input precision.
    setprecision(BigFloat, p)
    setprecision(ArbReal, p + 8)

    # Return as BigFloat type with precision p.
    return BigFloat("$total_sum")
end


# m = BesselI_minus_StruveL(10, BigFloat("2062.23523"))
# @show(m)

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
    R = BigFloat("$((v^2 - w^2) / (w^2 * v))")
    a = BigFloat("$(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))")
    b = BigFloat("$(R * β / sinh(β * v / 2))")

    # Coefficient independent of summation variables n & k.
     coefficient = BigFloat("$(2 * α * v^3 * β^(3 / 2) * exp(Ω * β / 2) / (3 * sqrt(π) * sinh(β / 2) * w^3))")

    # Initialise total value of double summation as a BigFloat.
    total_sum = BigFloat("0.0")
    next_current = BigFloat("1.0")
    n = BigInt(0)

    while abs(next_current) > abs(total_sum) * 1e-20 # Infinite summation from the Binomial expansion of (x^2 + a^2 - b * cos(v * x))^(-3/2). |b * cos(v * x) / (x^2 + a^2)| < 1 for v > 0 and β > 0. Here just limit summation to when next term adds negiglible amount to total sum.
        current = BigFloat("0.0")

        # Coefficient that depends on n summation variable.
        n_coefficient = BigFloat("$(-π * (-b)^n * (1 - exp(-Ω * β))  / (4 * a)^(n + 1) / 2)")

        # Finite sum over k for even n from (cos(v * x))^n expansion.
        for k = BigInt(-1):BigInt(floor((n - 1) / 2))

            # Coefficients dependent upon k.
            k_coeff = vcat(
                repeat([BigFloat("$(1 / (SpecialFunctions.gamma(k + 1) * SpecialFunctions.gamma(n - k + 1)))")], 4),
                repeat([BigFloat("$((1 + mod(n, 2)) / SpecialFunctions.gamma(n / 2 + 1)^2)")], 2)
            )

            # Arguments of the BesselI - StruveL functions.
            z_args = [
                BigFloat("$(Ω + 1 + v * (n - 2 * k))"),
                BigFloat("$(Ω - 1 + v * (n - 2 * k))"),
                BigFloat("$(Ω + 1 - v * (n - 2 * k))"),
                BigFloat("$(Ω - 1 - v * (n - 2 * k))"),
                BigFloat("$(Ω + 1)"),
                BigFloat("$(Ω - 1)"),
            ]

            # Sum over all arguments with their respective k-dependent coefficients. Function M(n, x) gives 'BesselI[n, x] - StruveL[-n, x]'.
            for (z, k_c) in zip(z_args, k_coeff)
                current = BigFloat("$(current + BesselI_minus_StruveL(BigInt(n), BigFloat("$(a * abs(z))")) * BigFloat("$(sign(z) * abs(z)^(n + 1) * k_c)"))")
            end
        end

        # Add nth contribution to the total sum.
        total_sum = BigFloat("$(total_sum + n_coefficient * current)")

        # Remember nth contribution for convergence comparison.
        next_current = BigFloat("$(n_coefficient * current)")

        # Move to next nth term.
        n = BigInt(n + 1)
    end

    # Hyperbolic integral that makes up the second term of ℜχ.
    I1 = quadgk(
        x -> BigFloat("$((1 - cosh(Ω * x)) * cosh(x - β / 2)/ exp(Ω * β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2))"),
        BigFloat("0.0"), BigFloat("$(β / 2)"), atol = 1e-8)[1]

    # Return final value obtained from double summation and hyperbolic integral.
    println("Frequency = $Ω \n Total sum = $(coefficient * total_sum) \n Integral = $(coefficient * I1) \n Difference = $(coefficient * total_sum + coefficient * I1) \n\n")
    return BigFloat("$(coefficient * total_sum)"), BigFloat("$(coefficient * I1)"), BigFloat("$(coefficient * (total_sum + I1))")
end

using ArbNumerics
using SpecialFunctions
using QuadGK
using Plots
plotly()
Plots.PlotlyBackend()
setprecision(256)
Ω_range = 3.01:0.1:5
ReX = [ℜχ(Ω, 30, 7, 5.8, 1.6) for Ω in Ω_range]
ReX1 = [log10(abs(x[1] + x[2])) for x in ReX]
# ReX2 = [log10(abs(log10(abs(x[2])))) for x in ReX]
ReX3 = [exp(x[3]) for x in ReX]
# t1 = [ℑχ(Ω, 3, 7, 5.8, 1.6) for Ω in Ω_range]
# p = plot(Ω_range, abs.(ReX1), yaxis=:log, label = "Struve expansion", xlabel = "Ω")
# plot!(p, Ω_range, abs.(ReX2), label = "Hyperbolic int")
# plot!(p, Ω_range, abs.(ReX3), label = "ReX")
# display(p)
q = plot(Ω_range,ReX1, label = "T=0", xlabel = "Ω", ylabel = "ReX")
# plot!(Ω_range, [exp(ℜχ(Ω, 7, 5.8, 1.6)) for Ω in Ω_range], label="Struve expansion")
# plot!(Ω_range, ReX3, label="ReX")
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
