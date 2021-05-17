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
function ReX(v=20,w=20)
    R=(v^2-w^2)/(w^2*v)

    Reχintegrand(u,Ω) =  (1-cos(Ω*u)*exp(im*u))/(R*(1-exp(im*v*u))-im*u)^(3/2)

    Ω=1
    #Reχ = quadgk(u->Reχintegrand(u,Ω), 0.0, Inf)
# OK, problematic as this is a complex (multi valued integration!)

    [  Reχintegrand(u,Ω) for u=0:20 ]
    # Oh yikes, and it explodes as u->0 takes it to ->0 in the denominator
end

"""
NOTE: This is old code that has been replaced by the full χ function.

ℜχ(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the real part of χ(Ω) in a zero temperature approximation (equation (16) in Devreese's et al.) for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
# function ℜχ(Ω, α, v, w)
#
#     # Initialise constants.
#     R = (v^2 - w^2) / (w^2 * v)
#
#     # Integrand from equation (16).
#     integrand(x, n) =
#         (
#             (n + 1 / 2) * x^(n - 1 / 2) * exp(-R * x) -
#             R * x^(n + 1 / 2) * exp(-R * x)
#         ) * log(abs((1 + n * v + x)^2 / (Ω^2 - (1 + n * v + x)^2))^(1 / 2))
#
#     # Initialise total sum, current sum (for convergence comparison) and zeroth term counter.
#     total_sum = 0.0
#     next_current = 0.0
#     n = 0
#
#     # Infinite summation. Here just limit summation to when next term adds negiglible amount to total sum.
#     while abs(next_current) >= abs(total_sum) * 1e-3
#         coef = R^n
#         inte = QuadGK.quadgk(x -> integrand(x, n), 0.0, Inf)[1]
#         # Calculate nth term of expansion from equation (16).
#         current = (-1)^n * (1 / (SpecialFunctions.gamma(n + 3 / 2) * SpecialFunctions.gamma(-n-1/2) * SpecialFunctions.gamma(n + 1))) * inte
#
#         # Remember nth term and add it to the total sum.
#         next_current = current
#         total_sum += next_current * coef
#
#         # Go to next term.
#         n += 1
#     end
#
#     # Return the total sum.
#     return 2 * α * (v / w)^3 / 3 * total_sum
# end
#
# include("hypergeometric_expansion.jl")
# import .hypergeom_exp
#
# include("bessel_minus_struve.jl")
# import .BesselI_minus_StruveL
#
# function arb_binomial(x, y; prec = 64)
#     setextrabits(0)
#     setprecision(ArbReal, prec)
#     x = ArbReal("$x")
#     y = ArbReal("$y")
#     one = ArbReal("1")
#     ArbNumerics.gamma(x + one) / (ArbNumerics.gamma(y + one) * ArbNumerics.gamma(x - y + one))
# end
#
# function ℜχ(Ω, β, α, v, w; prec = 24)
#
#     # Initialise precision of ArbReal to prec.
#     p = prec
#     setextrabits(0)
#     setprecision(ArbReal, prec)
#
#     Ω = ArbReal("$Ω")
#     β = ArbReal("$β")
#     α = ArbReal("$α")
#     v = ArbReal("$v")
#     w = ArbReal("$w")
#
#     R = ArbReal((v^2 - w^2) / (w^2 * v))
#     a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
#     b = ArbReal(R * β / sinh(β * v / 2))
#
#     coefficient = ArbReal(α * β^(3/2) * v^3 / (3 * a * √π * w^3 * sinh(β / 2)))
#
#     M_c(z, n) = hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 0; prec = p)
#     M_s(z, n) = hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 1; prec = p)
#     θ(z, n) = bessel_minus_struve.BesselI_minus_StruveL(n, z, a; prec = p)
#
#     c = ArbReal(cosh(Ω * β / 2))
#     s = ArbReal(sinh(Ω * β / 2))
#     two = ArbReal("2")
#
#     n = ArbReal("0")
#     result = ArbReal("0.0")
#     term  = ArbReal("1")
#     err = 1e-3  # Machine accuracy of specified precision prec.
#
#     while abs(ArbFloat(term)) > err * abs(ArbFloat(result))
#
#         if mod(n, 2) == 0
#             even_term = M_c(1, n)
#             for z in [Ω + 1, Ω - 1]
#                 even_term += (s * (M_s(z, n) + θ(z, n)) - c * M_c(z, n)) / two
#             end
#         else
#             even_term = ArbReal("0.0")
#         end
#
#         k_term = ArbReal("0.0")
#         for k in 0:floor(n / 2 - 1 / 2)
#             k_coeff = arb_binomial(n, k; prec = p)
#             for z in [1 + v * (n - 2 * k), 1 - v * (n - 2 * k)]
#                 k_term += k_coeff * M_c(z, n)
#             end
#             for z in [Ω + 1 + v * (n - 2 * k), Ω - 1 + v * (n - 2 * k), Ω + 1 - v * (n - 2 * k), Ω - 1 - v * (n - 2 * k)]
#                 k_term += k_coeff * (s * (M_s(z, n) + θ(z, n)) - c * M_c(z, n)) / two
#             end
#         end
#
#         term = arb_binomial(-3/2, n; prec = p) * (-b / (two * a))^n * (arb_binomial(n, n/2; prec = p) * even_term + k_term)
#         result += term
#
#         # println("term: n = ", n, "\nterm value: ", ball(coefficient * term), "\ncumulant result: ", ball(coefficient * result), "\n")
#         n += 1
#
#         # Double precision if rounding error in result exceeds accuracy specified by prec.
#         if radius(result) > err * abs(midpoint(result))
#             p *= 2
#             setprecision(ArbReal, p)
#
#             # println("Not precise enough. Error = ", radius(result), " > ", err * abs(midpoint(result)), ". Increasing precision to ", p, " bits.\n")
#
#             Ω = ArbReal("$Ω")
#             β = ArbReal("$β")
#             α = ArbReal("$α")
#             v = ArbReal("$v")
#             w = ArbReal("$w")
#
#             R = ArbReal((v^2 - w^2) / (w^2 * v))
#             a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
#             b = ArbReal(R * β / sinh(β * v / 2))
#
#             n = ArbReal("0")
#
#             M_c(z, n) = hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 0; prec = p)
#             M_s(z, n) = hypergeometric_expansion.hypergeom_exp(β * z / 2, n, β, a, 1; prec = p)
#             θ(z, n) = bessel_minus_struve.BesselI_minus_StruveL(n, z, a; prec = p)
#
#             result = ArbReal("0.0")
#             term  = ArbReal("1")
#
#             c = ArbReal(cosh(Ω * β / 2))
#             s = ArbReal(sinh(Ω * β / 2))
#             two = ArbReal("2")
#         end
#     end
#     println("Frequency: ", ArbReal(Ω, bits = prec), ". Final result: ", ball(ArbReal(coefficient * result, bits = prec)))
#     ArbReal(coefficient * result, bits = prec)
# end


"""
----------------------------------------------------------------------
Polaron absorption coefficient Γ(Ω).
----------------------------------------------------------------------
"""

"""
optical_absorption(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the absorption coefficient Γ(Ω) for the polaron at at finite temperatures (equation (11a) in Devreese's et al.) for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function optical_absorption(Ω, β, α, v, w)
    X = χ(Ω, β, α, v, w)[1]
    ℜχ = real(X)
    ℑχ = imag(X)
    result = Ω * ℑχ / (Ω^4 - 2 * Ω^2 * ℜχ + abs2(X))
    return result
end

"""
optical_absorption(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the absorption coefficient Γ(Ω) for the polaron at finite temperatures (equation (11b) in Devreese's et al.) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
# function optical_absorption(Ω, β, α, v, w)
#     ((Ω^2 - w^2)^2 / (Ω^5 - Ω^3 * v^2)^2) * imag(χ(Ω, β, α, v, w)[1])
# end
