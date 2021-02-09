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

        # Calculate nth term of expansion from equation (16).
        current =
            -(1 / gamma(n + 3 / 2)) *
            QuadGK.quadgk(x -> integrand(x, n), 0.0, Inf)[1]

        # Remember nth term and add it to the total sum.
        next_current = current
        total_sum += next_current

        # Go to next term.
        n += 1
    end

    # Return the total sum.
    return total_sum
end

"""
M(n::Int, x::Float64)

    Calculates BesselI[n, x] - StruveL[-n, x] using a series expansion of Gamma functions. n is the order and x the argument to the functions. Adaptive precision is used to produce the correct value, starting at the precision of the input argument x. Used for evaluating ℜχ at finite temperatures.
"""
function M(n, x)

    # Get starting precision.
    p = Int64(precision(x))
    setextrabits(0) # Sets the number of extra bits used to zero.
    setworkingprecision(ArbReal, bits = p)

    # Initialise n and x as ArbReal types for arbitrary precision.
    n = ArbReal(n, bits = p)
    x = ArbReal(Float64(x), bits = p)

    # Initialise the total sum, current sum and k values as ArbReal types.
    total_sum = ArbReal(0.0, bits = p)
    next_current = ArbReal(0.0, bits = p) # Used for comparisons for convergence.
    k = ArbReal(0, bits = p) # Is summed over.

    # Infinite summation. Here just limit summation to when next term adds negiglible amount to total sum.
    while abs(next_current) >= abs(total_sum) * 1e-8

        # Calculate kth term in expansion in terms of gamma functions at set precision p.
        current =
            sign(ArbReal(x, bits = p)) *
            ArbReal(x / 2, bits = p)^ArbReal(2 * k + n, bits = p) / (
                ArbNumerics.gamma(ArbReal(k + 1, bits = p)) *
                ArbNumerics.gamma(ArbReal(n + k + 1, bits = p))
            ) -
            ArbNumerics.ArbReal(
                x / 2,
                bits = p,
            )^ArbReal(2 * k - n + 1, bits = p) / (
                ArbNumerics.gamma(ArbReal(k + 3 / 2, bits = p)) *
                ArbNumerics.gamma(ArbReal(-n + k + 3 / 2, bits = p))
            )

        # Remember kth term for convergence comparison and add to total sum.
        next_current = ArbReal(current, bits = p)
        total_sum += ArbReal(next_current, bits = p)

        # Check the accuracy of kth term. If the order of the error is within 30 orders of magnitude of the estimated value, start entire summation again with incrementally doubled precision until the error is more than 10^30 smaller than the estimated value. (This could be massively optimised).
        if floor(log10(abs(ball(total_sum)[1]))) ==
           floor(log10(abs(ball(total_sum)[2]))) + 30
            p *= 2
            setworkingprecision(ArbReal, bits = p)
            total_sum = ArbReal(0.0, bits = p)
            next_current = ArbReal(0.0, bits = p)
            k = ArbReal(-1)
        end
        k += ArbReal(1)
    end

    # Return the total value.
    return total_sum
end

"""
ℜχ(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64, N::Int)

    Calculate the real part of χ(Ω) at finite temperatures for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function ℜχ(Ω, β, α, v, w)

    # Set arguments to BigFloat precision. Without this the calculations break down due to large values of hyperbolic functions.
    Ω = BigFloat(Ω)
    β = BigFloat(β)
    α = BigFloat(α)
    v = BigFloat(v)
    w = BigFloat(w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    # Coefficient independent of summation variables n & k.
    coefficient = 2 * α * v^3 * β^(3 / 2) / (3 * sqrt(π) * sinh(β / 2) * w^3)

    # Initialise total value of double summation as a BigFloat.
    total_sum = BigFloat(0.0)
    next_current = BigFloat(0.0)
    n = 0

    while abs(next_current) >= abs(total_sum) * 1e-3 # Infinite summation from the Binomial expansion of (x^2 + a^2 - b * cos(v * x))^(-3/2). |b * cos(v * x) / (x^2 + a^2)| < 1 for v > 0 and β > 0. Here just limit summation to when next term adds negiglible amount to total sum.

        current = BigFloat(0.0)

        # Coefficient that depends on n summation variable.
        n_coefficient = -π * sinh(Ω * β / 2) * (-b)^n / (4 * a)^(n + 1)

        # Finite sum over k for even n from (cos(v * x))^n expansion.
        for k = -1:floor((n - 1) / 2)

            # Coefficients dependent upon k.
            k_coeff = vcat(
                repeat(
                    [1 / (SpecialFunctions.gamma(k + 1) * SpecialFunctions.gamma(n - k + 1)),],
                    4,
                ),
                repeat(
                    [(1 + (-1)^n) * 2^n * gamma(n / 2 + 1 / 2) / (2 * sqrt(π) * SpecialFunctions.gamma(n / 2 + 1)),],
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
                current += M(n + 1, a * abs(z)) * sign(z) * abs(z)^(n + 1) * k_c
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
    I1 = quadgk(
        x -> (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2),
        BigFloat(0),
        BigFloat(β / 2),
        maxevals = 10^4,
        order = 7,
    )[1]

    # Return final value obtained from double summation and hyperbolic integral.
    return coefficient * (total_sum + I1)
end

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
