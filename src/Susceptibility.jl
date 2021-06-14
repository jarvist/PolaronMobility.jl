# Susceptibility.jl
#  Frequency-dependent susceptibility (and thus mobility) of the Feynman polaron state
#    A WORK IN PROGRESS.

# Data structure to store the results
struct susceptibility
    nu
    ImX
    μ
end
Susceptibility()=susceptibility([],[],[])

"""
function ImX(nurange,v,w,βred,α,ω,mb)

Impedance in (47a) from Feynman1962, directly solving freq dep without taking
Hellwarth1999 limit of v->0 .

Calculates a frequency dependent (over range of nu) susceptibility which can be linked back to mobility.

HERE BE DRAGONS!
Not well tested or validated code; the central numeric integration is highly
oscillatory and becomes intractable for large nu.
"""
function ImX(nurange,v,w,βred,α,ω,mb)
        @printf("\nAin't no Imaginary Impedance like a Feynman 1962 ImX...\n")
        println("ImX Call signature: v: $v w: $w βred: $βred α: $α ω: $ω mb: $mb")
    # Feynman, I love you - but using Nu, v; Omega, w in the same paper + formulas, for similar objects?!
    s=Susceptibility()
    for nu in nurange
        R=(v^2-w^2)/(w^2*v)     # FHIP1962, page 1011, eqn (47c). Note this is wrong in some textbooks / 1990s PRB.
        b=R*βred/sinh(βred*v/2) # FHIP1962, page 1010, eqn (47b)
        a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2)) # FHIP1962, page 1010, eqn (47b)
        k(u,a,b,v,nu) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u)*cos(nu*u) # integrand with cos(vu) term, as (47a)

        @printf("Numerical integration of FHIP1962(42a): nu=%.2f ",nu)
        println("   a: $(a) b:$(b) ")

        # This is a simple approximation of part of the
        # https://www.gnu.org/software/gsl/doc/html/integration.html#qawf-adaptive-integration-for-fourier-integrals
        # GSL adaptive method for Fourier integrals
        # We split the integral into a load of integrals, balanced at some of
        # the roots.
        c=(2*floor(nu)+1)*π/nu
        println("Fourier integral c: $(c)")
        fourier_range = [ c*i for i in 0:2501 ]

        # Catch Inf c range
        if nu==0
            fourier_range = [0, Inf]
        end

        # These params tweaked to get the best behaviour at reproducing
        # Mishchenko-Fig4 in a timely manner.
        @time n=quadgk(u->k(u,a,b,v,nu),fourier_range... ,
                       maxevals=10^6,rtol=0.0001, atol=1e-15, order=7) # numerical quadrature integration of (2)
        K=n[1]
        err=n[2]
        @printf(" quadgk: K=%g err=%g\n",K,err)

        if K<err # we've lost control of our errors, due to losing the oscillatory war with nu
            break # --> so give up.
        end

        # Full 47a constructed here
        ImX= 2*α/(3*sqrt(π)) * βred^(3/2) * (sinh(βred*nu/2))/sinh(βred/2) * (v^3/w^3) * K

        μ=ImX^-1 * (q)/(ω*mb)

        @printf(" %.3f %g %g\n",nu,ImX,μ)

        append!(s.nu,nu)
        append!(s.ImX,ImX)
        append!(s.μ,μ)
    end

    @printf("\n\n")
    return(s)
end

"""
----------------------------------------------------------------------
Memory Function.
----------------------------------------------------------------------
"""

"""
χ(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) of the polaron at finite temperatures (equation (35a) in FHIP 1962) for a given frequency Ω. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function χ(Ω, β, α, v, w)

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v)

    # FHIP1962, page 1009, eqn (35c).
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) * coth(β * v / 2) + x^2 / β - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36).
    S(x) = 2 * α / (3 * √π) * (exp(1im * x) + 2 * cos(x) / (exp(β) - 1)) / (D(x))^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(x) = (1 - exp(-1im * Ω * x)) * imag(S(x))
    quadgk(x -> integrand(x), 0.0, Inf)[1]
end

"""
χ(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) of the polaron at zero-temperatures (equation (35a) in FHIP 1962) for a given frequency Ω. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function χ(Ω, α, v, w)

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v)

    # FHIP1962, page 1009, eqn (35c) with β → ∞.
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36) with β → ∞.
    S(x) = 2 * α / (3 * √π) * exp(1im * x) / (D(x))^(3 / 2)

    # FHIP1962, page 1009, eqn (35a). Set upper limit < ∞ but very large so cos argument finite.
    integrand(x) = (1 - exp(-1im * Ω * x)) * imag(S(x))
    QuadGK.quadgk(x -> integrand(x), 0.0, 1e200)[1]
end

"""
χ_dc(β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function lim(Ω → 0){χ(Ω) / Ω} of the polaron at finite temperatures (equation (35a) in FHIP 1962) at zero frequency. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function χ_dc(β, α, v, w)

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v)

    # FHIP1962, page 1009, eqn (35c).
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) * coth(β * v / 2) + x^2 / β - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36).
    S(x) = 2 * α / (3 * √π) * (exp(1im * x) + 2 * cos(x) / (exp(β) - 1)) / (D(x))^(3 / 2)

    # Set frequency small enough to mimic Ω = 0 without generating numerical instabilities in integral.
    Ω = 1e-200

    # FHIP1962, page 1009, eqn (35a). Readily divided by Ω to get sinc function in integrand rather than sine which would give 0 all the time.
    integrand(x) = (1 - exp(-1im * Ω * x)) * imag(S(x)) / Ω
    QuadGK.quadgk(x -> integrand(x), 0.0, Inf)[1]
end

"""
----------------------------------------------------------------------
The Moblity of the Polaron.
----------------------------------------------------------------------
"""

"""
μ_ac(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the ac moblity μ(Ω) of the polaron at finite temperatues (equation (46) in FHIP 1962) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function μ_ac(Ω, β, α, v, w)
    Ω / imag(χ(Ω, β, α, v, w))
end

"""
μ_ac(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the ac moblity μ(Ω) of the polaron at zero-temperatures (equation (46) in FHIP 1962) for a given frequency Ω. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function μ_ac(Ω, α, v, w)
    Ω / imag(χ(Ω, α, v, w))
end

"""
μ_dc(β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the dc moblity μ(β) of the polaron at finite temperatures (equation (46) in FHIP 1962) at zero frequenc. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function μ_dc(β, α, v, w)
    1 / imag(χ_dc(β, α, v, w))
end

"""
NOTE: Old code, kept in-case the special functions expansion of Imχ has any future uses. This code, as well as the expansion for Reχ have been replaced with the χ function above.

ℑχ(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the imaginary part of χ(Ω) in a zero temperature approximation (equation (16) in Devreese's et al.) for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function ℑχ(Ω, α, v, w)

    if Ω < 1
        return 0.0
    else
        R = (v^2 - w^2) / (w^2 * v)
        coefficient = 2 / 3 * α * (v / w)^3

        total_sum = 0.0
        next_current = 0.0
        n = 0

        while next_current >= total_sum * 1e-8
            current =
                -coefficient *
                sqrt(π) *
                gamma(n + 1) *
                (-2 * R)^n *
                abs(Ω - 1 - n * v)^(n + 1 / 2) *
                exp(-R * abs(Ω - 1 - n * v)) *
                (1 + sign(Ω - 1 - n * v)) / (
                    gamma(-n - 1 / 2) *
                    BigCombinatorics.doublefactorial(2 * n + 1)
                )
            next_current = current
            total_sum += next_current
            n += 1
        end
        return total_sum
    end
end

"""
ℑχ(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the imaginary part of χ(Ω) at finite temperature using an infinite expansion of BesselK functions (see Appendix A in Devreese's et al.), for a given frequency Ω. β is the thermodynamic beta. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function ℑχ(Ω, β, α, v, w)

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
    coefficient =
        2 * α / (3 * sqrt(π)) * (v / w)^3 * β^(3 / 2) * sinh(Ω * β / 2) /
        sinh(β / 2)

    # Initialise total value of double summation as a BigFloat.
    total_sum = BigFloat(0.0)
    next_current = BigFloat(0.0)
    n = 0

    while next_current >= total_sum * 1e-3 # Infinite summation from the Binomial expansion of (x^2 + a^2 - b * cos(v * x))^(-3/2). |b * cos(v * x) / (x^2 + a^2)| < 1 for v > 0 and β > 0. Here just limit summation to when next term adds negiglible amount to total sum.

        current = BigFloat(0.0)

        # Coefficient that depends on n summation variable.
        n_coefficient = -2 * SpecialFunctions.beta(-1 / 2 - n, n + 1)^(-1) * (-1)^n * b^n * (4 * a)^(-n - 1) * sqrt(π) / SpecialFunctions.gamma(n + 3 / 2)

        # Summation over expansion of (cos(v * x))^n splits into sums over odd n and even n.

        # Finite sum over k for even n from (cos(v * x))^n expansion.
        for k = -1:floor((n - 1) / 2)

            # Coefficient that depends on k summation variable.
            k_coefficient_inverse = vcat(
                repeat([(n + 1) * SpecialFunctions.beta(n - k + 1, k + 1)], 4),
                repeat([(n + 1) * SpecialFunctions.beta(n / 2 + 1, n / 2 + 1)], 2,),
            )

            # Arguments of resultant cosines obtained from cos^n(vx) * cos(x) * cos(Ωx).
            y = [
                (Ω + 1 + v * (n - 2 * k)),
                (Ω + 1 - v * (n - 2 * k)),
                (Ω - 1 + v * (n - 2 * k)),
                (Ω - 1 - v * (n - 2 * k)),
                (Ω + 1),
                (Ω - 1),
            ]

            for (x, k_coeff) in zip(y, k_coefficient_inverse)  # Iterate over cosine arguments for brevity.

                """
                Integral over cos(yx) / (x^2 + a^2)^(n + 3/2) can be identified with BesselK functions.
                """

                # Sometimes SpecialFunctions.besselk overflow, so catch error if it does and instead use ArbNumerics.besselk.
                err = Nothing
                bessel = try
                    SpecialFunctions.besselk(Int(n + 1), abs(Float64(x)) * Float64(a))
                catch err
                    err
                end

                # If SpecialFunctions.besselk overflows or loses accuracy (i.e. equates to 0.0 which shouldn't happen in this range) use arbitrary precision.
                if err == SpecialFunctions.AmosException(2) || bessel == 0.0

                    setextrabits(0) # Do not use extra bits for evals.
                    p = Int64(precision(coefficient)) # Start at BigFloat prec.
                    bessel = ArbNumerics.besselk(n + 1, abs(ArbFloat(x * a, bits = p)))

                    # Increment precision until overflows / loss of accuracy ceases.
                    while bessel == NaN || bessel == 0.0
                        bessel = ArbNumerics.besselk(n + 1, abs(ArbFloat(x * a, bits = p)))
                        p += 1
                    end
                end

                # Calculate nth term in expansion in terms of BesselK functions.
                current += bessel * abs(x)^(n + 1) / k_coeff
            end
        end

        # Remember nth term for convergence comparison and add to total sum.
        next_current = n_coefficient * current
        total_sum += n_coefficient * current

        # Go to next n.
        n += 1
    end

    # Return final value obtained from double summation.
    return coefficient * total_sum
end
