# Susceptibility.jl
#  Frequency-dependent susceptibility (and thus mobility) of the Feynman polaron state
#    A WORK IN PROGRESS.

# Data structure to store the results
struct susceptibility
    nu
    ImX
    μ
end
Susceptibility() = susceptibility([], [], [])

"""
function ImX(nurange,v,w,βred,α,ω,mb)

Impedance in (47a) from Feynman1962, directly solving freq dep without taking
Hellwarth1999 limit of v->0 .

Calculates a frequency dependent (over range of nu) susceptibility which can be linked back to mobility.

HERE BE DRAGONS!
Not well tested or validated code; the central numeric integration is highly
oscillatory and becomes intractable for large nu.
"""
function ImX(nurange, v, w, βred, α, ω, mb)
    @printf("\nAin't no Imaginary Impedance like a Feynman 1962 ImX...\n")
    println("ImX Call signature: v: $v w: $w βred: $βred α: $α ω: $ω mb: $mb")
    # Feynman, I love you - but using Nu, v; Omega, w in the same paper + formulas, for similar objects?!
    s = Susceptibility()
    for nu in nurange
        R = (v^2 - w^2) / (w^2 * v)     # FHIP1962, page 1011, eqn (47c). Note this is wrong in some textbooks / 1990s PRB.
        b = R * βred / sinh(βred * v / 2) # FHIP1962, page 1010, eqn (47b)
        a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2)) # FHIP1962, page 1010, eqn (47b)
        k(u, a, b, v, nu) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) * cos(nu * u) # integrand with cos(vu) term, as (47a)

        @printf("Numerical integration of FHIP1962(42a): nu=%.2f ", nu)
        println("   a: $(a) b:$(b) ")

        # This is a simple approximation of part of the
        # https://www.gnu.org/software/gsl/doc/html/integration.html#qawf-adaptive-integration-for-fourier-integrals
        # GSL adaptive method for Fourier integrals
        # We split the integral into a load of integrals, balanced at some of
        # the roots.
        c = (2 * floor(nu) + 1) * π / nu
        println("Fourier integral c: $(c)")
        fourier_range = [c * i for i in 0:2501]

        # Catch Inf c range
        if nu == 0
            fourier_range = [0, Inf]
        end

        # These params tweaked to get the best behaviour at reproducing
        # Mishchenko-Fig4 in a timely manner.
        @time n = quadgk(u -> k(u, a, b, v, nu), fourier_range...,
            maxevals=10^6, rtol=0.0001, atol=1e-15, order=7) # numerical quadrature integration of (2)
        K = n[1]
        err = n[2]
        @printf(" quadgk: K=%g err=%g\n", K, err)

        if K < err # we've lost control of our errors, due to losing the oscillatory war with nu
            break # --> so give up.
        end

        # Full 47a constructed here
        ImX = 2 * α / (3 * sqrt(π)) * βred^(3 / 2) * (sinh(βred * nu / 2)) / sinh(βred / 2) * (v^3 / w^3) * K

        μ = ImX^-1 * (q) / (ω * mb)

        @printf(" %.3f %g %g\n", nu, ImX, μ)

        append!(s.nu, nu)
        append!(s.ImX, ImX)
        append!(s.μ, μ)
    end

    @printf("\n\n")
    return (s)
end
