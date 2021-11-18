# arb_hypgeom.jl

module arb_hypgeom
export one_f_two_fast, one_f_two

using ArbNumerics

"""
Implementation of the hypergeometric function 1F2 using the expansion:

    1F2(a1; {b1, b2}; z) = ∑_{k=0}^∞ (a1)_k z^k / ((b1)_k ⋅ (b2)_k ⋅ k!)

where (x)_k is the Pochhammer rising factorial = Γ(x+k)/Γ(x) with Γ(x) the Gamma function.
"""

function one_f_two(a, b, x; prec = 64) # z > 0 & n >= 0

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec + 8)

    a = ArbReal("$a")
    b = (ArbReal("$i") for i in b)
    x = ArbReal("$x")

    k = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    risingfact(w, n) = ArbReal(ArbNumerics.gamma(w + n) / ArbNumerics.gamma(w))

    while true
        term = ArbReal(risingfact(a, k) * x^k / (prod(risingfact.(b, k)) * ArbNumerics.gamma(k + 1)))

        # Break loop if term smaller than accuracy of result.
        if abs(term) < err
            break
        end

        result += term
        k += ArbReal("1")

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if ball(result)[2] > err
            setprecision(ArbReal, precision(result) * 2)
            k = ArbReal("0")
            a = ArbReal("$a")
            b = (ArbReal("$i") for i in b)
            x = ArbReal("$x")
            result = ArbReal("0.0")
        end
    end
    ArbReal(result, bits = prec + 8)
end

"""
Implementation of specific 1F2 hypergeometric function used for integrals:

    ∫_0^1 cosh(zx) x^{2m} dx = 1F2(m+1/2; {1/2, m+3/2}; z^2/4)/(2m+1)
                             = ∑_{t=0}^∞ z^{2t} / ((2t+2m+1)⋅(2t)!)

    ∫_0^1 sinh(zx) x^{2m} dx = z ⋅ 1F2(m+1; {3/2, m+2}; z^2/4)/(2m+2)
                             = ∑_{t=0}^∞ z^{2t+1} / ((2t+2m+2)⋅(2t+1)!)

which we can combine into a generic summation:

    ∑_{t=0}^∞ z^{2t+h} / ((2t+2m+1+h)⋅(2t+h)!)

which gives the cosh version for h=0 and the sinh version for h=1. This specialised 1F2 converges faster than generic 1F2 algorithm, one_f_two.
"""

function one_f_two_fast(x, m, h; prec = 64) # z > 0 & n >= 0

    # Initialise precision of ArbReal to prec.
    p = prec
    setextrabits(0)
    setprecision(ArbReal, p)

    x = ArbReal("$x")
    m = ArbReal("$m")
    h = ArbReal("$h")

    k = ArbReal("0")
    result = ArbReal("0.0")
    term = ArbReal("1.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    while abs(midpoint(term)) > err * abs(midpoint(result))
        term = ArbReal(x^(2 * k + h) / ((2 * k + 2 * m + 1 + h) * ArbNumerics.gamma(2 * k + 1 + h)))
        result += term
        # println("term: k = ", k, "\nterm value: ", ball(ArbReal(term, bits = prec)), "\ncumulant result: ", ball(ArbReal(result, bits = prec)), "\n")
        k += 1

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if radius(result) > err * abs(midpoint(result))
            p *= 2
            setprecision(ArbReal, p)
            # println("Not precise enough. Error = ", abs(radius(result)/midpoint(result)), " > ", err, ". Increasing precision to ", p, " bits.\n")
            x = ArbReal("$x")
            m = ArbReal("$m")
            h = ArbReal("$h")
            k = ArbReal("0")
            result = ArbReal("0.0")
            term = ArbReal("1.0")
        end
    end
    # println("x: ", ArbReal(x, bits = prec), ". Final result: ", ArbReal(result, bits = prec))
    ArbReal(result, bits = prec)
end

# β = ArbReal("2.0")
# α = ArbReal("7.0")
# v = ArbReal("5.8")
# w = ArbReal("1.6")
# R = ArbReal((v^2 - w^2) / (w^2 * v))
# a = ArbReal(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))
# z = ArbReal("90.61")
# m = ArbReal("1.0")
# @time c = one_f_two_fast(z, m, 0; prec = 64)
# @show(c)

end # end module

################################################################################

"""
Some plots and tests of the above functions to check that they produce the correct values to the specified precision.
"""

# using Plots
# using QuadGK
# plotly()
# Plots.PlotlyBackend()
#
# cosh_integral(z, m) = quadgk(t -> (m + 1) * cosh(z * t) * t^m, BigFloat(0), BigFloat(1))[1]
#
# # x_range = 0:200
# # m = 1
# # @time I = [cosh_integral(x, m) for x in x_range]
# # @time F = abs.([arb_hypgeom_1f2(m / 2 + 1 / 2, (1 / 2, m / 2 + 3/2), x^2 / 4.0; prec = 128) for x in x_range])
# # p = plot(x_range, F, yaxis=:log)
# # plot!(x_range, I)
# # display(p)
# # @show(I, F)
# # println(m / 2 + 1 / 2, (1 / 2, m / 2 + 3/2), x_range^2 / 4.0)
#
# """
# Conclusion: The integral is a bit faster than the hypergeometric function, but not as accurate. So, probably use 1F2 since the time difference doesn't seem to be major!
# """
#
# # @time F_test = arb_hypgeom_1f2(1/2, (1/2, 1/2), 1/2; prec = 512)
# # setprecision(BigFloat, 512)
# # t_test = parse(BigFloat, "2.17818355660857086398922206782012528343129403292165693281081574094992093930206091916837394310903919052930938811320102906774905335969974706202490049810447181898779259195640840870264051023127290510625691429186995752638681515944856002645204867940492548077249544")
# # @show(F_test, t_test)
# # @show(isequal("$F_test"[end-1], "$t_test"[end-1]))
#
# """
# The 1F2 function above can accurate produce high precision values as shown in the above test. But it can be simplified and optimized mathematically initially as shown below.
# """
#
# x_range = 0:200
# m = 1
# h = 0
# @time L = abs.([arb_hypergeom_1f2_fast(x, m, h; prec = 128) for x in x_range])
# # @time F = abs.([arb_hypgeom_1f2(m + 1 / 2, (1 / 2, m + 3/2), x^2 / 4.0; prec = 128) / (2 * m + 1) for x in x_range])
# p = plot(x_range, L, yaxis=:log)
# # plot!(x_range, F)
# display(p)
# @show(L)
