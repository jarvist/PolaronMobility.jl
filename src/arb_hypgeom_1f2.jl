# arb_hypgeom_1f2.jl

using ArbNumerics
using Plots
using QuadGK
plotly()
Plots.PlotlyBackend()

function arb_hypgeom_1f2(a, b, z; prec = 64) # z > 0 & n >= 0

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec + 8)

    a = ArbReal("$a")
    b = (ArbReal("$i") for i in b)
    z = ArbReal("$z")

    k = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    risingfact(x, n) = ArbReal(ArbNumerics.gamma(x + n) / ArbNumerics.gamma(x))

    while true
        term = ArbReal(risingfact(a, k) * z^k / (prod(risingfact.(b, k)) * ArbNumerics.gamma(k + 1)))

        # Break loop if term smaller than accuracy of result.
        if abs(term) < eps(result)
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
            z = ArbReal("$z")
            result = ArbReal("0.0")
        end
    end
    ArbReal(result, bits = prec + 8)
end

cosh_integral(z, m) = quadgk(t -> (m + 1) * cosh(z * t) * t^m, BigFloat(0), BigFloat(1))[1]

# x_range = 0:200
# m = 1
# @time I = [cosh_integral(x, m) for x in x_range]
# @time F = abs.([arb_hypgeom_1f2(m / 2 + 1 / 2, (1 / 2, m / 2 + 3/2), x^2 / 4.0; prec = 128) for x in x_range])
# p = plot(x_range, F, yaxis=:log)
# plot!(x_range, I)
# display(p)
# @show(I, F)
# println(m / 2 + 1 / 2, (1 / 2, m / 2 + 3/2), x_range^2 / 4.0)

"""
Conclusion: The integral is a bit faster than the hypergeometric function, but not as accurate. So, probably use 1F2 since the time difference doesn't seem to be major!
"""

F_test = arb_hypgeom_1f2(1/2, (1/2, 1/2), 1/2; prec = 512)
setprecision(BigFloat, 512)
t_test = parse(BigFloat, "2.17818355660857086398922206782012528343129403292165693281081574094992093930206091916837394310903919052930938811320102906774905335969974706202490049810447181898779259195640840870264051023127290510625691429186995752638681515944856002645204867940492548077249544")
@show(F_test, t_test)
@show(isequal("$F_test"[end-1], "$t_test"[end-1]))


"""
The 1F2 function above can accurate produce high precision values as shown in the above test.
"""
