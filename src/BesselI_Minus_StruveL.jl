using ArbNumerics
using Plots
plotly()
Plots.PlotlyBackend()

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
        if abs(term) < eps(result)^2
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

x = 1:200
@time S = abs.(BesselI_minus_StruveL.(1, x; prec = 256))
p = plot(x, S, yaxis=:log)
display(p)
@show(S)

# @time S_test = BesselI_minus_StruveL(1, 1/2; prec = 512)
# t_test = parse(BigFloat, "-0.4326676496002076900053886410785092156777939805674897110150085171427097841422348277548777662847685061084103473774743461230956021691253453509671770438054326932177452541623496380111820545155461149382516536713540154883355853639545855374918858504158554993")
# @show(S_test, t_test)
