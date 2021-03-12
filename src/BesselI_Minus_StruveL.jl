using ArbNumerics
using Plots
plotly()
Plots.PlotlyBackend()

function BesselI_minus_StruveL(n, z; prec = 64) # z > 0 & n >= 0

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec)

    k = ArbReal(0)
    result = ArbReal(0.0)
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
        k += ArbReal(1)

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if ball(result)[2] > err
            setprecision(ArbReal, precision(result) * 2)
            k = ArbReal(0)
            result = ArbReal(0.0)
        end
    end
    ArbReal(result, bits = prec)
end

x = 0.001:1:200
@time S = abs.(BesselI_minus_StruveL.(1, x; prec = 512))
p = plot(x, S, yaxis=:log)
display(p)
@show(S)
