
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

function arb_hypgeom_1f2(a, b, x; prec = 64) # z > 0 & n >= 0

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
            x = ArbReal("$x")
            result = ArbReal("0.0")
        end
    end
    ArbReal(result, bits = prec + 8)
end

function hypergeom_expansion(z, n, β, a, h; prec = 64)

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec + 8)

    z = ArbReal("$z")
    n = ArbReal("$n")
    β = ArbReal("$β")
    a = ArbReal("$a")

    m = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    m_binomial_coeff(m) = ArbReal("$(ArbNumerics.gamma(-n - 1//2) / (ArbNumerics.gamma(m + 1) * ArbNumerics.gamma(-n - m - 1//2)))")

    if h == 0
        one_f_two = function (m)
            ArbReal("$(arb_hypgeom_1f2(m + 1/2, (1/2, m + 3/2), (z / 2)^2; prec = prec) / (2 * m + 1))")
        end
    elseif h == 1
        one_f_two = function (m)
            ArbReal("$(z * arb_hypgeom_1f2(m + 1, (3/2, m + 2), (z / 2)^2; prec = prec) / (2 * m + 2))")
        end
    end

    while true
        term = ArbReal("$(m_binomial_coeff(m) * (-1)^m * (β / (2 * a))^(2 * n + 2 * m + 3) * one_f_two(m))")

        # Break loop if term smaller than accuracy of result.
        if abs(term) < eps(result)
            break
        end

        result += term
        m += ArbReal("1")

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if ball(result)[2] > err
            setprecision(ArbReal, precision(result) * 2)

            n = ArbReal("$n")
            z = ArbReal("$z")
            β = ArbReal("$β")
            a = ArbReal("$a")

            m = ArbReal("0")
            result = ArbReal("0.0")
        end
    end
    ArbReal(result, bits = prec + 8)
end

function hyperbolic_expansion(Ω, β, α, v, w; prec = 64)

    # Initialise precision of ArbReal to prec.
    setextrabits(0)
    setprecision(ArbReal, prec + 8)

    Ω = ArbReal("$Ω")
    β = ArbReal("$β")
    α = ArbReal("$α")
    v = ArbReal("$v")
    w = ArbReal("$w")

    R = ArbReal("$((v^2 - w^2) / (w^2 * v))")
    a = ArbReal("$(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))")
    b = ArbReal("$(R * β / sinh(β * v / 2))")

    coeff = ArbReal("$(-2 * α * β^(5 / 2) * v^3 / (3 * a^3 * w^3 * sinh(β / 2)))")

    n_coeff(n) = ArbReal("$(-b / (2 * a^2)^n)")

    k_coeff(k) = ArbReal("$(1 / (ArbNumerics.gamma(k + 1) * ArbNumerics.gamma(n - k + 1)))")

    M_c(n, z) = ArbReal("$(hypergeom_expansion(z, n, β, a, 0; prec = prec))")

    M_s(n, z) = ArbReal("$(hypergeom_expansion(z, n, β, a, 1; prec = prec))")


    z_1 = ArbReal("1")
    z_2 = [ArbReal("$(Ω + 1)"), ArbReal("$(Ω - 1)")]
    z_3(n, k) = [ArbReal("$(1 + v * (n - 2 * k))"), ArbReal("$(1 - v * (n - 2 * k))")]
    z_4(n, k) = [ArbReal("$(Ω + 1 + v * (n - 2 * k))"), ArbReal("$(Ω - 1 + v * (n - 2 * k))"), ArbReal("$(Ω + 1 - v * (n - 2 * k))"), ArbReal("$(Ω - 1 - v * (n - 2 * k))")]

    c = ArbReal("$(cosh(Ω * β / 2))")
    s = ArbReal("$(sinh(Ω * β / 2))")

    n = ArbReal("0")
    result = ArbReal("0.0")
    err = eps(result)  # Machine accuracy of specified precision prec.

    while n < 4

        if mod(n, 2) == 0
            even_term = ArbReal("$((M_c(n, 1) - sum(c .* M_c.(n, z_2) - s .* M_s.(n, z_2)) / 2) / (ArbNumerics.gamma(n / 2 + 1)^2))")
        else
            even_term = ArbReal("0.0")
        end

        k_term = ArbReal("0.0")
        for k in ArbReal("-1"):ArbReal("$(floor(n / 2 - 1 / 2))")
            k_coeff = ArbReal("$(ArbNumerics.gamma(k + 1) * ArbNumerics.gamma(n - k + 1))")
            if k_coeff == 0 || isnan(k_coeff)
                k_term += ArbReal("0.0")
            else
                k_term += ArbReal("$((sum(M_c.(n, z_3(n, k))) - sum(c .* M_c.(n, z_4(n, k)) - s .* M_s.(n, z_4(n, k))) / 2) / k_coeff)")
            end
        end

        term = ArbReal("$(n_coeff(n) * (even_term + k_term))")

        # Break loop if term smaller than accuracy of result.
        if abs(term) < eps(result)
            break
        end

        result += ArbReal("$term")
        @show(n, term, result * coeff)
        n += ArbReal("1")

        # Double precision if rounding error in result exceeds accuracy specified by prec.
        if ball(result)[2] > err
            setprecision(ArbReal, precision(result) * 2)

            Ω = ArbReal("$Ω")
            β = ArbReal("$β")
            α = ArbReal("$α")
            v = ArbReal("$v")
            w = ArbReal("$w")

            R = ArbReal("$((v^2 - w^2) / (w^2 * v))")
            a = ArbReal("$(sqrt(β^2 / 4 + R * β * coth(β * v / 2)))")
            b = ArbReal("$(R * β / sinh(β * v / 2))")

            coeff = ArbReal("$(-2 * α * β^(5 / 2) * v^3 / (3 * a^3 * w^3 * sinh(β / 2)))")

            n_coeff(n) = ArbReal("$(-b / (2 * a^2)^n)")

            k_coeff(k) = ArbReal("$(1 / (ArbNumerics.gamma(k + 1) * ArbNumerics.gamma(n - k + 1)))")

            M_c(n, z) = ArbReal("$(hypergeom_expansion(z, n, β, a, 0; prec = prec))")

            M_s(n, z) = ArbReal("$(hypergeom_expansion(z, n, β, a, 1; prec = prec))")

            z_1 = ArbReal("1")
            z_2 = [ArbReal("$(Ω + 1)"), ArbReal("$(Ω - 1)")]
            z_3(n, k) = [ArbReal("$(1 + v * (n - 2 * k))"), ArbReal("$(1 - v * (n - 2 * k))")]
            z_4(n, k) = [ArbReal("$(Ω + 1 + v * (n - 2 * k))"), ArbReal("$(Ω - 1 + v * (n - 2 * k))"), ArbReal("$(Ω + 1 - v * (n - 2 * k))"), ArbReal("$(Ω - 1 - v * (n - 2 * k))")]

            c = ArbReal("$(cosh(Ω * β / 2))")
            s = ArbReal("$(sinh(Ω * β / 2))")

            n = ArbReal("0")
            result = ArbReal("0.0")
        end
    end
    result = ArbReal("$(result * coeff)")
    @show(result)
    ArbReal(result, bits = prec + 8)
end

function hyperbolic_integral(Ω, β, α, v, w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 / 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    integrand(x) = (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2)

    integral = quadgk(x -> integrand(x), 0, β / 2)
    @show(coefficient * integral[1])
    return coefficient * integral[1]
end



Ω_range = 0.01:0.1:3
β = 3.0
α = 7.0
v = 5.8
w = 1.6

R_int = [hyperbolic_integral(Ω, β, α, v, w) for Ω in Ω_range]
R_exp = [hyperbolic_expansion(Ω, β, α, v, w; prec = 64) for Ω in Ω_range]
p = plot(Ω_range, R_exp)
plot!(Ω_range, R_int)
display(p)
