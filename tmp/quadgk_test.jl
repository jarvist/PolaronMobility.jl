
function cosh_integral(x)
    err = eps(x)
    integrand(t) = cosh(x * t)
    integral = quadgk(t -> integrand(t), BigFloat("0"), BigFloat("1"), atol = err)
end

function exact_form(x)
    sinh(x) / x
end

setprecision(BigFloat, 128)
# x_range = BigFloat("5")
# @time exact = [exact_form(x) for x in x_range]
# @show(exact)
# @time integral = [cosh_integral(x) for x in x_range]
# @show(integral)

function hyperbolic_integral_two(Ω, β, α, v, w)

    err = eps(Ω)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α * β^(3 // 2) * v^3 / (3 * sqrt(π) * sinh(β / 2) * w^3)

    integrand(x) = (1 - cosh(Ω * β * (1 - x) / 2)) * cosh(x * β / 2) / (a^2 - β^2 * x^2 / 4 - b * cosh(v * β * x / 2))^(3 // 2)

    integral = quadgk(x -> integrand(x), BigFloat("0.0"), BigFloat("1.0"), atol = err)

    return coefficient * β * integral[1] / 2
end

Ω_range = [BigFloat("0.01")]
β = BigFloat("3.0")
α = BigFloat("7.0")
v = BigFloat("5.8")
w = BigFloat("1.6")
# @code_warntype hyperbolic_integral_two(Ω_range, β, α, v, w)
@time hyp = hyperbolic_integral_two.(Ω_range, β, α, v, w)
@show(hyp)
