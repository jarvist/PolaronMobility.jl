
"""
This file checks whether QuadGK.quadgk() can integrate numerically to arbitrary precision using BigFloat types. To check, I look at the integral:

    ∫_0^1 cosh(xt) dt

which has a known closed analytic form:

    sinh(x)/x,  x!=0

which we can use for comparison. x has to be parsed as a BigFloat type using
BigFloat("\$x") (remove space in actual code) NOT BigFloat(x) which would just convert a Float64 type to a BigFloat (and so would only be accurate to 64 bits). To check that passing BigFloat("\$x") into sinh(x)/x gives an accurate arbitrary precision answer (and not plagued by rounding errors), I checked with wolframalpha answers which does produce arb. prec. answers and they matched and so passing BigFloat("\$x") into sinh(x)/x gives a good check of whether quadgk above works to arb. prec. too.

I confirmed that quadgk can give arb prec answers (I just using it incorrectly previously to give Float64 accurate results), so now using the correct method, I convert hyperbolic_integral_two() from check_hyperbolic_integral.jl:

    β/2 ∫_0^1 (1 - cosh(Ωβ(1-x)/2)) cosh(xβ/2) / (a^2 - β^2 x^2/4 - bcosh(vβx/2))^{3/2} dx

into arb prec version and see what it gives. Spoilers: it's super slow, further confirmation that it's now arb prec.
"""

"""
This function just evaluates the integral:

    ∫_0^1 cosh(xt) dt

using quadgk. To make it arbitrary precision, the limits have to be parsed properly into a BigFloat type using BigFloat("#") and the absolute tolerance of the integration has to be set to the machine accuracy of your required precision. This is obtained using eps(x::BigFloat) (e.g. eps(Float64) ~ 2.22e-16) and set with atol = eps(x::BigFloat). The precision of the BigFloats can then be changed using setprecision(BigFloat, precision::Int) externally.
"""
function cosh_integral(x)
    err = eps(x)
    integrand(t) = cosh(x * t)
    integral = quadgk(t -> integrand(t), BigFloat("0"), BigFloat("1"), atol = err)
end

"""
Just the exact closed analytic answer to the cosh integral used for comparison. x has to be parsed as a BigFloat("#") before passing into the function to evaluate it at higher precision properly.
"""
function exact_form(x)
    sinh(x) / x
end

"""
Check outputs of the two above functions.
"""

setprecision(BigFloat, 128) # set the precision of the BigFloat types.
x_range = BigFloat("5")

# for bits = 128 and x_range = BigFloat("5") if get:
# 0.005659 seconds (13.49 k allocations: 768.068 KiB)
@time exact = [exact_form(x) for x in x_range]
# exact = fill(14.84064211555775179540189439921291311993)
@show(exact)
# 0.374177 seconds (1.00 M allocations: 46.879 MiB, 3.29% gc time)
@time integral = [cosh_integral(x) for x in x_range]
# fill((14.84064211555775179540189439921291311988, 2.350988701644575015937473074444491355637e-38))
@show(integral)
# Can see that the integral uses large memory allocations. Last two digits differing is common for BigFloat precisions (natural rounding error from float arithmatic) but we can see that the answers match to 128 bit precision.

"""
Using the arbitrary precision quadgk method used above, apply this to the hyperbolic integral I have been trying to solve i.e. hyperbolic_integral_one() in check_hyperbolic_integrals.jl. I actually use hyperbolic_integral_two() here as the change of variables to the limits [0,1] makes it easier for quadgk to evaluate at higher βs.

This allocates a lot of memory to evaluating the integral. So I did use @code_typewarn to check that there are no type-stabilities, but it did not find any so I think this large memory allocation is just a result of small atol set by the machine precision. For 128 bits this is atol ~ 1.93e-34.
"""
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

"""
Check output and timing of the above function.
"""

setprecision(BigFloat, 128)
Ω_range = [BigFloat("12.01")]
β = BigFloat("4.0")
α = BigFloat("7.0")
v = BigFloat("5.8")
w = BigFloat("1.6")
# For Ω_range = [BigFloat("12.01")], β = BigFloat("4.0"), α = BigFloat("7.0"), v = BigFloat("5.8"), w = BigFloat("1.6") I get:
# 209.267907 seconds (1.02 G allocations: 37.016 GiB, 5.82% gc time)
@time hyp = hyperbolic_integral_two.(Ω_range, β, α, v, w)
# hyp = BigFloat[-7.214150399279151977978320029665100466133e+09]
@show(hyp)
# This takes a LONG TIME to evaluate, with a massive 1.02 G allocations to allocate 37.016 GiB memory!
