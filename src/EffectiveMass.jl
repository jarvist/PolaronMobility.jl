# EffectiveMass.jl

function polaron_effective_mass(v, w, α, ω, β)

    D(t) = 2 * (v^2 - w^2) * sin(v * t / 2) * sin(v * (t - im * ω * β) / 2) / sinh(v * ω * β / 2) / v^3 - im * w^2 * t * (1 + im * t / β / ω) / v^2

    # # FHIP1962, page 1009, eqn (36).
    S(t) = cos(t - im * β * ω / 2) / sinh(ω * β / 2) / D(t)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t) = imag(S(t)) * t^2

    integral = quadgk(t -> integrand(t/(1-t))/(1-t)^2, 0, 1-eps(Float64))[1]

    mass = α * integral / (3 * √π)

    return mass
end

polaron_effective_mass(v, w, α::Vector, ω::Vector, β) = sum(polaron_effective_mass.(v, w, α, ω, β))

