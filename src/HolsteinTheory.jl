# HolsteinTheory.jl

function holstein_B(v, w, α, γ, ω)

    d(τ) = D_imag(τ, v, w)

    integrand(τ) = exp(-τ) * erf(π * sqrt(d(τ) / γ / 2))^3 / d(τ)^(3/2)

    integral = quadgk(τ -> integrand(τ), 0, Inf)[1]

    return α * ω * integral * sqrt(γ / (2π)^3)
end

function holstein_energy(v, w, α, γ, ω)
    Ar = A(v, w, ω)
    Br = holstein_B(v, w, α, γ, ω)
    Cr = C(v, w, ω)
    # total_energy = 3 / 4 * (v - w)^2 / v - Br
    total_energy = -(Ar + Br + Cr)
    return total_energy, Ar, Br, Cr
end

function holsteinvw(v::Real, w::Real, αγωβ...; upper_limit=1e6)

    Δv = v .- w
    initial = [Δv + eps(Float64), w]

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [upper_limit, upper_limit]

    # The multiple phonon mode free energy function to minimise.
    f(x) = holstein_energy(x[1] .+ x[2], x[2], αγωβ...)[1]

    # Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff=:forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS())
    )

    # Extract the v and w parameters that minimised the free energy.
    Δv, w = Optim.minimizer(solution)
    E, A, B, C = holstein_energy(Δv .+ w, w, αγωβ...)

    # if Optim.converged(solution) == false
    #     @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    # end

    # Return the variational parameters that minimised the free energy.
    return Δv .+ w, w, E, A, B, C
end

# holsteinvw(α, γ, ω; upper_limit=1e6) = holsteinvw(3.4, 2.6, α, γ, ω; upper_limit=upper_limit)