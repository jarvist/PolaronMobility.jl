# HolsteinTheory.jl
#  A work in progress.

# Compared to Hellwarth's presentation of the Feynman variational approach for the Frohlich
# Hamiltonian, then only thing that changes with the Holstein Hamiltonian is in the
# calculation of the 'B' component (the coupled electron-phonon component).
function holstein_B(v, w, α, a, ω)
    d(τ) = D_imag(τ, v, w)

    integrand(τ) = exp(-τ) * erf(π / a * sqrt(d(τ) / 2))^3 / d(τ)^(3/2)

    integral = quadgk(τ -> integrand(τ), 0, Inf)[1]

    return α * ω * integral / √2 / π^(3/2)
end

function holstein_B(v, w, α, a, ω, β)
    d(τ) = D_imag(τ, v, w, ω, β)

    integrand(τ) = cosh(τ - ω * β / 2) / sinh(ω * β / 2) * erf(π / a * sqrt(d(τ) / 2))^3 / d(τ)^(3/2)

    integral = quadgk(τ -> integrand(τ * ω * β / 2), 0, 1)[1] * ω * β / 2

    return α * ω * integral / √2 / π^(3/2)
end

function holstein_energy(v, w, α, a, ωβ...)
    Ar = A(v, w, ωβ...)
    Br = holstein_B(v, w, α, a, ωβ...)
    Cr = C(v, w, ωβ...)
    total_energy = -(Ar + Br + Cr)
    return total_energy, Ar, Br, Cr
end

function holsteinvw(v::Real, w::Real, αaωβ...; upper_limit=Inf)
    Δv = v .- w
    initial = [Δv + eps(Float64), w]

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [upper_limit, upper_limit]

    # The multiple phonon mode free energy function to minimise.
    f(x) = holstein_energy(x[1] .+ x[2], x[2], αaωβ...)[1]

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
    E, A, B, C = holstein_energy(Δv .+ w, w, αaωβ...)

    # if Optim.converged(solution) == false
    #     @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    # end

    # Return the variational parameters that minimised the free energy.
    return Δv .+ w, w, E, A, B, C
end

# holsteinvw(α, γ, ω; upper_limit=1e6) = holsteinvw(3.4, 2.6, α, γ, ω; upper_limit=upper_limit)
