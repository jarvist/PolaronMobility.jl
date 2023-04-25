# HolsteinTheory.jl
#  A work in progress.

# Compared to Hellwarth's presentation of the Feynman variational approach for the Frohlich
# Hamiltonian, then only thing that changes with the Holstein Hamiltonian is in the
# calculation of the 'B' component (the coupled electron-phonon component).
function holstein_B(v, w, α, ω; dims=3)
    d(τ) = D_imag(τ, v, w) * dims / ω

    integrand(τ) = exp(-τ) * (erf(π * sqrt(d(τ))) / sqrt(d(τ)))^(dims)

    integral = quadgk(τ -> integrand(τ), 0, Inf)[1]

    return 2 * α * dims * integral / (2π)^(dims)
end

function holstein_B(v, w, α, ω, β; dims=3)
    d(τ) = D_imag(τ, v, w, ω, β) / dims * 2

    integrand(τ) = cosh(τ - ω * β / 2) / sinh(ω * β / 2) * (erf(π * sqrt(d(τ))) / sqrt(d(τ)))^(dims)

    integral = quadgk(τ -> integrand(τ), 0, β / 2)[1] 

    return 2 * α * dims * ω / (4π)^(dims/2) * integral
end

function holstein_energy(v, w, α, ω; dims=3)
    f = (v - w)^2 / v * ω / 4
    Br = holstein_B(v, w, α, ω; dims=dims)
    total_energy = f - Br
    return total_energy, f, Br
end

function holsteinvw(v::Real, w::Real, α, ω; dims=3, upper_limit=Inf)
    Δv = v .- w
    initial = [Δv + eps(Float64), w]

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [upper_limit, upper_limit]

    # The multiple phonon mode free energy function to minimise.
    f(x) = holstein_energy(x[1] .+ x[2], x[2], α, ω; dims=dims)[1]

    # Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff=:forward),
        lower,
        upper,
        initial,
        Fminbox(LBFGS())
    )

    # Extract the v and w parameters that minimised the free energy.
    Δv, w = Optim.minimizer(solution)
    E, f, B = holstein_energy(Δv .+ w, w, α, ω; dims=dims)

    # if Optim.converged(solution) == false
    #     @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    # end

    # Return the variational parameters that minimised the free energy.
    return Δv .+ w, w, E, f, B
end

# holsteinvw(α, γ, ω; upper_limit=1e6) = holsteinvw(3.4, 2.6, α, γ, ω; upper_limit=upper_limit)
