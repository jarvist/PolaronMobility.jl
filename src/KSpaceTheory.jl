# KSpaceTheory.jl

function frohlich_matrix(k; α=1, ω=1)
	ω^(3/2) * 2√2 * π * α / k^2
end

function holstein_matrix(k, α, ω, z)
    z * α * ω
end

function spherical_k_integral(τ, elph_matrix, v, w, α, ω; cutoff = [0, Inf])
	integral, _ = quadgk(k -> k^2 * elph_matrix(k; α = α, ω = ω) * exp(-k^2 * D_imag(τ * ω, v, w) / 2), cutoff[1], cutoff[2]) 
    return integral / 2π^2
end

function cartesian_k_integral(τ, elph_matrix, v, w, α, ω; dims = 3, cutoff = [-π, π])
    k_integrand(k) = elph_matrix(k, α, ω, 2 * dims) * exp(-k^2 * D_imag(τ * ω, v, w) / 2)
    integral, _ = quadgk(k -> k_integrand(k), cutoff[1], cutoff[2])
    return (integral / 2π)^dims
end

function frohlich_B(elph_matrix, v, w, α, ω, ; cutoff = [0, Inf])
    integral, _ = quadgk(τ -> exp(-τ * ω) * spherical_k_integral(τ, elph_matrix, v, w, α, ω; cutoff = cutoff), 0, Inf)
    return integral
end

function holstein_B_k(elph_matrix, v, w, α, ω; dims = 3, cutoff = [-π, π])
    integral, _ = quadgk(τ -> exp(-τ * ω) * cartesian_k_integral(τ, elph_matrix, v, w, α, ω; dims = dims, cutoff = cutoff), 0, Inf)
    return integral
end

function frohlich_energy(elph_matrix, v, w, α, ω; cutoff = [0, Inf])
    -(A(v, w, ω) + C(v, w, ω)) - frohlich_B(elph_matrix, v, w, α, ω; cutoff = cutoff)
end

function holstein_energy_k(elph_matrix, v, w, α, ω; dims = 3, cutoff = [-π, π])
    -dims / 3 * (A(v, w, ω) + C(v, w, ω)) - holstein_B_k(elph_matrix, v, w, α, ω; dims = dims, cutoff = cutoff)
end

function frohlich_vw_k(elph_matrix, v, w, α, ω; cutoff = [0, Inf], lower = [0, 0], upper = [Inf, Inf])

    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + eps(Float64), w]

    # Feynman 1955 athermal action 
    f(x) = frohlich_energy(elph_matrix, x[1] + x[2], x[2], α, ω; cutoff = cutoff)

    # Use Optim to optimise v and w to minimise enthalpy.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff=:forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
    )

    # Get v and w values that minimise the free energy.
    Δv, w = Optim.minimizer(solution) # v=Δv+w
	E = Optim.minimum(solution)

	if Optim.converged(solution) == false
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    end

    # Return variational parameters that minimise the free energy.
    return Δv + w, w, E
end

function holstein_vw_k(elph_matrix, v, w, α, ω; dims = 3, cutoff = [-π, π], lower = [0, 0], upper = [Inf, Inf])

    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + eps(Float64), w]

    # Feynman 1955 athermal action 
    f(x) = holstein_energy_k(elph_matrix, x[1] + x[2], x[2], α, ω; dims = dims, cutoff = cutoff)

    # Use Optim to optimise v and w to minimise enthalpy.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff=:forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
    )

    # Get v and w values that minimise the free energy.
    Δv, w = Optim.minimizer(solution) # v=Δv+w
	E = Optim.minimum(solution)

	if Optim.converged(solution) == false
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    end

    # Return variational parameters that minimise the free energy.
    return Δv + w, w, E
end