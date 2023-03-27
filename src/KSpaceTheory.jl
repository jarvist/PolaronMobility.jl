# KSpaceTheory.jl

function frohlich_matrix(k; α=1, ω=1)
	ω^(3/2) * 2√2 * π * α / k^2
end

function spherical_k_integral(τ, elph_matrix, α, ω, v, w; cutoff = Inf)
	return quadgk(k -> k^2 * elph_matrix(k; α=α, ω=ω) * exp(-k^2 * D_imag(τ, v * ω, w * ω) / 2), 0, cutoff)[1]
end

function frohlich_B(elph_matrix, α, ω, v, w; cutoff = Inf)
    quadgk(τ -> exp(-τ * ω) * spherical_k_integral(τ, elph_matrix, α, ω, v, w; cutoff = cutoff), 0, Inf)[1] / 2π^2
end

function frohlich_energy(elph_matrix, α, ω, v, w; cutoff = Inf)
    -A(v, w, ω) - C(v, w, ω) - frohlich_B(elph_matrix, α, ω, v, w; cutoff = cutoff)
end

function frohlich_vw_k(elph_matrix, α, ω; v = 3, w = 2, cutoff = Inf)
    # Limits of the optimisation.
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + eps(Float64), w]

    # Feynman 1955 athermal action 
    f(x) = frohlich_energy(elph_matrix, α, ω, x[1] + x[2], x[2]; cutoff = cutoff)

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
