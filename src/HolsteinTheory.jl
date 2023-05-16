# HolsteinTheory.jl
# A work in progress.

# Compared to Hellwarth's presentation of the Feynman variational approach for the Frohlich
# Hamiltonian, then only thing that changes with the Holstein Hamiltonian is in the
# calculation of the 'B' component (the coupled electron-phonon component).

function holstein_B(v, w, α, ω; dims=3)

    # Trial polaron green's function
	G(τ) = abs(D_imag(τ, v, w) / ω * dims)

    integral = quadgk(τ -> exp(-τ) * (erf(π * √G(τ)) / √G(τ))^dims, 0, Inf)[1]

    return 2 * dims * α / (π)^(dims/2) * integral
end

function holstein_B(v, w, α, ω, β; dims=3)

    # Trial polaron green's function
	G(τ) = D_imag(τ, v, w, ω, β) / ω * dims

    integral = quadgk(τ -> cosh(τ - ω * β / 2) / sinh(ω * β / 2) * (erf(π * √G(τ)) / √G(τ))^dims, 0, β / 2)[1]

    return 2 * dims * α / (π)^(dims/2) * integral
end

function holstein_energy(v, w, α, ωβ...; dims=3)

	kinetic_energy = -2 * dims - dims * (A(v, w, ωβ...) + C(v, w, ωβ...)) / 4

	potential_energy = -holstein_B(v, w, α, ωβ...; dims=dims)

    total_energy = kinetic_energy + potential_energy

	return total_energy, kinetic_energy, potential_energy
end

function holsteinvw(v, w, α, ωβ...; dims=3, upperlimit=Inf)

    # Limits of the optimisation.
    lower = [0.0, 1.0]
    upper = [upperlimit, 20.0]

    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + eps(Float64), w]

    # Feynman 1955 athermal action 
    f(x) = holstein_energy((x[1] + x[2]), x[2], α, ωβ...; dims=dims)[1]

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
	total_energy, kinetic_energy, potential_energy = holstein_energy(Δv + w, w, α, ωβ...; dims=dims)

	if Optim.converged(solution) == false
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    end

    # Return variational parameters that minimise the free energy.
    return Δv + w, w, total_energy, kinetic_energy, potential_energy
end
