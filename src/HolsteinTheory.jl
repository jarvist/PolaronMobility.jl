# HolsteinTheory.jl

"""
	Imaginary time polaron Green's function with temperature dependence.
"""
function polaron_propagator(τ, v, w, ω, β)
	(v^2 - w^2) / v^3 / ω * (1 - exp(-v * τ * ω)) * (1 - exp(-v * ω * (β - τ))) / (1 - exp(-v * β * ω)) + w^2 / v^2 * τ * ω * (1 - τ / β)
end

"""
	Imaginary time polaron Green's function at zero temperature.
"""
function polaron_propagator(τ, v, w, ω)
	w^2 * τ * ω / v^2 + (v^2 - w^2) / v^3 * (1 - exp(-v * τ * ω)) / ω
end

"""
	Imaginary time phonon Green's function at zero temperature.
"""
function phonon_propagator(τ, ω)
	exp(-ω * τ)
end

"""
	Imaginary time phonon Green's function with temperature dependence.
"""
function phonon_propagator(τ, ω, β)
	n = 1 / (exp(β * ω) - 1)
	n * exp(ω * τ) + (1 + n) * exp(-ω * τ)
end

"""
	Integrand for imaginary time integral for the holstein interaction energy at finite temperature. Here the k-space integral is evaluated analytically.
"""
function holstein_interaction_energy_integrand(τ, v, w, α, ω, β; dims = 3)
	coupling = holstein_coupling(1, α, ω; dims = dims)
	propagator = polaron_propagator(τ, v, w, ω, β)
	phonon_propagator(τ, ω, β) * (coupling * erf(π * sqrt(propagator * dims)) / sqrt(2π * propagator))^dims
end

"""
	Integrand for imaginary time integral for the holstein interaction energy at zero temperature. Here the k-space integral is evaluated analytically.
"""
function holstein_interaction_energy_integrand(τ, v, w, α, ω; dims = 3)
	coupling = holstein_coupling(1, α, ω; dims = dims)
	propagator = polaron_propagator(τ, v, w, ω)
	phonon_propagator(τ, ω) * (coupling * erf(π * sqrt(propagator * dims / ω)) / sqrt(2π * propagator))^dims
end

"""
	Electron-phonon interaction energy for the Holstein mode at finite temperature. Here the k-space integral is evaluated analytically.
"""
function holstein_interaction_energy(v, w, α, ω, β; dims = 3)
	integral, _ = quadgk(τ -> holstein_interaction_energy_integrand(τ, v, w, α, ω, β; dims = dims), 0, β/2)
	return integral
end

"""
	Electron-phonon interaction energy for the Holstein mode at finite temperature. Here the k-space integral is evaluated analytically.
"""
function holstein_interaction_energy(v, w, α, ω; dims = 3)
	integral, _ = quadgk(τ -> holstein_interaction_energy_integrand(τ, v, w, α, ω; dims = dims), 0, Inf)
	return integral
end

"""
	Total free energy for the Holstein model. Here the k-space integral is evaluated analytically.
"""
function holstein_energy(v, w, α, ωβ...; dims = 3)
	kinetic_energy = -2 * dims - electron_energy(v, w, ωβ...; dims = dims) / dims
	potential_energy = -holstein_interaction_energy(v, w, α, ωβ...; dims = dims) / 2.0^dims
	return (kinetic_energy + potential_energy), kinetic_energy, potential_energy
end

"""
	Optimization code to find the variational apramaters v and w that give the lowest upper-bound to the free energy.
"""
function vw_variation(energy, v, w; lower = [0, 0], upper = [Inf, Inf])

    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + eps(Float64), w]

    # Feynman 1955 athermal action 
    f(x) = energy(x[1] + x[2], x[2])[1]

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
	E, K, P = energy(Δv + w, w)

	if Optim.converged(solution) == false
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    end

    # Return variational parameters that minimise the free energy.
    return Δv + w, w, E, K, P
end

"""
	General memory function for a given structure factor with frequency dependence.
"""
function general_memory_function(Ω, structure_factor; limits = [0, Inf])
	integral, _ = quadgk(t -> (1 - exp(im * Ω * t)) / Ω * imag(structure_factor(t)), limits[1], limits[2])
	return integral
end

"""
	General memory function for a given structure factor in the DC limit.
"""
function general_memory_function(structure_factor; limits = [0, Inf])
	integral, _ = quadgk(t -> -im * t * imag(structure_factor(t)), limits[1], limits[2])
	return integral
end

"""
	The Dynamical Structure Factor for the Holstein model with temperature dependence. Here the k-space integral is evaulated analytically.
"""
function holstein_structure_factor(t, v, w, α, ω, β; dims = 3)
	
	coupling = holstein_coupling(1, α, ω; dims = dims)^dims

	propagator = polaron_propagator(im * t, v, w, ω, β)
	
	first_integral = erf(π * sqrt(propagator / 2)) / sqrt(2π) / propagator^(3/2) - exp(-π^2 * propagator / 2) / propagator
	second_integral = erf(π * sqrt(propagator / 2)) / sqrt(2π * propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω, β)

	prefactor * dims * coupling * first_integral * second_integral^(dims - 1)
end

"""
	The Dynamical Structure Factor for the Holstein model at zero dependence. Here the k-space integral is evaulated analytically.
"""
function holstein_structure_factor(t, v, w, α, ω; dims = 3)
	
	coupling = holstein_coupling(1, α, ω; dims = dims)^dims

	propagator = polaron_propagator(im * t, v, w, ω)
	
	first_integral = erf(π * sqrt(propagator / 2)) / sqrt(2π) / propagator^(3/2) - exp(-π^2 * propagator / 2) / propagator
	second_integral = erf(π * sqrt(propagator / 2)) / sqrt(2π * propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω)

	prefactor * dims * coupling * first_integral * second_integral^(dims - 1)
end

"""
	Memory function for the Holstein model with temperature dependence. Here the k-space integral is evaulated analytically.
"""
function holstein_memory_function(Ω, v, w, α, ω, β; dims = 3)
	 structure_factor(t) = holstein_structure_factor(t, v, w, α, ω, β; dims = dims)
	 return general_memory_function(Ω, structure_factor)
end

"""
	Memory function for the Holstein model at zero temperature. Here the k-space integral is evaulated analytically.
"""
function holstein_memory_function(Ω, v, w, α, ω; dims = 3)
	 structure_factor(t) = holstein_structure_factor(t, v, w, α, ω; dims = dims)
	 return general_memory_function(Ω, structure_factor, limits = [0, 1e4])
end

"""
	DC Mobility for the Holstein model. Here the k-space integral is evaulated analytically.
"""
function holstein_mobility(v, w, α, ω, β; dims = 3)
	 structure_factor(t) = holstein_structure_factor(t, v, w, α, ω, β; dims = dims)
	 1 / imag(general_memory_function(structure_factor))
end