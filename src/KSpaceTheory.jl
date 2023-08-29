# KSpaceTheory.jl

"""
	Integrand for a k-space integral in cartesian coordinates.
"""
function cartesian_k_integrand(k, coupling, propagator; rₚ = 1)
	coupling(k) * exp(-k^2 * rₚ^2 * propagator / 2)
end

"""
	Integrand for a k-space integral in spherical coordinates.
"""
function spherical_k_integrand(k, coupling, propagator; rₚ = 1)
	4π * k^2 * coupling(k) * exp(-k^2 * rₚ^2 * propagator / 2)
end

"""
	Integral over first Brillioun Zone in cartesian coordinates.
"""
function cartesian_k_integral(coupling, propagator; rₚ = 1, a = 1, limits = [-π, π])
	integral, _ = quadgk(k -> cartesian_k_integrand(k, coupling, propagator; rₚ = rₚ), limits[1] / a, limits[2] / a)
	integral * a / 2π
end

"""
	Integral over all k-space in spherical coordinates.
"""
function spherical_k_integral(coupling, propagator; rₚ = 1, limits = [0, Inf])
	integral, _ = quadgk(k -> spherical_k_integrand(k, coupling, propagator; rₚ = rₚ), limits[1], limits[2])
	integral / 8π^3
end

"""
	Electron-phonon coupling matrix for the Holstein model.
"""
function holstein_coupling(k, α, ω; dims = 3)
	(dims * α * ω)^(1/dims)
end

"""
	Electron-phonon coupling matrix for the Frohlich model.
"""
function frohlich_coupling(k, α, ω)
	ω^(3/2) * 2√2 * π * α / k^2
end

"""
	Integrand for imaginary time integral for the holstein interaction energy at finite temperature.
"""
function holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; dims = 3)
	coupling(k) = holstein_coupling(k, α, ω; dims = dims)
	propagator = polaron_propagator(τ, v, w, ω, β)
	phonon_propagator(τ, ω, β) * cartesian_k_integral(coupling, propagator; rₚ = 1, a = 1, limits = [-π, π])^dims
end

"""
	Integrand for imaginary time integral for the holstein interaction energy at zero temperature.
"""
function holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω; dims = 3)
	coupling(k) = holstein_coupling(k, α, ω; dims = dims)
	propagator = polaron_propagator(τ, v, w, ω)
	phonon_propagator(τ, ω) * (cartesian_k_integral(coupling, propagator; rₚ = 1, a = 1.0, limits = [-π, π]))^dims
end

"""
	Integrand for imaginary time integral for the frohlich interaction energy at finite temperature.
"""
function frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω, β)
	coupling(k) = frohlich_coupling(k, α, ω)
	propagator = polaron_propagator(τ, v, w, ω, β)
	phonon_propagator(τ, ω, β) * spherical_k_integral(coupling, propagator; rₚ = 1.0, limits = [0, Inf])
end

"""
	Integrand for imaginary time integral for the frohlich interaction energy at zero temperature.
"""
function frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω)
	coupling(k) = frohlich_coupling(k, α, ω)
	propagator = polaron_propagator(τ, v, w, ω)
	phonon_propagator(τ, ω) * spherical_k_integral(coupling, propagator; rₚ = 1.0, limits = [0, Inf])
end

"""
	Electron-phonon interaction energy for the Holstein mode at finite temperature.
"""
function holstein_interaction_energy_k_space(v, w, α, ω, β; dims = 3)
	integral, _ = quadgk(τ -> holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; dims = dims), 0, β/2)
	return integral
end

"""
	Electron-phonon interaction energy for the Holstein mode at zero temperature.
"""
function holstein_interaction_energy_k_space(v, w, α, ω; dims = 3)
	integral, _ = quadgk(τ -> holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω; dims = dims), 0, Inf)
	return integral
end

"""
	Electron-phonon interaction energy for the Frohlich model at finite temperature.
"""
function frohlich_interaction_energy_k_space(v, w, α, ω, β)
	integral, _ = quadgk(τ -> frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω, β), 0, β/2)
	return integral
end

"""
	Electron-phonon interaction energy for the Frohlich mode at zero temperature.
"""
function frohlich_interaction_energy_k_space(v, w, α, ω)
	integral, _ = quadgk(τ -> frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω), 0, Inf)
	return integral
end

"""
	Energy associated with the free electron at finite temperature.
"""
function electron_energy(v, w, ω, β; dims = 3)
	# dims / β / ω * (log(v / w) - 1 / 2 * log(2π * ω * β) - log(sinh(v * ω * β / 2) / sinh(w * ω * β / 2))) + dims / 4 * (v^2 - w^2) / v * (coth(v * ω * β / 2) - 2 / (v * ω * β)) * ω
	(A(v, w, ω, β) + C(v, w, ω, β)) / 3
end

"""
	Energy associated with the free electron at zero temperature.
"""
function electron_energy(v, w, ω; dims = 3)
	-(v - w) / 2 * ω + (1 / (4 * v)) * (v^2 - w^2) * ω
end

"""
	Total free energy for the Holstein model.
"""
function holstein_energy_k_space(v, w, α, ωβ...; dims = 3)
	kinetic_energy = -2 * dims - electron_energy(v, w, ωβ...; dims = dims) 
	potential_energy = -holstein_interaction_energy_k_space(v, w, α, ωβ...; dims = dims)
	return kinetic_energy + potential_energy, kinetic_energy, potential_energy
end

"""
	Total free energy for the Frohlich model.
"""
function frohlich_energy_k_space(v, w, α, ωβ...)
	kinetic_energy = -electron_energy(v, w, ωβ...) 
	potential_energy = -frohlich_interaction_energy_k_space(v, w, α, ωβ...)
	return kinetic_energy + potential_energy, kinetic_energy, potential_energy
end

"""
	The Dynamical Structure Factor for the Holstein model with temperature dependence.
"""
function holstein_structure_factor_k_space(t, v, w, α, ω, β; dims = 3)
	
	coupling_one(k) = holstein_coupling(k, α, ω; dims = dims) * k^2
	coupling_two(k) = holstein_coupling(k, α, ω; dims = dims)

	propagator = polaron_propagator(im * t, v, w, ω, β)
	
	first_integral = cartesian_k_integral(coupling_one, propagator)
	second_integral = cartesian_k_integral(coupling_two, propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω, β)

	prefactor * dims * first_integral * second_integral^(dims - 1)
end

"""
	The Dynamical Structure Factor for the Holstein model at zero temperature.
"""
function holstein_structure_factor_k_space(t, v, w, α, ω; dims = 3)
	
	coupling_one(k) = holstein_coupling(k, α, ω; dims = dims) * k^2
	coupling_two(k) = holstein_coupling(k, α, ω; dims = dims)

	propagator = polaron_propagator(im * t, v, w, ω)
	
	first_integral = cartesian_k_integral(coupling_one, propagator)
	second_integral = cartesian_k_integral(coupling_two, propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω)

	prefactor * dims * first_integral * second_integral^(dims - 1)
end

"""
	The Dynamical Structure Factor for the Frohlich model with temperature dependence.
"""
function frohlich_structure_factor_k_space(t, v, w, α, ω, β)
	
	coupling(k) = frohlich_coupling(k, α, ω) * k^2

	propagator = polaron_propagator(im * t, v, w, ω, β)
	
	integral = spherical_k_integral(coupling, propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω, β)

	prefactor * integral
end

"""
	The Dynamical Structure Factor for the Frohlich model at zero temperature.
"""
function frohlich_structure_factor_k_space(t, v, w, α, ω)
	
	coupling(k) = frohlich_coupling(k, α, ω) * k^2

	propagator = polaron_propagator(im * t, v, w, ω)
	
	integral = spherical_k_integral(coupling, propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω)

	prefactor * integral
end

"""
	Memory function for the Holstein model with temperature dependence.
"""
function holstein_memory_function_k_space(Ω, v, w, α, ω, β; dims = 3)
	 structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, ω, β; dims = dims)
	 return general_memory_function(Ω, structure_factor)
end

"""
	Memory function for the Holstein model at zero temperature.
"""
function holstein_memory_function_k_space(Ω, v, w, α, ω; dims = 3)
	 structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, ω; dims = dims)
	 return general_memory_function(Ω, structure_factor, limits = [0, 1e3])
end

"""
	DC Mobility for the Holstein model.
"""
function holstein_mobility_k_space(v, w, α, ω, β; dims = 3)
	 structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, ω, β; dims = dims)
	 1 / imag(general_memory_function(structure_factor))
end

"""
	Memory function for the Frohlich model with temperature dependence.
"""
function frohlich_memory_function_k_space(Ω, v, w, α, ω, β)
	 structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ω, β)
	 return general_memory_function(Ω, structure_factor; limits = [0, Inf])
end

"""
	Memory function for the Frohlich model at zero temperature.
"""
function frohlich_memory_function_k_space(Ω, v, w, α, ω)
	 structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ω)
	 return general_memory_function(Ω, structure_factor; limits = [0, 1e3])
end

"""
	DC Mobility for the Frohlich model.
"""
function frohlich_mobility_k_space(v, w, α, ω, β)
	 structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ω, β)
	 1 / imag(general_memory_function(structure_factor; limits = [0, Inf]))
end