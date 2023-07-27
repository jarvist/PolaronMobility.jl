### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ cd6cae06-27d1-11ee-1f4b-07833f940a86
using QuadGK, Optim, Plots, SpecialFunctions

# ╔═╡ 9ed21af8-f7c9-4cba-b172-e3b376665add
"""
	Unzips an array of tuples.
"""
unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

# ╔═╡ 1e3e504b-5389-4ef1-9d63-5ac4e10364ef
md"# Green's Functions"

# ╔═╡ e005016f-36e2-47a6-9172-8e32ac8bbff6
"""
	Imaginary time polaron Green's function with temperature dependence.
"""
function polaron_propagator(τ, v, w, ω, β)
	(v^2 - w^2) / v^3 / ω * (1 - exp(-v * τ * ω)) * (1 - exp(-v * ω * (β - τ))) / (1 - exp(-v * β * ω)) + w^2 / v^2 * τ * ω * (1 - τ / β)
end

# ╔═╡ 6c51fe56-6a20-48c2-b6e4-4075350b807e
"""
	Imaginary time polaron Green's function at zero temperature.
"""
function polaron_propagator(τ, v, w, ω)
	w^2 * τ * ω / v^2 + (v^2 - w^2) / v^3 * (1 - exp(-v * τ * ω)) / ω
end

# ╔═╡ 14b28d72-9bb9-47fb-a374-c0ec287d8440
"""
	Imaginary time phonon Green's function at zero temperature.
"""
function phonon_propagator(τ, ω)
	exp(-ω * τ)
end

# ╔═╡ 7bf51005-7e61-4f8d-b5ea-db37aae30988
"""
	Imaginary time phonon Green's function with temperature dependence.
"""
function phonon_propagator(τ, ω, β)
	n = 1 / (exp(β * ω) - 1)
	n * exp(ω * τ) + (1 + n) * exp(-ω * τ)
end

# ╔═╡ 3a7d73db-2d2b-43b8-ae39-adc01d9f5b30
md"# K-Space Integrals"

# ╔═╡ 2c04db4c-6366-48f5-96fb-78a579c727b1
"""
	Integrand for a k-space integral in cartesian coordinates.
"""
function cartesian_k_integrand(k, coupling, propagator; rₚ = 1)
	coupling(k) * exp(-k^2 * rₚ^2 * propagator / 2)
end

# ╔═╡ 50639eb7-3047-4a95-a668-f67a66775e4d
"""
	Integrand for a k-space integral in spherical coordinates.
"""
function spherical_k_integrand(k, coupling, propagator; rₚ = 1)
	4π * k^2 * coupling(k) * exp(-k^2 * rₚ^2 * propagator / 2)
end

# ╔═╡ 256f3d5d-f474-4151-a3db-bec8099615fa
"""
	Integral over first Brillioun Zone in cartesian coordinates.
"""
function cartesian_k_integral(coupling, propagator; rₚ = 1, a = 1, limits = [-π, π])
	integral, _ = quadgk(k -> cartesian_k_integrand(k, coupling, propagator; rₚ = rₚ), limits[1] / a, limits[2] / a)
	integral * a / 2π
end

# ╔═╡ 237bf7db-2e0e-4c82-8ec9-4ac65ddfdd5c
"""
	Integral over all k-space in spherical coordinates.
"""
function spherical_k_integral(coupling, propagator; rₚ = 1, limits = [0, Inf])
	integral, _ = quadgk(k -> spherical_k_integrand(k, coupling, propagator; rₚ = rₚ), limits[1], limits[2])
	integral / 8π^3
end

# ╔═╡ 65093555-b406-45c5-854a-fdea6dad4c9c
md"# Polaron Coupling Matrices"

# ╔═╡ d295e844-ff9f-497e-993a-101f5c2ea0a8
"""
	Electron-phonon coupling matrix for the Holstein model.
"""
function holstein_coupling(k, α, ω; dims = 3)
	(2 * dims * α * ω)^(1/dims)
end

# ╔═╡ e5bbf2cc-2650-4666-9650-5a91ecea304d
"""
	Electron-phonon coupling matrix for the Frohlich model.
"""
function frohlich_coupling(k, α, ω)
	ω^(3/2) * 2√2 * π * α / k^2
end

# ╔═╡ 577587db-25b0-463b-ab44-38522bcd0039
md"# Polaron Free Energies (K Space)"

# ╔═╡ f73c700b-2f00-4deb-85cc-abdaaa4118e5
"""
	Integrand for imaginary time integral for the holstein interaction energy at finite temperature.
"""
function holstein_interaction_energy_integrand(τ, v, w, α, ω, β; dims = 3)
	coupling(k) = holstein_coupling(k, α, ω; dims = dims)
	propagator = polaron_propagator(τ, v, w, ω, β)
	phonon_propagator(τ, ω, β) * cartesian_k_integral(coupling, propagator; rₚ = 1, a = 1, limits = [-π, π])^dims
end

# ╔═╡ 4143c668-699d-46ea-b825-bc4996fb26a4
"""
	Integrand for imaginary time integral for the holstein interaction energy at zero temperature.
"""
function holstein_interaction_energy_integrand(τ, v, w, α, ω; dims = 3)
	coupling(k) = holstein_coupling(k, α, ω; dims = dims)
	propagator = polaron_propagator(τ, v, w, ω)
	phonon_propagator(τ, ω) * (cartesian_k_integral(coupling, propagator; rₚ = 1, a = 1.0, limits = [-π, π]))^dims
end

# ╔═╡ 57ed148b-fef1-40f4-a113-a39a1aebb010
"""
	Integrand for imaginary time integral for the frohlich interaction energy at finite temperature.
"""
function frohlich_interaction_energy_integrand(τ, v, w, α, ω, β)
	coupling(k) = frohlich_coupling(k, α, ω)
	propagator = polaron_propagator(τ, v, w, ω, β)
	phonon_propagator(τ, ω, β) * spherical_k_integral(coupling, propagator; rₚ = 1.0, limits = [0, Inf])
end

# ╔═╡ f223d349-a0d7-4fc6-a6ed-f8b3c45c8dad
"""
	Integrand for imaginary time integral for the frohlich interaction energy at zero temperature.
"""
function frohlich_interaction_energy_integrand(τ, v, w, α, ω)
	coupling(k) = frohlich_coupling(k, α, ω)
	propagator = polaron_propagator(τ, v, w, ω)
	phonon_propagator(τ, ω) * spherical_k_integral(coupling, propagator; rₚ = 1.0, limits = [0, Inf])
end

# ╔═╡ ac60f707-8a87-47b8-9856-80194fb5161c
"""
	Electron-phonon interaction energy for the Holstein mode at finite temperature.
"""
function holstein_interaction_energy(v, w, α, ω, β; dims = 3)
	integral, _ = quadgk(τ -> holstein_interaction_energy_integrand(τ, v, w, α, ω, β; dims = dims), 0, β/2)
	return integral
end

# ╔═╡ 949b3f46-146c-4957-a14c-4c742d6d10ba
"""
	Electron-phonon interaction energy for the Holstein mode at zero temperature.
"""
function holstein_interaction_energy(v, w, α, ω; dims = 3)
	integral, _ = quadgk(τ -> holstein_interaction_energy_integrand(τ, v, w, α, ω; dims = dims), 0, Inf)
	return integral
end

# ╔═╡ 5fa89b72-fab5-45b3-831c-590021953079
"""
	Electron-phonon interaction energy for the Frohlich model at finite temperature.
"""
function frohlich_interaction_energy(v, w, α, ω, β)
	integral, _ = quadgk(τ -> frohlich_interaction_energy_integrand(τ, v, w, α, ω, β), 0, β/2)
	return integral
end

# ╔═╡ 303b4158-8a39-4733-84d2-ff1dbe6bb556
"""
	Electron-phonon interaction energy for the Frohlich mode at zero temperature.
"""
function frohlich_interaction_energy(v, w, α, ω)
	integral, _ = quadgk(τ -> frohlich_interaction_energy_integrand(τ, v, w, α, ω), 0, Inf)
	return integral
end

# ╔═╡ 0fef496e-6555-4f07-9586-c5d18314e386
"""
	Energy associated with the free electron at finite temperature.
"""
function electron_energy(v, w, ω, β; dims = 3)
	dims / β / ω * (log(v / w) - 1 / 2 * log(2π * ω * β) - log(sinh(v * ω * β / 2) / sinh(w * ω * β / 2))) + dims / 4 * (v^2 - w^2) / v * (coth(v * ω * β / 2) - 2 / (v * ω * β)) * ω
end

# ╔═╡ d94ab3bf-ef6e-4c1b-ba54-c025916e42e1
"""
	Energy associated with the free electron at zero temperature.
"""
function electron_energy(v, w, ω; dims = 3)
	-dims * (v - w) / 2 * ω + (dims / (4 * v)) * (v^2 - w^2) * ω
end

# ╔═╡ 62faa65d-d426-4d53-9816-ca1de8b24874
"""
	Total free energy for the Holstein model.
"""
function holstein_energy(v, w, α, ωβ...; dims = 3)
	kinetic_energy = -2 * dims - electron_energy(v, w, ωβ...; dims = dims) 
	potential_energy = -holstein_interaction_energy(v, w, α, ωβ...; dims = dims)
	return kinetic_energy + potential_energy, kinetic_energy, potential_energy
end

# ╔═╡ 293773c1-f0f1-4f74-8434-3b82f6cf8883
"""
	Total free energy for the Frohlich model.
"""
function frohlich_energy(v, w, α, ωβ...)
	kinetic_energy = -electron_energy(v, w, ωβ...) 
	potential_energy = -frohlich_interaction_energy(v, w, α, ωβ...)
	return kinetic_energy + potential_energy, kinetic_energy, potential_energy
end

# ╔═╡ b175dd29-6925-4352-9914-91d8e2f5fa9c
md"# Polaron Free Energy (Analytic)"

# ╔═╡ 0e99ecd9-f2ff-4536-8d8d-a04e42f6c98f
"""
	Integrand for imaginary time integral for the holstein interaction energy at finite temperature. Here the k-space integral is evaluated analytically.
"""
function holstein_interaction_energy_integrand_analytic(τ, v, w, α, ω, β; dims = 3)
	coupling = holstein_coupling(1, α, ω; dims = dims)
	propagator = polaron_propagator(τ, v, w, ω, β)
	phonon_propagator(τ, ω, β) * (coupling * erf(π * sqrt(propagator / 2)) / sqrt(2π * propagator))^dims
end

# ╔═╡ c1dd1475-4ad6-478f-9c18-7c20df32be81
"""
	Integrand for imaginary time integral for the holstein interaction energy at zero temperature. Here the k-space integral is evaluated analytically.
"""
function holstein_interaction_energy_integrand_analytic(τ, v, w, α, ω; dims = 3)
	coupling = holstein_coupling(1, α, ω; dims = dims)
	propagator = polaron_propagator(τ, v, w, ω)
	phonon_propagator(τ, ω) * (coupling * erf(π * sqrt(propagator / 2)) / sqrt(2π * propagator))^dims
end

# ╔═╡ fe21e743-2cf4-4e62-857a-aebbe539287e
"""
	Integrand for imaginary time integral for the frohlich interaction energy at finite temperature. Here the k-space integral is evaluated analytically.
"""
function frohlich_interaction_energy_integrand_analytic(τ, v, w, α, ω, β)
	coupling = frohlich_coupling(1, α, ω)
	propagator = polaron_propagator(τ, v, w, ω, β)
	phonon_propagator(τ, ω, β) * coupling * 4π / (2π)^3 * sqrt(π / 2 / propagator)
end

# ╔═╡ fe0642f4-09b9-4209-baf8-959d46d593ad
"""
	Integrand for imaginary time integral for the frohlich interaction energy at finite temperature. Here the k-space integral is evaluated analytically.
"""
function frohlich_interaction_energy_integrand_analytic(τ, v, w, α, ω)
	coupling = frohlich_coupling(1, α, ω)
	propagator = polaron_propagator(τ, v, w, ω)
	phonon_propagator(τ, ω) * coupling * 4π / (2π)^3 * sqrt(π / 2 / propagator)
end

# ╔═╡ b2622ccb-9c3d-40a3-9a4d-f248e90632dc
"""
	Electron-phonon interaction energy for the Holstein mode at finite temperature. Here the k-space integral is evaluated analytically.
"""
function holstein_interaction_energy_analytic(v, w, α, ω, β; dims = 3)
	integral, _ = quadgk(τ -> holstein_interaction_energy_integrand_analytic(τ, v, w, α, ω, β; dims = dims), 0, β/2)
	return integral
end

# ╔═╡ 7b893719-994e-4a8b-83f8-67b61f158189
"""
	Electron-phonon interaction energy for the Holstein mode at finite temperature. Here the k-space integral is evaluated analytically.
"""
function holstein_interaction_energy_analytic(v, w, α, ω; dims = 3)
	integral, _ = quadgk(τ -> holstein_interaction_energy_integrand_analytic(τ, v, w, α, ω; dims = dims), 0, Inf)
	return integral
end

# ╔═╡ 4bdfeceb-01a1-43f4-8634-54ed51f5ad82
"""
	Electron-phonon interaction energy for the Frohlich mode at finite temperature. Here the k-space integral is evaluated analytically.
"""
function frohlich_interaction_energy_analytic(v, w, α, ω, β)
	integral, _ = quadgk(τ -> frohlich_interaction_energy_integrand_analytic(τ, v, w, α, ω, β), 0, β/2)
	return integral
end

# ╔═╡ 11a19771-8cf0-4f4a-98d4-5a1060d500de
"""
	Electron-phonon interaction energy for the Frohlich mode at zero temperature. Here the k-space integral is evaluated analytically.
"""
function frohlich_interaction_energy_analytic(v, w, α, ω)
	integral, _ = quadgk(τ -> frohlich_interaction_energy_integrand_analytic(τ, v, w, α, ω), 0, Inf)
	return integral
end

# ╔═╡ 0c5fc742-029a-4f47-8897-c9c1483e90bd
"""
	Total free energy for the Holstein model. Here the k-space integral is evaluated analytically.
"""
function holstein_energy_analytic(v, w, α, ωβ...; dims = 3)
	kinetic_energy = -2 * dims - electron_energy(v, w, ωβ...; dims = dims) 
	potential_energy = -holstein_interaction_energy_analytic(v, w, α, ωβ...; dims = dims)
	return kinetic_energy + potential_energy, kinetic_energy, potential_energy
end

# ╔═╡ 9e8e3e8f-095a-4bd0-a44e-54e9a0be3ebe
"""
	Total free energy for the Frohlich model. Here the k-space integral is evaluated analytically.
"""
function frohlich_energy_analytic(v, w, α, ωβ...)
	kinetic_energy = -electron_energy(v, w, ωβ...) 
	potential_energy = -frohlich_interaction_energy_analytic(v, w, α, ωβ...)
	return kinetic_energy + potential_energy, kinetic_energy, potential_energy
end

# ╔═╡ fef541dc-568a-477f-9343-bee877a46b4a
md"# Optimization"

# ╔═╡ 3a1ade3b-c89e-4137-9a7d-2faa572bea7c
"""
	Optimization code to find the variational apramaters v and w that give the lowest upper-bound to the free energy.
"""
function optimize(energy, v, w; lower = [0, 0], upper = [Inf, Inf])

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

# ╔═╡ 9a43b38c-8769-426c-99a9-98cdc6fbebac
md"# Dynamical functions: Memory function, conductivity and mobility (K Space)"

# ╔═╡ f93fc7a5-f59d-4158-aa5b-eec93598cdf9
"""
	General memory function for a given structure factor with frequency dependence.
"""
function memory_function(Ω, structure_factor; limits = [0, Inf])
	integral, _ = quadgk(t -> (1 - exp(im * Ω * t)) / Ω * imag(structure_factor(t)), limits[1], limits[2], rtol=1e-4)
	return integral
end

# ╔═╡ f122d4aa-3d4b-43e7-9370-4fe64efaa961
"""
	General memory function for a given structure factor in the DC limit.
"""
function memory_function(structure_factor; limits = [0, Inf])
	integral, _ = quadgk(t -> -im * t * imag(structure_factor(t)), limits[1], limits[2], rtol=1e-4)
	return integral
end

# ╔═╡ fb7b8e2a-6d17-42ae-b264-5c3434643ec4
"""
	The Dynamical Structure Factor for the Holstein model with temperature dependence.
"""
function holstein_structure_factor(t, v, w, α, ω, β; dims = 3)
	
	coupling_one(k) = holstein_coupling(k, α, ω; dims = dims) * k^2
	coupling_two(k) = holstein_coupling(k, α, ω; dims = dims)

	propagator = polaron_propagator(im * t, v, w, ω, β)
	
	first_integral = cartesian_k_integral(coupling_one, propagator)
	second_integral = cartesian_k_integral(coupling_two, propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω, β)

	prefactor * dims * first_integral * second_integral^(dims - 1)
end

# ╔═╡ 5d2f12ce-5e37-4f31-91e7-e4efbd4442ce
"""
	The Dynamical Structure Factor for the Holstein model at zero temperature.
"""
function holstein_structure_factor(t, v, w, α, ω; dims = 3)
	
	coupling_one(k) = holstein_coupling(k, α, ω; dims = dims) * k^2
	coupling_two(k) = holstein_coupling(k, α, ω; dims = dims)

	propagator = polaron_propagator(im * t, v, w, ω)
	
	first_integral = cartesian_k_integral(coupling_one, propagator)
	second_integral = cartesian_k_integral(coupling_two, propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω)

	prefactor * dims * first_integral * second_integral^(dims - 1)
end

# ╔═╡ 72196d1e-a5cf-4ffd-9d6a-63ab54a83d54
"""
	The Dynamical Structure Factor for the Frohlich model with temperature dependence.
"""
function frohlich_structure_factor(t, v, w, α, ω, β)
	
	coupling(k) = frohlich_coupling(k, α, ω) * k^2

	propagator = polaron_propagator(im * t, v, w, ω, β)
	
	integral = spherical_k_integral(coupling, propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω, β)

	prefactor * integral
end

# ╔═╡ f3fe677b-5edd-4301-b35c-e707c2106649
"""
	The Dynamical Structure Factor for the Frohlich model at zero temperature.
"""
function frohlich_structure_factor(t, v, w, α, ω)
	
	coupling(k) = frohlich_coupling(k, α, ω) * k^2

	propagator = polaron_propagator(im * t, v, w, ω)
	
	integral = spherical_k_integral(coupling, propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω)

	prefactor * integral
end

# ╔═╡ 9684cc0e-aebf-4e51-86d6-7591266bc356
"""
	Memory function for the Holstein model with temperature dependence.
"""
function holstein_memory_function(Ω, v, w, α, ω, β; dims = 3)
	 structure_factor(t) = holstein_structure_factor(t, v, w, α, ω, β; dims = dims)
	 return memory_function(Ω, structure_factor)
end

# ╔═╡ 302b7385-a924-4849-8882-77d6a473e504
"""
	Memory function for the Holstein model at zero temperature.
"""
function holstein_memory_function(Ω, v, w, α, ω; dims = 3)
	 structure_factor(t) = holstein_structure_factor(t, v, w, α, ω; dims = dims)
	 return memory_function(Ω, structure_factor, limits = [0, 1e3])
end

# ╔═╡ f2a175aa-6854-4b0b-b811-2dd0783aca60
"""
	DC Mobility for the Holstein model.
"""
function holstein_mobility(v, w, α, ω, β; dims = 3)
	 structure_factor(t) = holstein_structure_factor(t, v, w, α, ω, β; dims = dims)
	 1 / imag(memory_function(structure_factor))
end

# ╔═╡ 9aa269c6-1df3-4553-ba61-c0b0fd2608c9
"""
	Memory function for the Frohlich model with temperature dependence.
"""
function frohlich_memory_function(Ω, v, w, α, ω, β)
	 structure_factor(t) = frohlich_structure_factor(t, v, w, α, ω, β)
	 return memory_function(Ω, structure_factor; limits = [0, Inf])
end

# ╔═╡ 24fec156-91bd-4761-9933-e7c5ade3693b
"""
	Memory function for the Frohlich model at zero temperature.
"""
function frohlich_memory_function(Ω, v, w, α, ω)
	 structure_factor(t) = frohlich_structure_factor(t, v, w, α, ω)
	 return memory_function(Ω, structure_factor; limits = [0, 1e3])
end

# ╔═╡ ae71295e-6fab-4360-b6ec-e2d4db10f9ae
"""
	DC Mobility for the Frohlich model.
"""
function frohlich_mobility(v, w, α, ω, β)
	 structure_factor(t) = frohlich_structure_factor(t, v, w, α, ω, β)
	 1 / imag(memory_function(structure_factor); limits = [0, Inf])
end

# ╔═╡ 7560cf27-f085-4f74-99f0-d33bfa542a9f
md"# Dynamical functions: Memory function, conductivity and mobility (Analytic)"

# ╔═╡ e8cf5096-e167-4a66-9fbd-ec1181602e77
"""
	The Dynamical Structure Factor for the Holstein model with temperature dependence. Here the k-space integral is evaulated analytically.
"""
function holstein_structure_factor_analytic(t, v, w, α, ω, β; dims = 3)
	
	coupling = holstein_coupling(1, α, ω; dims = dims)^dims

	propagator = polaron_propagator(im * t, v, w, ω, β)
	
	first_integral = erf(π * sqrt(propagator / 2)) / sqrt(2π) / propagator^(3/2) - exp(-π^2 * propagator / 2) / propagator
	second_integral = erf(π * sqrt(propagator / 2)) / sqrt(2π * propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω, β)

	prefactor * dims * coupling * first_integral * second_integral^(dims - 1)
end

# ╔═╡ ea7f5541-84b8-453a-aa0f-84f7df119b3f
"""
	The Dynamical Structure Factor for the Holstein model at zero dependence. Here the k-space integral is evaulated analytically.
"""
function holstein_structure_factor_analytic(t, v, w, α, ω; dims = 3)
	
	coupling = holstein_coupling(1, α, ω; dims = dims)^dims

	propagator = polaron_propagator(im * t, v, w, ω)
	
	first_integral = erf(π * sqrt(propagator / 2)) / sqrt(2π) / propagator^(3/2) - exp(-π^2 * propagator / 2) / propagator
	second_integral = erf(π * sqrt(propagator / 2)) / sqrt(2π * propagator)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω)

	prefactor * dims * coupling * first_integral * second_integral^(dims - 1)
end

# ╔═╡ 03b981b7-d257-4c89-a16d-613815486ac1
"""
	The Dynamical Structure Factor for the Frohlich model with temperature dependence. Here the k-space integral is evaulated analytically.
"""
function frohlich_structure_factor_analytic(t, v, w, α, ω, β)
	
	coupling = frohlich_coupling(1, α, ω)

	propagator = polaron_propagator(im * t, v, w, ω, β)
	
	integral = 4π / (2π)^3 * sqrt(π / 2) / propagator^(3/2)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω, β)

	prefactor * coupling * integral
end

# ╔═╡ 20afdee0-d3be-48fc-8f61-705680397432
"""
	The Dynamical Structure Factor for the Frohlich model at zero temperature. Here the k-space integral is evaulated analytically.
"""
function frohlich_structure_factor_analytic(t, v, w, α, ω)
	
	coupling = frohlich_coupling(1, α, ω)

	propagator = polaron_propagator(im * t, v, w, ω)
	
	integral = 4π / (2π)^3 * sqrt(π / 2) / propagator^(3/2)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω)

	prefactor * coupling * integral
end

# ╔═╡ 95ae687b-1cd3-43f5-bf28-d8a4f5cca022
"""
	Memory function for the Holstein model with temperature dependence. Here the k-space integral is evaulated analytically.
"""
function holstein_memory_function_analytic(Ω, v, w, α, ω, β; dims = 3)
	 structure_factor(t) = holstein_structure_factor_analytic(t, v, w, α, ω, β; dims = dims)
	 return memory_function(Ω, structure_factor)
end

# ╔═╡ 92964c87-3fa3-40de-8f08-4b773b229a6a
"""
	Memory function for the Holstein model at zero temperature. Here the k-space integral is evaulated analytically.
"""
function holstein_memory_function_analytic(Ω, v, w, α, ω; dims = 3)
	 structure_factor(t) = holstein_structure_factor_analytic(t, v, w, α, ω; dims = dims)
	 return memory_function(Ω, structure_factor, limits = [0, 1e3])
end

# ╔═╡ 7e00a70f-5a52-42f5-9056-efc47b4ab839
"""
	Memory function for the Frohlich model with temperature dependence. Here the k-space integral is evaulated analytically.
"""
function frohlich_memory_function_analytic(Ω, v, w, α, ω, β)
	 structure_factor(t) = frohlich_structure_factor_analytic(t, v, w, α, ω, β)
	 return memory_function(Ω, structure_factor; limits = [0, Inf])
end

# ╔═╡ b3818117-e954-4081-ad92-ce066d0fc1b6
"""
	Memory function for the Frohlich model at zero temperature. Here the k-space integral is evaulated analytically.
"""
function frohlich_memory_function_analytic(Ω, v, w, α, ω)
	 structure_factor(t) = frohlich_structure_factor_analytic(t, v, w, α, ω)
	 return memory_function(Ω, structure_factor; limits = [0, 1e3])
end

# ╔═╡ f2ede4d0-d51e-4e55-af43-cb22fffe30f4
"""
	DC Mobility for the Holstein model. Here the k-space integral is evaulated analytically.
"""
function holstein_mobility_analytic(v, w, α, ω, β; dims = 3)
	 structure_factor(t) = holstein_structure_factor_analytic(t, v, w, α, ω, β; dims = dims)
	 1 / imag(memory_function(structure_factor))
end

# ╔═╡ bb9fe3b2-c1d4-4806-8b14-d234e91ee8a9
"""
	DC Mobility for the Frohlich model. Here the k-space integral is evaulated analytically.
"""
function frohlich_mobility_analytic(v, w, α, ω, β)
	 structure_factor(t) = frohlich_structure_factor_analytic(t, v, w, α, ω, β)
	 1 / imag(memory_function(structure_factor; limits = [0, Inf]))
end

# ╔═╡ ed021864-b2de-4f4b-8aac-175896e7bd39
md"# Calculations"

# ╔═╡ b10da998-6267-49b1-88ce-09ee7774eebc
md"## Dependence on α (athermal)"

# ╔═╡ 9b079f65-9dc9-4578-bba0-c55dbe9563ef
begin
	αrange = 0.01:0.01:4 # Electron-phonon coupling (dimensionless)
	ω = 1 # Adiabaticity (units of J)
	d = 3 # Dimensionality of lattice for Holstein
	holstein_energies = []
	frohlich_energies = []
	holstein_energies_analytic = []
	frohlich_energies_analytic = []
	for α in αrange
		best_holstein_energy(v, w) = holstein_energy(v, w, α, ω; dims = d)
		best_frohlich_energy(v, w) = frohlich_energy(v, w, α .* 3, ω)
		best_holstein_energy_analytic(v, w) = holstein_energy_analytic(v, w, α, ω; dims = d)
		best_frohlich_energy_analytic(v, w) = frohlich_energy_analytic(v, w, α .* 3, ω)
		push!(holstein_energies, best_holstein_energy)
		push!(frohlich_energies, best_frohlich_energy)
		push!(holstein_energies_analytic, best_holstein_energy_analytic)
		push!(frohlich_energies_analytic, best_frohlich_energy_analytic)
	end
end

# ╔═╡ 1641996e-41f2-4e31-bb1d-bc7a1763274c
# Optimise the variational paramaters and energy for Holstein (analytic)
vhaα, whaα, Ehaα, Khaα, Phaα = unzip(optimize.(holstein_energies_analytic, 5, 2, upper=[100, 100]))

# ╔═╡ 0ae0ef1d-0c63-43c6-8e67-55b44dba6b67
# Optimise the variational paramaters and energy for Holstein (k-space)
vhα, whα, Ehα, Khα, Phα = unzip(optimize.(holstein_energies, 5, 2, upper=[100, 100]))

# ╔═╡ 1b2e2933-623f-4560-9527-e24fe55ca485
plot(αrange, [vhα whα vhaα whaα], title="GS Holstein Analytic and K-Space v,w", linestyle=:auto, linewidth = 2, xlabel = "α", ylabel = "v, w", labels = ["v k-space" "w k-space" "v analytic" "w analytic"], yticks = 1:9)

# ╔═╡ 37387a74-af02-4485-9edf-3a950fcb822a
plot(αrange, [Ehα Khα Phα Ehaα Khaα Phaα], title="GS Holstein Analytic and K-Space Energies", linestyle=:auto, linewidth = 2, xlabel = "α", ylabel = "Total Energy, Kinetic and Potential (J)", labels = ["E k-space" "K k-space" "P k-space" "E analytic" "K analytic" "P analytic"], yticks = -16:2:0)

# ╔═╡ ea012194-878d-433b-bf97-d2e8e90e7125
# Optimise the variational paramaters and energy for Frohlich (analytic)
vfaα, wfaα, Efaα, Kfaα, Pfaα = unzip(optimize.(frohlich_energies_analytic, 5, 2, upper=[100, 100]))

# ╔═╡ 2be6279b-f73e-4253-8f30-460d3444841d
# Optimise the variational paramaters and energy for Frohlich (k-space)
vfα, wfα, Efα, Kfα, Pfα = unzip(optimize.(frohlich_energies, 5, 2, upper=[100, 100]))

# ╔═╡ d11918c4-cb56-4566-9b4f-e7ddec05cdc5
plot(αrange, [vfα wfα vfaα wfaα], title="GS Frohlich Analytic and K-Space v,w", linestyle=:auto, linewidth = 2, xlabel = "α", ylabel = "v, w", labels = ["v k-space" "w k-space" "v analytic" "w analytic"], yticks = 1:2:19)

# ╔═╡ 52f64be4-46b8-4b5a-adc8-500e5d3a11bc
plot(αrange, [Efα Kfα Pfα Efaα Kfaα Pfaα], title="GS Frohlich Analytic and K-Space Energies", linestyle=:auto, linewidth = 2, legend = :bottomleft, xlabel = "α", ylabel = "Total Energy, Kinetic and Potential (J)", labels = ["E k-space" "K k-space" "P k-space" "E analytic" "K analytic" "P analytic"], yticks = -28:4:12)

# ╔═╡ 78878f4e-7c74-4a61-862e-b934cfedb10c
plot(αrange, [vhaα whaα vfaα wfaα], title="GS Comparison of Holstein and Frohlich v & w", linestyle=:auto, linewidth = 2, xlabel = "α", ylabel = "v, w", labels = ["v Holstein" "w Holstein" "v Frohlich" "w Frohlich"], yticks = 1:2:19)

# ╔═╡ a930fe34-d9e3-4dd7-bd9b-dd457e23587a
plot(αrange, [Ehaα .+ 6 Khaα .+ 6 Phaα Efaα Kfaα Pfaα], title="GS Comparison of Holstein and Frohlich Energies", linestyle=:auto, linewidth = 2, legend = :bottomleft, xlabel = "α", ylabel = "Total Energy, Kinetic and Potential", labels = ["E Holstein" "K Holstein" "P Holstein" "E Frohlich" "K Frohlich" "P Frohlich"], yticks = -28:4:12)

# ╔═╡ 180ca806-a36f-47e9-a287-07b5a9624a0c
md"## Dependence on Frequency Ω"

# ╔═╡ c1022521-a44f-4c4e-9fe7-92be7842651d
# Frequency Range
Ωrange = 0.1:0.1:20

# ╔═╡ 45b7b0d6-8ab8-4aa4-af9d-1da6dcdd7d26
begin
	α = 2
	β = 1 # Temperature T = 1
	best_holstein_energy_thermal(v, w) = holstein_energy_analytic(v, w, α, ω, β; dims = d)
	best_frohlich_energy_thermal(v, w) = frohlich_energy_analytic(v, w, α * 3, ω, β)
	best_holstein_energy_athermal(v, w) = holstein_energy_analytic(v, w, α, ω; dims = d)
	best_frohlich_energy_athermal(v, w) = frohlich_energy_analytic(v, w, α * 3, ω)
	vh, wh, Eh, Kh, Ph = optimize(best_holstein_energy_thermal, 5, 2, upper=[100, 100])
	vf, wf, Ef, Kf, Pf = optimize(best_frohlich_energy_thermal, 5, 2, upper=[100, 100])
	vh0, wh0, Eh0, Kh0, Ph0 = optimize(best_holstein_energy_athermal, 5, 2, upper=[100, 100])
	vf0, wf0, Ef0, Kf0, Pf0 = optimize(best_frohlich_energy_athermal, 5, 2, upper=[100, 100])
end

# ╔═╡ 36540645-f6b4-4415-95b2-4110cb8eab50
# ╠═╡ disabled = true
#=╠═╡
# Holstein mobility from k-space integral at T = 0 (Took 4000s)
χh0 = holstein_memory_function.(Ωrange, vh0, wh0, α, ω; dims = d)
  ╠═╡ =#

# ╔═╡ e31fd984-0a12-4c09-9249-3191d900aaec
# Saved Holstein Memory at zero temperature (takes a while to generate)
χh0_saved = [0.714308+0.000359172im
1.44904+0.000314577im
2.22692+0.000247171im
3.07586+0.000166002im
4.03345+8.24917e-5im
5.15511+6.33179e-6im
6.53083-5.38311e-5im
8.32657-9.29601e-5im
10.9193-0.000115636im
15.9211+0.0018468im
22.3336+7.56354im
22.2644+16.9262im
18.4343+24.7751im
12.4448+30.5718im
5.09489+34.4193im
-3.85212+36.6126im
-14.8963+32.0148im
-19.4273+25.0047im
-21.2523+18.8364im
-21.4796+13.5716im
-20.6724+9.18444im
-19.1477+5.60342im
-16.8872+2.88562im
-14.6384+1.50579im
-12.8168+0.74979im
-11.3388+0.333921im
-10.1371+0.118942im
-9.15834+0.0241481im
-8.36215-1.3158e-6im
-7.71908-1.9818e-6im
-7.17777-1.26125e-6im
-6.70973+1.30417e-5im
-6.29797+1.56332e-8im
-5.93061-1.6093e-5im
-5.59914+4.39882e-5im
-5.29697+2.60442e-6im
-5.01894+2.57487e-5im
-4.76102+3.13111e-6im
-4.51974-1.27429e-6im
-4.29228-2.31081e-5im
-4.076-1.96178e-7im
-3.86884+1.19525e-5im
-3.66874+1.41835e-5im
-3.47369-1.48267e-6im
-3.28181-2.04871e-8im
-3.09105+6.38322e-5im
-2.89929-1.71331e-6im
-2.70401+1.00471e-5im
-2.50229-6.32094e-5im
-2.29041+6.75686e-6im
-2.06349+5.40051e-6im
-1.81473+2.27459e-6im
-1.53407-2.05391e-6im
-1.20534-1.71706e-5im
-0.799332-1.02498e-5im
-0.248097-1.19029e-5im
0.656061+0.157971im
1.70032+0.89669im
2.49314+2.21208im
2.78457+3.95708im
2.39069+5.96159im
1.0439+8.06967im
-1.90589+8.92784im
-3.99312+8.36222im
-5.57254+7.41574im
-6.71421+6.20294im
-7.42395+4.82683im
-7.68806+3.37489im
-7.36461+1.9653im
-6.70259+1.13723im
-6.0984+0.629242im
-5.55941+0.312303im
-5.08888+0.125754im
-4.68775+0.0308771im
-4.35718+0.000357018im
-4.10323+3.4417e-6im
-3.89773-6.66468e-6im
-3.72294+1.70585e-6im
-3.5704-2.55643e-6im
-3.43476+1.37726e-5im
-3.31241+1.46601e-5im
-3.20075+2.19683e-5im
-3.09776-3.25089e-7im
-3.00188-7.6029e-6im
-2.91185+3.47485e-5im
-2.8269+4.90909e-6im
-2.74587-6.53358e-7im
-2.66817+6.96592e-6im
-2.59311-4.50089e-8im
-2.5201+1.44941e-7im
-2.44853+2.41606e-7im
-2.37783+2.84775e-7im
-2.3074+3.31745e-7im
-2.23655+4.10858e-7im
-2.16452+5.18092e-7im
-2.0904+6.01349e-7im
-2.01305-1.85874e-5im
-1.93093+3.61666e-7im
-1.84196+1.15586e-7im
-1.743-1.21993e-5im
-1.62911-1.47657e-6im
-1.49122-2.0946e-6im
-1.30609+0.00340405im
-1.04138+0.0556191im
-0.717158+0.230313im
-0.399561+0.581906im
-0.185613+1.13675im
-0.219172+1.89578im
-0.842331+2.5712im
-1.49047+2.8021im
-2.13768+2.84262im
-2.75369+2.68644im
-3.28387+2.34253im
-3.66409+1.83055im
-3.77025+1.18428im
-3.59515+0.742127im
-3.39959+0.446935im
-3.1964+0.242298im
-2.9969+0.107955im
-2.8113+0.0310045im
-2.65002+0.00116994im
-2.52918+8.75116e-6im
-2.43655+1.81599e-5im
-2.35938-4.82626e-6im
-2.29283+1.01454e-6im
-2.23407+1.06562e-6im
-2.18127+7.02253e-6im
-2.13319-3.89211e-6im
-2.0889+5.94968e-8im
-2.04772-4.91017e-6im
-2.00908-1.00413e-5im
-1.97261-8.23645e-7im
-1.9379-2.68213e-6im
-1.90472-1.87335e-6im
-1.87276+5.30589e-7im
-1.84183-5.14919e-6im
-1.81172+4.67904e-7im
-1.78222+5.19284e-7im
-1.75316+4.36296e-7im
-1.72432-1.07579e-6im
-1.69551+4.94391e-8im
-1.66647-1.3375e-7im
-1.63691-2.49325e-7im
-1.60648-2.84875e-7im
-1.57472-2.57578e-7im
-1.54096-2.02246e-7im
-1.50427-1.58817e-7im
-1.4631-1.51812e-7im
-1.41454+3.46395e-5im
-1.35208+0.00244588im
-1.26867+0.0179816im
-1.16501+0.0656179im
-1.05532+0.168015im
-0.97485+0.347109im
-1.0267+0.578631im
-1.14292+0.734356im
-1.30573+0.852048im
-1.50792+0.908942im
-1.72727+0.886664im
-1.93014+0.772324im
-2.05627+0.558844im
-2.03868+0.376117im
-1.99826+0.24584im
-1.94051+0.145045im
-1.87113+0.0711513im
-1.79666+0.023553im
-1.72514+0.00171474im
-1.67024-3.19806e-6im
-1.63027+3.04318e-7im
-1.5974-7.86344e-6im
-1.56911+2.39457e-7im
-1.5441+1.05948e-7im
-1.52146+4.42655e-7im
-1.50073+3.62264e-6im
-1.48151+4.02764e-7im
-1.46352-8.5689e-7im
-1.44654-1.40424e-6im
-1.43041-8.76726e-7im
-1.41498-1.67512e-6im
-1.40018-3.0136e-7im
-1.38589-9.81452e-6im
-1.37204+1.6597e-8im
-1.35856+1.77867e-7im
-1.34539+2.82861e-7im
-1.33246+3.06497e-7im
-1.31972+2.48798e-7im
-1.30711-7.83314e-6im
-1.29458-1.56239e-8im
-1.28203+1.38553e-5im
-1.26941-2.29906e-7im
-1.25658-3.34723e-7im
-1.24343-2.01715e-7im
-1.22976+3.52006e-5im
-1.21529+1.10967e-6im
-1.19959-1.47843e-5im
-1.18172+7.16373e-5im
-1.16006+0.00103244im
-1.13273+0.00568186im
-1.09947+0.0193611im
-1.06446+0.0500102im]

# ╔═╡ 7411b9e8-e515-4e4a-9f24-0d3c36844fb5
# Hosltein mobility from analytic at T = 0
χha0 = holstein_memory_function_analytic.(Ωrange, vh0, wh0, α, ω; dims = d)

# ╔═╡ 3ee52dc5-3665-4d0a-a9a2-c50c9ff8219c
plot(Ωrange, [real.(χh0_saved) imag.(χh0_saved) real.(χha0) imag.(χha0)], title="Holstein Memory for K-Space and Analytic T = 0", linestyle=:auto, linewidth = 2, xlabel = "Frequency Ω (J/ħ)", ylabel = "Memory Function χ(Ω)", labels = ["ℜχ k-space" "ℑχ k-space" "ℜχ analytic" "ℑχ analytic"])

# ╔═╡ d18c413d-7262-4041-9ca9-fd1c9e84b44a
# Holstein mobility from k-space integral at T = 1
χh = holstein_memory_function.(Ωrange, vh, wh, α, ω, β; dims = d)

# ╔═╡ e1df0094-3bfc-4cd4-8c28-d58804bd0209
# Hosltein mobility from analytic at T = 1
χha = holstein_memory_function_analytic.(Ωrange, vh, wh, α, ω, β; dims = d)

# ╔═╡ 5bbffee3-d9d1-4001-b17a-e951747271ac
plot(Ωrange, [real.(χh) imag.(χh) real.(χha) imag.(χha)], title="Holstein Memory for K-Space and Analytic T = 1", linestyle=:auto, linewidth = 2, xlabel = "Frequency Ω (J/ħ)", ylabel = "Memory Function χ(Ω)", labels = ["ℜχ k-space" "ℑχ k-space" "ℜχ analytic" "ℑχ analytic"])

# ╔═╡ b50f0096-ec99-4a7a-b076-be02dbb6f66d
plot(Ωrange, [real.(χha0) imag.(χha0) real.(χha) imag.(χha)], title="Holstein Memory for T = 0 and 1", linestyle=:auto, linewidth = 2, xlabel = "Frequency Ω (J/ħ)", ylabel = "Memory Function χ(Ω)", labels = ["ℜχ T=0" "ℑχ T=0" "ℜχ T = 1" "ℑχ T = 1"])

# ╔═╡ 49094531-c2fe-4eec-9740-c4b8fe22e578
# Frohlich mobility from analytic T = 1
χfa = frohlich_memory_function_analytic.(Ωrange, vf, wf, α * 3, ω, β)

# ╔═╡ b01bd192-90d9-486b-9847-f18ccc516b16
# Frohlich mobility from analytic T = 0
χfa0 = frohlich_memory_function_analytic.(Ωrange, vf0, wf0, α * 3, ω)

# ╔═╡ 892ef1e9-716a-4118-9415-863acc3fd9a6
plot(Ωrange, [real.(χha0) imag.(χha0) real.(χfa0) imag.(χfa0)], title="Comparison of Holstein and Frohlich Memory T = 0", linestyle=:auto, linewidth = 2, xlabel = "Frequency Ω (J/ħ)", ylabel = "Memory Function χ(Ω)", labels = ["ℜχ Holstein" "ℑχ Holstein" "ℜχ Frohlich" "ℑχ Frohlich"])

# ╔═╡ df2af8d4-9356-4bf9-a4ab-423eca2a8a83
plot(Ωrange, [real.(χha) imag.(χha) real.(χfa) imag.(χfa)], title="Comparison of Holstein and Frohlich Memory T = 1", linestyle=:auto, linewidth = 2, xlabel = "Frequency Ω (J/ħ)", ylabel = "Memory Function χ(Ω)", labels = ["ℜχ Holstein" "ℑχ Holstein" "ℜχ Frohlich" "ℑχ Frohlich"])

# ╔═╡ 63e75af5-bdfb-4152-9bcb-0110aa13a193
md"## Dependence on Temperature T = 1/β"

# ╔═╡ 785fd2a7-5db6-4242-9ed0-f3c7fd1b4955
# Temperature range T = 0.1 to 100
βrange = [1/x for x in 0.1:0.01:100]

# ╔═╡ b304bc21-d50e-4c84-8aa9-0daba2f98de9
# Holstein DC mobility (k-space)
μh = holstein_mobility.(vh, wh, α, ω, βrange; dims = d)

# ╔═╡ 3ad7f9d9-1620-4e9d-b577-2ba723220f4a
# Holstein DC mobility (analytic)
μha = holstein_mobility_analytic.(vh, wh, α, ω, βrange; dims = d)

# ╔═╡ d3f94402-1cb6-4aa5-8eb2-5d94040af4eb
plot(1 ./ βrange, [μh μha], yaxis = :log, xaxis = :log, title="Holstein K-Space and Analytic Mobility", linestyle=:auto, linewidth = 2, xlabel = "Temperature T (J/kB)", ylabel = "DC Mobility μ(T)", labels = ["μ k-space" "μ analytic"], yticks = ([2.0^x for x in -4:5], ["2^$x" for x in -4:5]), xticks = ([2.0^x for x in -3:6], ["2^$x" for x in -3:6]))

# ╔═╡ ef7889fd-1cf8-4424-ab0f-eeb32e2b6a30
# Frohlich DC mobility (analytic)
μf = frohlich_mobility_analytic.(vf, wf, α * 3, ω, βrange)

# ╔═╡ b4c58d45-f857-4e3b-a579-1b2ca930c778
plot(1 ./ βrange, [μf μh], yaxis = :log, xaxis = :log, title="Comparison of Holstein and Frohlich Mobility", linestyle=:auto, linewidth = 2, xlabel = "Temperature T (J/kB)", ylabel = "DC Mobility μ(T)", labels = ["μ Holstein" "μ Frohlich"], yticks = ([2.0^x for x in -6:5], ["2^$x" for x in -6:5]), xticks = ([2.0^x for x in -3:6], ["2^$x" for x in -3:6]))

# ╔═╡ 3740c483-035d-4e54-9588-50ec52eb4cbd
md"## Dependence on Adiabaticity / Phonon Frequency ω"

# ╔═╡ dc8ffa26-381a-4647-bccf-2fceb3dd5d69
# Adiabaticity / Phonon Frequency range ω = 0.1 to 2
ωrange = 0.1:0.01:2

# ╔═╡ 75ed1295-a710-4fa8-80ab-7be1af3fc54a
begin
	holstein_energies_ω = []
	frohlich_energies_ω = []
	holstein_energies_analytic_ω = []
	frohlich_energies_analytic_ω = []
	for ω in ωrange
		best_holstein_energy_ω(v, w) = holstein_energy(v, w, α, ω; dims = d)
		best_frohlich_energy_ω(v, w) = frohlich_energy(v, w, α .* 3, ω)
		best_holstein_energy_analytic_ω(v, w) = holstein_energy_analytic(v, w, α, ω; dims = d)
		best_frohlich_energy_analytic_ω(v, w) = frohlich_energy_analytic(v, w, α .* 3, ω)
		push!(holstein_energies_ω, best_holstein_energy_ω)
		push!(frohlich_energies_ω, best_frohlich_energy_ω)
		push!(holstein_energies_analytic_ω, best_holstein_energy_analytic_ω)
		push!(frohlich_energies_analytic_ω, best_frohlich_energy_analytic_ω)
	end
end

# ╔═╡ c4567ab7-7a04-4a66-8ffc-78992122ad9a
# Optimise the variational paramaters and energy for Holstein (k-space)
vhω, whω, Ehω, Khω, Phω = unzip(optimize.(holstein_energies_ω, 5, 2, upper=[200, 200]))

# ╔═╡ a11b5f7e-b252-4735-b539-95dbd37dadfd
# Optimise the variational paramaters and energy for Frohlich (k-space)
vfω, wfω, Efω, Kfω, Pfω = unzip(optimize.(frohlich_energies_ω, 5, 2, upper=[200, 200]))

# ╔═╡ d05e8814-360c-419d-8b2b-16f8cde823ad
# Optimise the variational paramaters and energy for Holstein (analytic)
vhaω, whaω, Ehaω, Khaω, Phaω = unzip(optimize.(holstein_energies_analytic_ω, 5, 2, upper=[200, 200]))

# ╔═╡ 89c437e9-b13e-4d71-a789-a2bd1a33c63b
# Optimise the variational paramaters and energy for Frohlich (analytic)
vfaω, wfaω, Efaω, Kfaω, Pfaω = unzip(optimize.(frohlich_energies_analytic_ω, 5, 2, upper=[200, 200]))

# ╔═╡ 4275c478-b59d-4574-805e-b1fed99e292d
plot(ωrange, [vhω whω vhaω whaω], title="GS Holstein Analytic and K-Space v,w", linestyle=:auto, linewidth = 2, xlabel = "ω", ylabel = "v, w", labels = ["v k-space" "w k-space" "v analytic" "w analytic"])

# ╔═╡ 3e029c3e-cb3a-42af-90a1-40ccf33f096a
plot(ωrange, [Ehω Khω Phω Ehaω Khaω Phaω], title="GS Holstein Analytic and K-Space Energies", linestyle=:auto, linewidth = 2, xlabel = "ω", ylabel = "Total Energy, Kinetic and Potential (J)", labels = ["E k-space" "K k-space" "P k-space" "E analytic" "K analytic" "P analytic"])

# ╔═╡ 8ea0f8d6-4e9c-4278-ac8e-cb9e149c3092
plot(ωrange, [vfω wfω vfaω wfaω], title="GS Frohlich Analytic and K-Space v,w", linestyle=:auto, linewidth = 2, xlabel = "ω", ylabel = "v, w", labels = ["v k-space" "w k-space" "v analytic" "w analytic"])

# ╔═╡ b0fd79fc-6da4-4fa6-b20d-b214258969fa
plot(ωrange, [Efω Kfω Pfω Efaω Kfaω Pfaω], title="GS Frohlich Analytic and K-Space Energies", linestyle=:auto, linewidth = 2, legend = :bottomleft, xlabel = "ω", ylabel = "Total Energy, Kinetic and Potential (J)", labels = ["E k-space" "K k-space" "P k-space" "E analytic" "K analytic" "P analytic"])

# ╔═╡ d849bfe9-2eaa-48e2-b0d1-7b7cfb0686e8
plot(ωrange, [vhaω whaω vfaω wfaω], title="GS Comparison of Holstein and Frohlich v & w", linestyle=:auto, linewidth = 2, xlabel = "ω", ylabel = "v, w", labels = ["v Holstein" "w Holstein" "v Frohlich" "w Frohlich"])

# ╔═╡ a0f847b4-ed1e-45cf-bf96-974f4ae5d358
plot(ωrange, [Ehaω .+ 6 Khaω .+ 6 Phaω Efaω Kfaω Pfaω], title="GS Comparison of Holstein and Frohlich Energies", linestyle=:auto, linewidth = 2, legend = :bottomleft, xlabel = "ω", ylabel = "Total Energy, Kinetic and Potential", labels = ["E Holstein" "K Holstein" "P Holstein" "E Frohlich" "K Frohlich" "P Frohlich"])

# ╔═╡ 999fb552-2430-4956-98e7-1faddd6e9a3e
md"# Self Consistent Extension"

# ╔═╡ ede8cba2-cb70-4c6f-94ef-2474ad7e4823
function polaron_self_energy(τ, impedence, β)
	integrand(Ω) = (
	(1 - exp(-Ω * τ)) / (1 - exp(-β * Ω)) + (1 - exp(Ω * τ)) / (exp(β * Ω) - 1)
	) / impedence(Ω)

	integral, _ = quadgk(Ω -> integrand(Ω), -Inf, Inf)

	return integral / (2π * im)
end

# ╔═╡ 538f0863-b572-4cf9-97e5-ceafac10b5be
function polaron_impedence(Ω, self_energy, β)
	integrand(τ) = (1 - cos(Ω * τ))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[compat]
Optim = "~1.7.6"
Plots = "~1.38.16"
QuadGK = "~2.8.2"
SpecialFunctions = "~2.3.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.1"
manifest_format = "2.0"
project_hash = "d6c55742439aa5a81796cef27d56a91be9ba8b7e"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f83ec24f76d4c8f525099b2ac475fc098138ec31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.11"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.ChangesOfVariables]]
deps = ["InverseFunctions", "LinearAlgebra", "Test"]
git-tree-sha1 = "2fba81a302a7be671aefe194f0525ef231104e7f"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.8"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "dd3000d954d483c1aad05fe1eb9e6a715c97013e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.22.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "4e88377ae7ebeaf29a047aa1ee40826e0b708a5d"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.7.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fe2838a593b5f776e1597e086dcd47560d94e816"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.3"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "cf25ccb972fec4e4817764d01c82386ae94f77b4"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.14"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "f0af9b12329a637e8fba7d6543f915fff6ba0090"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.4.2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "d73afa4a2bb9de56077242d98cf763074ab9a970"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.9"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f61f768bf090d97c532d24b64e07b237e9bb7b6b"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.9+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "cb56ccdd481c0dd7f975ad2b3b62d9eda088f7e2"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.14"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "eabe3125edba5c9c10b60a160b1779a000dc8b29"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.11"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "c3ce8e7420b3a6e071e0fe4745f5d4300e37b13f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.24"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1aa4b74f80b01c6bc2b89992b861b5f210e665b5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.21+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "e3a6546c1577bfd701771b477b794a52949e7594"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.6"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "4b2e829ee66d4218e0cef22c0a64ee37cf258c29"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "75ca67b2c6512ad2d0c767a7cfc55e75075f8bbc"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.16"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "364898e8f13f7eaaceec55fd3d08680498c0aa6e"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.4.2+3"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "7beb031cf8145577fbccacd94b8a8f4ce78428d3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "9cabadf6e7cd2349b6cf49f1915ad2028d65e881"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "Random"]
git-tree-sha1 = "c4d2a349259c8eba66a00a540d550f122a3ab228"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.15.0"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2222b751598bd9f4885c9ce9cd23e83404baa8ce"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.3+1"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═cd6cae06-27d1-11ee-1f4b-07833f940a86
# ╠═9ed21af8-f7c9-4cba-b172-e3b376665add
# ╟─1e3e504b-5389-4ef1-9d63-5ac4e10364ef
# ╠═e005016f-36e2-47a6-9172-8e32ac8bbff6
# ╠═6c51fe56-6a20-48c2-b6e4-4075350b807e
# ╠═14b28d72-9bb9-47fb-a374-c0ec287d8440
# ╠═7bf51005-7e61-4f8d-b5ea-db37aae30988
# ╟─3a7d73db-2d2b-43b8-ae39-adc01d9f5b30
# ╠═2c04db4c-6366-48f5-96fb-78a579c727b1
# ╠═50639eb7-3047-4a95-a668-f67a66775e4d
# ╠═256f3d5d-f474-4151-a3db-bec8099615fa
# ╠═237bf7db-2e0e-4c82-8ec9-4ac65ddfdd5c
# ╟─65093555-b406-45c5-854a-fdea6dad4c9c
# ╠═d295e844-ff9f-497e-993a-101f5c2ea0a8
# ╠═e5bbf2cc-2650-4666-9650-5a91ecea304d
# ╟─577587db-25b0-463b-ab44-38522bcd0039
# ╠═f73c700b-2f00-4deb-85cc-abdaaa4118e5
# ╠═4143c668-699d-46ea-b825-bc4996fb26a4
# ╠═57ed148b-fef1-40f4-a113-a39a1aebb010
# ╠═f223d349-a0d7-4fc6-a6ed-f8b3c45c8dad
# ╠═ac60f707-8a87-47b8-9856-80194fb5161c
# ╠═949b3f46-146c-4957-a14c-4c742d6d10ba
# ╠═5fa89b72-fab5-45b3-831c-590021953079
# ╠═303b4158-8a39-4733-84d2-ff1dbe6bb556
# ╠═0fef496e-6555-4f07-9586-c5d18314e386
# ╠═d94ab3bf-ef6e-4c1b-ba54-c025916e42e1
# ╠═62faa65d-d426-4d53-9816-ca1de8b24874
# ╠═293773c1-f0f1-4f74-8434-3b82f6cf8883
# ╟─b175dd29-6925-4352-9914-91d8e2f5fa9c
# ╠═0e99ecd9-f2ff-4536-8d8d-a04e42f6c98f
# ╠═c1dd1475-4ad6-478f-9c18-7c20df32be81
# ╠═fe21e743-2cf4-4e62-857a-aebbe539287e
# ╠═fe0642f4-09b9-4209-baf8-959d46d593ad
# ╠═b2622ccb-9c3d-40a3-9a4d-f248e90632dc
# ╠═7b893719-994e-4a8b-83f8-67b61f158189
# ╠═4bdfeceb-01a1-43f4-8634-54ed51f5ad82
# ╠═11a19771-8cf0-4f4a-98d4-5a1060d500de
# ╠═0c5fc742-029a-4f47-8897-c9c1483e90bd
# ╠═9e8e3e8f-095a-4bd0-a44e-54e9a0be3ebe
# ╠═fef541dc-568a-477f-9343-bee877a46b4a
# ╠═3a1ade3b-c89e-4137-9a7d-2faa572bea7c
# ╟─9a43b38c-8769-426c-99a9-98cdc6fbebac
# ╠═f93fc7a5-f59d-4158-aa5b-eec93598cdf9
# ╠═f122d4aa-3d4b-43e7-9370-4fe64efaa961
# ╠═fb7b8e2a-6d17-42ae-b264-5c3434643ec4
# ╠═5d2f12ce-5e37-4f31-91e7-e4efbd4442ce
# ╠═72196d1e-a5cf-4ffd-9d6a-63ab54a83d54
# ╠═f3fe677b-5edd-4301-b35c-e707c2106649
# ╠═9684cc0e-aebf-4e51-86d6-7591266bc356
# ╠═302b7385-a924-4849-8882-77d6a473e504
# ╠═f2a175aa-6854-4b0b-b811-2dd0783aca60
# ╠═9aa269c6-1df3-4553-ba61-c0b0fd2608c9
# ╠═24fec156-91bd-4761-9933-e7c5ade3693b
# ╠═ae71295e-6fab-4360-b6ec-e2d4db10f9ae
# ╟─7560cf27-f085-4f74-99f0-d33bfa542a9f
# ╠═e8cf5096-e167-4a66-9fbd-ec1181602e77
# ╠═ea7f5541-84b8-453a-aa0f-84f7df119b3f
# ╠═03b981b7-d257-4c89-a16d-613815486ac1
# ╠═20afdee0-d3be-48fc-8f61-705680397432
# ╠═95ae687b-1cd3-43f5-bf28-d8a4f5cca022
# ╠═92964c87-3fa3-40de-8f08-4b773b229a6a
# ╠═7e00a70f-5a52-42f5-9056-efc47b4ab839
# ╠═b3818117-e954-4081-ad92-ce066d0fc1b6
# ╠═f2ede4d0-d51e-4e55-af43-cb22fffe30f4
# ╠═bb9fe3b2-c1d4-4806-8b14-d234e91ee8a9
# ╟─ed021864-b2de-4f4b-8aac-175896e7bd39
# ╟─b10da998-6267-49b1-88ce-09ee7774eebc
# ╠═9b079f65-9dc9-4578-bba0-c55dbe9563ef
# ╠═1641996e-41f2-4e31-bb1d-bc7a1763274c
# ╠═0ae0ef1d-0c63-43c6-8e67-55b44dba6b67
# ╠═1b2e2933-623f-4560-9527-e24fe55ca485
# ╠═37387a74-af02-4485-9edf-3a950fcb822a
# ╠═ea012194-878d-433b-bf97-d2e8e90e7125
# ╠═2be6279b-f73e-4253-8f30-460d3444841d
# ╠═d11918c4-cb56-4566-9b4f-e7ddec05cdc5
# ╠═52f64be4-46b8-4b5a-adc8-500e5d3a11bc
# ╠═78878f4e-7c74-4a61-862e-b934cfedb10c
# ╠═a930fe34-d9e3-4dd7-bd9b-dd457e23587a
# ╟─180ca806-a36f-47e9-a287-07b5a9624a0c
# ╠═c1022521-a44f-4c4e-9fe7-92be7842651d
# ╠═45b7b0d6-8ab8-4aa4-af9d-1da6dcdd7d26
# ╠═36540645-f6b4-4415-95b2-4110cb8eab50
# ╟─e31fd984-0a12-4c09-9249-3191d900aaec
# ╠═7411b9e8-e515-4e4a-9f24-0d3c36844fb5
# ╠═3ee52dc5-3665-4d0a-a9a2-c50c9ff8219c
# ╠═d18c413d-7262-4041-9ca9-fd1c9e84b44a
# ╠═e1df0094-3bfc-4cd4-8c28-d58804bd0209
# ╠═5bbffee3-d9d1-4001-b17a-e951747271ac
# ╠═b50f0096-ec99-4a7a-b076-be02dbb6f66d
# ╠═49094531-c2fe-4eec-9740-c4b8fe22e578
# ╠═b01bd192-90d9-486b-9847-f18ccc516b16
# ╠═892ef1e9-716a-4118-9415-863acc3fd9a6
# ╠═df2af8d4-9356-4bf9-a4ab-423eca2a8a83
# ╟─63e75af5-bdfb-4152-9bcb-0110aa13a193
# ╠═785fd2a7-5db6-4242-9ed0-f3c7fd1b4955
# ╠═b304bc21-d50e-4c84-8aa9-0daba2f98de9
# ╠═3ad7f9d9-1620-4e9d-b577-2ba723220f4a
# ╠═d3f94402-1cb6-4aa5-8eb2-5d94040af4eb
# ╠═ef7889fd-1cf8-4424-ab0f-eeb32e2b6a30
# ╠═b4c58d45-f857-4e3b-a579-1b2ca930c778
# ╟─3740c483-035d-4e54-9588-50ec52eb4cbd
# ╠═dc8ffa26-381a-4647-bccf-2fceb3dd5d69
# ╠═75ed1295-a710-4fa8-80ab-7be1af3fc54a
# ╠═c4567ab7-7a04-4a66-8ffc-78992122ad9a
# ╠═a11b5f7e-b252-4735-b539-95dbd37dadfd
# ╠═d05e8814-360c-419d-8b2b-16f8cde823ad
# ╠═89c437e9-b13e-4d71-a789-a2bd1a33c63b
# ╠═4275c478-b59d-4574-805e-b1fed99e292d
# ╠═3e029c3e-cb3a-42af-90a1-40ccf33f096a
# ╠═8ea0f8d6-4e9c-4278-ac8e-cb9e149c3092
# ╠═b0fd79fc-6da4-4fa6-b20d-b214258969fa
# ╠═d849bfe9-2eaa-48e2-b0d1-7b7cfb0686e8
# ╠═a0f847b4-ed1e-45cf-bf96-974f4ae5d358
# ╠═999fb552-2430-4956-98e7-1faddd6e9a3e
# ╠═ede8cba2-cb70-4c6f-94ef-2474ad7e4823
# ╠═538f0863-b572-4cf9-97e5-ceafac10b5be
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
