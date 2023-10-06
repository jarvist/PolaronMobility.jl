# KSpaceTheory.jl

"""
    cartesian_k_integrand(k, coupling, propagator; rₚ = 1)

Calculate the integrand for a k-space integral in Cartesian coordinates.

# Arguments
- `k`: a scalar value representing the k-coordinate in k-space
- `coupling`: a function that takes a scalar k value and returns a scalar value representing the coupling strength
- `propagator`: a scalar value representing the propagator
- `rₚ`: an optional scalar value representing the radius

# Example Usage
```julia
coupling(k) = k^2  # define the coupling function
propagator = 0.5  # define the propagator value
result = cartesian_k_integrand(2.0, coupling, propagator; rₚ = 1)  # calculate the integrand for k = 2.0
println(result)  # print the result
```
Expected Output:
0.5
"""
function cartesian_k_integrand(k, coupling, propagator; rₚ = 1)
    coupling(k) * exp(-k^2 * rₚ^2 * propagator / 2)
end


"""
    spherical_k_integrand(k, coupling, propagator; rₚ = 1)

Calculate the integrand for a k-space integral in spherical coordinates.

## Arguments
- `k`: a scalar value representing the k-coordinate in k-space
- `coupling`: a function that takes a scalar k value and returns a scalar value representing the coupling strength
- `propagator`: a scalar value representing the propagator
- `rₚ`: an optional scalar value representing the radius

## Returns
The calculated integrand as a scalar value.

## Example
```julia
coupling(k) = k^2  # define the coupling function
propagator = 0.5  # define the propagator value
result = spherical_k_integrand(2.0, coupling, propagator; rₚ = 1)  # calculate the integrand for k = 2.0
println(result)  # print the result
```
Expected Output:
`0.5`
"""
function spherical_k_integrand(k, coupling, propagator; rₚ = 1)
    4π * k^2 * coupling(k) * exp(-k^2 * rₚ^2 * propagator / 2)
end


"""
    cartesian_k_integral(coupling, propagator; rₚ = 1, a = 1, limits = [-π, π])

Calculate the k-space integral in cartesian coordinates of the integrand `cartesian_k_integrand` over a specified range in k-space.

# Arguments
- `coupling`: A function that takes a scalar k value and returns a scalar value representing the coupling strength.
- `propagator`: A scalar value representing the propagator.
- `rₚ`: An optional scalar value representing the charactersitc polaron radius. Default value is 1.
- `a`: An optional scalar value representing ther lattice constant. Default value is 1.
- `limits`: An array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is 1D Brillouin zone.

# Returns
A scalar value representing the calculated k-space integral over the specified range in cartesian coordinates.

# Example
```julia
coupling(k) = k^2  # define the coupling function
propagator = 0.5  # define the propagator value
result = cartesian_k_integral(coupling, propagator; rₚ = 1, a = 1, limits = [-π, π])  # calculate the integral
println(result)  # print the result
```
Expected Output:
A scalar value representing the calculated k-space integral over the specified range in cartesian coordinates.
"""
function cartesian_k_integral(coupling, propagator; rₚ = 1, a = 1, limits = [-π, π])
    integral, _ = quadgk(k -> cartesian_k_integrand(k, coupling, propagator; rₚ = rₚ), limits[1] / a, limits[2] / a)
    integral * a / 2π
end


"""
    spherical_k_integral(coupling, propagator; rₚ = 1, limits = [0, Inf])

Calculate the k-space integral in spherical coordinates of the integrand `spherical_k_integrand` over a specified radius in k-space.

# Arguments
- `coupling`: A function that takes a scalar k value and returns a scalar value representing the coupling strength.
- `propagator`: A scalar value representing the propagator.
- `rₚ`: An optional scalar value representing the characteristic polaron radius. Default value is 1.
- `limits`: An array of two scalar values representing the lower and upper limits of the integration range in spherical coordinates. Default is [0, Inf] for infinite sphere.

# Returns
A scalar value representing the calculated integral over the specified range in spherical coordinates.

# Example
```julia
coupling(k) = k^2  # define the coupling function
propagator = 0.5  # define the propagator value
result = spherical_k_integral(coupling, propagator; rₚ = 1, limits = [0, Inf])  # calculate the integral
println(result)  # print the result
```
Expected Output:
A scalar value representing the calculated k-space integral over the specified range in spherical coordinates.
"""
function spherical_k_integral(coupling, propagator; rₚ = 1, limits = [0, Inf])
    integral, _ = quadgk(k -> spherical_k_integrand(k, coupling, propagator; rₚ = rₚ), limits[1], limits[2])
    integral / 8π^3
end


"""
    holstein_coupling(k, α, ω; dims = 1)

Calculate the coupling strength for the Holstein lattice polaron model.

# Arguments
- `k`: a scalar value representing the k-coordinate in k-space
- `α`: a scalar value representing the coupling constant
- `ω`: a scalar value representing the phonon frequency
- `dims`: an optional scalar value representing the dimensionality of the system. Default value is 3.

# Returns
The coupling strength for the Holstein model.

# Example
```julia
result = holstein_coupling(2.0, 0.5, 1.0; dims = 3)
println(result)
```
Expected Output:
6.0
"""
function holstein_coupling(k, α, ω; dims = 3)
    2 * α * ω * dims
end


"""
    frohlich_coupling(k, α, ω)

Calculate the coupling strength for the Frohlich continuum polaron model.

# Arguments
- `k`: a scalar value representing the k-coordinate in k-space
- `α`: a scalar value representing the coupling constant
- `ω`: a scalar value representing the phonon frequency

# Returns
The coupling strength for the Frohlich continuum polaron model.

# Example
```julia
result = frohlich_coupling(2.0, 0.5, 1.0)
println(result)
```
Expected Output:
`6.0`
"""
function frohlich_coupling(k, α, ω)
    ω^(3/2) * 2√2 * π * α / k^2
end

"""
    holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])

Calculate the integrand for the Holstein interaction energy in k-space at finite temperature.

# Arguments
- `τ`: a scalar value representing the imaginary time.
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.
- `dims`: an optional scalar value representing the dimensionality of the system. Default value is 3.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an optional array of two scalar values representing the lower and upper limits of the integration range in k-space. Default value is [-π, π].

# Returns
A scalar value representing the integrand for the Holstein interaction energy in k-space at finite temperature.
"""
function holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])
    coupling(k) = holstein_coupling(k, α, ω; dims = dims)
    propagator = polaron_propagator(τ, v, w, β)
    phonon_propagator(τ, ω, β) * cartesian_k_integral(coupling, propagator; rₚ = rₚ, a = a, limits = limits)^dims
end


"""
    holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω; dims = 3, rₚ = 1, a = 1, limits = [-π, π])

Calculate the integrand for the Holstein interaction energy in k-space at zero temperature.

# Arguments
- `τ`: a scalar value representing the imaginary time.
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `dims`: an optional scalar value representing the dimensionality of the system. Default value is 3.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an optional array of two scalar values representing the lower and upper limits of the integration range in k-space. Default value is [-π, π].

# Returns
A scalar value representing the integrand for the Holstein interaction energy in k-space at zero temperature.
"""
function holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω; dims = 3, rₚ = 1, a = 1, limits = [-π, π])
    coupling(k) = holstein_coupling(k, α, ω; dims = dims)
    propagator = polaron_propagator(τ, v, w)
    phonon_propagator(τ, ω) * carteisan_k_integral(coupling, propagator; rₚ = rₚ, a = a, limits = limits)^dims
end


"""
    frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω, β)

Calculate the integrand for the Frohlich interaction energy in k-space at finite temperaure.

## Arguments
- `τ`: a scalar value representing the imaginary time.
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `limits`: an optional array of two scalar values representing the lower and upper limits of the integration range in k-space. Default value is [0, Inf].

## Returns
A scalar value representing the integrand for the Frohlich interaction energy in k-space at finite temperature.
"""
function frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; rₚ = 1, limits = [0, Inf])
    coupling(k) = frohlich_coupling(k, α, ω)
    propagator = polaron_propagator(τ, v, w, β)
    phonon_propagator(τ, ω, β) * spherical_k_integral(coupling, propagator; rₚ = rₚ, limits = limits)
end

"""
    frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω; rₚ = 1, limits = [0, Inf])

Calculate the integrand for the Frohlich interaction energy in k-space at finite zero temperaure.

## Arguments
- `τ`: a scalar value representing the imaginary time.
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `limits`: an optional array of two scalar values representing the lower and upper limits of the integration range in k-space. Default value is [0, Inf].

## Returns
A scalar value representing the integrand for the Frohlich interaction energy in k-space at zero temperature.
"""
function frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω; rₚ = 1, limits = [0, Inf])
	coupling(k) = frohlich_coupling(k, α, ω)
	propagator = polaron_propagator(τ, v, w, ω)
	phonon_propagator(τ, ω) * spherical_k_integral(coupling, propagator; rₚ = rₚ, limits = limits)
end


"""
    holstein_interaction_energy_k_space(v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])

Calculate the Holstein polaron interaction energy in k-space at finite temperaure.

## Arguments
- `τ`: a scalar value representing the imaginary time.
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.
- `dims`: an optional scalar value representing the dimensionality of the system. Default value is 3.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an optional array of two scalar values representing the lower and upper limits of the integration range in k-space. Default value is [-π, π].

## Returns
A scalar value representing the Holstein polaron interaction energy in k-space at finite temperature.
"""
function holstein_interaction_energy_k_space(v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])
    interation_energy, _ = quadgk(τ -> holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; dims = dims, rₚ = rₚ, a = a, limits = limits), 0, β/2)
    return interation_energy
end


"""
    holstein_interaction_energy_k_space(v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])

Calculate the Holstein polaron interaction energy in k-space at zero temperaure.

## Arguments
- `τ`: a scalar value representing the imaginary time.
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `dims`: an optional scalar value representing the dimensionality of the system. Default value is 3.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an optional array of two scalar values representing the lower and upper limits of the integration range in k-space. Default value is [-π, π].

## Returns
A scalar value representing the Holstein polaron interaction energy in k-space at zero temperature.
"""
function holstein_interaction_energy_k_space(v, w, α, ω; dims = 3, rₚ = 1, a = 1, limits = [-π, π])
	interation_energy, _ = quadgk(τ -> holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω; dims = dims, rₚ = rₚ, a = a, limits = limits), 0, Inf)
	return interation_energy
end

"""
    frohlich_interaction_energy_k_space(v, w, α, ω, β; rₚ = 1, limits = [0, Inf])

Calculate the Frohlich polaron interaction energy in k-space at finite temperaure.

## Arguments
- `τ`: a scalar value representing the imaginary time.
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an optional array of two scalar values representing the lower and upper limits of the integration range in k-space. Default value is [-π, π].

## Returns
A scalar value representing the Frohlich polaron interaction energy in k-space at finite temperature.
"""
function frohlich_interaction_energy_k_space(v, w, α, ω, β; rₚ = 1, limits = [0, Inf])
	interation_energy, _ = quadgk(τ -> frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; rₚ = rₚ, limits = limits), 0, β/2)
	return interation_energy
end

"""
    frohlich_interaction_energy_k_space(v, w, α, ω, β; rₚ = 1, limits = [0, Inf])

Calculate the Frohlich polaron interaction energy in k-space at zero temperaure.

## Arguments
- `τ`: a scalar value representing the imaginary time.
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an optional array of two scalar values representing the lower and upper limits of the integration range in k-space. Default value is [-π, π].

## Returns
A scalar value representing the Frohlich polaron interaction energy in k-space at zero temperature.
"""
function frohlich_interaction_energy_k_space(v, w, α, ω; rₚ = 1, limits = [0, Inf])
	interation_energy, _ = quadgk(τ -> frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω; rₚ = rₚ, limits = limits), 0, Inf)
	return interation_energy
end


"""
    electron_energy(v, w, ω, β)

Calculate the free electron energy at finite temperature.

# Arguments
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.

# Returns
A scalar value representing the calculated free electron energy at finite temperature.

# Example
```julia
v = 0.5
w = 1.0
ω = 2.0
β = 0.2
result = electron_energy(v, w, ω, β; dims = 3)
println(result)
```
Expected Output:
A scalar value representing the calculated free electron energy at finite temperature.
"""
function electron_energy(v, w, ω, β)
    -(A(v, w, ω, β) + C(v, w, ω, β)) / 3
end


"""
    electron_energy(v, w, ω)

Calculate the free electron energy at zero temperature.

# Arguments
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `ω`: a scalar value representing the phonon frequency.

# Returns
A scalar value representing the calculated free electron energy at zero temperature.

# Example
```julia
v = 0.5
w = 1.0
ω = 2.0
β = 0.2
result = electron_energy(v, w, ω)
println(result)
```
Expected Output:
A scalar value representing the calculated free electron energy at finite temperature.
"""
function electron_energy(v, w, ω)
	(v - w)^2 / (4 * v) * ω
end


"""
    holstein_energy_k_space(v, w, α, ωβ...; dims = 3, rₚ = 1, a = 1, limits = [-π, π])

Calculate the total energy, kinetic energy, and interaction energy of the Holstein lattice polaron.

## Arguments
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: (optional) a scalar value representing the inverse temperature.
- `dims`: The number of dimensions of the system (default is 3).
- `rₚ`: The characteristic polaron radius (default is 1).
- `a`: The lattice constant (default is 1).
- `limits`: The limits of integration for the interaction energy calculation (default is [-π, π]).

## Returns
- `total_energy`: The calculated total polaron energy.
- `kinetic_energy`: The calculated polaron kinetic energy.
- `interaction_energy`: The calculated polaron interaction energy.
"""
function holstein_energy_k_space(v, w, α, ωβ...; dims = 3, rₚ = 1, a = 1, limits = [-π, π])
    kinetic_energy = -2 * dims + electron_energy(v, w, ωβ...) 
    interaction_energy = -holstein_interaction_energy_k_space(v, w, α, ωβ...; dims = dims, rₚ = rₚ, a = a, limits = limits)
    return kinetic_energy + interaction_energy, kinetic_energy, interaction_energy
end


"""
    frohlich_energy_k_space(v, w, α, ωβ...; rₚ = 1, limits = [0, Inf])

Calculate the total energy, kinetic energy, and interaction energy of the Frohlich lattice polaron.

## Arguments
- `v`: a scalar value representing a variational paramater.
- `w`: a scalar value representing a variational paramater.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: (optional) a scalar value representing the inverse temperature.
- `rₚ`: The characteristic polaron radius (default is 1).
- `limits`: The limits of integration for the interaction energy calculation (default is [0, Inf]).

## Returns
- `total_energy`: The calculated total polaron energy.
- `kinetic_energy`: The calculated polaron kinetic energy.
- `interaction_energy`: The calculated polaron interaction energy.
"""
function frohlich_energy_k_space(v, w, α, ωβ...; rₚ = 1, limits = [0, Inf])
	kinetic_energy = electron_energy(v, w, ωβ...) 
	interaction_energy = -frohlich_interaction_energy_k_space(v, w, α, ωβ...; rₚ = rₚ, limits = limits)
	return kinetic_energy + interaction_energy, kinetic_energy, interaction_energy
end


"""
    holstein_structure_factor_k_space(t, v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])

Calculate the structure factor in k-space for the Holstein lattice polaron model at finite temperature.

## Arguments
- `t`: a scalar value representing the real time.
- `v`: a scalar value representing a variational parameter.
- `w`: a scalar value representing a variational parameter.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.
- `dims`: an optional scalar value representing the dimensionality of the system. Default value is 3.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is 1D Brillouin zone.

## Returns
A scalar value representing the calculated structure factor in k-space for the Holstein lattice polaron model at finite temperature.
"""
function holstein_structure_factor_k_space(t, v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])
	
	coupling_one(k) = holstein_coupling(k, α, ω; dims = dims) * k^2
	coupling_two(k) = holstein_coupling(k, α, ω; dims = dims)

	propagator = polaron_propagator(im * t, v, w, β)
	
	first_integral = cartesian_k_integral(coupling_one, propagator; rₚ = rₚ, a = a, limits = limits)
	second_integral = cartesian_k_integral(coupling_two, propagator; rₚ = rₚ, a = a, limits = limits)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω, β)

	prefactor * dims * first_integral * second_integral^(dims - 1)
end


"""
    holstein_structure_factor_k_space(t, v, w, α, ω; dims = 3, rₚ = 1, a = 1, limits = [-π, π])

Calculate the structure factor in k-space for the Holstein lattice polaron model at zero temperature.

## Arguments
- `t`: a scalar value representing the real time.
- `v`: a scalar value representing a variational parameter.
- `w`: a scalar value representing a variational parameter.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `dims`: an optional scalar value representing the dimensionality of the system. Default value is 3.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is 1D Brillouin zone.

## Returns
A scalar value representing the calculated structure factor in k-space for the Holstein lattice polaron model at zero temperature.
"""
function holstein_structure_factor_k_space(t, v, w, α, ω; dims = 3, rₚ = 1, a = 1, limits = [-π, π])
	
	coupling_one(k) = holstein_coupling(k, α, ω; dims = dims) * k^2
	coupling_two(k) = holstein_coupling(k, α, ω; dims = dims)

	propagator = polaron_propagator(im * t, v, w)
	
	first_integral = cartesian_k_integral(coupling_one, propagator; rₚ = rₚ, a = a, limits = limits)
	second_integral = cartesian_k_integral(coupling_two, propagator; rₚ = rₚ, a = a, limits = limits)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω)

	prefactor * dims * first_integral * second_integral^(dims - 1)
end


"""
	frohlich_structure_factor_k_space(t, v, w, α, ω, β; rₚ = 1, limits = [0, Inf])

Calculate the structure factor in k-space for the Frohlich continuum polaron model at finite temperature.

## Arguments
- `t`: a scalar value representing the real time.
- `v`: a scalar value representing a variational parameter.
- `w`: a scalar value representing a variational parameter.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `limits`: an array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is an infinite sphere.

## Returns
A scalar value representing the calculated structure factor in k-space for the Frohlich continuum polaron model at finite temperature.
"""
function frohlich_structure_factor_k_space(t, v, w, α, ω, β; rₚ = 1, limits = [0, Inf])
	
	coupling(k) = frohlich_coupling(k, α, ω) * k^2

	propagator = polaron_propagator(im * t, v, w, β)
	
	integral = spherical_k_integral(coupling, propagator; rₚ = rₚ, limits = limits)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω, β)

	prefactor * integral
end


"""
	frohlich_structure_factor_k_space(t, v, w, α, ω; rₚ = 1, limits = [0, Inf])

Calculate the structure factor in k-space for the Frohlich continuum polaron model at zero temperature.

## Arguments
- `t`: a scalar value representing the real time.
- `v`: a scalar value representing the optimal variational parameter.
- `w`: a scalar value representing the optimal variational parameter.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `limits`: an array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is an infinite sphere.

## Returns
A scalar value representing the calculated structure factor in k-space for the Frohlich continuum polaron model at zero temperature.
"""
function frohlich_structure_factor_k_space(t, v, w, α, ω; rₚ = 1, limits = [0, Inf])
	
	coupling(k) = frohlich_coupling(k, α, ω) * k^2

	propagator = polaron_propagator(im * t, v, w)
	
	integral = spherical_k_integral(coupling, propagator; rₚ = rₚ, limits = limits)
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω)

	prefactor * integral
end


"""
    holstein_memory_function_k_space(Ω, v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])

Calculate the memory function for the Holstein model in k-space at finite temperature and frequency.

## Arguments
- `Ω`: a scalar value representing the frequency at which the memory function is evaluated.
- `v`: a scalar value representing the optimal variational parameter.
- `w`: a scalar value representing the optimal variational parameter.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.
- `dims`: an optional scalar value representing the dimensionality of the system. Default value is 3.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is 1D Brillouin zone.

## Returns
A scalar value representing the memory function of the Holstein model in k-space at finite temperature evaluated at the frequency `Ω`.
"""
function holstein_memory_function_k_space(Ω, v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])
	 structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, ω, β; dims = dims, rₚ = rₚ, a = a, limits = limits)
	 return general_memory_function(Ω, structure_factor)
end


"""
    holstein_memory_function_k_space(Ω, v, w, α, ω; dims = 3, rₚ = 1, a = 1, limits = [-π, π])

Calculate the memory function for the Holstein model in k-space at zero temperature and finite frequency.

## Arguments
- `Ω`: a scalar value representing the frequency at which the memory function is evaluated.
- `v`: a scalar value representing the optimal variational parameter.
- `w`: a scalar value representing the optimal variational parameter.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `dims`: an optional scalar value representing the dimensionality of the system. Default value is 3.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is 1D Brillouin zone.

## Returns
A scalar value representing the memory function of the Holstein model in k-space at zero temperature evaluated at the frequency `Ω`.
"""
function holstein_memory_function_k_space(Ω, v, w, α, ω; dims = 3, rₚ = 1, a = 1, limits = [-π, π])
	 structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, ω; dims = dims, rₚ = rₚ, a = a, limits = limits)
	 return general_memory_function(Ω, structure_factor, limits = [0, 1e4])
end


"""
    holstein_mobility_k_space(v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])

Calculate the DC mobility in k-space for a Holstein polaron system at finite temperature.

## Arguments
- `v`: a scalar value representing the optimal variational parameter.
- `w`: a scalar value representing the optimal variational parameter.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.
- `dims`: an optional scalar value representing the dimensionality of the system. Default value is 3.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `a`: an optional scalar value representing the lattice constant. Default value is 1.
- `limits`: an array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is 1D Brillouin zone.


## Returns
The DC mobility in k-space for the Holstein polaron system at finite temperature.
"""
function holstein_mobility_k_space(v, w, α, ω, β; dims = 3, rₚ = 1, a = 1, limits = [-π, π])
    structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, ω, β; dims = dims, rₚ = rₚ, a = a, limits = limits)
    1 / imag(general_memory_function(structure_factor))
end


"""
    frohlich_memory_function_k_space(Ω, v, w, α, ω, β; rₚ = 1, limits = [0, Inf])

Calculate the memory function for the Frohlich model in k-space at finite temperature and frequency.

## Arguments
- `Ω`: a scalar value representing the frequency at which the memory function is evaluated.
- `v`: a scalar value representing the optimal variational parameter.
- `w`: a scalar value representing the optimal variational parameter.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `limits`: an array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is an infinite sphere.

## Returns
A scalar value representing the memory function of the Frohlich model in k-space at finite temperature evaluated at the frequency `Ω`.
"""
function frohlich_memory_function_k_space(Ω, v, w, α, ω, β; rₚ = 1, limits = [0, Inf])
	 structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ω, β; rₚ = rₚ, limits = limits)
	 return general_memory_function(Ω, structure_factor; limits = [0, Inf])
end


"""
    frohlich_memory_function_k_space(Ω, v, w, α, ω; rₚ = 1, limits = [0, Inf])

Calculate the memory function for the Frohlich model in k-space at zero temperature and finite frequency.

## Arguments
- `Ω`: a scalar value representing the frequency at which the memory function is evaluated.
- `v`: a scalar value representing the optimal variational parameter.
- `w`: a scalar value representing the optimal variational parameter.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `limits`: an array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is an infinite sphere.

## Returns
A scalar value representing the memory function of the Frohlich model in k-space at finite temperature evaluated at the frequency `Ω`.
"""
function frohlich_memory_function_k_space(Ω, v, w, α, ω; rₚ = 1, limits = [0, Inf])
	 structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ω; rₚ = rₚ, limits = limits)
	 return general_memory_function(Ω, structure_factor; limits = [0, 1e4])
end


"""
    frohlich_mobility_k_space(v, w, α, ω, β; rₚ = 1, limits = [0, Inf])

Calculate the DC mobility in k-space for a Frohlich polaron system at finite temperature.

## Arguments
- `v`: a scalar value representing the optimal variational parameter.
- `w`: a scalar value representing the optimal variational parameter.
- `α`: a scalar value representing the coupling constant.
- `ω`: a scalar value representing the phonon frequency.
- `β`: a scalar value representing the inverse temperature.
- `rₚ`: an optional scalar value representing the characteristic polaron radius. Default value is 1.
- `limits`: an array of two scalar values representing the lower and upper limits of the integration range in k-space. Default is an infinite sphere.


## Returns
The DC mobility in k-space for the Frohlich polaron system at finite temperature.
"""
function frohlich_mobility_k_space(v, w, α, ω, β; rₚ = 1, limits = [0, Inf])
	 structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ω, β; rₚ = rₚ, limits = limits)
	 1 / imag(general_memory_function(structure_factor; limits = [0, Inf]))
end

