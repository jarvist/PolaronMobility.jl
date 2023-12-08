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
function cartesian_k_integrand(k, coupling, propagator)
    coupling(k) * exp(-k^2 * propagator / 2)
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
function spherical_k_integrand(k, coupling, propagator; dims = 3)
    2 * π^(dims/2) / gamma(dims/2) * k^(dims-1) * coupling(k) * exp(-k^2 * propagator / 2)
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
function cartesian_k_integral(coupling, propagator; limits = [-π, π])
    integral, _ = quadgk(k -> cartesian_k_integrand(k, coupling, propagator), limits[1], limits[2])
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
function spherical_k_integral(coupling, propagator; dims = 3, limits = [0, Inf])
    integral, _ = quadgk(k -> spherical_k_integrand(k, coupling, propagator; dims = dims), limits[1], limits[2])
    integral / (2π)^dims
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
function frohlich_coupling(k, α, ω; mb = 1, dims = 3)
    r_p = 1 / sqrt(2)
    ω^2 * α * r_p * gamma((dims - 1) / 2) * (2√π / k)^(dims - 1)
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
function holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; dims = 3)
    momentum_cutoff = gamma(dims / 2 + 1)^(1 / dims) * (2√π)
    coupling(k) = holstein_coupling(k, α, ω; dims = dims)
    propagator = polaron_propagator(τ, v, w, β) 
    phonon_propagator(τ, ω, β) * spherical_k_integral(coupling, propagator; dims = dims, limits = [0, momentum_cutoff]) 
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
function holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω; dims = 3)
    momentum_cutoff = gamma(dims / 2 + 1)^(1/ dims) * (2√π)
    coupling(k) = holstein_coupling(k, α, ω; dims = dims) 
    propagator = polaron_propagator(τ, v, w) 
    phonon_propagator(τ, ω) * spherical_k_integral(coupling, propagator; dims = dims, limits = [0, momentum_cutoff]) 
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
function frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; limits = [0, Inf], dims = 3)
    coupling(k) = frohlich_coupling(k, α, ω; dims = dims)
    propagator = polaron_propagator(τ, v, w, β)
    phonon_propagator(τ, ω, β) * spherical_k_integral(coupling, propagator; limits = limits, dims = dims)
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
function frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω; limits = [0, Inf], dims = 3)
	coupling(k) = frohlich_coupling(k, α, ω; dims = dims)
	propagator = polaron_propagator(τ, v, w) 
	phonon_propagator(τ, ω) * spherical_k_integral(coupling, propagator; limits = limits, dims = dims)
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
function holstein_interaction_energy_k_space(v, w, α, ω, β; dims = 3)
    interaction_energy, _ = quadgk(τ -> holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; dims = dims), 0, β/2)
    return interaction_energy
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
function holstein_interaction_energy_k_space(v, w, α, ω; dims = 3)
	interaction_energy, _ = quadgk(τ -> holstein_interaction_energy_integrand_k_space(τ, v, w, α, ω; dims = dims), 0, Inf)
	return interaction_energy 
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
function frohlich_interaction_energy_k_space(v, w, α, ω, β; limits = [0, Inf], dims = 3)
	interaction_energy, _ = quadgk(τ -> frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω, β; limits = limits, dims = dims), 0, β/2)
	return interaction_energy
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
function frohlich_interaction_energy_k_space(v, w, α, ω; limits = [0, Inf], dims = 3)
	interaction_energy, _ = quadgk(τ -> frohlich_interaction_energy_integrand_k_space(τ, v, w, α, ω; limits = limits, dims = dims), 0, Inf)
	return interaction_energy
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
function electron_energy(v, w, ω, β; dims = 3)
    if β == Inf
        return electron_energy(v, w, ω; dims = dims)
    end
    -(A(v, w, ω, β) + C(v, w, ω, β)) * dims / 3
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
function electron_energy(v, w, ω; dims = 3)
	(v - w)^2 / (4 * v) * ω * dims
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
function holstein_energy_k_space(v, w, α, ωβ...; dims = 3)
    kinetic_energy = electron_energy(v, w, ωβ...; dims = dims) 
    interaction_energy = -holstein_interaction_energy_k_space(v, w, α, ωβ...; dims = dims) 
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
function frohlich_energy_k_space(v, w, α, ωβ...; limits = [0, Inf], dims = 3)
	kinetic_energy = electron_energy(v, w, ωβ...; dims = dims)
	interaction_energy = -frohlich_interaction_energy_k_space(v, w, α, ωβ...; limits = limits, dims = dims)
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
function holstein_structure_factor_k_space(t, v, w, α, ω, β; dims = 3)
	
    momentum_cutoff = gamma(dims / 2 + 1)^(1/ dims) * (2√π)

	coupling(k) = holstein_coupling(k, α, ω; dims = dims) * k^2

	propagator = polaron_propagator(im * t, v, w, β)
	
	integral = spherical_k_integral(coupling_one, propagator; limits = [0, momentum_cutoff])
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω, β)

	prefactor * integral 
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
function holstein_structure_factor_k_space(t, v, w, α, ω; dims = 3)

    momentum_cutoff = gamma(dims / 2 + 1)^(1/ dims) * (2√π)
	
	coupling(k) = holstein_coupling(k, α, ω; dims = dims) * k^2

	propagator = polaron_propagator(im * t, v, w)
	
	integral = spherical_k_integral(coupling, propagator; limits = [0, momentum_cutoff])
	
	prefactor = 2 / 3 * phonon_propagator(im * t, ω)

	prefactor * integral
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
function frohlich_structure_factor_k_space(t, v, w, α, ω, β; limits = [0, Inf], dims = 3)
	
	coupling(k) = frohlich_coupling(k, α, ω; dims = dims) * k^2

	propagator = polaron_propagator(im * t, v, w, β)
	
	integral = spherical_k_integral(coupling, propagator; limits = limits, dims = dims)
	
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
function frohlich_structure_factor_k_space(t, v, w, α, ω; limits = [0, Inf], dims = 3)
	
	coupling(k) = frohlich_coupling(k, α, ω; dims = dims) * k^2

	propagator = polaron_propagator(im * t, v, w)
	
	integral = spherical_k_integral(coupling, propagator; limits = limits, dims = dims)
	
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
function holstein_memory_function_k_space(Ω, v, w, α, ω, β; dims = 3)
	 structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, ω, β; dims = dims)
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
function holstein_memory_function_k_space(Ω, v, w, α, ω; dims = 3)
	 structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, ω; dims = dims)
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
function holstein_mobility_k_space(v, w, α, ω, β; dims = 3)
    structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, ω, β; dims = dims)
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
function frohlich_memory_function_k_space(Ω, v, w, α, ω, β; limits = [0, Inf])
	 structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ω, β; limits = limits)
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
function frohlich_memory_function_k_space(Ω, v, w, α, ω; limits = [0, Inf])
	 structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ω; limits = limits)
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
function frohlich_mobility_k_space(v, w, α, ω, β; limits = [0, Inf])
	 structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ω, β; limits = limits)
	 1 / imag(general_memory_function(structure_factor; limits = [0, 1e5]))
end

