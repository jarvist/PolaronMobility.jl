# HolsteinTheory.jl
# Everything specifically associated with calculating properties of the Holstein polaron model.

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
function holstein_coupling(k, α, J, ω; dims = 3)
    α * ω * J * dims
end

holstein_coupling(k, α::Vector, ω::Vector; dims = 3) = sum(holstein_coupling(k, α[j], ω[j]; dims = dims) for j in eachindex(α))

"""
    holstein_interaction_energy(v, w, α, ωβ...; dims = 3)

Electron-phonon interaction energy for the Holstein mode at finite temperature. Here the k-space integral is evaluated analytically.

# Arguments
- `v`: A scalar representing a variational parameter.
- `w`: A scalar representing a variational parameter.
- `α`: A scalar representing the electron-phonon coupling.
- `ω`: A scalar representing the adiabaticity.
- `β`: A scalar representing the inverse temperature.
- `dims`: An optional argument representing the number of dimensions. Default is 3.

# Returns
- `integral`: The electron-phonon interaction energy for the Holstein mode at finite temperature.

# Example
```julia
v = 0.2
w = 0.1
α = 0.5
ω = 0.3
β = 1.0
result = holstein_interaction_energy(v, w, α, ω, β)
println(result)
```
This example calculates the electron-phonon interaction energy for given values of `v`, `w`, `α`, `ω`, and `β`. The result is then printed.
"""
function holstein_interaction_energy(v, w, α, J, ωβ...; a = 1, dims = 3)
    momentum_cutoff = gamma(dims / 2 + 1)^(1 / dims) * (2√π)
    coupling = holstein_coupling(1, α, J, ωβ[1]; dims = dims) 
    propagator(τ) = length(ωβ) == 1 ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, ωβ[2])
    integrand(τ) = phonon_propagator(τ, ωβ...) * P(dims, momentum_cutoff^2 * propagator(τ)) / (propagator(τ))^(dims/2)
    upper_limit = length(ωβ) == 1 ? Inf : ωβ[2] / 2
    integral, _ = quadgk(τ -> integrand(τ), 0, upper_limit)
    return integral * coupling / (4π)^(dims / 2) 
end

holstein_interaction_energy(v, w, α::Vector, ω::Vector, β; dims = 3) = sum(holstein_interaction_energy(v, w, α[j], ω[j], β; dims = dims) for j in eachindex(α))
holstein_interaction_energy(v, w, α::Vector, ω::Vector; dims = 3) = sum(holstein_interaction_energy(v, w, α[j], ω[j]; dims = dims) for j in eachindex(α))

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

## Returns
A scalar value representing the Holstein polaron interaction energy in k-space at finite temperature.
"""

function holstein_interaction_energy_k_space(v, w, α, J, ωβ...; a = 1, dims = 3)
    momentum_cutoff = gamma(dims / 2 + 1)^(1 / dims) * (2√π) / a
    polaron_radius = a * sqrt(J / ωβ[1])
    coupling(k) = holstein_coupling(k, α, J, ωβ[1]; dims = dims)
    propagator(τ) = length(ωβ) == 1 ? polaron_propagator(τ, v, w) * ωβ[1] : polaron_propagator(τ, v, w, ωβ[2]) * ωβ[1]
    integrand(τ) = phonon_propagator(τ, ωβ...) * spherical_k_integral(coupling, propagator(τ); dims = dims, limits = [0, momentum_cutoff], radius = polaron_radius)
    upper_limit = length(ωβ) == 1 ? Inf : ωβ[2] / 2
    integral, _ = quadgk(τ -> integrand(τ), 0, upper_limit)
	return integral * a^3
end

holstein_interaction_energy_k_space(v, w, α::Vector, ω::Vector, β; dims = 3) = sum(holstein_interaction_energy_k_space(v, w, α[j], J, ω[j], β; dims = dims) for j in eachindex(α))
holstein_interaction_energy_k_space(v, w, α::Vector, ω::Vector; dims = 3) = sum(holstein_interaction_energy_k_space(v, w, α[j], J, ω[j]; dims = dims) for j in eachindex(α))

"""
	Total free energy for the Holstein model. Here the k-space integral is evaluated analytically.
"""
function holstein_energy(v, w, α, J, ωβ...; dims = 3)
	A, C = length(ωβ) == 1 ? trial_energy(v, w; dims = 1) : trial_energy(v, w, ωβ[2]; dims = 1)
	B = holstein_interaction_energy(v, w, α, J, ωβ...; dims = dims)
    return -(A + B + C), A, B, C
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
function holstein_energy_k_space(v, w, α, J, ωβ...; a = 1, dims = 3)
    A, C = length(ωβ) == 1 ? trial_energy(v, w; dims = 1) : trial_energy(v, w, ωβ[2]; dims = 1)
    B = holstein_interaction_energy_k_space(v, w, α, J, ωβ...; a = 1, dims = dims)
    return -(A + B + C), A, B, C
end

"""
    holstein_structure_factor(t, v, w, α, ω, β; dims = 3)

Calculate the structure factor for the Holstein model.

## Arguments
- `t`: A scalar representing the time.
- `v`: A scalar representing a parameter.
- `w`: A scalar representing a parameter.
- `α`: A scalar representing a coupling constant.
- `ω`: A scalar representing a coupling constant.
- `β`: A scalar representing the inverse temperature.
- `dims`: An optional argument representing the number of dimensions. Default is 3.

## Returns
The value of the structure factor for the given inputs.

## Example
```julia
t = 0.5
v = 0.2
w = 0.1
α = 0.3
ω = 0.4
β = 1.0
result = holstein_structure_factor(t, v, w, α, ω, β)
println(result)
```
This example demonstrates how to use the `holstein_structure_factor` function to calculate the structure factor for a given set of parameters. The function is called with the values of `t`, `v`, `w`, `α`, `ω`, and `β` as arguments, and the result is then printed.
"""
function holstein_structure_factor(t, v, w, α, J, ωβ...; dims = 3)
    momentum_cutoff = gamma(dims / 2 + 1)^(1 / dims) * (2√π)
	coupling = holstein_coupling(1, α, J, ωβ[1]; dims = dims)
	propagator = length(ωβ) == 1 ? polaron_propagator(im * t, v, w) : polaron_propagator(im * t, v, w, ωβ[2])
	integral = P(dims + 2, propagator * momentum_cutoff^2 / 2) * (propagator)^(-dims/2 - 1) * 2^(1 - dims / 2) * π^(-dims / 2) * dims / 2
    return 2 / dims * coupling * integral * phonon_propagator(im * t, ωβ...) 
end

"""
    holstein_structure_factor_k_space(t, v, w, α, ωβ...; dims = 3)

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
function holstein_structure_factor_k_space(t, v, w, α, J, ωβ...; a = 1, dims = 3)
    momentum_cutoff = gamma(dims / 2 + 1)^(1/ dims) * (2√π) / a
    polaron_radius = a * sqrt(J / ωβ[1])
	coupling(k) = holstein_coupling(1, α, J, ωβ[1]; dims = dims) * k^2 * ωβ[1]
	propagator = length(ωβ) == 1 ? polaron_propagator(im * t, v, w) * ωβ[1] : polaron_propagator(im * t, v, w, ωβ[2]) * ωβ[1]
	integral = spherical_k_integral(coupling, propagator; dims = dims, limits = [0, momentum_cutoff], radius = polaron_radius)
	return 2 / dims * integral * phonon_propagator(im * t, ωβ...) * a^dims
end

"""
    holstein_memory_function(Ω, v, w, α, ω, β; dims = 3)

Calculate the memory function using the `general_memory_function` function.

## Arguments
- `Ω`: a parameter representing a frequency
- `v`, `w`, `α`, `ω`, `β`: parameters used to calculate the structure factor
- `dims`: an optional parameter representing the dimensions (default value is 3)

## Returns
The result of the `general_memory_function` function, which represents the calculated memory function.

## Example
```julia
result = holstein_memory_function(Ω, v, w, α, ω, β; dims = 3)
```
In this example, the `holstein_memory_function` is called with the input parameters `Ω`, `v`, `w`, `α`, `ω`, and `β`, and the optional parameter `dims` set to 3. The function calculates the memory function using the `general_memory_function` function and returns the result.
"""
function holstein_memory_function(Ω, v, w, α, J, ωβ...; dims = 3)
	structure_factor(t) = holstein_structure_factor(t, v, w, α, J, ωβ...; dims = dims)
	return polaron_memory_function(Ω, structure_factor)
end

holstein_memory_function(Ω, v, w, α::Vector, J, ω::Vector, β; dims = 3) = sum(holstein_memory_function(Ω, v, w, α[j], J, ω[j], β; dims = dims) for j in eachindex(α))
holstein_memory_function(Ω, v, w, α::Vector, J,  ω::Vector; dims = 3) = sum(holstein_memory_function(Ω, v, w, α[j], J, ω[j]; dims = dims) for j in eachindex(α))

"""
    holstein_memory_function_k_space(Ω, v, w, α, ωβ...; dims = 3)

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
function holstein_memory_function_k_space(Ω, v, w, α, J, ωβ...; a = 1, dims = 3)
	 structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, J, ωβ...; a = a, dims = dims)
	 return polaron_memory_function(Ω, structure_factor)
end

# holstein_memory_function_k_space(Ω, v, w, α::Vector, J, ω::Vector, β; dims = 3) = sum(holstein_memory_function_k_space(Ω, v, w, α[j], J, ω[j], β; dims = dims) for j in eachindex(α))
# holstein_memory_function_k_space(Ω, v, w, α::Vector, J, ω::Vector; dims = 3) = sum(holstein_memory_function_k_space(Ω, v, w, α[j], J, ω[j]; dims = dims) for j in eachindex(α))

"""
    holstein_mobility(v, w, α, ω, β; dims = 3)

Calculate the mobility using the `general_memory_function` and `holstein_structure_factor` functions.

## Arguments
- `v`: Parameter used in the `holstein_structure_factor` function.
- `w`: Parameter used in the `holstein_structure_factor` function.
- `α`: Parameter used in the `holstein_structure_factor` function.
- `ω`: Parameter used in the `holstein_structure_factor` function.
- `β`: Parameter used in the `holstein_structure_factor` function.
- `dims`: The dimensionality of the system. Default is 3.

## Returns
The calculated mobility.

## Example
```julia
v = 1.0
w = 2.0
α = 0.5
ω = 1.0
β = 0.2
dims = 3

result = holstein_mobility(v, w, α, ω, β, dims=dims)
println(result)
```
This code calculates the mobility using the given parameters and prints the result.
"""
function inverse_holstein_mobility(v, w, α, J, ω, β; dims = 3)
    structure_factor(t) = holstein_structure_factor(t, v, w, α, J, ω, β; dims = dims)
    return abs(imag(polaron_memory_function(structure_factor)))
end

inverse_holstein_mobility(v, w, α::Vector, J, ω::Vector, β; dims = 3) = sum(inverse_holstein_mobility(v, w, α[j], J,  ω[j], β; dims = dims) for j in eachindex(α))

"""
    polaron_mobility(v, w, α, ω, β)

The polaron mobility.

See also [`inverse_polaron_mobility`](@ref)
"""
holstein_mobility(v, w, α, J, ω, β; dims = 3) = 1 / inverse_holstein_mobility(v, w, α, J, ω, β; dims = dims)
holstein_mobility(v, w, α, J, ω; dims = 3) = reduce_array(repeat([Inf], length(α)))


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
function inverse_holstein_mobility_k_space(v, w, α, J, ω, β; a = 1, dims = 3)
    structure_factor(t) = holstein_structure_factor_k_space(t, v, w, α, J, ω, β; a = a, dims = dims)
    return abs(imag(polaron_memory_function(structure_factor)))
end

# inverse_holstein_mobility_k_space(v, w, α::Vector, J, ω::Vector, β; dims = 3) = sum(inverse_holstein_mobility_k_space(v, w, α[j], J, ω[j], β; dims = dims) for j in eachindex(α))

holstein_mobility_k_space(v, w, α, J, ω, β; a = 1, dims = 3) = 1 / inverse_holstein_mobility_k_space(v, w, α, J, ω, β; a = a, dims = dims)

holstein_mobility_k_space(v, w, α, J, ω; a = 1, dims = 3) = reduce_array(repeat([Inf], length(α)))

function holstein_complex_impedence(Ω, v, w, α, J, ωβ...; dims = 3)
    -im * (Ω - holstein_memory_function(v, w, α, J, ωβ...; dims = dims))
end

function holstein_complex_conductivity(Ω, v, w, α, J, ωβ...; dims = 3)
    return 1 / holstein_complex_impedence(Ω, v, w, α, J, ωβ...; dims = dims)
end

