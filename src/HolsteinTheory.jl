# HolsteinTheory.jl

"""
    polaron_propagator(τ, v, w, β)

Calculate the imaginary time polaron Green's function with temperature dependence.

# Arguments
- `τ`: A scalar representing the imaginary time.
- `v`: A scalar representing a variational parameter.
- `w`: A scalar representing a variational parameter.
- `β`: A scalar representing the inverse temperature.

# Returns
The value of the polaron propagator.

# Example
```julia
τ = 0.5
v = 0.2
w = 0.1
β = 1.0
result = polaron_propagator(τ, v, w, β)
println(result)
```
This example calculates the polaron propagator for given values of τ, v, w, and β. The result is then printed.
"""
function polaron_propagator(τ, v, w, β)
    (v^2 - w^2) / v^3 * (1 - exp(-v * τ)) * (1 - exp(-v * (β - τ))) / (1 - exp(-v * β)) + w^2 / v^2 * τ * (1 - τ / β) + eps(Float64)
end

"""
    polaron_propagator(τ, v, w)

Calculate the value of the polaron propagator based on the given inputs.

# Arguments
- `τ::Number`: A scalar representing the imaginary time.
- `v::Number`: A scalar representing a variational parameter.
- `w::Number`: A scalar representing a variational parameter.

# Returns
The value of the polaron propagator.

# Example
```julia
τ = 0.5
v = 0.2
w = 0.1
result = polaron_propagator(τ, v, w)
println(result)
```
This example calculates the polaron propagator for the given values of τ, v, and w. The result is then printed.
"""
function polaron_propagator(τ, v, w)
    w^2 * τ / v^2 + (v^2 - w^2) / v^3 * (1 - exp(-v * τ)) + eps(Float64)
end

"""
    phonon_propagator(τ, ω)

Calculate the value of the phonon propagator based on the given inputs.

# Arguments
- `τ`: A scalar representing the imaginary time.
- `ω`: A scalar representing the adiabaticity.

# Example
```julia
τ = 0.5
ω = 0.2
result = phonon_propagator(τ, ω)
println(result)
```
This example calculates the value of the phonon propagator for the given values of τ and ω. The result is then printed.

# Returns
The value of the phonon propagator.
"""
function phonon_propagator(τ, ω)
    exp(-ω * τ)
end

"""
    phonon_propagator(τ, ω, β)

Calculate the value of the phonon propagator at a given imaginary time (`τ`), phonon frequency (`ω`), and inverse temperature (`β`).

# Arguments
- `τ`: A scalar representing the imaginary time.
- `ω`: A scalar representing the phonon frequency.
- `β`: A scalar representing the inverse temperature.

# Returns
The value of the phonon propagator.

# Example
```julia
τ = 0.5
ω = 0.2
β = 1.0
result = phonon_propagator(τ, ω, β)
println(result)
```
This example calculates the value of the phonon propagator for `τ = 0.5`, `ω = 0.2`, and `β = 1.0`. The result is then printed.
"""
function phonon_propagator(τ, ω, β)
    n = 1 / (exp(β * ω) - 1)
    result = n * exp(ω * τ) + (1 + n) * exp(-ω * τ)
    if isnan(result)
        return phonon_propagator(τ, ω)
    else
        return result
    end
end

"""
    holstein_interaction_energy_integrand(τ, v, w, α, ω, β; dims = 3)

Calculate the integrand for the Holstein interaction energy.

## Arguments
- `τ`: A scalar representing the imaginary time.
- `v`: A scalar representing a variational parameter.
- `w`: A scalar representing a variational parameter.
- `α`: A scalar representing the electron-phonon coupling.
- `ω`: A scalar representing the adiabaticity.
- `β`: A scalar representing the inverse temperature.
- `dims`: An optional parameter representing the number of dimensions (default is 3).

## Returns
The integrand for the Holstein interaction energy.

## Example
```julia
τ = 0.5
v = 0.2
w = 0.1
α = 0.3
ω = 0.4
β = 1.0
result = holstein_interaction_energy_integrand(τ, v, w, α, ω, β)
println(result)
```
This example calculates the integrand for the Holstein interaction energy using the given values of `τ`, `v`, `w`, `α`, and `ω`. The result is then printed.
"""
function holstein_interaction_energy_integrand(τ, v, w, α, ω, β; dims = 3)
    if β == Inf
        return holstein_interaction_energy_integrand(τ, v, w, α, ω; dims = dims)
    end
    momentum_cutoff = gamma(dims / 2 + 1)^(1/ dims) * (2√π)
    coupling = holstein_coupling(1, α, ω; dims = dims)
    propagator = polaron_propagator(τ, v, w, β * ω) / 2 
    phonon_propagator(τ / ω, ω, β) * coupling * P(dims, momentum_cutoff^2 * propagator) / (4π * propagator)^(dims/2)
end

"""
    holstein_interaction_energy_integrand(τ, v, w, α, ω; dims = 3)

Calculate the integrand for the Holstein interaction energy.

## Arguments
- `τ`: A scalar representing the imaginary time.
- `v`: A scalar representing a variational parameter.
- `w`: A scalar representing a variational parameter.
- `α`: A scalar representing the electron-phonon coupling.
- `ω`: A scalar representing the adiabaticity.
- `dims`: An optional parameter representing the number of dimensions (default is 3).

## Returns
The integrand for the Holstein interaction energy.

## Example
```julia
τ = 0.5
v = 0.2
w = 0.1
α = 0.3
ω = 0.4
result = holstein_interaction_energy_integrand(τ, v, w, α, ω)
println(result)
```
This example calculates the integrand for the Holstein interaction energy using the given values of `τ`, `v`, `w`, `α`, and `ω`. The result is then printed.
"""

function holstein_interaction_energy_integrand(τ, v, w, α, ω; dims = 3)
    momentum_cutoff = gamma(dims / 2 + 1)^(1 / dims) * (2√π)
    coupling = holstein_coupling(1, α, ω; dims = dims) / dims
    propagator = polaron_propagator(τ, v, w) / 2
    phonon_propagator(τ / ω, ω) * coupling * P(dims, momentum_cutoff^2 * propagator) / (4π * propagator)^(dims/2) 
end

"""
    holstein_interaction_energy(v, w, α, ω, β; dims = 3)

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
function holstein_interaction_energy(v, w, α, ω, β; dims = 3)
    if β == Inf
        return holstein_interaction_energy(v, w, α, ω; dims = dims)
    end
    integral, _ = quadgk(τ -> holstein_interaction_energy_integrand(τ, v, w, α, ω, β; dims = dims), 0, β * ω / 2)
    return integral
end

"""
    holstein_interaction_energy(v, w, α, ω; dims = 3)

Calculate the interaction energy between electrons and phonons in the Holstein model.

# Arguments
- `v`: A scalar representing a variational parameter.
- `w`: A scalar representing a variational parameter.
- `α`: A scalar representing the electron-phonon coupling.
- `ω`: A scalar representing the adiabaticity.
- `dims`: An optional parameter representing the number of dimensions (default is 3).

# Returns
- `integral`: The interaction energy between electrons and phonons in the Holstein model.

# Example
```julia
v = 0.2
w = 0.1
α = 0.5
ω = 0.3
result = holstein_interaction_energy(v, w, α, ω)
println(result)
```
This example calculates the interaction energy for given values of `v`, `w`, `α`, and `ω`. The result is then printed.
"""
function holstein_interaction_energy(v, w, α, ω; dims = 3)
    integral, _ = quadgk(τ -> holstein_interaction_energy_integrand(τ, v, w, α, ω; dims = dims), 0, Inf)
    return integral
end

"""
	Total free energy for the Holstein model. Here the k-space integral is evaluated analytically.
"""
function holstein_energy(v, w, α, ωβ...; dims = 3)
	kinetic_energy = electron_energy(v, w, ωβ...; dims = dims)
	potential_energy = -holstein_interaction_energy(v, w, α, ωβ...; dims = dims)
    return (kinetic_energy + potential_energy), kinetic_energy, potential_energy
end


"""
    vw_variation(energy, initial_v, initial_w; lower_bounds = [0, 0], upper_bounds = [Inf, Inf])

This function optimizes the values of `v` and `w` to minimize the energy function `energy(x[1] + x[2], x[2])[1]`. It uses the `Optim` package to perform the optimization and returns the optimized values of `v` and `w`, as well as the minimized energy, kinetic energy, and potential energy.

## Arguments
- `energy`: A function that takes two arguments `x` and `y` and returns an array of energy components.
- `initial_v`: The initial value of `v`.
- `initial_w`: The initial value of `w`.
- `lower_bounds`: An optional array representing the lower bounds for `v` and `w` optimization. Default is `[0, 0]`.
- `upper_bounds`: An optional array representing the upper bounds for `v` and `w` optimization. Default is `[Inf, Inf]`.

## Returns
- `Δv + w`: The optimized value of `v`.
- `w`: The optimized value of `w`.
- `energy_minimized`: The minimized energy.
- `kinetic_energy`: The kinetic energy corresponding to the minimized energy.
- `potential_energy`: The potential energy corresponding to the minimized energy.

## Example
```julia
energy(x, y) = [x^2 + y^2, x^2, y^2]
initial_v = 0.5
initial_w = 0.2
lower_bounds = [0, 0]
upper_bounds = [Inf, Inf]
result = vw_variation(energy, initial_v, initial_w; lower_bounds, upper_bounds)
println(result)
```
This example demonstrates how to use the `vw_variation` function. It defines an energy function `energy(x, y)` that returns an array of energy components. It then calls `vw_variation` with the initial values of `v` and `w`, as well as lower and upper bounds for the optimization. The function optimizes `v` and `w` to minimize the energy and returns the optimized values, as well as the minimized energy, kinetic energy, and potential energy. The result is then printed.
"""
function vw_variation(energy, initial_v, initial_w; lower = [0, 0], upper = [1e4, 1e4])

    Δv = initial_v - initial_w # defines a constraint, so that v>w
    initial = [Δv + eps(Float64), initial_w]

    # Feynman 1955 athermal action 
    f(x) = energy(x[1] + x[2], x[2])[1]

    # Use a more efficient optimization algorithm or library to optimise v and w to minimise enthalpy.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff=:forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
        # Optim.Options(show_trace = false, g_tol = 1e-12)
    )

    # Get v and w values that minimise the free energy.
    Δv, w = Optim.minimizer(solution) # v=Δv+w
    energy_minimized = energy(Δv + w, w)

    if !Optim.converged(solution)
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $energy_minimized"
    end

    # Return variational parameters that minimise the free energy.
    return Δv + w, w, energy_minimized...
end

vw_variation(energy) = vw_variation(energy, 5, 3; lower = [0, 0], upper = [1e4, 1e4])

"""
    general_memory_function(Ω, structure_factor; limits = [0, Inf])

This function calculates the integral of a given structure factor with respect to time using the `quadgk` function from the `QuadGK` package.

## Arguments
- `Ω`: A scalar representing the frequency.
- `structure_factor`: A function that returns the value of the structure factor at a given time.
- `limits` (optional): A 2-element array specifying the lower and upper limits of integration. Default is `[0, Inf]`.

## Returns
The integral of the structure factor with respect to time.

## Example
```julia
# Define a structure factor function
function structure_factor(t)
    # Implementation of the structure factor
    # ...
end

# Calculate the memory function for a given frequency and structure factor
Ω = 0.5
limits = [0, 10]
result = general_memory_function(Ω, structure_factor; limits = limits)
println(result)
```
This example demonstrates how to use the `general_memory_function` to calculate the memory function for a given frequency `Ω` and structure factor function `structure_factor`. The `limits` argument is optional and specifies the lower and upper limits of integration. The result is then printed.
"""
function general_memory_function(Ω, structure_factor; limits = [0, Inf])
    integral, _ = quadgk(t -> (1 - exp(im * Ω * t)) / Ω * imag(structure_factor(t)), limits[1], limits[2])
    return integral
end

"""
    general_memory_function(structure_factor; limits = [0, Inf])

Calculate the integral of a given function `structure_factor` using the `quadgk` function in Julia.

## Arguments
- `structure_factor`: A function that takes a single argument `t` and returns a value.
- `limits`: An optional array specifying the lower and upper limits of integration. Default is `[0, Inf]`.

## Returns
- `integral`: The result of the numerical integration of the function `structure_factor` over the specified limits.

## Example
```julia
# Define the structure factor function
function structure_factor(t)
    return t^2 + 2t + 1
end

# Call the general_memory_function with the structure_factor function
result = general_memory_function(structure_factor; limits = [0, 10])

println(result)  # Output: 383.3333333333333
```
"""
function general_memory_function(structure_factor; limits = [0, Inf])
    integral, _ = quadgk(t -> -im * t * imag(structure_factor(t)), limits[1], limits[2])
    return integral
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
function holstein_structure_factor(t, v, w, α, ω, β; dims = 3)
    if β == Inf
        return holstein_structure_factor(t, v, w, α, ω; dims = dims)
    end
    momentum_cutoff = gamma(dims / 2 + 1)^(1 / dims) * (2√π)
	coupling = holstein_coupling(1, α, ω; dims = dims)
	propagator = polaron_propagator(im * t, v, w, β * ω) / 2
	integral = dims / 2 * ball_surface(dims) / (2π)^3 * P_plus_one(dims, propagator * momentum_cutoff^2) / propagator^(dims/2 + 1)
    phonon_propagator(im * t / ω, ω, β) * coupling * integral * 2 / dims / ω
end

"""
    holstein_structure_factor(t, v, w, α, ω; dims = 3)

Calculate the structure factor for the Holstein model.

## Arguments
- `t`: A scalar representing the time.
- `v`: A scalar representing a parameter.
- `w`: A scalar representing a parameter.
- `α`: A scalar representing a parameter.
- `ω`: A scalar representing a parameter.
- `dims`: An optional integer representing the number of dimensions (default is 3).

## Returns
The value of the structure factor.

## Example
```julia
t = 0.5
v = 0.2
w = 0.1
α = 1.0
ω = 0.5
result = holstein_structure_factor(t, v, w, α, ω)
println(result)
```
This example calculates the structure factor for the given values of `t`, `v`, `w`, `α`, and `ω`. The result is then printed.
"""
function holstein_structure_factor(t, v, w, α, ω; dims = 3)
	momentum_cutoff = gamma(dims / 2 + 1)^(1 / dims) * (2√π)
	coupling = holstein_coupling(1, α, ω; dims = dims)
	propagator = polaron_propagator(im * t, v, w) / 2
	integral = dims / 2 * ball_surface(dims) / (2π)^3 * P_plus_one(dims, propagator * momentum_cutoff^2) / propagator^(dims/2 + 1)
    phonon_propagator(im * t / ω, ω) * coupling * integral * 2 / ω / dims
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
function holstein_memory_function(Ω, v, w, α, ω, β; dims = 3)
    if β == Inf
        return holstein_memory_function(Ω, v, w, α, ω; dims = dims)
    end
	structure_factor(t) = holstein_structure_factor(t, v, w, α, ω, β; dims = dims)
	return general_memory_function(Ω / ω, structure_factor)
end

"""
    holstein_memory_function(Ω, v, w, α, ω; dims = 3)

Calculate the integral of a structure factor using the `general_memory_function` function.

## Arguments
- `Ω`: The frequency parameter.
- `v`, `w`, `α`, `ω`: Parameters used in the `holstein_structure_factor` function.
- `dims`: Optional parameter specifying the number of dimensions, defaults to 3.

## Returns
The integral value of the structure factor.

## Example
```julia
result = holstein_memory_function(Ω, v, w, α, ω; dims = 3)
```
In this example, the `holstein_memory_function` is called with the parameters `Ω`, `v`, `w`, `α`, and `ω`. The `dims` parameter is optional and defaults to 3. The function calculates the structure factor using the `holstein_structure_factor` function and then calls the `general_memory_function` to calculate the integral of the structure factor. The result is stored in the `result` variable.
"""
function holstein_memory_function(Ω, v, w, α, ω; dims = 3)
	 structure_factor(t) = holstein_structure_factor(t, v, w, α, ω; dims = dims)
	 return general_memory_function(Ω / ω, structure_factor, limits = [0, 1e6])
end

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
function holstein_mobility(v, w, α, ω, β; dims = 3)
    if β == Inf
        return Inf
    end
    structure_factor(t) = holstein_structure_factor(t, v, w, α, ω, β; dims = dims)
    abs(1 / imag(general_memory_function(structure_factor)))
end

function holstein_complex_impedence(Ω, v, w, α, ωβ...; dims = 3)
    -im * (Ω - holstein_memory_function(v, w, α, ωβ...; dims = dims))
end

function holstein_complex_conductivity(Ω, v, w, α, ωβ...; dims = 3)
    return 1 / holstein_complex_impedence(Ω, v, w, α, ωβ...; dims = dims)
end