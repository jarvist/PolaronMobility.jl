# PolaronFunctions.jl

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
    result = n * exp(τ * ω) + (1 + n) * exp(-τ * ω)
    if isnan(result)
        return phonon_propagator(τ, ω)
    else
        return result
    end
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
    exp(-τ * ω)
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
function vw_variation(energy, initial_v, initial_w; lower = [0, 0], upper = [Inf, Inf], warn = false)

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

    if !Optim.converged(solution) && warn
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $energy_minimized"
    end

    # Return variational parameters that minimise the free energy.
    return Δv + w, w, energy_minimized...
end

vw_variation(energy) = vw_variation(energy, 3.1, 3; lower = [0, 0], upper = [Inf, Inf])

function vw_variation(energy, initial_v::Vector, initial_w::Vector; lower = [0, 0], upper = [Inf, Inf], warn = false)

  if length(initial_v) != length(initial_w)
      return error("The number of variational parameters v & w must be equal.")
  end

  N_params = length(initial_v)

  Δv = initial_v .- initial_w
  initial = vcat(Δv .+ repeat([eps(Float64)], N_params), initial_w)

  # Limits of the optimisation.
  lower = vcat(fill(lower, N_params)...)
  upper = vcat(fill(upper, N_params)...)

  # The multiple phonon mode free energy function to minimise.
  f(x) = energy([x[2*n-1] for n in 1:N_params] .+ [x[2*n] for n in 1:N_params], [x[2*n] for n in 1:N_params])[1]

  # Use a more efficient optimization algorithm or library to optimise v and w to minimise enthalpy.
  solution = Optim.optimize(
    Optim.OnceDifferentiable(f, initial; autodiff=:forward),
    lower,
    upper,
    initial,
    Fminbox(BFGS()),
    # Optim.Options(show_trace = false, g_tol = 1e-12)
  )

  # Extract the v and w parameters that minimised the free energy.
  var_params = Optim.minimizer(solution)

  # Separate the v and w parameters into one-dimensional arrays (vectors).
  Δv = [var_params[2*n-1] for n in 1:N_params]
  w = [var_params[2*n] for n in 1:N_params]
  energy_minimized = energy(Δv .+ w, w)

  if !Optim.converged(solution) && warn
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $energy_minimized"
  end

  # Return the variational parameters that minimised the free energy.
  return Δv .+ w, w, energy_minimized...
end

"""
    polaron_memory_function(Ω, structure_factor; limits = [0, Inf])

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
function polaron_memory_function(Ω, structure_factor; limits = [0, Inf])
    if iszero(Ω)
      return polaron_memory_function(structure_factor; limits = limits)
    end
    integral, _ = quadgk(t -> (1 - exp(im * Ω * t)) / Ω * imag(structure_factor(t)), 0, Inf)
    return integral
end


"""
    general_memory_function(structure_factor; limits = [0, Inf])

Calculate the integral of a given function `structure_factor` using the `quadgk` function in Julia.

## Arguments
- `structure_factor`: A function that takes a single argument `t` and returns a value.
- `limits`: An optional array specifying the lower and upper limits of integration. Default is `[0, Inf]`.

## Returnstyp
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
function polaron_memory_function(structure_factor; limits = [0, Inf])
    integral, _ = quadgk(t -> -im * t * imag(structure_factor(t)), eps(Float64), Inf)
    return integral
end

# Utility / general use functions

# QOL function for removing singleton dimensions.
reduce_array(a) = length(a) == 1 ? only(a) : Array(dropdims(a, dims=tuple(findall(size(a) .== 1)...)))

# Regularized lower gamma function for half dimensions n/2.
function P(n, z)
  if n == 1
    erf(sqrt(z))
  elseif n == 2
    1 - exp(-z)
  elseif n > 2
    P_plus_one(n, z)
  end
end

# Recursion relations for higher half dimensions n/2 + m where m ∈ Integer.
function P_plus_one(n, z)
  if n<3
    return P(n, z)
  else
    return P(n-2, z) - exp(-z) * z^(n/2-1) / gamma(n/2) 
  end
end

# Surface area of an n-dimensional sphere.
function ball_surface(n)
  2 * π^(n/2) / gamma(n/2)
end

# n-dimensional spherical integral for polaron theory. 
function spherical_k_integral(coupling, propagator; dims = 3, radius = 1/sqrt(2), limits = [0, Inf])
  integrand(k) = k^(dims-1) * coupling(k) * exp(-k^2 * radius^2 * propagator)
  integral, _ = quadgk(k -> integrand(k), limits[1], limits[2])
  return integral * ball_surface(dims) / (2π)^dims
end