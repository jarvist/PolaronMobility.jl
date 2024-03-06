# TrialPolaron.jl
# This file contains all functions associated with calculating properties of the trial polaron system.

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
    c = (v^2 - w^2) / v^3
    result = c * (phonon_propagator(0, v, β) - phonon_propagator(τ, v, β)) + (1 - c * v) * τ * (1 - τ / β)
    return result
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
    c = (v^2 - w^2) / v^3
    result = (1 - v * c) * τ + c * (phonon_propagator(0, v) - phonon_propagator(τ, v))
    return result
end

"""
    D(τ, v, w, β)

Calculates the recoil function (a generalisation of D(u) in Eqn. (35c) in FHIP 1962) that approximates the exact influence (recoil effects) of the phonon bath on the electron with the influence of the fictitious masses attached by springs to the electron. It appears in the exponent of the intermediate scattering function.

# Arguments
- `τ::Float64`: is the imaginary time variable.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `β::Union{Float64, Vector{Float64}}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.

See FHIP 1962: https://doi.org/10.1103/PhysRev.127.1004.
"""
function polaron_propagator(τ, v::Vector, w::Vector, β)
    return τ * (1 - τ / β) + sum((h_i(i, v, w) / v[i]^2) * ((1 + exp(-v[i] * β) - exp(-v[i] * τ) - exp(v[i] * (τ - β))) / (v[i] * (1 - exp(-v[i] * β))) - τ * (1 - τ / β)) for i in eachindex(v))
end

"""
    D(τ, v, w)

Calculates the recoil function at zero-temperature.

# Arguments
- `τ::Float64`: is the imaginary time variable.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
"""
function polaron_propagator(τ, v::Vector, w::Vector)
    return τ + sum((h_i(i, v, w) / v[i]^2) * ((1 - exp(-v[i] * τ)) / v[i] - τ) for i in eachindex(v))
end

# Hellwarth et al. 1999 PRB - Part IV; T-dep of the Feynman variation parameter

# In Julia we have 'Multiple dispatch', so let's just construct the free
# energies (temperature-dependent) with the same name as before, but with the thermodynamic beta where required.

# Define Osaka's free-energies (Hellwarth1999 version) as Julia functions.

"""
    A(v, w, ω, β)

Hellwarth's A expression from Eqn. (62b) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
A(v, w, β) = β == Inf ? A(v, w) : 3 / β * (log(v / w) - (v - w) * β / 2 - log(1 - exp(-v * β)) + log(1 - exp(-w * β))) 

A(v, w) = -3 * (v - w) / 2


"""
    C(v, w, ω, β)

Hellwarth's C expression from Eqn. (62e) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
C(v, w, β) = 3 / 4 * (v^2 - w^2) / v * (coth(v * β / 2) - 2 / (v * β)) 

C(v, w) = (3 / (4 * v)) * (v^2 - w^2)

# Extending the Feynman theory to multiple phonon branches

#  Extending the Feynman theory to multiple variational parameters
# Multiple Parameter Polaron Free Energy
# Calculate the polaron free energy, generalised from Osaka's expression to the case where multiple phonon modes are present in the material.

"""
    κ_i(i, v, w)

Calculates the spring-constant coupling the electron to the 'ith' fictitious mass that approximates the exact electron-phonon interaction with a harmonic coupling to a massive fictitious particle. 

Required for calculating the polaron free energy.

Note: Not to be confused with the number of physical phonon branches; many phonon branches could be approximated with one or two etc. fictitious masses for example. The number of fictitious mass does not necessarily need to match the number of phonon branches.

# Arguments
- `i::Integer`: enumerates the current fictitious mass.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
"""
function κ_i(i, v::Vector, w::Vector)
    κ = v[i]^2 - w[i]^2
    κ *= prod(j != i ? (v[j]^2 - w[i]^2) / (w[j]^2 - w[i]^2) : 1 for j in eachindex(v))
    return κ
end

"""
    h_i(i, v, w)

Calculates the normal-mode (the eigenmodes) frequency of the coupling between the electron and the `ith' fictitious mass that approximates the exact electron-phonon interaction with a harmonic coupling to a massive fictitious particle. 

Required for calculating the polaron free energy.

Note: Not to be confused with the number of physical phonon branches; many phonon branches could be approximated with one or two etc. fictitious masses for example. The number of fictitious mass does not necessarily need to match the number of phonon branches.

# Arguments
- `i::Integer`: enumerates the current fictitious mass.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
"""
function h_i(i, v::Vector, w::Vector)
    h = v[i]^2 - w[i]^2
    h *= prod(j != i ? (w[j]^2 - v[i]^2) / (v[j]^2 - v[i]^2) : 1 for j in eachindex(v))
    return h
end

"""
    C_ij(i, j, v, w)

Calculates the element to the coupling matrix C_ij (a generalisation of Feynman's `C` coupling variational parameter in Feynman 1955) between the electron and the `ith' and `jth' fictitious masses that approximates the exact electron-phonon interaction with a harmonic coupling to a massive fictitious particle. 

Required for calculating the polaron free energy.

Note: Not to be confused with the number of physical phonon branches; many phonon branches could be approximated with one or two etc. fictitious masses for example. The number of fictitious mass does not necessarily need to match the number of phonon branches.

# Arguments
- `i::Integer, j::Integer`: enumerate the current fictitious masses under focus (also the index of the element in the coupling matrix C)
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
function C_ij(i, j, v::Vector, w::Vector)
    return w[i] * κ_i(i, v, w) * h_i(j, v, w) / (4 * (v[j]^2 - w[i]^2))
end

"""
    C(v, w, β)

Generalisation of the C function from Eqn. (62e) in Hellwarth et al. 1999. This is the expected value of the trial action <S_0> taken w.r.t trial action.

Required for calculating the polaron free energy.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `β::Union{Float64, Vector{Float64}}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.


See Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.
"""
function C(v::Vector, w::Vector, β)
    # Sum over the contributions from each fictitious mass.
    s = sum(C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2) - 2 / (β * v[j])) for i in eachindex(v), j in eachindex(w))

    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    return 3 * s
end

"""
    C(v, w)

Calculates `C(v, w, β)` but at zero-temperature, `β = Inf`.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
"""
function C(v::Vector, w::Vector)
    # Sum over the contributions from each fictitious mass.
    s = sum(C_ij(i, j, v, w) / (v[j] * w[i]) for i in eachindex(v), j in eachindex(w))
    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    return 3 * s
end

"""
    A(v, w, β)

Generalisation of the A function from Eqn. (62b) in Hellwarth et al. 1999. This is the Helmholtz free energy of the trial model.

Required for calculating the polaron free energy.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `β::Union{Float64, Vector{Float64}}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.

See Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.
"""
function A(v::Vector, w::Vector, β)
    # Sum over the contributions from each fictitious mass.
    s = -log(2π * β) / 2 + sum(v[i] == w[i] ? 0 :
                               log(v[i]) - log(w[i]) - β / 2 * (v[i] - w[i]) - log(1 - exp(-v[i] * β)) + log(1 - exp(-w[i] * β))
                               for i in eachindex(v))
    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    3 / β * s
end

"""
    A(v, w, n)

Calculates `A(v, w, β)` but at zero-temperature, `β = Inf`.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
"""
function A(v::Vector, w::Vector)
    s = sum(v .- w)
    return -3 * s / 2 
end

"""
    trial_energy(v, w, ω, β)

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
function trial_energy(v, w, β...; dims = 3)
    Ar = A(v, w, β...) * dims / 3
    Cr = C(v, w, β...) * dims / 3 
    return Ar, Cr
end