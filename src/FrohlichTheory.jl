# FrohlichTheory.jl
# Everything specifically associated with calculating properties of the Frohlich polaron model.

"""
    frohlichalpha(ε_Inf, ε_S, freq, m_eff)

Calculates the Frohlich alpha parameter, for a given dielectric constant, frequency (f) of phonon in Hertz, and effective mass (in units of the bare electron mass).

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
function frohlichalpha(ϵ_optic, ϵ_static, freq, m_eff)
    ω = freq * 2π * 1e12 # frequency to angular velocity
    # Note: we need to add a 4*pi factor to the permitivity of freespace.
    # This gives numeric agreement with literature values.  This is required as
    # the contemporary 1950s and 1960s literature implicitly used atomic units,
    # where the electric constant ^-1 has this factor baked in, k_e=1/(4πϵ_0).
    α = 1 / 2 / (4 * π * ϵ_0) *           # Units: m/F
        (1 / ϵ_optic - 1 / ϵ_static) *   # Units: none
        (eV^2 / (ħ * ω)) *               # Units: F
        sqrt(2 * m_eff * me * ω / ħ)    # Units: 1/m
    return α
end

"""
    frohlichalpha(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff)

Calculates the partial dielectric electron-phonon coupling parameter for a given longitudinal optical phonon mode. 

This decomposes the original Frohlich alpha coupling parameter (defined for a single phonon branch) into contributions from multiple phonon
branches.

# Arguments
- `ϵ_optic::Float64`: is the optical dielectric constant of the material.
- `ϵ_ionic::Float64`: is the ionic dielectric contribution from the phonon mode.
- `ϵ_total::Float64`: is the total ionic dielectric contribution from all phonon modes of the material.
- `phonon_mode_freq::Float64`: is the frequency of the phonon mode (THz).
- `m_eff::Float64` is the band mass of the electron (in units of electron mass m_e).
"""
function frohlichalpha(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff)

    # The Rydberg energy unit
    Ry = eV^4 * me / (2 * ħ^2)

    # Angular phonon frequency for the phonon mode (rad Hz).
    ω = phonon_mode_freq * 2π * 1e12

    # The static dielectric constant. Calculated here instead of inputted so that ionic modes are properly normalised.
    ϵ_static = ϵ_total + ϵ_optic

    # The contribution to the electron-phonon parameter from the currrent phonon mode. 1 / (4π ϵ_0) is the dielectric normalisation.
    α_j = (m_eff * Ry / (ħ * ω))^(1 / 2) * ϵ_ionic / (4π * ϵ_0) / (ϵ_optic * ϵ_static)

    return α_j
end

# Partial dielectric electron-phonon coupling parameter, decomposed into ionic dielectric contributions from each phonon mode of the material.  

"""
    ϵ_ionic_mode(phonon_mode_freq, ir_activity, volume)

Calculate the ionic contribution to the dielectric function for a given phonon mode.

# Arguments
- `phonon_mode_freq::Float64`: is the frequency of the mode in THz.
- `ir_activity::Float64`: is the infra-red activity of the mode in e²amu⁻¹.
- `volume::Float64`: is the volume of the unit cell of the material in m³.
"""
function ϵ_ionic_mode(phonon_mode_freq, ir_activity, volume) # single ionic mode

    # Angular phonon frequency for the phonon mode (rad Hz)
    ω_j = phonon_mode_freq * 2π * 1e12

    # Dielectric contribution from a single ionic phonon mode
    ϵ_mode = eV^2 * ir_activity / (3 * volume * ω_j^2 * amu) 

    # Normalise ionic dielectric contribution with 1 / (4π ϵ_0) (NB: the 4π has been pre-cancelled)
    return ϵ_mode / ϵ_0
end

"""
    ϵ_total(freqs_and_ir_activity, volume)

Calculate the total ionic contribution to the dielectric function from all phonon modes.

# Arguments
- `freqs_and_ir_activity::Matrix{Float64}`: is a matrix containeing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e²amu⁻¹) in the second column.
- `volume::Float64`: is the volume of the unit cell of the material in m^3.
"""
function ϵ_total(freqs_and_ir_activity, volume) # total ionic contribution to dielectric

    # Extract phonon frequencies (THz)
    phonon_freqs = freqs_and_ir_activity[:, 1]

    # Extra infra-red activities (e^2 amu^-1)
    ir_activity = freqs_and_ir_activity[:, 2]

    # Sum over all ionic contribution from each phonon mode
    total_ionic = 0

    for t in eachindex(phonon_freqs)
        total_ionic += ϵ_ionic_mode(phonon_freqs[t], ir_activity[t], volume)
    end

    return total_ionic
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
function frohlich_coupling(k, α, ω; dims = 3)
    r_p = sqrt(1 / 2)
    ω^2 * α * r_p * gamma((dims - 1) / 2) * (2√π / k)^(dims - 1)
end

frohlich_coupling(k, α::Vector, ω::Vector; mb = 1, dims = 3) = sum(frohlich_coupling(k, α[j], ω[j]; mb = mb, dims = dims) for j in eachindex(α))


"""
    frohlich_interaction_energy(v, w, α, ωβ...)

Integral of Eqn. (31) in Feynman 1955. Part of the overall ground-state energy expression.

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
function frohlich_interaction_energy(v, w, α, ωβ...; dims = 3)
    coupling = frohlich_coupling(1, α, ωβ[1]; dims = dims)
    propagator(τ) = length(ωβ) == 1 ? polaron_propagator(τ, v, w) * ωβ[1] : polaron_propagator(τ, v, w, ωβ[2]) * ωβ[1]
    integrand(τ) = phonon_propagator(τ, ωβ...) / sqrt(propagator(τ))
    upper_limit = length(ωβ) == 1 ? Inf : ωβ[2] / 2
    integral, _ = quadgk(τ -> integrand(τ), 0, upper_limit)
    return coupling * ball_surface(dims) / (2π)^dims * sqrt(π / 2) * integral
end

frohlich_interaction_energy(v, w, α::Vector, ω::Vector; dims = 3) = sum(frohlich_interaction_energy(v, w, α[j], ω[j]; dims = dims) for j in eachindex(α))

frohlich_interaction_energy(v, w, α::Vector, ω::Vector, β; dims = 3) = sum(frohlich_interaction_energy(v, w, α[j], ω[j], β; dims = dims) for j in eachindex(α))

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
function frohlich_interaction_energy_k_space(v, w, α, ωβ...; limits = [0, Inf], dims = 3)
    coupling(k) = frohlich_coupling(k, α, ωβ[1]; dims = dims)
    propagator(τ) = length(ωβ) == 1 ? polaron_propagator(τ, v, w) * ωβ[1] : polaron_propagator(τ, v, w, ωβ[2]) * ωβ[1]
    integrand(τ) = phonon_propagator(τ, ωβ...) * spherical_k_integral(coupling, propagator(τ); dims = dims, limits = limits)
    upper_limit = length(ωβ) == 1 ? Inf : ωβ[2] / 2
	integral, _ = quadgk(τ -> integrand(τ), 0, upper_limit)
	return integral
end

"""
    F(v, w, α, ω, β)

Calculates the Helmholtz free energy of the polaron for a material with multiple phonon branches. 
    
This generalises the Osaka 1959 (below Eqn. (22)) and Hellwarth. et al 1999 (Eqn. (62a)) free energy expressions.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `α::Union{Float64, Vector{Float64}}`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz).
- `β::Union{Float64, Vector{Float64}}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.

See Osaka, Y. (1959): https://doi.org/10.1143/ptp.22.437 and Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.
"""
function frohlich_energy(v, w, α, ωβ...; dims = 3, mb = 1)
    A, C = length(ωβ) == 1 ? trial_energy(v, w; dims = dims) : trial_energy(v, w, ωβ[2]; dims = dims)
    B = frohlich_interaction_energy(v, w, α, ωβ...; dims = dims) 
    return -(A + B + C), A, B, C
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
	A, C = length(ωβ) == 1 ? trial_energy(v, w; dims = dims) : trial_energy(v, w, ωβ[2]; dims = dims)
	B = frohlich_interaction_energy_k_space(v, w, α, ωβ...; limits = limits, dims = dims)
	return -(A + B + C), A, B, C
end

"""
----------------------------------------------------------------------
The Complex Impedence and Conductivity of the Polaron.
----------------------------------------------------------------------
"""

# - Polaron optical absorption, complex impedence and complex conductivity
# References:
# [1] R. P. Feynman, R. W. Hellwarth, C. K. Iddings, and P. M. Platzman, Mobility of slow electrons in a polar crystal, PhysicalReview127, 1004 (1962)
# [2] Devreese, J. De Sitter, and M. Goovaerts, “Optical absorption of polarons in thefeynman-hellwarth-iddings-platzman approximation,”Phys. Rev. B, vol. 5, pp. 2367–2381, Mar 1972.
# [3] F. Peeters and J. Devreese, “Theory of polaron mobility,” inSolid State Physics,pp. 81–133, Elsevier, 1984.

function frohlich_structure_factor(t, v, w, α, ωβ...; dims = 3)
	coupling = frohlich_coupling(1, α, ωβ[1]; dims = dims) * ωβ[1]
    propagator = length(ωβ) == 1 ? polaron_propagator(im * t, v, w) * ωβ[1] / 2 : polaron_propagator(im * t, v, w, ωβ[2]) * ωβ[1] / 2
	integral = ball_surface(dims) / (2π)^dims * √π / 4 / propagator^(3/2)
	return 2 / dims * coupling * integral * phonon_propagator(im * t, ωβ...) 
end

"""
	frohlich_structure_factor_k_space(t, v, w, α, ωβ...; limits = [0, Inf])

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
function frohlich_structure_factor_k_space(t, v, w, α, ωβ...; limits = [0, Inf], dims = 3)
	coupling(k) = frohlich_coupling(k, α, ωβ[1]; dims = dims) * k^2 * ωβ[1]
	propagator = length(ωβ) == 1 ? polaron_propagator(im * t, v, w) * ωβ[1] : polaron_propagator(im * t, v, w, ωβ[2]) * ωβ[1]
	integral = spherical_k_integral(coupling, propagator; limits = limits, dims = dims)
	return 2 / dims * phonon_propagator(im * t, ωβ...) * integral
end

frohlich_structure_factor_k_space(t, v, w, α::Vector, ω::Vector; limits = [0, Inf], dims = 3) = sum(frohlich_structure_factor_k_space(t, v, w, α[j], ω[j]; limits = limits, dims = dims) for j in eachindex(α))

frohlich_structure_factor_k_space(t, v, w, α::Vector, ω::Vector, β; limits = [0, Inf], dims = 3) = sum(frohlich_structure_factor_k_space(t, v, w, α[j], ω[j], β; limits = limits, dims = dims) for j in eachindex(α))

function frohlich_memory_function(Ω, v, w, α, ωβ...; dims = 3)
    structure_factor(t) = frohlich_structure_factor(t, v, w, α, ωβ...; dims = dims)
	return polaron_memory_function(Ω, structure_factor)
end

frohlich_memory_function(Ω, v, w, α::Vector, ω::Vector; dims = 3) = sum(frohlich_memory_function(Ω, v, w, α[j], ω[j]; dims = dims) for j in eachindex(α))

frohlich_memory_function(Ω, v, w, α::Vector, ω::Vector, β; dims = 3) = sum(frohlich_memory_function(Ω, v, w, α[j], ω[j], β; dims = dims) for j in eachindex(α))

"""
    frohlich_memory_function_k_space(Ω, v, w, α, ω, β; rₚ = 1, limits = [0, Inf])

Calculate the memory function for the Frohlich model in k-space at finite temperature and frequency.

## Arguments
- `Ω`: a scalar value representing the frequency at which the memorpy function is evaluated.
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
function frohlich_memory_function_k_space(Ω, v, w, α, ωβ...; dims = 3, limits = [0, Inf])
	structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ωβ...; dims = dims, limits = limits)
	return polaron_memory_function(Ω, structure_factor)
end


"""
    frohlich_complex_impedence(Ω, β, α, v, w; rtol = 1e-3, T = nothing, verbose = false)

Calculate the complex impedence Z(Ω) of the polaron at finite temperatures for a given frequency Ω (Eqn. (41) in FHIP 1962). 

# Arguments
- `Ω::Float64`: is the frequency (2π THz) of applied electric field.
- `β::Float64`: is the reduced thermodynamic betas. 
- `α::Float64`: is the Frohlich alpha coupling parameter.
- `v::Float64`: is the 'v' variational parameter.
- `w::Float64`: is the 'w' variational parameter. 
- `rtol`: relative tolerance for the accuracy of any quadrature integrations.
- `T`: is a token used by `make_polaron()` to keep track of the temperature for printing during a calculation. Do not alter.
- `verbose`: is used by `make_polaron()` to specify whether or not to print. Ignore.

See FHIP 1962: https://doi.org/10.1103/PhysRev.127.1004.

"""
function frohlich_complex_impedence(Ω, v, w, α, ωβ...; dims = 3)
    -im * (Ω + frohlich_memory_function(Ω, v, w, α, ωβ...; dims = dims))
end

"""
   frohlich_complex_conductivity(Ω, β, α, v, w; rtol = 1e-3)

Calculate the complex conductivity σ(Ω) of the polaron at finite temperatures for a given frequency Ω (equal to 1 / Z(Ω) with Z the complex impedence). 

# Arguments
- `Ω::Float64`: is the frequency (2π THz) of applied electric field.
- `β::Float64`: is the reduced thermodynamic betas. 
- `α::Float64`: is the Frohlich alpha coupling parameter.
- `v::Float64`: is the 'v' variational parameter.
- `w::Float64`: is the 'w' variational parameter. 
- `rtol`: relative tolerance for the accuracy of any quadrature integrations.

See also [`polaron_complex_impedence`](@ref)
"""
function frohlich_complex_conductivity(Ω, v, w, α, ωβ...; dims = 3)
    return 1 / frohlich_complex_impedence(Ω, v, w, α, ωβ...; dims = dims)
end

"""
    inverse_polaron_mobility(v, w, α, ω, β)

Calculate the inverse of the dc mobility μ of the polaron at finite temperatues (Eqn. (11.5) in F. Peeters and J. Devreese 1984) for a given frequency Ω.

# Arguments
- `v::Float64`: is the 'v' variational parameter.
- `w::Float64`: is the 'w' variational parameter. 
- `α::Float64`: is the Frohlich alpha coupling parameter.
- `ω::Float64`: is the angular phonon frequency.
- `β::Float64`: is the reduced thermodynamic beta. 
- `verbose`: is used by `make_polaron()` to specify whether or not to print. Ignore.

See F. Peeters and J. Devreese 1984: https://doi.org/10.1016/S0081-1947(08)60312-4.

See also [`polaron_mobility`](@ref), [`polaron_complex_conductivity`](@ref)
"""
function inverse_frohlich_mobility(v, w, α, ω, β; dims = 3)
    structure_factor(t) = frohlich_structure_factor(t, v, w, α, ω, β; dims = dims)
    return abs(imag(polaron_memory_function(structure_factor)))
end

"""
    inverse_polaron_mobility(v, w, α::Vector, ω::Vector, β::Vector)

inverse of the polaron mobility, but for multiple phonon modes.
"""
inverse_frohlich_mobility(v, w, α::Vector, ω::Vector, β; dims = 3) = sum(inverse_frohlich_mobility(v, w, α[j], ω[j], β; dims = dims) for j in eachindex(α))

"""
    polaron_mobility(v, w, α, ω, β)

The polaron mobility.

See also [`inverse_polaron_mobility`](@ref)
"""
frohlich_mobility(v, w, α, ω, β; dims = 3) = 1 / inverse_frohlich_mobility(v, w, α, ω, β; dims = dims)
frohlich_mobility(v, w, α, ω; dims = 3) = reduce_array(repeat([Inf], length(α)))


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
function frohlich_mobility_k_space(v, w, α, ω, β; dims = 3, limits = [0, Inf])
	structure_factor(t) = frohlich_structure_factor_k_space(t, v, w, α, ω, β; dims = dims, limits = limits)
	1 / imag(polaron_memory_function(structure_factor))
end

"""
    inverse_FHIP_mobility_lowT(v, w, α, ω, β)

FHIP low-temperature mobility, final result of Feynman1962.
[1.60] in Devreese2016 page 36; 6th Edition of Frohlich polaron notes (ArXiv).

See also [`FHIP_mobility_lowT`](@ref)
"""
function inverse_FHIP_mobility_lowT(v, w, α, ω, β)
    if β == Inf
        return 0
    else
        μ = (w / v)^3 * 3 / (4 * ω^2 * α * β) * exp(ω * β) * exp((v^2 - w^2) * ω / (w^2 * v))
        return 1 / μ
    end
end

"""
    inverse_FHIP_mobility_lowT(v, w, α::Vector, ω::Vector, β::Vector)

Inverse FHIP mobility for multiple phonon modes.

See also [`FHIP_mobility_lowT`](@ref)
"""
inverse_FHIP_mobility_lowT(v, w, α::Vector, ω::Vector, β) = sum(inverse_FHIP_mobility_lowT.(v, w, α, ω, β))

"""
    FHIP_mobility_lowT(v, w, α, ω, β)

FHIP mobility in the low-temperature approximation.

See also [`inverse_FHIP_mobility_lowT`](@ref)
"""
FHIP_mobility_lowT(v, w, α, ω, β) = 1 ./ inverse_FHIP_mobility_lowT(v, w, α, ω, β)

"""
    inverse_Kadanoff_mobility_lowT(v, w, α, ω, β)

Kadanoff low-temperaure mobility, constructed around Boltzmann equation.
Adds factor of 3 / (2 * β) c.f. FHIP, correcting phonon emission behaviour.

See also [`Kadanoff_mobility_lowT`](@ref)
"""
function inverse_Kadanoff_mobility_lowT(v, w, α, ω, β)
    
    if β == Inf
        return 0, 0, 0
    else

        # [1.61] in Devreese2016 - Kadanoff's Boltzmann eqn derived mobility.
        μ_Devreese2016 = (w / v)^3 / (2 * ω^2 * α) * exp(ω * β) * exp((v^2 - w^2) * ω / (w^2 * v))

        # From 1963 Kadanoff (particularly on right-hand-side of page 1367), Eqn. (9), 
        # we define equilibrium number of phonons (just from temperature T and phonon ω):
        # N̄ = (exp(β) - 1)^-1.
        # But! We find that:
        N̄ = exp(-β * ω) 

        # is only way to get Kadanoff1963 to be self-consistent with
        # FHIP, and later statements (Devreese) of the Kadanoff mobility.
        # It suggests that Kadanoff used the wrong identy for Nbar in Eqn. (23b) for
        # the Γ₀ function, and should have used a version with the -1 to
        # account for Bose / phonon statistics!

        # Between Eqns. (23) and (24) in Kadanoff1963, for small momenta skip intergration
        # and sing the fictitious mass M:
        M = (v^2 - w^2) / w^2

        # we get for Γ₀:
        Γ₀ = 2 * α * N̄ * (M + 1)^(1 / 2) * exp(-M * ω / v) * ω^2

        # NB: Kadanoff1963 uses ħ = ω = mb = 1 units. 
        # Factor of omega to get it as a rate relative to phonon frequency
        # Factor of omega*hbar to get it as a rate per energy window.
        # Hence, for the mobility from Eqn. (25) in Kadanoff1963 (SI effective mass q / mb):
        μ_Kadanoff1963 = 1 / ((M + 1) * Γ₀) 

        # # Energy loss is simply energy * rate:
        # energy_loss = ω * Γ₀ / 2π 

        # # Lifetime is:
        # τ = 1 / Γ₀

        return 1 / μ_Devreese2016, 1 / μ_Kadanoff1963, 1/Γ₀
    end
end

"""
    inverse_Kadanoff_mobility_lowT(v, w, α::Vector, ω::Vector, β::Vector)

Inverse Kadanoff mobility for multiple phonon modes.

See also [`Kadanoff_mobility_lowT`](@ref)
"""
inverse_Kadanoff_mobility_lowT(v, w, α::Vector, ω::Vector, β) = map(+, map(x -> 1 ./ x, inverse_Kadanoff_mobility_lowT.(v, w, α, ω, β))...)

"""
    Kadanoff_mobility_lowT(v, w, α, ω, β)

Kadanoff mobility in the low-temperature approximation.

See also [`inverse_Kadanoff_mobility_lowT`](@ref)
"""
Kadanoff_mobility_lowT(v, w, α, ω, β) = 1 ./ inverse_Kadanoff_mobility_lowT(v, w, α, ω, β)

"""
    inverse_Hellwarth_mobility(v, w, α, ω, β)

Calculates the DC mobility using Hellwarth et al. 1999 Eqn. (2).
Directly performs contour integration in Feynman1962, for finite temperature DC mobility.
Eqns. (2) and (1) are going back to the general (pre low-T limit) formulas in Feynman1962.  
To evaluate these, you need to do the explicit contour integration to get the polaron self-energy.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.

See also [`Hellwarth_mobility`](@ref)
"""
function inverse_Hellwarth_mobility(v, w, α, ω, β)

    if β == Inf
        return 0, 0
    else

        R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)
        b = R * β / sinh(β * v / 2) # Feynman1962 version; page 1010, Eqn (47b)
        a = sqrt((β / 2)^2 + R * β * coth(β * v / 2))
        k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u) + eps(Float64))^(-3 / 2) * cos(ω * u) # integrand in (2)
        K = quadgk(u -> k(u, a, b, v), 0, Inf, rtol=1e-3)[1] # numerical quadrature integration of (2)

        # Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS = α / (3 * sqrt(π)) * (β)^(5 / 2) / sinh(ω * β / 2) * (v^3 / w^3) * K
        μ = RHS * ω^(3/2)

        # Hellwarth1999/Biaggio1997, b=0 version... 'Setting b=0 makes less than 0.1% error'
        # So let's test this:
        K_0 = quadgk(u -> k(u, a, 0, v), 0, Inf, rtol=1e-3)[1] # Inserted b=0 into k(u, a, b, v).
        RHS_0 = α / (3 * sqrt(π)) * (β)^(5 / 2) / sinh(ω * β / 2) * (v^3 / w^3) * K_0
        μ_0 = RHS_0 * ω^(3/2)

        return μ, μ_0
    end
end

"""
    inverse_Hellwarth_mobility(v, w, α::Vector, ω::Vector, β::Vector)

Inverse Hellwarth mobility for multiple phonon modes.

See also [`Hellwarth_mobility`](@ref)
"""
inverse_Hellwarth_mobility(v, w, α::Vector, ω::Vector, β) = map(+, map(x -> 1 ./ x, inverse_Hellwarth_mobility.(v, w, α, ω, β))...)

"""
    Hellwarth_mobility(v, w, α, ω, β)

The Hellwarth polaron mobility.

See also [`inverse_Hellwarth_mobility`](@ref)
"""
Hellwarth_mobility(v, w, α, ω, β) = 1 ./ inverse_Hellwarth_mobility(v, w, α, ω, β)

function polaron_effective_mass(v, w, α, ω, β, Ω)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t) = (1 - exp(im * Ω * t / ω)) * imag(S(t, v, w, ω, β))

    integral = real(quadgk(t -> integrand(t), 0, Inf)[1])

    mass = 1 - 2 * α * ω^2 * integral / (3 * √π * Ω^2)

    return mass
end

polaron_effective_mass(v, w, α::Vector, ω::Vector, β) = sum(polaron_effective_mass.(v, w, α, ω, β))