# MemoryFunction.jl
# References:
# [1] R. P. Feynman, R. W. Hellwarth, C. K. Iddings, and P. M. Platzman, Mobility of slow electrons in a polar crystal, PhysicalReview127, 1004 (1962) https://doi.org/10.1103/PhysRev.127.1004.
# [2] Devreese, J. De Sitter, and M. Goovaerts, “Optical absorption of polarons in thefeynman-hellwarth-iddings-platzman approximation,”Phys. Rev. B, vol. 5, pp. 2367–2381, Mar 1972.
# [3] F. Peeters and J. Devreese, “Theory of polaron mobility,” inSolid State Physics,pp. 81–133, Elsevier, 1984.

"""
    polaron_memory_function(Ω, β, α, v, w; ω = 1.0, rtol = 1e-3)

Calculate the memory function χ(Ω) of the polaron at finite temperatures (equation (35) in FHIP 1962) for a given frequency Ω. 
    
Includes conditionals for zero-temperature and DC limits.  
    
# Arguments
- `Ω::Float64`: electric field frequency (2π THz).
- `β::Union{Float64, Vector{Float64}` reduced thermodynamic temperatures (unitless). Either a single value for one phonon frequency or a vector of values for multiple phonon frequencies.
- `α::Union{Float64, Vector{Float64}}`: Frohlich alpha coupling parameter (unitless). Either a single value for one phonon frequency or a vector of values for multiple phonon frequencies.
- `v::Union{Float64, Vector{Float64}}`: variational 'v' parameters that minimise the polaron free energy (unitless).
- `w::Union{Float64, Vector{Float64}}`: variational 'w' parameters that minimise the polaron free energy (unitless).
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
- `rtol = 1e-3`: specifies the relative error tolerance for the QuadGK integral.

See R. P. Feynman, R. W. Hellwarth, C. K. Iddings, and P. M. Platzman, Mobility of slow electrons in a polar crystal, PhysicalReview127, 1004 (1962). https://doi.org/10.1103/PhysRev.127.1004.

See also [`polaron_memory_function_thermal`](@ref), [`polaron_memory_function_athermal`](@ref), [`polaron_memory_function_dc`](@ref).
"""
function polaron_memory_function(Ω, β, α, v, w; ω=1.0)

    # Zero temperature and frequency is just zero.
    if Ω == 0.0 && any(x -> x == Inf64, β)
        return 0.0 + im * Inf64

        # DC zero frequency limit at finite temperatures.
    elseif Ω == 0 && any(x -> x != Inf64, β)
        return polaron_memory_function_dc(β, α, v, w; ω=ω)

        # Zero temperature limit at AC finite frequencies.
    elseif Ω != 0 && any(x -> x == Inf64, β)
        return polaron_memory_function_athermal(Ω, α, v, w; ω=ω)

        # Finite temperatures and frequencies away from zero limits.
    elseif Ω != 0 && any(x -> x != Inf64, β)
        return polaron_memory_function_thermal(Ω, β, α, v, w; ω=ω)

        # Any other frequencies or temperatures (e.g. negative or complex) prints error message.
    else
        println("Photon frequency Ω and thermodynamic temperature β must be ≥ 0.0.")
    end
end

"""
----------------------------------------------------------------------
Multiple Branch Polaron Memory Function and Complex Conductivity
----------------------------------------------------------------------

This section of the code is dedicated to calculating the polaron memory function and complex conductivity, generalised from FHIP's expression to the case where multiple phonon modes are present in the material.  
"""

"""
    function polaron_memory_function_thermal(Ω, β, α, v, w; ω = 1.0, rtol = 1e-3)

Calculate polaron complex memory function at finite temperature inclusive of multiple phonon branches j, each with angular frequency ω[j] (2π THz). 

# Arguments
- `Ω::Float64`: is the frequency (THz) of applied electric field.
- `β::Vector{Float64}`: is an array of reduced thermodynamic betas, one for each phonon frequency ω[j]. 
- `α::Vector{Float64}`: is an array of decomposed Frohlich alphas, one for each phonon frequency ω[j]. 
- `v::Union{Float64, Vector{Float64}}`: is a vector or single value of the 'v' variational parameter(s).
- `w::Union{Float64, Vector{Float64}}`: is a vector or single value of the 'w' variational parameter(s).
- `m_eff::Float64`: is the bare-band mass of the particle (typically electron / hole, in units of electron mass m_e).
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
- `rtol = 1e-3`: specifies the relative error tolerance for the QuadGK integral.

See also [`polaron_memory_function`](@ref).
"""
function polaron_memory_function_thermal(Ω, β, α, v, w; ω=1.0)

    v, w = v[1], w[1]
    D(τ, β) = τ * (1 - τ / β) + (v^2 - w^2) * ((1 + exp(-v * β) - exp(-v * τ) - exp(v * (τ - β))) / (v * (1 - exp(-v * β))) - τ * (1 - τ / β)) / v^2

    # FHIP1962, page 1009, eqn (36).
    S(t, β) = (exp(1im * t) + exp(-1im * t - β)) / (1 - exp(-β)) / D(-1im * t, β)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β, Ω) = (1 - exp(1im * Ω * 2π * t)) * imag(S(t, β))

    integral = quadgk.(t -> integrand.(t, β, Ω ./ ω), 0.0, Inf64, rtol=1e-4)[1]

    memory = sum(2 .* α .* ω .^ 2 .* integral ./ (3 * √π * Ω * 2π))

    return memory
end

"""
    function polaron_memory_function_athermal(Ω, α, v, w; ω = 1.0, rtol = 1e-3)

Calculate polaron complex memory function at zero-temperature inclusive of multiple phonon branches j, each with angular frequency ω[j] (2π THz). 

# Arguments
- `Ω::Float64`: is the frequency (THz) of applied electric field.
- `α::Vector{Float64}`: is an array of decomposed Frohlich alphas, one for each phonon frequency ω[j]. 
- `v::Union{Float64, Vector{Float64}}`: is a vector or single value of the 'v' variational parameter(s).
- `w::Union{Float64, Vector{Float64}}`: is a vector or single value of the 'w' variational parameter(s).
- `m_eff::Float64`: is the bare-band mass of the particle (typically electron / hole, in units of electron mass m_e).
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
- `rtol = 1e-3`: specifies the relative error tolerance for the QuadGK integral.

See also [`polaron_memory_function`](@ref).
"""
function polaron_memory_function_athermal(Ω, α, v, w; ω=1.0)

    v, w = v[1], w[1]
    D(τ) = τ + ((v^2 - w^2) / v^2) * ((1 - exp(-v * τ)) / v - τ)

    # FHIP1962, page 1009, eqn (36).
    S(t) = exp(1im * t) / D(-1im * t)^(3 / 2)

    # R = (v^2 - w^2) / w^2 / v

    # # FHIP1962, page 1009, eqn (36).
    # S(t) = exp(1im * t) / (R * (1 - exp(1im * v * t)) - 1im * t)^(3/2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, Ω) = (1 - exp(1im * 2π * Ω * t)) * imag(S(t))

    memory = sum(2 .* α .* ω .^2 .* quadgk.(t -> integrand.(t, Ω ./ ω), 0.0, 1e3, rtol=1e-4)[1] ./ (3 * √π * Ω * 2π))

    return memory
end

"""
    function polaron_memory_function_dc(β, α, v, w; ω = 1.0, rtol = 1e-3)

Calculate zero-frequency polaron complex memory function at finite temperature inclusive of multiple phonon branches j, each with angular frequency ω[j] (2π THz). 

# Arguments
- `β::Vector{Float64}`: is an array of reduced thermodynamic betas, one for each phonon frequency ω[j]. 
- `α::Vector{Float64}`: is an array of decomposed Frohlich alphas, one for each phonon frequency ω[j]. 
- `v::Union{Float64, Vector{Float64}}`: is a vector or single value of the 'v' variational parameter(s).
- `w::Union{Float64, Vector{Float64}}`: is a vector or single value of the 'w' variational parameter(s).
- `m_eff::Float64`: is the bare-band mass of the particle (typically electron / hole, in units of electron mass m_e).
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
- `rtol = 1e-3`: specifies the relative error tolerance for the QuadGK integral.

See also [`polaron_memory_function`](@ref).
"""
function polaron_memory_function_dc(β, α, v, w; ω=1.0)

    v, w = v[1], w[1]
    D(τ, β) = τ * (1 - τ / β) + (v^2 - w^2) * ((1 + exp(-v * β) - exp(-v * τ) - exp(v * (τ - β))) / (v * (1 - exp(-v * β))) - τ * (1 - τ / β)) / v^2

    # FHIP1962, page 1009, eqn (36).
    S(t, β) = (exp(1im * t) + exp(-1im * t - β)) / (1 - exp(-β)) / D(-1im * t, β)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β) = -im * t * imag(S(t, β))

    integral = quadgk.(t -> integrand.(t, β), 0.0, Inf64, rtol=1e-4)[1]

    memory = sum(2 .* α .* ω .* integral ./ (3 * √π))

    return memory
end

