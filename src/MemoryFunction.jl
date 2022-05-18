# MemoryFunction.jl
# References:
# [1] R. P. Feynman, R. W. Hellwarth, C. K. Iddings, and P. M. Platzman, Mobility of slow electrons in a polar crystal, PhysicalReview127, 1004 (1962)
# [2] Devreese, J. De Sitter, and M. Goovaerts, “Optical absorption of polarons in thefeynman-hellwarth-iddings-platzman approximation,”Phys. Rev. B, vol. 5, pp. 2367–2381, Mar 1972.
# [3] F. Peeters and J. Devreese, “Theory of polaron mobility,” inSolid State Physics,pp. 81–133, Elsevier, 1984.

"""
polaron_memory_function(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) of the polaron at finite temperatures (equation (35)
    in FHIP 1962 [1]) for a given frequency Ω. This function includes the zero-temperature
    and DC limits. β is the thermodynamic beta. v and w are the variational polaron
    parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol
    specifies the relative error tolerance for the QuadGK integral and corresponds to the
    error of the entire function. 
        
    Finite temperature and finite frequency memory function, including the limits to zero
    frequency Ω → 0 or zero temperature β → ∞.  

    [1] R. P. Feynman, R. W. Hellwarth, C. K. Iddings, and P. M. Platzman, Mobility of slow
    electrons in a polar crystal, PhysicalReview127, 1004 (1962)

"""
function polaron_memory_function(Ω, β, α, v, w; ω = 1.0, rtol = 1e-3)

    # Zero temperature and frequency is just zero.
	if Ω == 0 && any(x -> x == Inf64, β)
		return 0.0 + im * 0.0

    # DC zero frequency limit at finite temperatures.
	elseif Ω == 0 && any(x -> x != Inf64, β)
		return polaron_memory_function_dc(β, α, v, w; ω = ω, rtol = rtol)

    # Zero temperature limit at AC finite frequencies.
	elseif Ω != 0 && any(x -> x == Inf64, β)
		return polaron_memory_function_athermal(Ω, α, v, w; ω = ω, rtol = rtol)

    # Finite temperatures and frequencies away from zero limits.
	elseif Ω != 0 && any(x -> x != Inf64, β)
		return polaron_memory_function_thermal(Ω, β, α, v, w; ω = ω, rtol = rtol)

    # Any other frequencies or temperatures (e.g. negative or complex) prints error message.
	else
		println("Photon frequency Ω and thermodynamic temperature β must be ≥ 0.0.")
	end
end

"""
----------------------------------------------------------------------
Multiple Branch Polaron Memory Function and Complex Conductivity
----------------------------------------------------------------------

This section of the code is dedicated to calculating the polaron memory function and complex
conductivity, generalised from FHIP's expression to the case where multiple phonon modes are
present in the material.  
"""

"""
function multi_memory_function(Ω::Float64, β::Array{Float64}(undef, 1), α::Array{Float64}(undef, 1), v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), ω::Array{Float64}(undef, 1), m_eff::Float64)

    Calculate polaron complex memory function inclusive of multiple phonon branches j, each
    with angular frequency ω[j] (rad THz).

     - Ω is the frequency (THz) of applied electric field.
     - β is an array of reduced thermodynamic betas, one for each phonon frequency ω[j]. 
     - α is an array of decomposed Frohlich alphas, one for each phonon frequency ω[j]. 
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - m_eff is the is the conduction band mass of the particle (typically electron / hole, in units of electron mass m_e).
"""
function polaron_memory_function_thermal(Ω, β, α, v, w; ω = 1.0, rtol = 1e-3)

    # FHIP1962, page 1009, eqn (36).
    S(t, β) = (exp(1im * t) + exp(-1im * t - β)) / (1 - exp(-β)) / D_j(-1im * t, β, v, w)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β, Ω) = (1 - exp(1im * Ω * 2π * t)) * imag(S(t, β))

    integral = map((x, y) -> quadgk(t -> integrand(t, x, Ω / y), 0.0, Inf64, rtol = rtol)[1], β, ω)

     memory = sum(2 .* α .* ω .^2 .* integral ./ (3 * √π * Ω * 2π))

    return memory
end

function polaron_memory_function_athermal(Ω, α, v, w; ω = 1.0, rtol = 1e-3)
    # FHIP1962, page 1009, eqn (36).
    S(t) = exp(1im * t) / D_j(-1im * t, v, w)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, Ω) = (1 - exp(1im * 2π * Ω * t)) * imag(S(t))

    memory = 0.0

    for j in 1:length(ω)

        integral = quadgk(t -> integrand(t, Ω / ω[j]), 0.0, 1 / rtol, rtol = rtol)[1]
        memory += 2 * α[j] * ω[j]^2 * integral / (3 * √π * Ω * 2π)

    end

    return memory
end

function polaron_memory_function_dc(β, α, v, w; ω = 1.0, rtol = 1e-3)
    # FHIP1962, page 1009, eqn (36).
    S(t, β) = (exp(1im * t) + exp(-1im * t - β)) / (1 - exp(-β)) / D_j(-1im * t, β, v, w)^(3 / 2)
    
    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β) = -im * t * imag(S(t, β))

    integral = map(x -> quadgk(t -> integrand(t, x), 0.0, 1 / rtol, rtol = rtol)[1], β)

    memory = sum(2 .* α .* ω .* integral ./ (3 * √π))

    return memory
end

