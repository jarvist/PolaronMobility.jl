# MemoryFunction.jl
# References:
# [1] R. P. Feynman, R. W. Hellwarth, C. K. Iddings, and P. M. Platzman, Mobility of slow electrons in a polar crystal, PhysicalReview127, 1004 (1962)
# [2] Devreese, J. De Sitter, and M. Goovaerts, “Optical absorption of polarons in thefeynman-hellwarth-iddings-platzman approximation,”Phys. Rev. B, vol. 5, pp. 2367–2381, Mar 1972.
# [3] F. Peeters and J. Devreese, “Theory of polaron mobility,” inSolid State Physics,pp. 81–133, Elsevier, 1984.

"""
polaron_memory_function_thermal(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64; rtol = 1e-3)

    Calculates the memory function χ(Ω) of the polaron at finite temperatures (equation (35) in FHIP 1962 [1]) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral and corresponds to the error of the entire function. 

    Finite temperature and finite frequency memory function. Does not explicitly include the limits to zero frequency Ω → 0 or zero temperature β → ∞. Although, those limits can be approximated by taking Ω ⪅ O(1e-3) or β ⪆ O(100), explict and exact functions are provided.
"""
function polaron_memory_function_thermal(Ω, β, α, v, w; ω = 1.0, rtol = 1e-3)

    v = v[1]
    w = w[1]

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v)  

    # FHIP1962, page 1009, eqn (35c).
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) * coth(β[1] * v / 2) + x^2 / β[1] - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36).
    S(x) = 2 * α / (3 * √π * 2π) * (exp(1im * x) + 2 * cos(x) / (exp(β[1]) - 1)) / (D(x))^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(x) = (1 - exp(1im * Ω * 2π * x / ω)) * imag(S(x)) / Ω

    # Integrate using adapative quadrature algorithm with relative error tolerance rtol.
    integral, error = ω^2 .* quadgk(x -> integrand(x), 0.0, Inf, rtol = rtol)

    return integral
end

"""
polaron_memory_function_athermal(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) of the polaron at zero temperatures (equations (12, 13) in DSG 1972 [2]) for a given frequency Ω. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral and corresponds to the error of the entire function. 
        
    Zero temperature and finite frequency memory function. Explicitly includes the limit to zero temperature β → ∞.
"""
function polaron_memory_function_athermal(Ω, α, v, w; ω = 1.0, rtol = 1e-3)

    v = v[1]
    w = w[1]

	# FHIP1962, page 1011, eqn (47c).
	R = (v^2 - w^2) / (w^2 * v)

	# FHIP1962, page 1009, eqn (35c) taken to athermal limit.
	D(x) = w^2 / v^2 * (R * (1 - exp(im * v * x)) - 1im * x)

	# FHIP1962, page 1009, eqn (36) taken to athermal limit.
	S(x) = 2 * α / (3 * √π * 2π) * (exp(1im * x)) / (D(x))^(3 / 2)

	# FHIP1962, page 1009, eqn (35a) taken to athermal limit.
	integrand(x) = (1 - exp(1im * Ω * 2π * x / ω)) * imag(S(x)) / Ω

    # Integrand is an exponentially decaying oscillation. Provide an upper cut-off to the integral to ensure convergence. Using Inf results in a NaN. This mimics having a finite sum in the correspond hypergeometric function expansion for the integral.
	integral, error = ω^2 .* quadgk(x -> integrand(x), 0.0, 1e205, rtol = rtol) # 

    # When the frequency Ω ≤ 0 the imaginary part of the memory function is zero. This corresponds to no phonon emission below the phonon frequency of the longitudinal optical mode ω_LO.
	# if Ω <= 1
	# 	return real(integral) + im * 0.0
	# else
	# 	return integral
	# end

    return integral
end

"""
polaron_memory_function_dc(β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function lim(Ω → 0){χ(Ω) / Ω} of the polaron at finite temperatures (where χ(Ω) is from equation (35) in FHIP 1962 [1]) at zero frequency. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral and corresponds to the error of the entire function. 
        
    Zero frequency memory function. Explicitly includes the limit to zero frequency Ω → 0.
"""
function polaron_memory_function_dc(β, α, v, w; ω = 1.0, rtol = 1e-3)

    v = v[1]
    w = w[1]

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v)

    # FHIP1962, page 1009, eqn (35c).
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) * coth(β[1] * v / 2) + x^2 / β[1] - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36).
    S(x) = 2 * α / (3 * √π) * (exp(1im * x) + 2 * cos(x) / (exp(β[1]) - 1)) / (D(x))^(3 / 2)

    # FHIP1962, page 1009, eqn (35a) taken to dc zero frequency limit.
    integrand(x) = -im * x * imag(S(x))

    # Integrate using adapative quadrature algorithm with relative error tolerance rtol.
    integral, error = ω .* quadgk(x -> integrand(x), 0.0, Inf, rtol = rtol)

    return integral
end

"""
polaron_memory_function(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) of the polaron at finite temperatures (equation (35) in FHIP 1962 [1]) for a given frequency Ω. This function includes the zero-temperature and DC limits. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral and corresponds to the error of the entire function. 
        
    Finite temperature and finite frequency memory function, including the limits to zero frequency Ω → 0 or zero temperature β → ∞.
"""
function polaron_memory_function(Ω, β, α, v, w; ω = 1.0, rtol = 1e-3)

    # Zero temperature and frequency is just zero.
	if Ω == 0 && any(x -> x == Inf, β)
		return 0.0 + im * 0.0

    # DC zero frequency limit at finite temperatures.
	elseif Ω == 0 && any(x -> x != Inf, β)
		return polaron_memory_function_dc(β, α, v, w; ω = ω, rtol = rtol)

    # Zero temperature limit at AC finite frequencies.
	elseif Ω != 0 && any(x -> x == Inf, β)
		return polaron_memory_function_athermal(Ω, α, v, w; ω = ω, rtol = rtol)

    # Finite temperatures and frequencies away from zero limits.
	elseif Ω != 0 && any(x -> x != Inf, β)
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

This section of the code is dedicated to calculating the polaron memory function and complex conductivity, generalised from FHIP's expression to the case where multiple phonon modes are present in the material.
"""

"""
function multi_memory_function(Ω::Float64, β::Array{Float64}(undef, 1), α::Array{Float64}(undef, 1), v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), ω::Array{Float64}(undef, 1), m_eff::Float64)

    Calculate polaron complex memory function inclusive of multiple phonon branches j, each with angular frequency ω[j] (rad THz).

     - Ω is the frequency (THz) of applied electric field.
     - β is an array of reduced thermodynamic betas, one for each phonon frequency ω[j]. 
     - α is an array of decomposed Frohlich alphas, one for each phonon frequency ω[j]. 
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - m_eff is the is the conduction band mass of the particle (typically electron / hole, in units of electron mass m_e).
"""
@noinline function polaron_memory_function_thermal(Ω, β::Array, α::Array, v, w; ω = 1.0, rtol = 1e-3)

    # FHIP1962, page 1009, eqn (36).
    S(t, β) = cos(t - 1im * β / 2) / sinh(β / 2) / D_j(-1im * t, β, v, w)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β, Ω) = (1 - exp(1im * Ω * 2π * t)) * imag(S(t, β))

    memory = 0.0

    # Sum over the phonon modes.
    @simd for j in 1:length(ω) 
        
        # print out the current photon frequency and phonon mode frequency (THz).
        # println("Photon frequency = $ν, Phonon mode frequency = $(ω[j] / 2π)")

        # Add the contribution to the memory function from the `jth` phonon mode.
        @inbounds memory += 2 * α[j] * ω[j]^2 * quadgk(t -> integrand(t, β[j], Ω / ω[j]), 0.0, Inf, rtol = rtol)[1] / (3 * √π * Ω * 2π)
    end

    # Print out the value of the memory function.
    # println("Memory function: ", memory)

    return memory
end

@noinline function polaron_memory_function_athermal(Ω, α::Array, v, w; ω = 1.0, rtol = 1e-3)

    # FHIP1962, page 1009, eqn (36).
    S(t) = exp(im * t) / D_j(-1im * t, v, w)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, Ω) = (1 - exp(1im * 2π * Ω * t)) * imag(S(t))

    memory = 0.0

    # Sum over the phonon modes.
    @simd for j in 1:length(ω) 
        
        # print out the current photon frequency and phonon mode frequency (THz).
        # println("Photon frequency = $ν, Phonon mode frequency = $(ω[j] / 2π)")

        # Add the contribution to the memory function from the `jth` phonon mode.
        @inbounds memory += 2 * α[j] * ω[j]^2 * quadgk(t -> integrand(t, Ω / ω[j]), 0.0, 1e205, rtol = rtol)[1] / (3 * √π * Ω * 2π)
    end

    # Print out the value of the memory function.
    # println("Memory function: ", memory)

    return memory
end

@noinline function polaron_memory_function_dc(β::Array, α::Array, v, w; ω = 1.0, rtol = 1e-3)

    # FHIP1962, page 1009, eqn (36).
    S(t, β) = cos(t - 1im * β / 2) / sinh(β / 2) / D_j(-1im * t, β, v, w)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β) = -im * t * imag(S(t, β))

    memory = 0.0

    # Sum over the phonon modes.
    @simd for j in 1:length(ω) 
        
        # print out the current photon frequency and phonon mode frequency (THz).
        # println("Photon frequency = $ν, Phonon mode frequency = $(ω[j] / 2π)")

        # Add the contribution to the memory function from the `jth` phonon mode.
        @inbounds memory += 2 * α[j] * ω[j] * quadgk(t -> integrand(t, β[j]), 0.0, Inf, rtol = rtol)[1] / (3 * √π)
    end

    # Print out the value of the memory function.
    # println("Memory function: ", memory)

    return memory
end