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
function polaron_memory_function_thermal(Ω, β, α, v, w; rtol = 1e-3)

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v) \omega   

    # FHIP1962, page 1009, eqn (35c).
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) * coth(β * v / 2) + x^2 / β - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36).
    S(x) = 2 * α / (3 * √π) * (exp(1im * x) + 2 * cos(x) / (exp(β) - 1)) / (D(x))^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(x) = (1 - exp(1im * Ω * x)) * imag(S(x)) / Ω

    # Integrate using adapative quadrature algorithm with relative error tolerance rtol.
    integral, error = quadgk(x -> integrand(x), 0.0, Inf, rtol = rtol)

    return integral
end

"""
polaron_memory_function_athermal(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) of the polaron at zero temperatures (equations (12, 13) in DSG 1972 [2]) for a given frequency Ω. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral and corresponds to the error of the entire function. 
        
    Zero temperature and finite frequency memory function. Explicitly includes the limit to zero temperature β → ∞.
"""
function polaron_memory_function_athermal(Ω, α, v, w; rtol = 1e-3)

	# FHIP1962, page 1011, eqn (47c).
	R = (v^2 - w^2) / (w^2 * v)

	# FHIP1962, page 1009, eqn (35c) taken to athermal limit.
	D(x) = w^2 / v^2 * (R * (1 - exp(im * v * x)) - 1im * x)

	# FHIP1962, page 1009, eqn (36) taken to athermal limit.
	S(x) = 2 * α / (3 * √π) * (exp(1im * x)) / (D(x))^(3 / 2)

	# FHIP1962, page 1009, eqn (35a) taken to athermal limit.
	integrand(x) = (1 - exp(1im * Ω * x)) * imag(S(x)) / Ω

    # Integrand is an exponentially decaying oscillation. Provide an upper cut-off to the integral to ensure convergence. Using Inf results in a NaN. This mimics having a finite sum in the correspond hypergeometric function expansion for the integral.
	integral, error = quadgk(x -> integrand(x), 0.0, 1e205, rtol = rtol) # 

    # When the frequency Ω ≤ 0 the imaginary part of the memory function is zero. This corresponds to no phonon emission below the phonon frequency of the longitudinal optical mode ω_LO.
	if Ω <= 1
		return real(integral) + im * 0.0
	else
		return integral
	end
end

"""
polaron_memory_function_dc(β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function lim(Ω → 0){χ(Ω) / Ω} of the polaron at finite temperatures (where χ(Ω) is from equation (35) in FHIP 1962 [1]) at zero frequency. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral and corresponds to the error of the entire function. 
        
    Zero frequency memory function. Explicitly includes the limit to zero frequency Ω → 0.
"""
function polaron_memory_function_dc(β, α, v, w; rtol = 1e-3)

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v)

    # FHIP1962, page 1009, eqn (35c).
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) * coth(β * v / 2) + x^2 / β - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36).
    S(x) = 2 * α / (3 * √π) * (exp(1im * x) + 2 * cos(x) / (exp(β) - 1)) / (D(x))^(3 / 2)

    # FHIP1962, page 1009, eqn (35a) taken to dc zero frequency limit.
    integrand(x) = im * x * imag(S(x))

    # Integrate using adapative quadrature algorithm with relative error tolerance rtol.
    integral, error = quadgk(x -> integrand(x), 0.0, Inf, rtol = rtol)

    return integral
end

"""
polaron_memory_function(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) of the polaron at finite temperatures (equation (35) in FHIP 1962 [1]) for a given frequency Ω. This function includes the zero-temperature and DC limits. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral and corresponds to the error of the entire function. 
        
    Finite temperature and finite frequency memory function, including the limits to zero frequency Ω → 0 or zero temperature β → ∞.
"""
function polaron_memory_function(Ω, β, α, v, w; rtol = 1e-3)

    # Zero temperature and frequency is just zero.
	if Ω == 0 && β == Inf
		return 0.0 + im * 0.0

    # DC zero frequency limit at finite temperatures.
	elseif Ω == 0 && β != Inf
		return polaron_memory_function_dc(β, α, v, w; rtol = rtol)

    # Zero temperature limit at AC finite frequencies.
	elseif Ω != 0 && β == Inf
		return polaron_memory_function_athermal(Ω, α, v, w; rtol = rtol)

    # Finite temperatures and frequencies away from zero limits.
	elseif Ω != 0 && β != Inf
		return polaron_memory_function_thermal(Ω, β, α, v, w; rtol = rtol)

    # Any other frequencies or temperatures (e.g. negative or complex) prints error message.
	else
		println("Photon frequency Ω and thermodynamic temperature β must be ≥ 0.0.")
	end
end

"""
----------------------------------------------------------------------
The Complex Impedence and Conductivity of the Polaron.
----------------------------------------------------------------------
"""

"""
polaron_complex_impedence(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the complex impedence Z(Ω) of the polaron at finite temperatures for a given frequency Ω (equation (41) in FHIP 1962 [1]). β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral in the memory function. 
"""
function polaron_complex_impedence(Ω, β, α, v, w; rtol = 1e-3)
	return -im * Ω + im * polaron_memory_function(Ω, β, α, v, w; rtol = rtol)
end

"""
polaron_complex_conductivity(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the complex conductivity σ(Ω) of the polaron at finite temperatures for a given frequency Ω (equal to 1 / Z(Ω) with Z the complex impedence). β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral in the memory function. 
"""
function polaron_complex_conductivity(Ω, β, α, v, w; rtol = 1e-3)
	return 1 / polaron_complex_impedence(Ω, β, α, v, w; rtol = rtol)
end

"""
----------------------------------------------------------------------
The Moblity of the Polaron.
----------------------------------------------------------------------
"""

"""
polaron_mobility(β::Float64, α::Float64, v::Float64, w::Float64; rtol = 1e-3)

    Calculate the dc mobility μ of the polaron at finite temperatues (equation (11.5) in [3]) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error for the integral to reach.
"""

function polaron_mobility(β, α, v, w; rtol = 1e-3)
	return 1 / imag(polaron_memory_function_dc(β, α, v, w; rtol = rtol))
end