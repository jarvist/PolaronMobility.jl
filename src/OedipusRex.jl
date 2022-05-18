# OedipusRex.jl
# - Polaron optical absorption, complex impedence and complex conductivity
# References:
# [1] R. P. Feynman, R. W. Hellwarth, C. K. Iddings, and P. M. Platzman, Mobility of slow electrons in a polar crystal, PhysicalReview127, 1004 (1962)
# [2] Devreese, J. De Sitter, and M. Goovaerts, “Optical absorption of polarons in thefeynman-hellwarth-iddings-platzman approximation,”Phys. Rev. B, vol. 5, pp. 2367–2381, Mar 1972.
# [3] F. Peeters and J. Devreese, “Theory of polaron mobility,” inSolid State Physics,pp. 81–133, Elsevier, 1984.

"""
----------------------------------------------------------------------
Polaron absorption coefficient Γ(Ω).
----------------------------------------------------------------------
"""

"""
optical_absorption(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64; rtol = 1e-3)

    Calculate the absorption coefficient Γ(Ω) for the polaron at at finite temperatures
    (equation (11a) in [1]) for a given frequency Ω.  β is thermodynamic beta. v and w are
    the variational Polaron parameters that minimise the free energy, for the supplied
    α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral
    in the memory function.  

    [1] Devreese, J., De Sitter, J., & Goovaerts, M. (1972). Optical Absorption of Polarons
    in the Feynman-Hellwarth-Iddings-Platzman Approximation. Physical Review B, 5(6),
    2367–2381. doi:10.1103/physrevb.5.2367 

"""
function optical_absorption(Ω, β, α, v, w; rtol = 1e-3)
    real(complex_conductivity(Ω, β, α, v, w; rtol = rtol))
end

"""
----------------------------------------------------------------------
The Complex Impedence and Conductivity of the Polaron.
----------------------------------------------------------------------
"""

"""
polaron_complex_impedence(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the complex impedence Z(Ω) of the polaron at finite temperatures for a given
    frequency Ω (equation (41) in FHIP 1962 [1]). β is the thermodynamic beta. v and w are
    the variational polaron parameters that minimise the free energy, for the supplied
    α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK
    integral in the memory function.  
"""
function polaron_complex_impedence(Ω, β, α, v, w; ω = 1.0, rtol = 1e-3, T = nothing, verbose = false)
	impedance = -im * Ω * 2π + im * polaron_memory_function(Ω, β, α, v, w; ω = ω, rtol = rtol)

    if verbose
        println("\e[2K", "Process: $(count) / $processes ($(round.(count / processes * 100, digits = 1)) %) | T = $(round.(T, digits = 3)) | Ω = $(round.(Ω, digits = 3)) | Z = $(round.(impedance, digits = 3))")
        print("\033[F")

        global count += 1
    end

    return impedance
end

"""
polaron_complex_conductivity(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the complex conductivity σ(Ω) of the polaron at finite temperatures for
    a given frequency Ω (equal to 1 / Z(Ω) with Z the complex impedence). β is the
    thermodynamic beta. v and w are the variational polaron parameters that minimise the
    free energy, for the supplied α Frohlich coupling. rtol specifies the relative error
    tolerance for the QuadGK integral in the memory function.  
"""
function polaron_complex_conductivity(Ω, β, α, v, w; ω = 0.0, rtol = 1e-3)
	return 1 / polaron_complex_impedence(Ω, β, α, v, w; ω = ω, rtol = rtol)
end

"""
polaron_mobility(β::Float64, α::Float64, v::Float64, w::Float64; rtol = 1e-3)

    Calculate the dc mobility μ of the polaron at finite temperatues (equation (11.5) in
    [3]) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational
    polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
    rtol specifies the relative error for the integral to reach.  
"""
function polaron_mobility(β, α, v, w; ω = 1.0, rtol = 1e-3, T = nothing, verbose = false)
    
	mobility = abs(1 / imag(polaron_memory_function_dc(β, α, v, w; ω = ω, rtol = rtol)))

    if verbose
        println("\e[2K", "Process: $(count) / $processes ($(round.(count / processes * 100, digits = 1)) %) | T = $(round.(T, digits = 3)) | μ = $(round.(mobility, digits = 3))")
        print("\033[F")

        global count += 1
    end

    return mobility
end

