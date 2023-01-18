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
    optical_absorption(Ω, β, α, v, w; rtol = 1e-3)

Calculate the absorption coefficient Γ(Ω) for the polaron at at finite temperatures (Eqn. (11a) in DSG 1972) for a given frequency Ω.  

# Arguments
- `Ω::Float64`: is the frequency (2π THz) of applied electric field.
- `β::Float64`: is the reduced thermodynamic betas. 
- `α::Float64`: is the Frohlich alpha coupling parameter.
- `v::Float64`: is the 'v' variational parameter.
- `w::Float64`: is the 'w' variational parameter. 
- `rtol`: relative tolerance for the accuracy of any quadrature integrations.

See DSG 1972: https://doi.org/10.1103/PhysRevB.5.2367.
"""
function optical_absorption(Ω, β, α, v, w)
    real(complex_conductivity(Ω, β, α, v, w))
end

"""
----------------------------------------------------------------------
The Complex Impedence and Conductivity of the Polaron.
----------------------------------------------------------------------
"""

"""
    polaron_complex_impedence(Ω, β, α, v, w; rtol = 1e-3, T = nothing, verbose = false)

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
function polaron_complex_impedence(Ω, β, α, v, w; ω=1.0)
    impedance = -im * Ω * 2π + im * polaron_memory_function(Ω, β, α, v, w; ω=ω)

    return impedance
end

"""
    polaron_complex_conductivity(Ω, β, α, v, w; rtol = 1e-3)

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
function polaron_complex_conductivity(Ω, β, α, v, w; ω=1.0)
    return 1 / polaron_complex_impedence(Ω, β, α, v, w; ω=ω)
end

"""
    polaron_mobility(β, α, v, w; rtol = 1e-3, T = nothing, verbose = false)

Calculate the dc mobility μ of the polaron at finite temperatues (Eqn. (11.5) in F. Peeters and J. Devreese 1984) for a given frequency Ω.

# Arguments
- `β::Float64`: is the reduced thermodynamic betas. 
- `α::Float64`: is the Frohlich alpha coupling parameter.
- `v::Float64`: is the 'v' variational parameter.
- `w::Float64`: is the 'w' variational parameter. 
- `rtol`: relative tolerance for the accuracy of any quadrature integrations.
- `T`: is a token used by `make_polaron()` to keep track of the temperature for printing during a calculation. Do not alter.
- `verbose`: is used by `make_polaron()` to specify whether or not to print. Ignore.

See F. Peeters and J. Devreese 1984: https://doi.org/10.1016/S0081-1947(08)60312-4.

See also [`polaron_complex_conductivity`](@ref)
"""
function polaron_mobility(β, α, v, w; ω=1.0)
    if any(x -> x == Inf, β)
        return Inf
    else
        return abs(1 / imag(polaron_memory_function_dc(β, α, v, w; ω=ω)))
    end
end

