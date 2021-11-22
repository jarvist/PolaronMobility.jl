# OedipusRex.jl
# - Polaron optical absorption
#     A WORK IN PROGRESS!

"""
----------------------------------------------------------------------
Polaron absorption coefficient Γ(Ω).
----------------------------------------------------------------------
"""

"""
optical_absorption(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the absorption coefficient Γ(Ω) for the polaron at at finite temperatures
    (equation (11a) in [1]) for a given frequency Ω.  β is thermodynamic beta.  v and w are
    the variational Polaron parameters that minimise the free energy, for the supplied
    α Frohlich coupling.  

    [1] Devreese, J., De Sitter, J., & Goovaerts, M. (1972). Optical Absorption of Polarons
    in the Feynman-Hellwarth-Iddings-Platzman Approximation. Physical Review B, 5(6),
    2367–2381.  doi:10.1103/physrevb.5.2367 

"""
function optical_absorption(Ω, β, α, v, w)
    X = χ(Ω, β, α, v, w)[1]
    ℜχ = real(X)
    ℑχ = imag(X)
    result = Ω * ℑχ / (Ω^4 - 2 * Ω^2 * ℜχ + abs2(X))
    return result
end

