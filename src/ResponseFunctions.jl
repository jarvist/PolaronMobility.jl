# ReponseFunctions.jl
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
function frohlich_complex_impedence(Ω, v, w, α, ωβ...)
    -im * (Ω + frohlich_memory_function(Ω, v, w, α, ωβ...))
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
function frohlich_complex_conductivity(Ω, v, w, α, ωβ...)
    return 1 / frohlich_complex_impedence(Ω, v, w, α, ωβ...)
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
function inverse_frohlich_mobility(v, w, α, ω, β)
    structure_factor(t) = frohlich_structure_factor(t, v, w, α, ω, β)
    return abs(imag(general_memory_function(structure_factor))) / ω
end

"""
    inverse_polaron_mobility(v, w, α::Vector, ω::Vector, β::Vector)

inverse of the polaron mobility, but for multiple phonon modes.
"""
inverse_frohlich_mobility(v, w, α::Vector, ω::Vector, β) = sum(inverse_frohlich_mobility.(v, w, α, ω, β))

"""
    polaron_mobility(v, w, α, ω, β)

The polaron mobility.

See also [`inverse_polaron_mobility`](@ref)
"""
frohlich_mobility(v, w, α, ω, β) = 1 ./ inverse_frohlich_mobility(v, w, α, ω, β)

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
        μ = (w / v)^3 * 3 / (4 * ω^2 * α * β) * exp(ω * β) * exp((v^2 - w^2) / (w^2 * v))
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
        μ_Devreese2016 = (w / v)^3 / (2 * ω * α) * exp(ω * β) * exp((v^2 - w^2) / (w^2 * v))

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
        Γ₀ = 2 * α * N̄ * (M + 1)^(1 / 2) * exp(-M / v) * ω 

        # NB: Kadanoff1963 uses ħ = ω = mb = 1 units. 
        # Factor of omega to get it as a rate relative to phonon frequency
        # Factor of omega*hbar to get it as a rate per energy window.
        # Hence, for the mobility from Eqn. (25) in Kadanoff1963 (SI effective mass q / mb):
        μ_Kadanoff1963 = 1 / ((M + 1) * Γ₀) 

        # # Energy loss is simply energy * rate:
        # energy_loss = ω * Γ₀ / 2π 

        # # Lifetime is:
        # τ = 1 / Γ₀

        return 1 / μ_Devreese2016, 1 / μ_Kadanoff1963, Γ₀
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
        b = R * ω * β / sinh(ω * β * v / 2) # Feynman1962 version; page 1010, Eqn (47b)
        a = sqrt((ω * β / 2)^2 + R * ω * β * coth(ω * β * v / 2))
        k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)
        K = quadgk(u -> k(u, a, b, v), 0, 1e5)[1] # numerical quadrature integration of (2)

        # Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS = α / (3 * sqrt(π)) * (ω * β)^(5 / 2) / sinh(ω * β / 2) * (v^3 / w^3) * K
        μ = RHS

        # Hellwarth1999/Biaggio1997, b=0 version... 'Setting b=0 makes less than 0.1% error'
        # So let's test this:
        K_0 = quadgk(u -> k(u, a, 0, v), 0, 1e5)[1] # Inserted b=0 into k(u, a, b, v).
        RHS_0 = α / (3 * sqrt(π)) * (ω * β)^(5 / 2) / sinh(ω * β / 2) * (v^3 / w^3) * K_0
        μ_0 = RHS_0     

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

