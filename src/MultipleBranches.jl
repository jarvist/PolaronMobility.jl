# MultipleBranches.jl
#  Extending the Feynman theory to multiple phonon branches

"""
    frohlichPartial((f, ϵ_mode); ϵ_o, ϵ_s, meff)

Calculate a (partial) dielectric electron-phonon coupling element.

# Arguments
- `f`: frequency of mode in THz.
- `ϵ_mode`: this mode's contribution to dielectric.
- `ϵ_o`: optical dielectric.
-  `ϵ_s`: total static dielectric contribution.
"""
function frohlichPartial((f, ϵ_mode); ϵ_o, ϵ_s, meff)
    ω = f * 1e12 * 2π
    α = 1 / (4π * ϵ_0) * ϵ_mode / (ϵ_o * ϵ_s) * (q^2 / ħ) * sqrt(meff * me / (2 * ω * ħ))
    return α
end

# deprecated signature wrapped via multiple dispatch
frohlichPartial(ϵ_o, ϵ_s, ϵ_mode, f, meff) = frohlichPartial((f, ϵ_mode), ϵ_o = ϵ_o, ϵ_s = ϵ_s, meff = meff)

rows(M::Matrix) = map(x -> reshape(getindex(M, x, :), :, size(M)[2]), 1:size(M)[1])

"""
    IRtoDielectric(IRmodes, volume)

From absolute value of IR activities of phonon modes, generate a per-mode
contribution to the low-frequency dielectric constant.

IRmodes are tuples f, S with Frequency in THz; InfraRed activity in e^2 amu^-1.
"""
function IRtoDielectric(IRmodes,volume)
    ϵ = 0.0 #* q^2 * THz^-2 * amu^-1 * m^-3
    for r in rows(IRmodes)
        f, S = r # frequency in THz; activity in e^2 amu^-1
        f = f * 1E12 #* THz
        ω = 2π * f
        S = S * q^2 / amu
        ϵ_mode = S / ω^2 / volume
        ϵ_mode /= 3 # take isotropic average = divide by 3
        ϵ += ϵ_mode
        println("Mode f= $f S= $S ϵ_mode = $(ϵ_mode / ϵ_0)")
    end
    println("Raw ionic dielectric contribution: $ϵ absolute $(ϵ / ϵ_0) relative")
    return ϵ / ϵ_0
end

"""
    IRtoalpha(IR, volume)

Calculates contribution to dielectric constant for each polar phonon mode, and thereby the Frohlich alpha contribution for this mode.
"""
function IRtoalpha(IR; volume, ϵ_o, ϵ_s, meff)
    ϵ = 0.0 #* q^2 * THz^-2 * amu^-1 * m^-3
    α_sum = 0.0
    for r in rows(IR)
        f, S = r # frequency in THz; activity in e^2 amu^-1
        ω = 2π * f * 1e12
        S = S * q^2 / amu
        ϵ_mode = S / ω^2 / volume
        ϵ_mode /= 3 # take isotropic average = divide by 3
        ϵ_mode /= ϵ_0 # reduced dielectric constant
        ϵ += ϵ_mode
        #println("Mode f= $f S= $S ϵ_mode = $(upreferred(ϵ_mode/u"ϵ0"))")
        α = frohlichPartial(ϵ_o, ϵ_s, ϵ_mode, f, meff)
        α_sum += α
        if (α > 0.1)
            println("Notable Mode f = $f, α_partial = $α.")
        end
    end
    println("Sum alpha: $(α_sum)")
    return α_sum
end

"""
    DieletricFromIRmode(IRmode)

Calculate dielectric from an individual mode.

IRmode is a tuple f, S with Frequency in THz; InfraRed activity in e^2 amu^-1.
"""
function DielectricFromIRmode(IRmode; volume)
    f, S = IRmode # Assumes Frequency in THz; InfraRed activity in e^2 amu^-1
    ω = 2π * f * 1e12
    S = S * q^2 / amu
    ϵ_mode = S / ω^2 / volume
    ϵ_mode /= 3 # take isotropic average = divide by 3
    ϵ_mode /= ϵ_0 # reduced dielectric
    return ϵ_mode
end

"""
    Hellwarth1999mobilityRHS((α, (v, w) ,f), effectivemass, T)

Calculates the DC mobility using Hellwarth et al. 1999 Eqn. (2).

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
function Hellwarth1999mobilityRHS((α, (v, w) ,f), effectivemass, T)
    mb = effectivemass * MassElectron
    ω = f * 1e12 * 2π
    βred = ħ * ω / (kB * T)

    R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)
    b = R * βred / sinh(βred * v / 2) # Feynman1962 version; page 1010, Eqn (47b)
    a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))
    k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3/2) * cos(u) # integrand in (2)
    K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

    # Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
    RHS = α / (3 * sqrt(π)) * βred^(5/2) / sinh(βred / 2) * (v^3 / w^3) * K
    μ = RHS^(-1) * q / (ω * mb)

    return 1 / μ
end

"""
----------------------------------------------------------------------
Multiple Branch Frohlich Alpha
----------------------------------------------------------------------

Partial dielectric electron-phonon coupling parameter, decomposed into ionic dielectric contributions from each phonon mode of the matieral.  
"""

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
    ω_j = 2π * phonon_mode_freq * 1e12 

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
    total_ionic = 0.0
    
    for t in eachindex(phonon_freqs)
        total_ionic += ϵ_ionic_mode(phonon_freqs[t], ir_activity[t], volume)
    end

    return total_ionic
end

"""
    effective_freqs(freqs_and_ir_activity, num_var_params)

Generates a matrix of effective phonon modes with frequencies and infra-red activities derived from a larger matrix using the Principal Component Analysis (PCA) method.

# Arguments
- `freqs_and_ir_activity::Matrix{Float64}`: is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e²amu⁻¹) in the second column.
- `num_var_params::Integer`: is the number of effective modes required (which needs to be less than the number of modes in `freqs_and_ir_activity``).

*** POSSIBLY REDUNDANT ***
"""
function effective_freqs(freqs_and_ir_activity, num_var_params) # PCA Algorithm

    # Check that the number of effective modes is less than the number of actual phonon modes.
    if num_var_params >= size(freqs_and_ir_activity)[1]

        println("The number of effective phonon modes has to be less than the total number of phonon modes.")

    else

        # Centralise data by subtracting the columnwise mean
        standardized_matrix = freqs_and_ir_activity' .- mean(freqs_and_ir_activity', dims = 2) 

        # Calculate the covariance matrix S' * S. Matrix size is (n - 1) x (n - 1) for number of params (here n = 2)
        covariance_matrix = standardized_matrix' * standardized_matrix 

        # Extract eigenvectors of the covariance matrix
        eigenvectors = eigvecs(covariance_matrix) 

        # Project the original data along the covariance matrix eigenvectors and undo the centralisation
        reduced_matrix = standardized_matrix[:, 1:num_var_params] * eigenvectors[1:num_var_params, 1:num_var_params] * 
        eigenvectors[1:num_var_params, 1:num_var_params]' .+ mean(freqs_and_ir_activity', dims = 2)

        # Resultant matrix is positive definite and transposed.
        return abs.(reduced_matrix')
    end
end

"""
    multi_frohlichalpha(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff)

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
function multi_frohlichalpha(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff) 

    # The Rydberg energy unit
    Ry = eV^4 * me / (2 * ħ^2)

    # Angular phonon frequency for the phonon mode (rad Hz).
    ω = 2π * 1e12 * phonon_mode_freq 

    # The static dielectric constant. Calculated here instead of inputted so that ionic modes are properly normalised.
    ϵ_static = ϵ_total + ϵ_optic

    # The contribution to the electron-phonon parameter from the currrent phonon mode. 1 / (4π ϵ_0) is the dielectric normalisation.
    α_j  = (m_eff * Ry / (ħ * ω))^(1 / 2) * ϵ_ionic / (4π * ϵ_0) / (ϵ_optic * ϵ_static)

    return α_j
end

"""
----------------------------------------------------------------------
Multiple Branch Polaron Free Energy
----------------------------------------------------------------------

Calculate the polaron free energy, generalised from Osaka's expression to the case where multiple phonon modes are present in the material.
"""

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
function κ_i(i, v, w)
    κ = v[i]^2 - w[i]^2
    κ *= prod(j != i ? (v[j]^2 - w[i]^2) / (w[j]^2 - w[i]^2) : 1.0 for j in eachindex(v))
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
function h_i(i, v, w)
    h = v[i]^2 - w[i]^2
    h *= prod(j != i ? (w[j]^2 - v[i]^2) / (v[j]^2 - v[i]^2) : 1.0 for j in eachindex(v))
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
function C_ij(i, j, v, w)
    C = w[i] * κ_i(i, v, w) * h_i(j, v, w) / (4 * (v[j]^2 - w[i]^2))
    return C
end

"""
    D_j(τ, β, v, w)

Calculates the recoil function (a generalisation of D(u) in Eqn. (35c) in FHIP 1962) that approximates the exact influence (recoil effects) of the phonon bath on the electron with the influence of the fictitious masses attached by springs to the electron. It appears in the exponent of the intermediate scattering function.

# Arguments
- `τ::Float64`: is the imaginary time variable.
- `β::Float64`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.

See FHIP 1962: https://doi.org/10.1103/PhysRev.127.1004.
"""
function D_j(τ, β, v, w)
    D = τ * (1 - τ / β) + sum((h_i(i, v, w) / v[i]^2) * ((1 + exp(-v[i] * β) - exp(-v[i] * τ) - exp(v[i] * (τ - β))) / (v[i] * (1 - exp(-v[i] * β))) - τ * (1 - τ / β)) for i in eachindex(v))
    return D
end

"""
    D_j(τ, v, w)

Calculates the recoil function at zero-temperature.

# Arguments
- `τ::Float64`: is the imaginary time variable.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.

See also ['D_j'](@ref).
"""
function D_j(τ, v, w)
    D = τ + sum((h_i(i, v, w) / v[i]^2) * ((1 - exp(-v[i] * τ)) / v[i] - τ) for i in eachindex(v))
    return D
end

"""
    B_j(α, β, v, w)

Generalisation of the B function from Eqn. (62c) in Hellwarth et al. 1999. This is the expected value of the exact action <S_j> taken w.r.t trial action, given for the 'jth' phonon mode.

Required for calculating the polaron free energy.

# Arguments
- `α::Float64`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
- `β::Float64`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.

See Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.

See also ['B'](@ref).
"""
function B_j(α, β, v, w)
    B_integrand(τ) = cosh(τ - β / 2) / sqrt(abs(D_j(τ, β, v, w)))
    B = α / (√π * sinh(β / 2)) * quadgk(τ -> B_integrand(τ), 0.0, β / 2.0)[1]
    return B
end

"""
    B_j(α, v, w; rtol = 1e-3)

Calculates `B_j(α, β, v, w)` but at zero-temperature, `β = Inf`.

# Arguments
- `α::Float64`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.

See also [`B_j`](@ref).
"""
function B_j(α, v, w)
    B_integrand(τ) = exp(-abs(τ)) / sqrt(abs(D_j(abs(τ), v, w)))
    B = α / √π * quadgk(τ -> B_integrand(τ), 0.0, Inf64)[1]
    return B
end

"""
    C_j(β, v, w, n)

Generalisation of the C function from Eqn. (62e) in Hellwarth et al. 1999. This is the expected value of the trial action <S_0> taken w.r.t trial action.

Required for calculating the polaron free energy.

# Arguments
- `β::Float64`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `n::Integer`: is the number of phonon modes.

See Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.

See also ['C'](@ref).
"""
function C_j(β, v, w, n)
    # Sum over the contributions from each fictitious mass.
    s = sum(C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2)  - 2 / (β * v[j])) for i in eachindex(v), j in eachindex(w))

    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    return 3 * s / n
end

"""
    C_j(v, w, n)

Calculates `C_j(β, v, w, n)` but at zero-temperature, `β = Inf`.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `n::Integer`: is the number of phonon modes.

See also [`C_j`](@ref).
"""
function C_j(v, w, n)
    # Sum over the contributions from each fictitious mass.
    s = sum(C_ij(i, j, v, w) / (v[j] * w[i]) for i in eachindex(v), j in eachindex(w))
    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    return 3 * s / n
end

"""
    A_j(β, v, w, n)

Generalisation of the A function from Eqn. (62b) in Hellwarth et al. 1999. This is the Helmholtz free energy of the trial model.

Required for calculating the polaron free energy.

# Arguments
- `β::Float64`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `n::Integer`: is the number of phonon modes.

See Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.

See also ['A'](@ref).
"""
function A_j(β, v, w, n)
    # Sum over the contributions from each fictitious mass.
    s = -log(2π * β) / 2 + sum(v[i] == w[i] ? 0.0 : 
    log(v[i]) -  log(w[i]) - β / 2 * (v[i] - w[i]) - log(1 - exp(-v[i] * β)) + log(1 - exp(-w[i] * β))
    for i in eachindex(v))
    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    3 / β * s / n
end

"""
    A_j(v, w, n)

Calculates `A_j(β, v, w, n)` but at zero-temperature, `β = Inf`.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `n::Integer`: is the number of phonon modes.

See also [`A_j`](@ref).
"""
function A_j(v, w, n)
    s = sum(v .- w)
    return -3 * s / (2 * n)
end

"""
    multi_F(v, w, α, β; ω = 1.0, T = nothing, verbose = false)

Calculates the Helmholtz free energy of the polaron for a material with multiple phonon branches. 
    
This generalises the Osaka 1959 (below Eqn. (22)) and Hellwarth. et al 1999 (Eqn. (62a)) free energy expressions.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `α::Float64`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
- `β::Float64`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
- `T`: is a token used by `make_polaron()` to keep track of the temperature for printing during a calculation. Do not alter.
- `verbose`: is used by `make_polaron()` to specify whether or not to print. Ignore.

See Osaka, Y. (1959): https://doi.org/10.1143/ptp.22.437 and Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.

See also [`F`](@ref).
"""
function multi_F(v, w, α, β; ω = 1.0, T = nothing, verbose = false)

    # Add contribution to the total free energy from the phonon mode.
    F = sum(-(B_j(α[j], β[j], v, w) + C_j(β[j], v, w, length(ω)) + A_j(β[j], v, w, length(ω))) * ω[j] for j in eachindex(ω))

    # Print the free energy.
    if verbose
        println("\e[2K", "Process: $(count) / $processes ($(round.(count / processes * 100, digits = 1)) %) | T = $(round.(T, digits = 3)) | F = $(round.(F, digits = 3))")
        print("\033[F") 
        
        global count += 1
    end
        
    # Free energy in units of meV
    return F
end

"""
    multi_F(v, w, α; ω = 1.0, verbose = false)

Calculates the zero-temperature ground-state energy of the polaron for a material with multiple phonon branches. Similar to `multi_F(v, w, α, β)` but with `β = Inf`. Generalises Eqn. (33) in Feynman 1955.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `α::Float64`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
- `verbose`: is used by `make_polaron()` to specify whether or not to print. Ignore.  

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.

See also [`multi_F`](@ref).
"""
function multi_F(v, w, α; ω = 1.0, verbose = false)

    # Add contribution to the total free energy from the phonon mode.
	F = sum(-(B_j(α[j], v, w) + C_j(v, w, length(ω)) + A_j(v, w, length(ω))) * ω[j] for j in eachindex(ω))

    # Print the free energy.
    if verbose
        println("\e[2K", "Process: $(count) / $processes ($(round.(count / processes * 100, digits = 1)) %) | T = 0.0 | F = $(round.(F, digits = 3))")
        print("\033[F") 
        global count += 1  
    end
        
    # Free energy in units of meV
    return F
end

"""
    var_params(α, β; v = 0.0, w = 0.0, ω = 1.0, N = 1, show_trace = false, T = nothing, verbose = false)

Minimises the multiple phonon mode free energy function for a set of vₚ and wₚ variational parameters. The variational parameters follow the inequality: v₁ > w₁ > v₂ > w₂ > ... > vₙ > wₙ. Generalises `feynmanvw` to multiple variational parameters.

# Arguments
- `α::Vector{Float64}`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
- `β::Vector{Float64}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.
- `v::Float64, w::Float64`: determines if the function should start with a random initial set of variational parameters (v, w = 0.0) or a given set of variational parameter values.
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
- `N::Integer`: specifies the number of variational parameter pairs, v_p and w_p, to use in minimising the free energy.
- `show_trace::Bool`: shows the optimsation trace from `Optim.jl`.
- `T`: is a token used by `make_polaron()` to keep track of the temperature for printing during a calculation. Do not alter.
- `verbose`: is used by `make_polaron()` to specify whether or not to print. Ignore.

See also [`multi_F`](@ref), [`feynmanvw`](@ref).
"""
function var_params(α, β; v = 0.0, w = 0.0, ω = 1.0, N = 1, show_trace = false, T = nothing, verbose = false) # N number of v and w params

    if N != length(v) != length(w)
        return error("The number of variational parameters v & w must be equal to N.")
    end

    # Use a random set of N initial v and w values.
    if v == 0.0 || w == 0.0
		# Intial guess for v and w parameters.
    	initial = [x for x in 1.0:(2.0 * N)] # initial guess around 4 and ≥ 1.
	else
        Δv = v .- w
        initial = vcat(Δv, w)
    end

    # Limits of the optimisation.
    lower = fill(0.0, 2 * N)
    upper = fill(Inf64, 2 * N)

	# The multiple phonon mode free energy function to minimise.
	f(x) = multi_F([x[2 * n - 1] for n in 1:N] .+ [x[2 * n] for n in 1:N], [x[2 * n] for n in 1:N], α, β; ω = ω)

	# Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
	solution = Optim.optimize(
		Optim.OnceDifferentiable(f, initial; autodiff = :forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
		Optim.Options(show_trace = show_trace), # Set time limit for asymptotic convergence if needed.
	)

	# Extract the v and w parameters that minimised the free energy.
	var_params = Optim.minimizer(solution)

	# Separate the v and w parameters into one-dimensional arrays (vectors).
	Δv = [var_params[2 * n - 1] for n in 1:N]
	w = [var_params[2 * n] for n in 1:N]

    # if Optim.converged(solution) == false
    #     @warn "Failed to converge T = $T K variational solution. v = $(Δv .+ w), w = $w."
    # end

    # Print the variational parameters that minimised the free energy.
	if verbose
        println("\e[2K", "Process: $(count) / $processes ($(round.(count / processes * 100, digits = 1)) %) | T = $(round.(T, digits = 3)) | v = $(round.(Δv .+ w, digits = 3)) | w = $(round.(w, digits = 3))")
        print("\033[F")   

        global count += 1
    end

    # Return the variational parameters that minimised the free energy.
    return (Δv .+ w, w)
end

"""
    var_params(α; v = 0.0, w = 0.0, ω = 1.0, N = 1, show_trace = false, verbose = false)

Minimises the multiple phonon mode free energy function for a set of vₚ and wₚ variational parameters at zero-temperature. Similar to `var_params(α, β)` but with `β = Inf`.

# Arguments
- `α::Vector{Float64}`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
- `v::Float64, w::Float64`: determines if the function should start with a random initial set of variational parameters (v, w = 0.0) or a given set of variational parameter values.
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
- `N::Integer`: specifies the number of variational parameter pairs, v_p and w_p, to use in minimising the free energy.
- `show_trace::Bool`: shows the optimsation trace from `Optim.jl`.
- `verbose`: is used by `make_polaron()` to specify whether or not to print. Ignore.

See also [`multi_F`](@ref), [`feynmanvw`](@ref), [`var_param`](@ref).
"""
function var_params(α; v = 0.0, w = 0.0, ω = 1.0, N = 1, show_trace = false, verbose = false) # N number of v and w params
 
    if N != length(v) != length(w)
        return error("The number of variational parameters v & w must be equal to N.")
    end

    # Use a random set of N initial v and w values.
    if v == 0.0 || w == 0.0
		# Intial guess for v and w parameters.
    	initial = [x for x in 1.0:(2.0 * N)] # initial guess around 4 and ≥ 1.
	else
        Δv = v .- w
        initial = vcat(Δv .+ rtol, w)
    end

    # Limits of the optimisation.
    lower = fill(0.0, 2 * N)
    upper = fill(100.0, 2 * N)

	# The multiple phonon mode free energy function to minimise.
	f(x) = multi_F([x[2 * n - 1] for n in 1:N] .+ [x[2 * n] for n in 1:N], [x[2 * n] for n in 1:N], α; ω = ω)

	# Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
	solution = Optim.optimize(
		Optim.OnceDifferentiable(f, initial; autodiff = :forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
		Optim.Options(show_trace = show_trace), # Set time limit for asymptotic convergence if needed.
	)

	# Extract the v and w parameters that minimised the free energy.
	var_params = Optim.minimizer(solution)

	# Separate the v and w parameters into one-dimensional arrays (vectors).
	Δv = [var_params[2 * n - 1] for n in 1:N]
	w = [var_params[2 * n] for n in 1:N]

    # if Optim.converged(solution) == false
    #     @warn "Failed to converge T = 0 K variational solution. v = $(Δv .+ w), w = $w."
    # end

	# Print the variational parameters that minimised the free energy.
	if verbose
        println("\e[2K", "Process: $(count) / $processes ($(round.(count / processes * 100, digits = 1)) %) | T = 0.0 | v = $(round.(Δv .+ w, digits = 3)) | w = $(round.(w, digits = 3))")
        print("\033[F")    

        global count += 1
    end

    # Return the variational parameters that minimised the free energy.
    return (Δv .+ w, w)
end

