# FeynmanTheory.jl

"""
    frohlichalpha(ε_Inf, ε_S, freq, m_eff)

Calculates the Frohlich alpha parameter, for a given dielectric constant, frequency (f) of phonon in Hertz, and effective mass (in units of the bare electron mass).

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
function frohlichalpha(ϵ_optic, ϵ_static, freq, m_eff)
    ω = 2π * freq * 1e12 # frequency to angular velocity
    # Note: we need to add a 4*pi factor to the permitivity of freespace.
    # This gives numeric agreement with literature values.  This is required as
    # the contemporary 1950s and 1960s literature implicitly used atomic units,
    # where the electric constant ^-1 has this factor baked in, k_e=1/(4πϵ_0).
    α = 0.5 / (4 * π * ϵ_0) *           # Units: m/F
        (1 / ϵ_optic - 1 / ϵ_static) *   # Units: none
        (eV^2 / (ħ * ω)) *               # Units: F
        sqrt(2 * me * m_eff * ω / ħ)    # Units: 1/m
    return α
end

# Athermal (Feynman 1955) model.
# Set up equations for the polaron free energy, which we will variationally improve upon.

"""
    fF(τ, v, w)

Integrand of Eqn. (31) in Feynman 1955. Part of the overall ground-state energy expression.

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
fF(τ, v, w) = (abs(w^2 * τ + (v^2 - w^2) / v * (1 - exp(-v * τ))))^(-0.5) * exp(-τ)


"""
    AF(v, w, α)

Integral of Eqn. (31) in Feynman 1955. Part of the overall ground-state energy expression.

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
AF(v, w, α) = π^(-0.5) * α * v * quadgk(τ -> fF(τ, v, w), 0.0, Inf64)[1]

"""
    F(τ, v, w)

Ground state energy expression. Eqn. (33) in Feynman 1955.

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
F(v, w, α) = (3 / (4 * v)) * (v - w)^2 - AF(v, w, α)

# Hellwarth et al. 1999 PRB - Part IV; T-dep of the Feynman variation parameter

# In Julia we have 'Multiple dispatch', so let's just construct the free
# energies (temperature-dependent) with the same name as before, but with the thermodynamic beta where required.

# Define Osaka's free-energies (Hellwarth1999 version) as Julia functions.

"""
    A(v, w, β)

Hellwarth's A expression from Eqn. (62b) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
A(v, w, β) = 3 / β * (log(v / w) - 1 / 2 * log(2π * β) - log(sinh(v * β / 2) / sinh(w * β / 2)))

"""
    Y(x, v, β)

Hellwarth's Y expression from Eqn. (62d) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression. Contained in denominator of the integrand of Eqn. (62c).

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
Y(x, v, β) = 1 / (1 - exp(-v * β)) * (1 + exp(-v * β) - exp(-v * x) - exp(v * (x - β)))

"""
    f(x, v, w, β)

Integrand of Hellwarth's B expression from Eqn. (62c) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
f(x, v, w, β) = (exp(β - x) + exp(x)) / sqrt(abs(w^2 * x * (1 - x / β) + Y(x, v, β) * (v^2 - w^2) / v))

"""
    B(v, w, β, α)

Hellwarth's B expression from Eqn. (62c) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
B(v, w, α, β) = α * v / (sqrt(π) * (exp(β) - 1)) * quadgk(x -> f(x, v, w, β), 0.0, β / 2)[1]

"""
    C(v, w, β)

Hellwarth's C expression from Eqn. (62e) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
C(v, w, β) = 3 / 4 * (v^2 - w^2) / v * (coth(v * β / 2) - 2 / (v * β))

# 62a
"""
    F(v, w, β, α)

Hellwarth's total *free energy* expression from Eqn. (62a) in Hellwarth et al. 1999 PRB.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
F(v, w, α, β) = -(A(v, w, β) + B(v, w, α, β) + C(v, w, β))

# Let's wrap the Feynman athermal variation approximation in a simple function
# Variational optimsation methods for finding the polaron ground state at finite or zero temperatures.

"""
    feynmanvw(αβ...; v = 3.0, w = 3.0)

Calculate v and w variational polaron parameters, for the supplied α Frohlich coupling, and inverse reduced temperature β. Returns (v, w, E).

This version uses the Osaka thermal action symmetrised for computation.

See Hellwarth et al. 1999: https://doi.org/10.1103/PhysRevB.60.299.

# Examples
```jldoctest
v, w, E = feynmanvw(2.39, 0.36, v = 3.0, w = 1.0)
```
"""
function feynmanvw(αβ...; v = 3.0, w = 3.0) # v, w defaults

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [Inf64, Inf64]
    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + eps(Float64), w]

    # Thermal action 
    f(x) = F(x[1] + x[2], x[2], αβ...)

    # Use Optim to optimise v and w to minimise enthalpy.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff = :forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
    )

    # Get v and w values that minimise the free energy.
    Δv, w = Optim.minimizer(solution) # v=Δv+w
    E = Optim.minimum(solution)

    # Return variational parameters that minimise the free energy.
    return Δv + w, w, E
end

# Extending the Feynman theory to multiple phonon branches

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
        standardized_matrix = freqs_and_ir_activity' .- mean(freqs_and_ir_activity', dims=2)

        # Calculate the covariance matrix S' * S. Matrix size is (n - 1) x (n - 1) for number of params (here n = 2)
        covariance_matrix = standardized_matrix' * standardized_matrix

        # Extract eigenvectors of the covariance matrix
        eigenvectors = eigvecs(covariance_matrix)

        # Project the original data along the covariance matrix eigenvectors and undo the centralisation
        reduced_matrix = standardized_matrix[:, 1:num_var_params] * eigenvectors[1:num_var_params, 1:num_var_params] *
                         eigenvectors[1:num_var_params, 1:num_var_params]' .+ mean(freqs_and_ir_activity', dims=2)

        # Resultant matrix is positive definite and transposed.
        return abs.(reduced_matrix')
    end
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
    ω = 2π * 1e12 * phonon_mode_freq

    # The static dielectric constant. Calculated here instead of inputted so that ionic modes are properly normalised.
    ϵ_static = ϵ_total + ϵ_optic

    # The contribution to the electron-phonon parameter from the currrent phonon mode. 1 / (4π ϵ_0) is the dielectric normalisation.
    α_j = (m_eff * Ry / (ħ * ω))^(1 / 2) * ϵ_ionic / (4π * ϵ_0) / (ϵ_optic * ϵ_static)

    return α_j
end

# Multiple Branch Polaron Free Energy

"""
    multi_F(v, w, α, β, ω)

Calculates the Helmholtz free energy of the polaron for a material with multiple phonon branches. 
This generalises the Osaka 1959 (below Eqn. (22)) and Hellwarth. et al 1999 (Eqn. (62a)) free energy expressions.

# Arguments
- `v::Float64`: the v variational parameter.
- `w::Float64`: the w variational parameter.
- `α::Vector{Float64}`: is a vector of the partial dielectric electron-phonon coupling parameter for each phonon mode.  
- `β::Vector{Float64}`: is a vector of the reduced thermodynamic temperature ħωⱼ/(kT) for each phonon mode.
- `ω::Vector{Float64}`: is a vector of phonon mode frequencies (units 2π THz). 

See Osaka, Y. (1959): https://doi.org/10.1143/ptp.22.437 and Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.

See also [`F`](@ref).
"""
function F(v, w, α::Vector, ω::Vector, β)

    N_modes = length(ω)

    # Add contribution to the total free energy from the phonon mode.
    F = -sum((B.(v, w, α, β) .+ C.(v, w, β) ./ N_modes .+ A.(v, w, β) ./ N_modes) .* ω)

    # Free energy in units of meV
    return F
end

"""
    F(v, w, α, ω)

Calculates the zero-temperature ground-state energy of the polaron for a material with multiple phonon branches. Similar to `F(v, w, α, ω, β)` but with `β = Inf`. Generalises Eqn. (33) in Feynman 1955.

# Arguments
- `v::Float64`: the v variational parameter.
- `w::Float64`: the w variational parameter.
- `α::Vector{Float64}`: is a vector of the partial dielectric electron-phonon coupling parameter for each phonon mode. 
- `ω::Vector{Float64}`: is a vector of phonon mode frequencies (units 2π THz).  

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.

See also [`multi_F`](@ref).
"""
function F(v, w, α::Vector, ω::Vector)

    N_modes = length(ω)

    # Add contribution to the total free energy from the phonon mode.
	F = sum(((3 ./ (4 .* v)) .* (v .- w) .^2 / N_modes .- AF.(v, w, α)) .* ω) 

    # Free energy in units of meV
    return F
end

"""
    feynmanvw(α, β, ω; v = 3.0, w = 3.0)

Minimises the multiple phonon mode free energy function for a single v and w variational parameter. Generalises `feynmanvw` to multiple phonon modes.

# Arguments
- `α::Vector{Float64}`: is a vector of the partial dielectric electron-phonon coupling parameter for each phonon mode.  
- `β::Vector{Float64}`: is a vector of the reduced thermodynamic temperature ħωⱼ/(kT) for each phonon mode.
- `ω::Vector{Float64}`: is a vector of phonon mode frequencies (units 2π THz). 
- `v::Float64, w::Float64`: specifies initial guess for the variational parameters (v, w = 3.0, the zero-coupling solution, by default).

See also [`F`](@ref), [`feynmanvw`](@ref).
"""
function feynmanvw(αβω::Vector{Float64}...; v = 3.0, w = 3.0) # N number of v and w params

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [Inf64, Inf64]
    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + eps(Float64), w]

	# The multiple phonon mode free energy function to minimise.
	f(x) = F(x[1] + x[2], x[2], αβω...)

    # Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff=:forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS())
    )

    # Extract the v and w parameters that minimised the free energy.
    Δv, w = Optim.minimizer(solution)

    # Extract the free energy E.
    E = Optim.minimum(solution)

    # Check if the minimisation converged.
    if Optim.converged(solution) == false
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E."
    end

    # Return the variational parameters `v` and `w` and the minimised the free energy `E`.
    return (Δv + w, w, E)
end


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
frohlichPartial(ϵ_o, ϵ_s, ϵ_mode, f, meff) = frohlichPartial((f, ϵ_mode), ϵ_o=ϵ_o, ϵ_s=ϵ_s, meff=meff)

rows(M::Matrix) = map(x -> reshape(getindex(M, x, :), :, size(M)[2]), 1:size(M)[1])

"""
    IRtoDielectric(IRmodes, volume)

From absolute value of IR activities of phonon modes, generate a per-mode
contribution to the low-frequency dielectric constant.

IRmodes are tuples f, S with Frequency in THz; InfraRed activity in e^2 amu^-1.
"""
function IRtoDielectric(IRmodes, volume)
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
