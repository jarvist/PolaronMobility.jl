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
    α = 1 / 2 / (4 * π * ϵ_0) *           # Units: m/F
        (1 / ϵ_optic - 1 / ϵ_static) *   # Units: none
        (eV^2 / (ħ * ω)) *               # Units: F
        sqrt(2 * me * m_eff * ω / ħ)    # Units: 1/m
    return α
end

# Athermal (Feynman 1955) model.
# Set up equations for the polaron free energy, which we will variationally improve upon.

"""
    B(τ, v, w)

Integrand of Eqn. (31) in Feynman 1955. Part of the overall ground-state energy expression.

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
B_integrand(τ, v, w) = (abs(w^2 * τ + (v^2 - w^2) / v * (1 - exp(-v * τ))))^(-1 / 2) * exp(-τ)

"""
    B(v, w, α)

Integral of Eqn. (31) in Feynman 1955. Part of the overall ground-state energy expression.

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
B(v, w, α) = π^(-1 / 2) * α * v * quadgk(τ -> B_integrand(τ, v, w), 0, Inf)[1]

A(v, w) = -3 * (v - w) / 2

C(v, w) = (3 / (4 * v)) * (v^2 - w^2)

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
Y(x, v, β) = 1 / (1 - exp(-v * β)) * (1 + exp(-v * β) - exp(-v * x) - exp(v * (x - β))) + eps(x)

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
B(v, w, α, β) = α * v / (sqrt(π) * (exp(β) - 1)) * quadgk(x -> f(x, v, w, β), 0, β / 2)[1]

"""
    C(v, w, β)

Hellwarth's C expression from Eqn. (62e) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
C(v, w, β) = 3 / 4 * (v^2 - w^2) / v * (coth(v * β / 2) - 2 / (v * β))

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
    total_ionic = 0

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

#  Extending the Feynman theory to multiple variational parameters
# Multiple Parameter Polaron Free Energy
# Calculate the polaron free energy, generalised from Osaka's expression to the case where multiple phonon modes are present in the material.

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
function κ_i(i, v::Vector, w::Vector)
    κ = v[i]^2 - w[i]^2
    κ *= prod(j != i ? (v[j]^2 - w[i]^2) / (w[j]^2 - w[i]^2) : 1 for j in eachindex(v))
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
function h_i(i, v::Vector, w::Vector)
    h = v[i]^2 - w[i]^2
    h *= prod(j != i ? (w[j]^2 - v[i]^2) / (v[j]^2 - v[i]^2) : 1 for j in eachindex(v))
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
function C_ij(i, j, v::Vector, w::Vector)
    C = w[i] * κ_i(i, v, w) * h_i(j, v, w) / (4 * (v[j]^2 - w[i]^2))
    return C
end

"""
    D(τ, v, w, β)

Calculates the recoil function (a generalisation of D(u) in Eqn. (35c) in FHIP 1962) that approximates the exact influence (recoil effects) of the phonon bath on the electron with the influence of the fictitious masses attached by springs to the electron. It appears in the exponent of the intermediate scattering function.

# Arguments
- `τ::Float64`: is the imaginary time variable.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `β::Union{Float64, Vector{Float64}}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.

See FHIP 1962: https://doi.org/10.1103/PhysRev.127.1004.
"""
function D(τ, v::Vector, w::Vector, β)
    return τ * (1 - τ / β) + sum((h_i(i, v, w) / v[i]^2) * ((1 + exp(-v[i] * β) - exp(-v[i] * τ) - exp(v[i] * (τ - β))) / (v[i] * (1 - exp(-v[i] * β))) - τ * (1 - τ / β)) for i in eachindex(v))
end

"""
    D(τ, v, w)

Calculates the recoil function at zero-temperature.

# Arguments
- `τ::Float64`: is the imaginary time variable.
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
"""
function D(τ, v::Vector, w::Vector)
    return τ + sum((h_i(i, v, w) / v[i]^2) * ((1 - exp(-v[i] * τ)) / v[i] - τ) for i in eachindex(v))
end

"""
    B(v, w, α, β)

Generalisation of the B function from Eqn. (62c) in Hellwarth et al. 1999. This is the expected value of the exact action <S_j> taken w.r.t trial action, given for the 'jth' phonon mode.

Required for calculating the polaron free energy.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `α::Union{Float64, Vector{Float64}}`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
- `β::Union{Float64, Vector{Float64}}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.

See Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.
"""
function B(v::Vector, w::Vector, α, β)
    B_integrand(τ) = cosh(τ - β / 2) / sqrt(abs(D(τ, v, w, β)))
    return α / (√π * sinh(β / 2)) * quadgk(τ -> B_integrand(τ * β / 2), 0.0, 1.0)[1] * β / 2
end

"""
    B(v, w, α; rtol = 1e-3)

Calculates `B(v, w, α, β)` but at zero-temperature, `β = Inf`.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `α::Union{Float64, Vector{Float64}}`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
"""
function B(v::Vector, w::Vector, α)
    B_integrand(τ) = exp(-abs(τ)) / sqrt(abs(D(abs(τ), v, w)))
    return α / √π * quadgk(τ -> B_integrand(τ), 0, Inf)[1]
end

"""
    C(v, w, β)

Generalisation of the C function from Eqn. (62e) in Hellwarth et al. 1999. This is the expected value of the trial action <S_0> taken w.r.t trial action.

Required for calculating the polaron free energy.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `β::Union{Float64, Vector{Float64}}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.


See Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.
"""
function C(v::Vector, w::Vector, β)
    # Sum over the contributions from each fictitious mass.
    s = sum(C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2) - 2 / (β * v[j])) for i in eachindex(v), j in eachindex(w))

    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    return 3 * s
end

"""
    C(v, w)

Calculates `C(v, w, β)` but at zero-temperature, `β = Inf`.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
"""
function C(v::Vector, w::Vector)
    # Sum over the contributions from each fictitious mass.
    s = sum(C_ij(i, j, v, w) / (v[j] * w[i]) for i in eachindex(v), j in eachindex(w))
    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    return 3 * s
end

"""
    A(v, w, β)

Generalisation of the A function from Eqn. (62b) in Hellwarth et al. 1999. This is the Helmholtz free energy of the trial model.

Required for calculating the polaron free energy.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `β::Union{Float64, Vector{Float64}}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.

See Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.
"""
function A(v::Vector, w::Vector, β)
    # Sum over the contributions from each fictitious mass.
    s = -log(2π * β) / 2 + sum(v[i] == w[i] ? 0 :
                               log(v[i]) - log(w[i]) - β / 2 * (v[i] - w[i]) - log(1 - exp(-v[i] * β)) + log(1 - exp(-w[i] * β))
                               for i in eachindex(v))
    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    3 / β * s
end

"""
    A(v, w, n)

Calculates `A(v, w, β)` but at zero-temperature, `β = Inf`.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
"""
function A(v::Vector, w::Vector)
    s = sum(v .- w)
    return -3 * s / 2
end

"""
    F(v, w, α, ω, β)

Calculates the Helmholtz free energy of the polaron for a material with multiple phonon branches. 
    
This generalises the Osaka 1959 (below Eqn. (22)) and Hellwarth. et al 1999 (Eqn. (62a)) free energy expressions.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `α::Union{Float64, Vector{Float64}}`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz).
- `β::Union{Float64, Vector{Float64}}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.

See Osaka, Y. (1959): https://doi.org/10.1143/ptp.22.437 and Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.
"""
function F(v, w, α, ω, β)

    # Insurance to avoid breaking the integrals with Infinite beta.
    if any(x -> x == Inf, β)
        return F(v, w, α, ω)
    end

    N_modes = length(ω)

    Ar = Vector{Any}(undef, N_modes)
    Br = Vector{Any}(undef, N_modes)
    Cr = Vector{Any}(undef, N_modes)
    energy = 0.0

    # Add contribution to the total free energy from the phonon mode.
    for j in eachindex(ω)
        Ar[j] = A(v, w, β[j]) / N_modes * ω[j]
        Br[j] = B(v, w, α[j], β[j]) * ω[j]
        Cr[j] = C(v, w, β[j]) / N_modes * ω[j]
        energy -= (Ar[j] + Br[j] + Cr[j])
    end

    # Free energy in units of meV
    return energy, sum(Ar), sum(Br), sum(Cr)
end

"""
    F(v, w, α, ω)

Calculates the zero-temperature ground-state energy of the polaron for a material with multiple phonon branches. Similar to `F(v, w, α, ω, β)` but with `β = Inf`. Generalises Eqn. (33) in Feynman 1955.

# Arguments
- `v::Vector{Float64}`: is a vector of the v variational parameters.
- `w::Vector{Float64}`: is a vector of the w variational parameters.
- `α::Union{Float64, Vector{Float64}}`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). 

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
function F(v, w, α, ω)

    N_modes = length(ω)

    Ar = A(v, w)
    Br = Vector{Any}(undef, N_modes)
    Cr = C(v, w)
    energy = -(Ar + Cr) * sum(ω) / N_modes

    # Add contribution to the total free energy from the phonon mode.
    for j in eachindex(ω)
        Br[j] = B(v, w, α[j]) * ω[j]
        energy -= Br[j]
    end

    # Free energy in units of meV
    return energy, Ar, sum(Br), Cr
end

"""
    feynmanvw(v, w, αωβ...; upper_limit = Inf64)

Minimises the multiple phonon mode free energy function for a set of vₚ and wₚ variational parameters. The variational parameters follow the inequality: v₁ > w₁ > v₂ > w₂ > ... > vₙ > wₙ. Generalises `feynmanvw` to multiple variational parameters.

# Arguments
- `v::Vector{Float64}`: vector of initial v parameters.
- `w::Vector{Float64}`: vector of initial w parameters.
- `α::Union{Float64, Vector{Float64}}`: is the partial dielectric electron-phonon coupling parameter for one or many phonon modes.  
- `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz) for one or many phonon modes.
- `β::Union{Float64, Vector{Float64}}`: is the reduced thermodynamic temperature ħωⱼ/(kT) for one or many phonon modes.

See also [`F`](@ref).
"""
function feynmanvw(v::Vector, w::Vector, αωβ...; upper_limit=1e6)

    if length(v) != length(w)
        return error("The number of variational parameters v & w must be equal.")
    end

    N_params = length(v)

    Δv = v .- w
    initial = vcat(Δv .+ repeat([eps(Float64)], N_params), w)

    # Limits of the optimisation.
    lower = fill(0.0, 2 * N_params)
    upper = fill(upper_limit, 2 * N_params)

    # The multiple phonon mode free energy function to minimise.
    f(x) = F([x[2*n-1] for n in 1:N_params] .+ [x[2*n] for n in 1:N_params], [x[2*n] for n in 1:N_params], αωβ...)[1]

    # Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff=:forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS())
    )

    # Extract the v and w parameters that minimised the free energy.
    var_params = Optim.minimizer(solution)

    # Separate the v and w parameters into one-dimensional arrays (vectors).
    Δv = [var_params[2*n-1] for n in 1:N_params]
    w = [var_params[2*n] for n in 1:N_params]
    E, A, B, C = F(Δv .+ w, w, αωβ...)

    # if Optim.converged(solution) == false
    #     @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    # end

    # Return the variational parameters that minimised the free energy.
    return Δv .+ w, w, E, A, B, C
end

function feynmanvw(v::Real, w::Real, αωβ...; upper_limit=1e6)

    Δv = v .- w
    initial = [Δv + eps(Float64), w]

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [upper_limit, upper_limit]

    # The multiple phonon mode free energy function to minimise.
    f(x) = F(x[1] .+ x[2], x[2], αωβ...)[1]

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
    E, A, B, C = F(Δv .+ w, w, αωβ...)

    # if Optim.converged(solution) == false
    #     @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    # end

    # Return the variational parameters that minimised the free energy.
    return Δv .+ w, w, E, A, B, C
end
