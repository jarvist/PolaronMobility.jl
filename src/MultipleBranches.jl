# MultipleBranches.jl
# #  Extending the Feynman theory to multiple phonon branches

# """
#     frohlichPartial((f, ϵ_mode); ϵ_o, ϵ_s, meff)

# Calculate a (partial) dielectric electron-phonon coupling element.

# # Arguments
# - `f`: frequency of mode in THz.
# - `ϵ_mode`: this mode's contribution to dielectric.
# - `ϵ_o`: optical dielectric.
# -  `ϵ_s`: total static dielectric contribution.
# """
# function frohlichPartial((f, ϵ_mode); ϵ_o, ϵ_s, meff)
#     ω = f * 1e12 * 2π
#     α = 1 / (4π * ϵ_0) * ϵ_mode / (ϵ_o * ϵ_s) * (q^2 / ħ) * sqrt(meff * me / (2 * ω * ħ))
#     return α
# end

# # deprecated signature wrapped via multiple dispatch
# frohlichPartial(ϵ_o, ϵ_s, ϵ_mode, f, meff) = frohlichPartial((f, ϵ_mode), ϵ_o=ϵ_o, ϵ_s=ϵ_s, meff=meff)

# rows(M::Matrix) = map(x -> reshape(getindex(M, x, :), :, size(M)[2]), 1:size(M)[1])

# """
#     IRtoDielectric(IRmodes, volume)

# From absolute value of IR activities of phonon modes, generate a per-mode
# contribution to the low-frequency dielectric constant.

# IRmodes are tuples f, S with Frequency in THz; InfraRed activity in e^2 amu^-1.
# """
# function IRtoDielectric(IRmodes, volume)
#     ϵ = 0.0 #* q^2 * THz^-2 * amu^-1 * m^-3
#     for r in rows(IRmodes)
#         f, S = r # frequency in THz; activity in e^2 amu^-1
#         f = f * 1E12 #* THz
#         ω = 2π * f
#         S = S * q^2 / amu
#         ϵ_mode = S / ω^2 / volume
#         ϵ_mode /= 3 # take isotropic average = divide by 3
#         ϵ += ϵ_mode
#         println("Mode f= $f S= $S ϵ_mode = $(ϵ_mode / ϵ_0)")
#     end
#     println("Raw ionic dielectric contribution: $ϵ absolute $(ϵ / ϵ_0) relative")
#     return ϵ / ϵ_0
# end

# """
#     IRtoalpha(IR, volume)

# Calculates contribution to dielectric constant for each polar phonon mode, and thereby the Frohlich alpha contribution for this mode.
# """
# function IRtoalpha(IR; volume, ϵ_o, ϵ_s, meff)
#     ϵ = 0.0 #* q^2 * THz^-2 * amu^-1 * m^-3
#     α_sum = 0.0
#     for r in rows(IR)
#         f, S = r # frequency in THz; activity in e^2 amu^-1
#         ω = 2π * f * 1e12
#         S = S * q^2 / amu
#         ϵ_mode = S / ω^2 / volume
#         ϵ_mode /= 3 # take isotropic average = divide by 3
#         ϵ_mode /= ϵ_0 # reduced dielectric constant
#         ϵ += ϵ_mode
#         #println("Mode f= $f S= $S ϵ_mode = $(upreferred(ϵ_mode/u"ϵ0"))")
#         α = frohlichPartial(ϵ_o, ϵ_s, ϵ_mode, f, meff)
#         α_sum += α
#         if (α > 0.1)
#             println("Notable Mode f = $f, α_partial = $α.")
#         end
#     end
#     println("Sum alpha: $(α_sum)")
#     return α_sum
# end

# """
#     DieletricFromIRmode(IRmode)

# Calculate dielectric from an individual mode.

# IRmode is a tuple f, S with Frequency in THz; InfraRed activity in e^2 amu^-1.
# """
# function DielectricFromIRmode(IRmode; volume)
#     f, S = IRmode # Assumes Frequency in THz; InfraRed activity in e^2 amu^-1
#     ω = 2π * f * 1e12
#     S = S * q^2 / amu
#     ϵ_mode = S / ω^2 / volume
#     ϵ_mode /= 3 # take isotropic average = divide by 3
#     ϵ_mode /= ϵ_0 # reduced dielectric
#     return ϵ_mode
# end

# """
#     Hellwarth1999mobilityRHS((α, (v, w) ,f), effectivemass, T)

# Calculates the DC mobility using Hellwarth et al. 1999 Eqn. (2).

# See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
# """
# function Hellwarth1999mobilityRHS((α, (v, w), f), effectivemass, T)
#     mb = effectivemass * MassElectron
#     ω = f * 1e12 * 2π
#     βred = ħ * ω / (kB * T)

#     R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)
#     b = R * βred / sinh(βred * v / 2) # Feynman1962 version; page 1010, Eqn (47b)
#     a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))
#     k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)
#     K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

#     # Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
#     RHS = α / (3 * sqrt(π)) * βred^(5 / 2) / sinh(βred / 2) * (v^3 / w^3) * K
#     μ = RHS^(-1) * q / (ω * mb)

#     return 1 / μ
# end

# """
# ----------------------------------------------------------------------
# Multiple Branch Frohlich Alpha
# ----------------------------------------------------------------------

# Partial dielectric electron-phonon coupling parameter, decomposed into ionic dielectric contributions from each phonon mode of the matieral.  
# """

# """
#     ϵ_ionic_mode(phonon_mode_freq, ir_activity, volume)

# Calculate the ionic contribution to the dielectric function for a given phonon mode.

# # Arguments
# - `phonon_mode_freq::Float64`: is the frequency of the mode in THz.
# - `ir_activity::Float64`: is the infra-red activity of the mode in e²amu⁻¹.
# - `volume::Float64`: is the volume of the unit cell of the material in m³.
# """
# function ϵ_ionic_mode(phonon_mode_freq, ir_activity, volume) # single ionic mode

#     # Angular phonon frequency for the phonon mode (rad Hz)
#     ω_j = 2π * phonon_mode_freq * 1e12

#     # Dielectric contribution from a single ionic phonon mode
#     ϵ_mode = eV^2 * ir_activity / (3 * volume * ω_j^2 * amu)

#     # Normalise ionic dielectric contribution with 1 / (4π ϵ_0) (NB: the 4π has been pre-cancelled)
#     return ϵ_mode / ϵ_0
# end

# """
#     ϵ_total(freqs_and_ir_activity, volume)

# Calculate the total ionic contribution to the dielectric function from all phonon modes.

# # Arguments
# - `freqs_and_ir_activity::Matrix{Float64}`: is a matrix containeing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e²amu⁻¹) in the second column.
# - `volume::Float64`: is the volume of the unit cell of the material in m^3.
# """
# function ϵ_total(freqs_and_ir_activity, volume) # total ionic contribution to dielectric

#     # Extract phonon frequencies (THz)
#     phonon_freqs = freqs_and_ir_activity[:, 1]

#     # Extra infra-red activities (e^2 amu^-1)
#     ir_activity = freqs_and_ir_activity[:, 2]

#     # Sum over all ionic contribution from each phonon mode
#     total_ionic = 0.0

#     for t in eachindex(phonon_freqs)
#         total_ionic += ϵ_ionic_mode(phonon_freqs[t], ir_activity[t], volume)
#     end

#     return total_ionic
# end

# """
#     effective_freqs(freqs_and_ir_activity, num_var_params)

# Generates a matrix of effective phonon modes with frequencies and infra-red activities derived from a larger matrix using the Principal Component Analysis (PCA) method.

# # Arguments
# - `freqs_and_ir_activity::Matrix{Float64}`: is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e²amu⁻¹) in the second column.
# - `num_var_params::Integer`: is the number of effective modes required (which needs to be less than the number of modes in `freqs_and_ir_activity``).

# *** POSSIBLY REDUNDANT ***
# """
# function effective_freqs(freqs_and_ir_activity, num_var_params) # PCA Algorithm

#     # Check that the number of effective modes is less than the number of actual phonon modes.
#     if num_var_params >= size(freqs_and_ir_activity)[1]

#         println("The number of effective phonon modes has to be less than the total number of phonon modes.")

#     else

#         # Centralise data by subtracting the columnwise mean
#         standardized_matrix = freqs_and_ir_activity' .- mean(freqs_and_ir_activity', dims=2)

#         # Calculate the covariance matrix S' * S. Matrix size is (n - 1) x (n - 1) for number of params (here n = 2)
#         covariance_matrix = standardized_matrix' * standardized_matrix

#         # Extract eigenvectors of the covariance matrix
#         eigenvectors = eigvecs(covariance_matrix)

#         # Project the original data along the covariance matrix eigenvectors and undo the centralisation
#         reduced_matrix = standardized_matrix[:, 1:num_var_params] * eigenvectors[1:num_var_params, 1:num_var_params] *
#                          eigenvectors[1:num_var_params, 1:num_var_params]' .+ mean(freqs_and_ir_activity', dims=2)

#         # Resultant matrix is positive definite and transposed.
#         return abs.(reduced_matrix')
#     end
# end

# """
#     multi_frohlichalpha(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff)

# Calculates the partial dielectric electron-phonon coupling parameter for a given longitudinal optical phonon mode. 

# This decomposes the original Frohlich alpha coupling parameter (defined for a single phonon branch) into contributions from multiple phonon
# branches.

# # Arguments
# - `ϵ_optic::Float64`: is the optical dielectric constant of the material.
# - `ϵ_ionic::Float64`: is the ionic dielectric contribution from the phonon mode.
# - `ϵ_total::Float64`: is the total ionic dielectric contribution from all phonon modes of the material.
# - `phonon_mode_freq::Float64`: is the frequency of the phonon mode (THz).
# - `m_eff::Float64` is the band mass of the electron (in units of electron mass m_e).
# """
# function multi_frohlichalpha(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff)

#     # The Rydberg energy unit
#     Ry = eV^4 * me / (2 * ħ^2)

#     # Angular phonon frequency for the phonon mode (rad Hz).
#     ω = 2π * 1e12 * phonon_mode_freq

#     # The static dielectric constant. Calculated here instead of inputted so that ionic modes are properly normalised.
#     ϵ_static = ϵ_total + ϵ_optic

#     # The contribution to the electron-phonon parameter from the currrent phonon mode. 1 / (4π ϵ_0) is the dielectric normalisation.
#     α_j = (m_eff * Ry / (ħ * ω))^(1 / 2) * ϵ_ionic / (4π * ϵ_0) / (ϵ_optic * ϵ_static)

#     return α_j
# end

# """
# ----------------------------------------------------------------------
# Multiple Branch Polaron Free Energy
# ----------------------------------------------------------------------
# """

# """
#     multi_F(v, w, α, β; ω = 1.0, T = nothing, verbose = false)

# Calculates the Helmholtz free energy of the polaron for a material with multiple phonon branches. 
    
# This generalises the Osaka 1959 (below Eqn. (22)) and Hellwarth. et al 1999 (Eqn. (62a)) free energy expressions.

# # Arguments
# - `v::Vector{Float64}`: is a vector of the v variational parameters.
# - `w::Vector{Float64}`: is a vector of the w variational parameters.
# - `α::Float64`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
# - `β::Float64`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.
# - `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
# - `T`: is a token used by `make_polaron()` to keep track of the temperature for printing during a calculation. Do not alter.
# - `verbose`: is used by `make_polaron()` to specify whether or not to print. Ignore.

# See Osaka, Y. (1959): https://doi.org/10.1143/ptp.22.437 and Hellwarth, R. W., Biaggio, I. (1999): https://doi.org/10.1103/PhysRevB.60.299.

# See also [`F`](@ref).
# """
# function multi_F(v, w, α, β; ω = 1.0)

#     N_modes = length(ω)

#     # Add contribution to the total free energy from the phonon mode.
#     F = -sum((B.(v, w, β, α) .+ C.(v, w, β) ./ N_modes .+ A.(v, w, β) ./ N_modes) .* ω)

#     # Free energy in units of meV
#     return F
# end

# """
#     multi_F(v, w, α; ω = 1.0, verbose = false)

# Calculates the zero-temperature ground-state energy of the polaron for a material with multiple phonon branches. Similar to `multi_F(v, w, α, β)` but with `β = Inf`. Generalises Eqn. (33) in Feynman 1955.

# # Arguments
# - `v::Vector{Float64}`: is a vector of the v variational parameters.
# - `w::Vector{Float64}`: is a vector of the w variational parameters.
# - `α::Float64`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
# - `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
# - `verbose`: is used by `make_polaron()` to specify whether or not to print. Ignore.  

# See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.

# See also [`multi_F`](@ref).
# """
# function multi_F(v, w, α; ω = 1.0)

#     N_modes = length(ω)

#     # Add contribution to the total free energy from the phonon mode.
# 	F = sum(((3 ./ (4 .* v)) .* (v .- w) .^2 / N_modes .- AF.(v, w, α)) .* ω) 

#     # Free energy in units of meV
#     return F
# end

# """
#     extended_feynmanvw(α, β; v = 0.0, w = 0.0, ω = 1.0, N = 1, show_trace = false, T = nothing, verbose = false)

# Minimises the multiple phonon mode free energy function for a set of vₚ and wₚ variational parameters. The variational parameters follow the inequality: v₁ > w₁ > v₂ > w₂ > ... > vₙ > wₙ. Generalises `feynmanvw` to multiple variational parameters.

# # Arguments
# - `α::Vector{Float64}`: is the partial dielectric electron-phonon coupling parameter for the 'jth' phonon mode.  
# - `β::Vector{Float64}`: is the reduced thermodynamic temperature ħωⱼ/(kT) associated with the 'jth' phonon mode.
# - `v::Float64, w::Float64`: determines if the function should start with a random initial set of variational parameters (v, w = 0.0) or a given set of variational parameter values.
# - `ω::Union{Float64, Vector{Float64}}`: phonon mode frequencies (2π THz). Predefined as `ω = 1.0` for a single mode in polaron units.
# - `N::Integer`: specifies the number of variational parameter pairs, v_p and w_p, to use in minimising the free energy.
# - `show_trace::Bool`: shows the optimsation trace from `Optim.jl`.
# - `T`: is a token used by `make_polaron()` to keep track of the temperature for printing during a calculation. Do not alter.
# - `verbose`: is used by `make_polaron()` to specify whether or not to print. Ignore.

# See also [`multi_F`](@ref), [`feynmanvw`](@ref).
# """
# function extended_feynmanvw(α...; v = 4.0, w = 2.0, ω = 1.0, show_trace = false, N=1) # N number of v and w params

#     # Limits of the optimisation.
#     lower = [0.0, 0.0]
#     upper = [Inf64, Inf64]
#     Δv = v - w # defines a constraint, so that v>w
#     initial = [Δv + eps(Float64), w]

# 	# The multiple phonon mode free energy function to minimise.
# 	f(x) = multi_F(x[1] + x[2], x[2], α...; ω = ω)

#     # Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
#     solution = Optim.optimize(
#         Optim.OnceDifferentiable(f, initial; autodiff=:forward),
#         lower,
#         upper,
#         initial,
#         Fminbox(BFGS()),
#         Optim.Options(show_trace=show_trace), # Set time limit for asymptotic convergence if needed.
#     )

#     # Extract the v and w parameters that minimised the free energy.
#     Δv, w = Optim.minimizer(solution)
#     energy = Optim.minimum(solution)

#     # Separate the v and w parameters into one-dimensional arrays (vectors).

#     if Optim.converged(solution) == false
#         @warn "Failed to converge. v = $(Δv .+ w), w = $w."
#     end

#     # Return the variational parameters that minimised the free energy.
#     return (Δv + w, w, energy)
# end



