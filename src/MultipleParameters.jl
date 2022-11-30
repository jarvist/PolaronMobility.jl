# MultipleParameters.jl
#  Extending the Feynman theory to multiple variational parameters

"""
----------------------------------------------------------------------
Multiple Parameter Polaron Free Energy
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
    B = α / (√π * sinh(β / 2)) * quadgk(τ -> B_integrand(τ * β / 2), 0.0, 1.0)[1] * β / 2
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
    s = sum(C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2) - 2 / (β * v[j])) for i in eachindex(v), j in eachindex(w))

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
function multi_param_F(v, w, α, β; ω = 1.0)

    # Add contribution to the total free energy from the phonon mode.
    F = sum(-(B_j(α[j], β[j], v, w) + C_j(β[j], v, w, length(ω)) + A_j(β[j], v, w, length(ω))) * ω[j] for j in eachindex(ω))

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
function multi_param_F(v, w, α; ω = 1.0)

    # Add contribution to the total free energy from the phonon mode.
	F = sum(-(B_j(α[j], v, w) + C_j(v, w, length(ω)) + A_j(v, w, length(ω))) * ω[j] for j in eachindex(ω))

    # Free energy in units of meV
    return F
end

"""
    multi_param_feynmanvw(α, β; v = 0.0, w = 0.0, ω = 1.0, N = 1, show_trace = false, T = nothing, verbose = false)

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
function multi_param_feynmanvw(α...; v = 0.0, w = 0.0, ω = 1.0, N = 1, show_trace = false) # N number of v and w params

    if N != length(v) != length(w)
        return error("The number of variational parameters v & w must be equal to N.")
    end

    # Use a random set of N initial v and w values.
    if v == 0.0 || w == 0.0
        # Intial guess for v and w parameters.
        initial = [x for x in 1.0:(2.0*N)] # initial guess around 4 and ≥ 1.
    else
        Δv = v .- w
        initial = vcat(Δv, w)
    end

    # Limits of the optimisation.
    lower = fill(0.0, 2 * N)
    upper = fill(Inf64, 2 * N)

	# The multiple phonon mode free energy function to minimise.
	f(x) = multi_F([x[2 * n - 1] for n in 1:N] .+ [x[2 * n] for n in 1:N], [x[2 * n] for n in 1:N], α...; ω = ω)

    # Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff=:forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
        Optim.Options(show_trace=show_trace), # Set time limit for asymptotic convergence if needed.
    )

    # Extract the v and w parameters that minimised the free energy.
    var_params = Optim.minimizer(solution)
    energy = Optim.minimum(solution)

    # Separate the v and w parameters into one-dimensional arrays (vectors).
    Δv = [var_params[2*n-1] for n in 1:N]
    w = [var_params[2*n] for n in 1:N]

    if Optim.converged(solution) == false
        @warn "Failed to converge. v = $(Δv .+ w), w = $w."
    end

    # Return the variational parameters that minimised the free energy.
    return (Δv .+ w, w, energy)
end


