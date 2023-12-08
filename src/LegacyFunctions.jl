# LegacyFunctions.jl
# Some old functions that are still potentially useful and have use in some tests.

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
    f(x) = frohlich_energy([x[2*n-1] for n in 1:N_params] .+ [x[2*n] for n in 1:N_params], [x[2*n] for n in 1:N_params], αωβ...)[1]

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
    E, A, B, C = frohlich_energy(Δv .+ w, w, αωβ...)

    # if Optim.converged(solution) == false
    #     @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    # end

    # Return the variational parameters that minimised the free energy.
    return Δv .+ w, w, E, A, B, C
end

function feynmanvw(v::Real, w::Real, αωβ...; upper = [Inf, Inf], lower = [0, 0])
    
    Δv = v .- w
    initial = [Δv + eps(Float64), w] 

    # The multiple phonon mode free energy function to minimise.
    f(x) = frohlich_energy(x[1] .+ x[2], x[2], αωβ...)[1]

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
    E, A, B, C = frohlich_energy(Δv .+ w, w, αωβ...)

    if Optim.converged(solution) == false
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $E"
    end

    # Return the variational parameters that minimised the free energy.
    return Δv .+ w, w, E, A, B, C
end

feynmanvw(αωβ...) = feynmanvw(3.4, 2.6, αωβ...)

