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

feynmanvw(αωβ...) = vw_variation((v,w)->frohlich_energy(v, w, αωβ...))

feynmanvw(α) = feynmanvw(α, 1)

