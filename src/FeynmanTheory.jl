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
fF(τ, v, w) = (abs(w^2 * τ + (v^2 - w^2) / v * (1 - exp(-v*τ))))^(-0.5) * exp(-τ)


"""
    AF(v, w, α)

Integral of Eqn. (31) in Feynman 1955. Part of the overall ground-state energy expression.

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
AF(v, w, α) = π^(-0.5) * α * v * quadgk(τ -> fF(τ, v, w), 0, Inf)[1]

"""
    F(τ, v, w)

Ground state energy expression. Eqn. (33) in Feynman 1955.

See Feynman 1955: http://dx.doi.org/10.1103/PhysRev.97.660.
"""
F(v, w, α) = (3 / (4 * v)) * (v - w)^2 - AF(v, w, α)

# Let's wrap the Feynman athermal variation approximation in a simple function

"""
    feynmanvw(α; v = 0.0, w = 0.0)

Calculate v and w variational polaron parameters, for the supplied α Frohlich coupling. Returns (v, w).

This version uses the original athermal action (Feynman 1955).
"""
function feynmanvw(α; v = 3.0, w = 3.0) # v, w defaults

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [Inf, Inf]
    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + 0.01, w]

    # Feynman 1955 athermal action 
    f(x) = F(x[1] + x[2], x[2], α)

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

    # Return variational parameters that minimise the free energy.
    return Δv + w, w
end

# Hellwarth et al. 1999 PRB - Part IV; T-dep of the Feynman variation parameter

# In Julia we have 'Multiple dispatch', so let's just construct the free
# energies (temperature-dependent) with the same name as before, but with the thermodynamic beta where required.

# Define Osaka's free-energies (Hellwarth1999 version) as Julia functions.

"""
    A(v, w, β)

Hellwarth's A expression from Eqn. (62b) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
A(v, w, β) = 3 / β * (log(v / w) - 1/2 * log(2π * β) - log(sinh(v * β / 2) / sinh(w * β / 2)))

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
B(v, w, β, α) = α * v / (sqrt(π) * (exp(β) - 1)) * quadgk(x -> f(x, v, w,β), 0, β/2)[1]

"""
    C(v, w, β)

Hellwarth's C expression from Eqn. (62e) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
C(v, w, β) = 3/4 * (v^2 - w^2) / v * (coth(v * β / 2) - 2 / (v * β))
# 62a
"""
    F(v, w, β, α)

Hellwarth's total free energy expression from Eqn. (62a) in Hellwarth et al. 1999 PRB.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
F(v, w, β, α) = -(A(v, w, β) + B(v, w, β, α) + C(v, w, β))

# Variational optimsation methods for finding the polaron ground state at finite temperature.

"""
    feynmanvw(α, β; v = 0.0, w = 0.0)

Calculate v and w variational polaron parameters, for the supplied α Frohlich coupling, and inverse reduced temperature β. Returns (v, w).

This version uses the Osaka thermal action symmetrised for computation.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.

# Examples
```jldoctest
juila> v, w = feynmanvw(2.39, 0.36, v = 3.0, w = 1.0)
(19.86, 16.96)
```
"""
function feynmanvw(α, β; v = 3.0, w = 3.0) # v, w defaults

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [Inf, Inf]
    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + 0.01, w]

    # Feynman 1955 athermal action 
    f(x) = F(x[1] + x[2], x[2], β, α)

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

    # Return variational parameters that minimise the free energy.
    return Δv + w, w
end

