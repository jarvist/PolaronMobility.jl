# HolsteinTheory.jl
# A work in progress.

# Compared to Hellwarth's presentation of the Feynman variational approach for the Frohlich
# Hamiltonian, then only thing that changes with the Holstein Hamiltonian is in the
# calculation of the 'B' component (the coupled electron-phonon component).

function holstein_B(v, w, α, ω; dims=3)

    G_polaron(τ) = ((v^2 - w^2) / v^3 * (1 - exp(-v * ω * τ)) + w^2 / v^2 * ω * τ + eps(Float64)) * dims / ω

    G_phonon(τ) = exp(-ω * τ)

    integrand(τ) = G_phonon(τ) * erf(π * sqrt(G_polaron(τ)))^dims / G_polaron(τ)^(dims/2)

    integral = quadgk(τ -> integrand(τ), 0, Inf)[1]

    return 2 * dims * ω * α * (2π)^(-dims/2) * integral
end

function holstein_B(v, w, α, ω, β; dims=3)

    G_polaron(τ) = ((v^2 - w^2) / v^3 * (1 - exp(-v * ω * τ)) * (1 - exp(-v * ω * (β - τ))) / (1 - exp(-v * β * ω)) + w^2 / v^2 * ω * τ * (1 - τ / β) + eps(Float64)) * dims / ω

    G_phonon(τ) = exp(-ω * τ) / (1 - exp(-β * ω)) + exp(ω * τ) / (exp(β * ω) - 1)

    integrand(τ) = G_phonon(τ) * erf(π * sqrt(G_polaron(τ)))^dims / G_polaron(τ)^(dims/2)

    integral = quadgk(τ -> integrand(τ), 0, β / 2)[1]

    return 2 * dims * ω * α * (2π)^(-dims/2) * integral
end

function holstein_energy(v, w, α, ωβ...; dims=3)

	kinetic_energy = -dims * (A(v, w, ωβ...) + C(v, w, ωβ...)) / 4 - 2 * dims

	potential_energy = -holstein_B(v, w, α, ωβ...; dims=dims)

    total_energy = kinetic_energy + potential_energy

	return total_energy, kinetic_energy, potential_energy
end

function holsteinvw(v, w, α, ωβ...; dims=3, upperlimit=Inf)

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [upperlimit, upperlimit]

    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + eps(Float64), w]

    # Feynman 1955 athermal action 
    f(x) = holstein_energy((x[1] + x[2]), x[2], α, ωβ...; dims=dims)[1]

    # Use Optim to optimise v and w to minimise enthalpy.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff=:forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
    )

    # Get v and w values that minimise the free energy.
    Δv, w = Optim.minimizer(solution) # v=Δv+w
	total_energy, kinetic_energy, potential_energy = holstein_energy(Δv + w, w, α, ωβ...; dims=dims)

	if Optim.converged(solution) == false
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $total_energy"
    end

    # Return variational parameters that minimise the free energy.
    return Δv + w, w, total_energy, kinetic_energy, potential_energy
end

function holstein_memory_function(v, w, α, ω, β, Ω; dims=3)

    n(x) = 1 / (exp(x * β) - 1)

    dim_lattice = 2 * dims # 'z'. Linear / Square / Cubic etc.

    length_scale_squared = dim_lattice / ω

    g_squared = length_scale_squared * α

    fictitious_mass = v^2 / w^2 

    fictitious_radius = (v^2 - w^2) / (2 * v^3) 

    G_phonon(t) = (1 + n(ω)) * exp(-im * ω * t) + n(ω) * exp(im * ω * t)

    G_polaron(t) = ((2 * fictitious_mass)^(-1) * (-im * t * ω + t^2 * ω / β) + fictitious_radius * (1 - exp(im * v * t * ω) + 4 * n(v * ω) * sin(v * t * ω / 2)^2)) * length_scale_squared

    S(t) = dims * length_scale_squared / (2π)^dims * G_phonon(t) * (√π / 2 * erf(π * √G_polaron(t)) / (G_polaron(t))^(3/2) - π * exp(-π^2 * G_polaron(t)) / G_polaron(t)) * (√π * erf(π * √G_polaron(t)) / √G_polaron(t))^(dims-1)

    integrand(t) = (1 - exp(im * Ω * ω * t)) * imag(S(t)) / Ω

    integral = quadgk(t -> integrand(t), 0, 1e4)[1]

    return 2 / 3 * g_squared * integral
end

function holstein_memory_function(v, w, α, ω, β; dims=3)

    n(x) = 1 / (exp(x * β) - 1)

    dim_lattice = 2 * dims # 'z'. Linear / Square / Cubic etc.

    length_scale_squared = dim_lattice / ω

    g_squared = length_scale_squared * α

    fictitious_mass = v^2 / w^2 

    fictitious_radius = (v^2 - w^2) / (2 * v^3) 

    G_phonon(t) = (1 + n(ω)) * exp(-im * ω * t) + n(ω) * exp(im * ω * t)

    G_polaron(t) = ((2 * fictitious_mass)^(-1) * (-im * t * ω + t^2 * ω / β) + fictitious_radius * (1 - exp(im * v * t * ω) + 4 * n(v * ω) * sin(v * t * ω / 2)^2)) * length_scale_squared 

    S(t) = dims * length_scale_squared / (2π)^dims * G_phonon(t) * (√π / 2 * erf(π * √G_polaron(t)) / (G_polaron(t))^(3/2) - π * exp(-π^2 * G_polaron(t)) / G_polaron(t)) * (√π * erf(π * √G_polaron(t)) / √G_polaron(t))^(dims-1)

    integrand(t) = -im * t * imag(S(t))

    integral = quadgk(t -> integrand(t), 0, 1e4)[1]

    return 2 / 3 * g_squared * integral
end
