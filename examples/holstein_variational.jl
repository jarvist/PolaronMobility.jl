### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 66d2d0e0-ee56-11ec-21f9-5d6ae25cff48
begin
    using Optim
    using QuadGK
    using SpecialFunctions
    using Plots
end

# ╔═╡ 9f6474e5-41f1-4b96-aeee-d299404a3718
md"""
## Feynman's variational theory for the (continious) Frohlich model.
"""

# ╔═╡ 9429ae6b-55e1-4fb2-9a88-ae7636c35bb5
begin
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
    AF(v, w, α) = π^(-0.5) * α * v * quadgk(τ -> fF(τ, v, w), 0.0, Inf)[1]

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
    function feynmanvw_old(α; v=3.0, w=3.0) # v, w defaults

        # Limits of the optimisation.
        lower = [0.0, 0.0]
        upper = [Inf, Inf]
        Δv = v - w # defines a constraint, so that v>w
        initial = [Δv + 0.01, w]

        # Feynman 1955 athermal action 
        f(x) = F(x[1] + x[2], x[2], α)

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

        # Return variational parameters that minimise the free energy.
        return Δv + w, w
    end
end

# ╔═╡ 1cf1e5cc-301c-4779-86cb-5b6436ecebf0
v, w = feynmanvw_old(2.39)

# ╔═╡ 2a609c87-4579-490e-99e7-0d0c25e11333
F(v, w, 2.39)

# ╔═╡ 7bf82f19-e6fc-4c3f-b7e8-a02bebf0d0dd
md"""
# Explicit k-space integral version of Frohlich model.
"""

# ╔═╡ 6cd58deb-fb04-44fc-a854-dce5b5cb8565
function D(τ, v, w)
    abs(w^2 / v^2 * τ + (v^2 - w^2) / v^3 * (1 - exp(-v * τ)))
end

# ╔═╡ f7a07dc5-ba45-4870-b94c-44ae5d85436d
md"""
### Test first that it reproduces the Feynman analytic result. (R = 0).
"""

# ╔═╡ 82f23e07-599f-4372-8e82-e6086ab1b40a
v, w

# ╔═╡ cc237871-08fd-4658-8fa5-8658d9394b3b
md"""
### Now try Holstein model with k-space integral.
"""

# ╔═╡ 907b2475-0749-4d85-b09e-6eed75d35e68
md"""
### Now the k-space integral for the Holstein can be done and is expressed in terms of an erf function. So lets do it as it will be quicker!
"""

# ╔═╡ 63bf4a80-ed9c-46fc-9ede-6d0b05fc2519
t_range = LinRange(0.0, (6π^2)^(1/3), 100)

# ╔═╡ dde4d55d-702a-4ced-9e37-13f89f78d3d5
pi^2

# ╔═╡ e7a4f46e-8fcc-48ef-9c5c-5790b42d7648
"""
Let's now plot the α dependence of Frohlich vs Holstein!
"""

# ╔═╡ 78936a9d-9567-4def-bf89-85be9e62e12d
α = 0.1:0.1:12

# ╔═╡ 8d726e23-5a40-4a10-8405-7c9775a0c9d0
frohlich_params = feynmanvw_old.(α)

# ╔═╡ 2cdfdcbc-a6e7-444c-a9a0-47529ddb091a
begin
    vf = [i[1] for i in frohlich_params]
    wf = [i[2] for i in frohlich_params]
end

# ╔═╡ 35ce9514-521c-4985-a85e-f236f67417c9
frohlich_energies = [F(vf[i], wf[i], α[i]) for i in 1:length(α)]

# ╔═╡ 3f5069df-89e2-4ac7-96e4-f3c79936fdd9
md"""
## Thermal (free energy) version of the theory...
"""

# ╔═╡ 9c1284f7-f503-4c4d-be6a-e4e04043246a
"""
    A(v, w, β)

Hellwarth's A expression from Eqn. (62b) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
A(v, w, β) = 3 / β * (log(v / w) - 1 / 2 * log(2π * β) - log(sinh(v * β / 2) / sinh(w * β / 2)))

# ╔═╡ b36f7997-f29c-4a92-bb00-9e82ae7d2a05
"""
    Y(x, v, β)

Hellwarth's Y expression from Eqn. (62d) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression. Contained in denominator of the integrand of Eqn. (62c).

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
Y(x, v, β) = 1 / (1 - exp(-v * β)) * (1 + exp(-v * β) - exp(-v * x) - exp(v * (x - β)))

# ╔═╡ d4e1ac9e-84d7-4809-ac1c-c24a9263db84
"""
    f(x, v, w, β)

Integrand of Hellwarth's B expression from Eqn. (62c) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
f_Fro(x, v, w, β) = (exp(β - x) + exp(x)) / sqrt(abs(w^2 * x * (1 - x / β) + Y(x, v, β) * (v^2 - w^2) / v))


# ╔═╡ ca1d9705-b49e-43f8-bd17-ab4a7e28b578
D(x, v, w, β) = w^2 / v^2 * x * (1 - x / β) + Y(x, v, β) * (v^2 - w^2) / v^3

# ╔═╡ 00fc8b3d-c485-43e2-92f0-451539ecf1e6
"""
 k_integral

R is cutoff for Holstei in k-space integral (in k-space spherical volume, R)

if R==0, Frohlich model (power should be 0, as the K^2 cancels with the 1/K^2 dependency in the el-ph matrix element)

Power is the K^Power
"""
function k_integral(t, v, w; power=0, R=0)
    if R == 0 # Frohlich
        return quadgk(k -> k^power * exp(-k^2 * D(t, v, w) / 2), 0.0, Inf)[1]
    else # Holstein, note limits on integral, R=k-space volume
        return quadgk(k -> k^power * exp(-k^2 * D(t, v, w) / 2), 0.0, R)[1] / 2
    end
end

# ╔═╡ dcf7f24a-6ef4-4613-9778-e2637f5e6c16
function Frohlich_el_integral(α, v, w; power=0, R=0)
    α * quadgk(t -> exp(-t) * k_integral(t, v, w; power=power, R=R), 0.0, Inf)[1] * √2 / pi
end

# ╔═╡ 6714b04e-9b26-42d0-89eb-f6c510cd5f56
function Frohlich_energy(α, v, w; power=0, R=0)
    3 * (v - w)^2 / (4 * v) - Frohlich_el_integral(α, v, w; power=power, R=R)
end

# ╔═╡ 429d0098-9d78-47e5-915b-cf16c7fce660
function var_params_k(α; v=3.0, w=3.0, power=0, R=0)
    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [Inf, Inf]
    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + 0.01, w]

    # Feynman 1955 athermal action 
    f(x) = Frohlich_energy(α, x[1] + x[2], x[2]; power=power, R=R)

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

    # Return variational parameters that minimise the free energy.
    return Δv + w, w
end

# ╔═╡ 71118d82-37fb-4504-b944-c3461de91233
v_k, w_k = var_params_k(2.39; v=v, w=w)

# ╔═╡ 80cd6d37-92d3-4757-a62e-490cffb4b63c
vh_k, wh_k = var_params_k(2.39; v=v, w=w, power=2, R=(6 * π^2)^(1 / 3))

# ╔═╡ 2be860db-52b0-4514-872f-db79dc454323
Frohlich_energy(2.39, v_k, w_k)

# ╔═╡ 410d234c-25fe-4752-926e-f6c5d263ef05
Frohlich_energy(2.39, vh_k, wh_k; power=2, R=(6 * π^2)^(1 / 3))

# ╔═╡ 14f66830-c111-4466-bce9-965b2e48d2d4
k_int = k_integral.(t_range, v, w; power=2, R=(6π^2)^(1/3))

# ╔═╡ 0b9cb699-eace-4f78-a547-0a897aa6c361
function momentum_cutoff(τ, k_c, v, w; V_0 = (2π)^3, r_p = sqrt(1/2), a = 1.0)
	if τ > 0.0
		V_0 / (2π)^3 * π / r_p^3 / D(τ, v, w)^(3/2) * (
	√π * erf(√D(τ, v, w) * k_c * r_p) - 2 * √D(τ, v, w) * k_c * r_p * exp(-k_c^2 * r_p^2 * D(τ, v, w)) 
	) / 4pi
	elseif τ == 0.0
		2pi^2
	end
end

# ╔═╡ 9c89b255-acb4-495c-9cb4-b71b05c13708
mom_cut = momentum_cutoff.(t_range, (6π^2)^(1/3), v, w; r_p = sqrt(1/2), a = 1.0)/2

# ╔═╡ eb401f35-2fb2-4432-aa9a-28d393be71fb
begin
	plot(t_range, k_int)
	plot!(t_range, mom_cut)
end

# ╔═╡ 079d7328-c2e1-4681-a4e3-e19a53841d25
function holstein_elph_integral(α, v, w; R=(6 * π^2)^(1 / 3))
    integrand(τ) = exp(-τ) * momentum_cutoff(τ, R, v, w)
    integral = α * quadgk(τ -> integrand(τ), 0.0, Inf)[1] / pi / sqrt(2)
    return integral
end

# ╔═╡ bc20477e-2ce3-4377-9af0-5d479256175a
function holstein_energy(α, v, w)
    3 * (v - w)^2 / (4 * v) - holstein_elph_integral(α, v, w)
end

# ╔═╡ 380c5c1b-01e1-455b-a757-46c87fc45285
function holstein_variational_solution(α; v=3.0, w=3.0)
    # Limits of the optimisation.
    lower = [0.0, 0.3]
    upper = [Inf, Inf]
    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + 0.01, w]

    # Feynman 1955 athermal action 
    f(x) = holstein_energy(α, x[1] + x[2], x[2])

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

    # Return variational parameters that minimise the free energy.
    return Δv + w, w
end

# ╔═╡ 4e99725c-0895-4bc3-97ba-70eb675c7f47
f_Hol(x, v, w, β) = (exp(β - x) + exp(x)) * (erf(sqrt(D(x, v, w, β) / 2) * (6 * π^2)^(1 / 3)) / (2 * D(x, v, w, β)^(3 / 2)) - (6 * π^2)^(1 / 3) * exp(-D(x, v, w, β) * (6 * π^2)^(2 / 3) / 2) / (sqrt(2π) * D(x, v, w, β)))

# ╔═╡ e8b8e4de-e548-447e-82c5-979124ce6cf3
"""
    B(v, w, β, α)

Hellwarth's B expression from Eqn. (62c) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
B_Fro(v, w, β, α) = α * v / (sqrt(π) * (exp(β) - 1)) * quadgk(x -> f_Fro(x, v, w, β), 0, β / 2)[1]


# ╔═╡ f87703c9-fa55-420a-b3cc-c6b947ef456a
B_Hol(v, w, β, α) = α / (sqrt(2) * pi * (exp(β) - 1)) * quadgk(x -> f_Hol(x, v, w, β), 0, β / 2)[1]

# ╔═╡ c1358b12-4f68-4485-8ee7-73ef672f3ed8
"""
    C(v, w, β)

Hellwarth's C expression from Eqn. (62e) in Hellwarth et al. 1999 PRB. Part of the overall free energy expression.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
C(v, w, β) = 3 / 4 * (v^2 - w^2) / v * (coth(v * β / 2) - 2 / (v * β))

# ╔═╡ fc699405-e1d6-48dd-aff3-a7ae7fb0114b
"""
    F(v, w, β, α)

Hellwarth's total free energy expression from Eqn. (62a) in Hellwarth et al. 1999 PRB.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
F_Fro(v, w, β, α) = -(A(v, w, β) + B_Fro(v, w, β, α) + C(v, w, β))

# ╔═╡ 9ed1bad6-9ac5-444f-8d02-974bdd955b09
F_Hol(v, w, β, α) = -(A(v, w, β) + B_Hol(v, w, β, α) + C(v, w, β))

# ╔═╡ 5d5bf742-89c2-45e6-81f6-3b22e90b7bd0
"""
Check that it roughly gives the Holstein ground-state energy (high beta)
"""

# ╔═╡ 2455abad-abf2-4b0a-bb3f-30d97de9a9ec
function feynmanvw(α, β; v=3.0, w=3.0) # v, w defaults

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [Inf, Inf]
    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + 0.01, w]

    # Feynman 1955 athermal action 
    f(x) = F_Fro(x[1] + x[2], x[2], β, α)

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

    # Return variational parameters that minimise the free energy.
    return Δv + w, w
end

# ╔═╡ 4d5e8c8e-2225-43fd-a934-0c8dd9f21585
function holstein_variational_solution(α, β; v=3.0, w=3.0) # v, w defaults

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [Inf, Inf]
    Δv = v - w # defines a constraint, so that v>w
    initial = [Δv + 0.01, w]

    # Feynman 1955 athermal action 
    f(x) = F_Hol(x[1] + x[2], x[2], β, α)

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

    # Return variational parameters that minimise the free energy.
    return Δv + w, w
end

# ╔═╡ 466bc0e9-4cd2-42f5-8bbc-1317b5c1c7b6
vh, wh = holstein_variational_solution(2.39; v=v, w=w)

# ╔═╡ a1d89e74-7ab6-4c8f-8748-63ede869decb
holstein_energy(2.39, vh, wh)

# ╔═╡ 60254d3a-97b0-4c7e-a6c3-ad9ea0d93c97
F_Hol(vh, wh, 100, 2.39)

# ╔═╡ 5f5785c4-b9be-4246-820d-3e3ccf841cc6
holstein_params = holstein_variational_solution.(α)
# holstein_params = var_params_k.(α, power=2, R=(6 * π^2)^(1 / 3))

# ╔═╡ 2ce6cc8d-0582-4f08-9e2f-5317f2f5b2e0
begin
    vh_new = [i[1] for i in holstein_params]
    wh_new = [i[2] for i in holstein_params]
end

# ╔═╡ f5e2a904-5e2f-4dcd-8691-b7620d24a7f0
begin
    pv = plot(α, vf, label="Fro v", legend=:topleft, linewidth=2, xlabel="alpha, α", ylabel="v and w", minorgrid=true)
    plot!(α, wf, label="Fro w", linewidth=2)
    plot!(α, vh_new, label="Hol v", linewidth=2, linestyle=:dash)
    plot!(α, wh_new, label="Hol w", linewidth=2, linestyle=:dash)
end

# ╔═╡ a8852b06-3f60-4bac-a610-8b0a1b8cf560
savefig(pv, "vw_comparison.pdf")

# ╔═╡ 650e36ab-7ac8-4d0b-ad74-f5693c2c1586
holstein_energies = [holstein_energy(α[i], vh_new[i], wh_new[i]) for i in 1:length(α)]
# holstein_energies = [Frohlich_energy(α[i], vh_new[i], wh_new[i], power=2, R=(6 * π^2)^(1 / 3)) for i in 1:length(α)]

# ╔═╡ b2454e06-146c-4fc1-8928-70c7982922ec
begin
    pe = plot(α, frohlich_energies .* 1000 * 1.05e-34 / 1.6e-19 * 1e12 * 2.25 * 2π, label="Fro", linewidth=2, xlabel="alpha, α", ylabel="Energy (meV)")
    plot!(α, holstein_energies .* 1000 * 1.05e-34 / 1.6e-19 * 1e12 * 2.25 * 2π, label="Hol", linewidth=2, linestyle=:dash)
end

# ╔═╡ 67154f5a-e61e-4a5e-abf8-215bf8c8646e
savefig(pe, "energy_comparison.pdf")

# ╔═╡ af11d3d0-b93d-4916-baad-5c81cec45ab4
"""
Temperature dependence free-energy and v/w comparisons
"""

# ╔═╡ 28f3a13d-cee1-4415-be15-4cefc474207d
T_range = 1:400

# ╔═╡ 7cd14da3-799d-4a53-a446-10a6842e5214
hol_params = [holstein_variational_solution(2.39, 2π * 2.25 * 1e12 * 1.05e-34 / T / 1.38e-23) for T in T_range]

# ╔═╡ 2baddbea-169a-4558-9f47-1c43e1269f1f
fro_params = [feynmanvw(2.39, 2π * 2.25 * 1e12 * 1.05e-34 / T / 1.38e-23) for T in T_range]

# ╔═╡ 1d45213c-3921-4a67-a76c-f37db47a14ad
begin
    vh_β = [i[1] for i in hol_params]
    wh_β = [i[2] for i in hol_params]
end

# ╔═╡ 3d423a59-b295-4026-8666-e08c38f01747
begin
    vf_β = [i[1] for i in fro_params]
    wf_β = [i[2] for i in fro_params]
end

# ╔═╡ 0e771039-19b8-4cec-bb12-8ddcb4b6975e
hol_energy_β = [F_Hol(vh_β[T], wh_β[T], 2π * 2.25 * 1e12 * 1.05e-34 / T_range[T] / 1.38e-23, 2.39) for T in 1:length(T_range)]

# ╔═╡ 07e267b3-03e8-4b1f-af0a-6c5ea6c57684
fro_energy_β = [F_Fro(vf_β[T], wf_β[T], 2π * 2.25 * 1e12 * 1.05e-34 / T_range[T] / 1.38e-23, 2.39) for T in 1:length(T_range)]

# ╔═╡ 1d08e7de-026b-4b1c-ba29-d12ce9f2b4ef
begin
    pvt = plot(T_range, vf_β, label="Fro v", linewidth=2, legend=:topleft, xlabel="Temperature (K)", ylabel="v and w params")
    plot!(T_range, wf_β, label="Fro w", linewidth=2)
    plot!(T_range, vh_β, label="Hol v", linestyle=:dash, linewidth=2)
    plot!(T_range, wh_β, label="Hol w", linestyle=:dash, linewidth=2)
end

# ╔═╡ a6ff82d8-ac2f-4034-ad64-c3800d67ecfa
begin
    pft = plot(T_range, fro_energy_β .* 1000 * 1.05e-34 / 1.6e-19 * 1e12 * 2.25 * 2π, label="Fro", xlabel="Temperature (K)", ylabel="Free energy (ħω)", linewidth=2)
    plot!(T_range, hol_energy_β .* 1000 * 1.05e-34 / 1.6e-19 * 1e12 * 2.25 * 2π, label="Hol", linewidth=2, linestyle=:dash)
end

# ╔═╡ a0cd088c-0412-4b1d-ba32-441f59e20c3c
savefig(pvt, "vw_comparison_temp.pdf")

# ╔═╡ 35864bca-4813-45ab-b579-ba9b40492189
savefig(pft, "free_energy_comparison.pdf")

# ╔═╡ a835bb67-f4d0-4b42-8cda-de6b44936ca2
"""
Have a look at effective mass
"""

# ╔═╡ 5ac9e949-4b55-42f5-9dd4-1de517051112
function fro_eff_mass(v, w, α)
    integrand(τ) = τ^2 * exp(-τ) / (D(τ, v, w))^(3 / 2)
    integral = α * quadgk(τ -> integrand(τ), 0.0, Inf)[1] / (3 * sqrt(π))
    return 1 + integral
end

# ╔═╡ 6a848cac-dc22-4419-8833-7279a4940834
function hol_eff_mass(v, w, α)
    R = (6 * π^2)^(1 / 3)
    integrand(τ) = τ^2 * exp(-τ) * (3 * erf(sqrt(D(τ, v, w) / 2) * R) / (2 * D(τ, v, w)^(5 / 2)) - R * exp(-D(τ, v, w) * R^2 / 2) * (D(τ, v, w) * R^2 + 3) / (sqrt(2π) * D(τ, v, w)^2))
    integral = α * quadgk(τ -> integrand(τ), 0.0, Inf)[1] / (3 * sqrt(2) * pi)
    return 1 + integral
end

# ╔═╡ 2041244d-1a95-4e9f-846b-0f64d5fca357
α[1:30]

# ╔═╡ 4d9a6e4c-0b35-4a55-ac12-4c21e52603b7
fro_m = [fro_eff_mass(vf[i], wf[i], α[i]) for i in 1:length(α)]

# ╔═╡ 4a4d9cd9-b991-40c9-9e09-18ab10f0e7ba
hol_m = [hol_eff_mass(vh_new[i], wh_new[i], α[i]) for i in 1:length(α)]

# ╔═╡ 308b415e-8965-4d31-b374-078add023990
begin
    pm = plot(α, 1 ./ fro_m, label="Fro", xlabel="α", ylabel="m₀ / m*", linewidth=2)
    plot!(α, 1 ./ hol_m, label="Hol", linewidth=2, linestyle=:dash)
end

# ╔═╡ 4b001539-7208-4100-ae75-6c1232f540f8
savefig(pm, "mass_comparison_alpha.pdf")

# ╔═╡ 614fe5ba-f265-420d-b914-0ecadedfbacc
md"""
# FHIP mobility theory
"""

# ╔═╡ 750b1d30-79c4-4f5a-8df4-3b6fcffb7d49
function P(x)
    1 / (exp(x) - 1)
end

# ╔═╡ 31c30c6c-d68a-48f7-a9f7-c6a6572367cb
function D_c(u, β, v, w)
    w^2 / v^2 * ((v^2 - w^2) / (w^2 * v) * (1 - exp(im * v * u) + 4 * P(β * v) * (sin(v * u / 2))^2) - im * u + u^2 / β)
end

# ╔═╡ 4e17762c-461b-43b5-b037-da144af61e27
function fro_S(u, α, β, v, w)
    2 * α / (3 * √π) * (exp(im * u) + 2 * P(β) * cos(u)) / D_c(u, β, v, w)^(3 / 2)
end

# ╔═╡ babe6b05-bb06-42b2-84aa-27b2c51d7226
function hol_S(u, α, β, v, w)
    R = (6 * π^2)^(1 / 3)
    2 * α / (3 * √π) * (exp(im * u) + 2 * P(β) * cos(u)) * (erf(sqrt(D_c(u, β, v, w) / 2) * R) / (2 * D_c(u, β, v, w)^(3 / 2)) - R * exp(-D_c(u, β, v, w) * R^2 / 2) / (sqrt(2π) * D_c(u, β, v, w))) * 2π
end

# ╔═╡ 2042d4f1-b830-4d76-a38f-240506168618
function fro_chi(Ω, α, β, v, w)
    integrand(u) = (1 - exp(im * Ω * u)) * imag(fro_S(u, α, β, v, w))
    integral = quadgk(u -> integrand(u), 0.0, Inf, rtol=1e-4)[1]
    return integral
end

# ╔═╡ bab68f24-c8fb-4b8d-865f-5e9e760513ca
function hol_chi(Ω, α, β, v, w)
    integrand(u) = (1 - exp(im * Ω * u)) * imag(hol_S(u, α, β, v, w))
    integral = quadgk(u -> integrand(u), 0.0, Inf, rtol=1e-4)[1]
    return integral
end

# ╔═╡ 9655e230-37d7-4504-8440-34078996a22d
function fro_impedance(Ω, α, β, v, w)
    return (Ω^2 - fro_chi(Ω, α, β, v, w)) / (-im * Ω)
end

# ╔═╡ 8455cb9b-e0f7-4832-a880-84a31d7b5d11
function hol_impedance(Ω, α, β, v, w)
    return (Ω^2 - hol_chi(Ω, α, β, v, w)) / (-im * Ω)
end

# ╔═╡ 75830e58-af79-4f94-b6a9-92df688355ce
md"""# Study freq-dep impedance in rubrene"""

# ╔═╡ 9452286f-e801-4db0-b884-f0c7b2373879
Ω_range = 0.01:0.1:10.0

# ╔═╡ 36ce62de-bf89-4e9d-beb5-47da3cfebc3d
fro_z = fro_impedance.(Ω_range, 2.39, 100, v_k, w_k)

# ╔═╡ dae9f84a-f7b9-4ec1-8399-d39410df2ab9
hol_z = hol_impedance.(Ω_range, 2.39, 100, vh_k, wh_k)

# ╔═╡ fb473fb1-1e78-48c8-b20d-aa64f0d35ba2
begin
    rubrene_real_imped = plot(xlabel="Omega (THz)", ylabel="Re(Χ)")
    plot!(Ω_range, -real.(1 ./ fro_z), linewidth=2, label="Fro")
    plot!(Ω_range, -real.(1 ./ hol_z), linewidth=2, linestyle=:dash, label="Hol")
end

# ╔═╡ 907cc295-5315-467f-8dc1-bec11b3b284a
savefig(rubrene_real_imped, "rubrene_real_imped.pdf")

# ╔═╡ 68ab6be9-2f77-4bc3-900a-9087b45a2a33
md"# Imag Imped"

# ╔═╡ 86dd4b1d-51e9-4f41-b88e-fd51d639baa3
begin
    rubrene_imag_imped = plot(xlabel="Omega (THz)", ylabel="Im(Χ)")
	
    plot!(Ω_range, -imag.(1 ./ fro_z), ylims = (-1, 2), linewidth=2, label="Fro")
    plot!(Ω_range, -imag.(1 ./ hol_z), linewidth=2, linestyle=:dash, label="Hol")
end

# ╔═╡ 207db8eb-8de1-4e29-8a95-679fcbc902f2
savefig(rubrene_imag_imped, "rubrene_imag_imped.pdf")

# ╔═╡ 5b150eed-0d08-4edc-a7a4-8e421bca66d0
md""" # DC Chi functions, for DC mobility."""

# ╔═╡ 63eaadaa-9f2b-4de5-964e-7e1a8dfd0121
function fro_chi_dc(α, β, v, w)
    integrand(u) = -im * u * imag(fro_S(u, α, β, v, w))
    integral = quadgk(u -> integrand(u), 0.0, Inf, rtol=1e-4)[1]
    return integral
end

# ╔═╡ 2e6d1157-62db-4774-887d-4872a1b64d4d
function hol_chi_dc(α, β, v, w)
    integrand(u) = -im * u * imag(hol_S(u, α, β, v, w))
    integral = quadgk(u -> integrand(u), 0.0, Inf, rtol=1e-4)[1]
    return integral
end

# ╔═╡ 8e489020-256a-4014-9df3-99ef4bc35f6b
function fro_mobility(α, β, v, w)
    abs(1 / imag(fro_chi_dc(α, β, v, w)))
end

# ╔═╡ 71653a17-89dd-4cfd-8c5c-abd6c655c362
function hol_mobility(α, β, v, w)
    abs(1 / imag(hol_chi_dc(α, β, v, w)))
end

# ╔═╡ cec30175-2096-42a6-965a-5cf5aec4fcaf
fro_μ = [fro_mobility(2.39, 2π * 2.25 * 1e12 * 1.05e-34 / T / 1.38e-23, vf_β[T], wf_β[T]) for T in 15:length(T_range)]

# ╔═╡ 5c758d90-a0b3-4e4e-a492-d96536ccbe6f
hol_μ = [hol_mobility(2.39, 2π * 2.25 * 1e12 * 1.05e-34 / T / 1.38e-23, vh_β[T], wh_β[T]) for T in 15:length(T_range)]

# ╔═╡ 0459ac5b-48bf-4645-8d91-a4c61c07578e
begin
    MAPI_mob_plot = plot(T_range[15:end], fro_μ .* 1.6e-19 / (2π * 2.25 * 0.12 * 1e12 * 9.11e-31) * 100^2, label="Fro", minorgrid=true, ylims=(50, 200), linewidth=2, xlabel="Temperature (K)", ylabel="Mobility (cm^2/Vs)", title="MAPI Mobility")
    plot!(T_range[15:end], hol_μ .* 1.6e-19 / (2π * 2.25 * 0.12 * 1e12 * 9.11e-31) * 100^2, label="Hol", linewidth=2, linestyle=:dash)
end

# ╔═╡ 4bd4af25-59ef-4008-b039-52ada489f81c
savefig(MAPI_mob_plot, "MAPI_mob_plot.pdf")

# ╔═╡ 639bd491-ccaa-4913-b4e9-c8f20977bce5
"""
300K Rubrene mobility. Get ~ 13cm^2/Vs for Holstein and ~ 20.6 cm^2/Vs Frohlich
"""

# ╔═╡ 799c1983-d4c4-4882-9024-29a492c7b2ba
rubrene_params_hol = [holstein_variational_solution(2.21, 1e13 * 1.41 * 1.05e-34 / T / 1.38e-23) for T in T_range]

# ╔═╡ 2858d507-0ce0-4a98-ad35-fa12cdb5fa89
rubrene_params_fro = [feynmanvw(2.21, 1e13 * 1.41 * 1.05e-34 / T / 1.38e-23) for T in T_range]

# ╔═╡ d1b784a5-ba76-46f4-8613-53503997ec05
begin
    rubrene_vh = [i[1] for i in rubrene_params_hol]
    rubrene_wh = [i[2] for i in rubrene_params_hol]
end

# ╔═╡ d41fc336-e9d1-4f00-8a68-fe2f26103a82
begin
    rubrene_vf = [i[1] for i in rubrene_params_fro]
    rubrene_wf = [i[2] for i in rubrene_params_fro]
end

# ╔═╡ 0014b87c-e9e1-45a2-be18-ca587adcc1af
rubrene_mob_hol = [hol_mobility(2.21, 1e13 * 1.41 * 1.05e-34 / T_range[T] / 1.38e-23, rubrene_vh[T], rubrene_wh[T]) for T in 15:length(T_range)] .* 1.6e-19 / (1e13 * 1.41 * 9.11e-31 * 0.88) * 100^2

# ╔═╡ 560edf4e-fa8f-4c29-be6b-f51d17783f36
rubrene_mob_fro = [fro_mobility(2.21, 1e13 * 1.41 * 1.05e-34 / T_range[T] / 1.38e-23, rubrene_vf[T], rubrene_wf[T]) for T in 15:length(T_range)] .* 1.6e-19 / (1e13 * 1.41 * 9.11e-31 * 0.88) * 100^2

# ╔═╡ 0cc49e62-1191-43f6-976a-4b649b8348f6
begin
    rubrene_mob_plot = plot(T_range[15:end], rubrene_mob_fro, label="Fro", minorgrid=true, ylims=(0, 40), linewidth=2, xlabel="Temperature (K)", ylabel="Mobility (cm^2/Vs)", title="Rubrene Mobility")
    plot!(T_range[15:end], rubrene_mob_hol, label="Hol", linewidth=2, linestyle=:dash)
end

# ╔═╡ fe290c3a-7a5d-4d74-b5e3-495714598f8a
savefig(rubrene_mob_plot, "rubrene_mob_plot.pdf")

# ╔═╡ 3a8efacf-2936-4517-9fa3-538809d7c67f
vrh, wrh = holstein_variational_solution(2.21, 1e13 * 1.41 * 1.05e-34 / 300 / 1.38e-23)

# ╔═╡ 5e3e532c-10bc-4ba7-8a97-84c5b84e0b1b
vrf, wrf = feynmanvw(2.21, 1e13 * 1.41 * 1.05e-34 / 300 / 1.38e-23)

# ╔═╡ 30c5b617-5d93-40fc-9b65-91f5dd1d2a1b
hol_mobility(2.21, 1e13 * 1.41 * 1.05e-34 / 300 / 1.38e-23, vrh, wrh) * 1.6e-19 / (1e13 * 1.41 * 9.11e-31 * 0.88) * 100^2

# ╔═╡ f170dfc6-cdd2-4a13-9db7-3035f159d9b2
fro_mobility(2.21, 1e13 * 1.41 * 1.05e-34 / 300 / 1.38e-23, vrf, wrf) * 1.6e-19 / (1e13 * 1.41 * 9.11e-31 * 0.88) * 100^2

# ╔═╡ c29f44ce-759e-4a52-bbf5-2cdc4136114e


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[compat]
Optim = "~1.7.0"
Plots = "~1.30.1"
QuadGK = "~2.4.2"
SpecialFunctions = "~2.1.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0310e08cb19f5da31d08341c6120c047598f5b9c"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.5.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "844b061c104c408b24537482469400af6075aae4"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.5"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "d3ba08ab64bdfd27234d3f61956c966266757fe6"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.7"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "04ed1f0029b6b3af88343e439b995141cb0d0b8d"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.17.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "a69dd6db8a809f78846ff259298678f0d6212180"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.34"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "57f7cde02d7a53c9d1d28443b9f11ac5fbe7ebc9"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.3"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c98aea696662d09e215ef7cda5296024a9646c75"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.4"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "bc9f7725571ddb4ab2c4bc74fa397c1c5ad08943"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.69.1+0"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "e07a1b98ed72e3cdd02c6ceaab94b8a606faca40"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.2.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "fe9aea4ed3ec6afdfbeb5a4f39a2208909b162a6"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.5"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "45b288af6956e67e621c5cbb2d75a261ab58300b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.20"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "1903afc76b7d01719d9c30d3c7d501b61db96721"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.4"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "151d91d63d8d6c1a5789ecb7de51547e00480f1b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.4"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "c95373e73290cf50a8a22c3375e4625ded5c5280"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "d0a61518267b44a70427c0b690b5e993a4f5fe01"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.30.2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6954a456979f23d05085727adb17c4551c19ecd1"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.12"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "b03a3b745aa49b566f128977a7dd1be8711c5e71"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.14"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═66d2d0e0-ee56-11ec-21f9-5d6ae25cff48
# ╟─9f6474e5-41f1-4b96-aeee-d299404a3718
# ╠═9429ae6b-55e1-4fb2-9a88-ae7636c35bb5
# ╠═1cf1e5cc-301c-4779-86cb-5b6436ecebf0
# ╠═2a609c87-4579-490e-99e7-0d0c25e11333
# ╟─7bf82f19-e6fc-4c3f-b7e8-a02bebf0d0dd
# ╠═6cd58deb-fb04-44fc-a854-dce5b5cb8565
# ╠═00fc8b3d-c485-43e2-92f0-451539ecf1e6
# ╠═dcf7f24a-6ef4-4613-9778-e2637f5e6c16
# ╠═6714b04e-9b26-42d0-89eb-f6c510cd5f56
# ╠═429d0098-9d78-47e5-915b-cf16c7fce660
# ╟─f7a07dc5-ba45-4870-b94c-44ae5d85436d
# ╠═82f23e07-599f-4372-8e82-e6086ab1b40a
# ╠═71118d82-37fb-4504-b944-c3461de91233
# ╠═2be860db-52b0-4514-872f-db79dc454323
# ╟─cc237871-08fd-4658-8fa5-8658d9394b3b
# ╠═80cd6d37-92d3-4757-a62e-490cffb4b63c
# ╠═410d234c-25fe-4752-926e-f6c5d263ef05
# ╟─907b2475-0749-4d85-b09e-6eed75d35e68
# ╠═0b9cb699-eace-4f78-a547-0a897aa6c361
# ╠═63bf4a80-ed9c-46fc-9ede-6d0b05fc2519
# ╠═dde4d55d-702a-4ced-9e37-13f89f78d3d5
# ╠═14f66830-c111-4466-bce9-965b2e48d2d4
# ╠═9c89b255-acb4-495c-9cb4-b71b05c13708
# ╠═eb401f35-2fb2-4432-aa9a-28d393be71fb
# ╠═079d7328-c2e1-4681-a4e3-e19a53841d25
# ╠═bc20477e-2ce3-4377-9af0-5d479256175a
# ╠═380c5c1b-01e1-455b-a757-46c87fc45285
# ╠═466bc0e9-4cd2-42f5-8bbc-1317b5c1c7b6
# ╠═a1d89e74-7ab6-4c8f-8748-63ede869decb
# ╠═e7a4f46e-8fcc-48ef-9c5c-5790b42d7648
# ╠═78936a9d-9567-4def-bf89-85be9e62e12d
# ╠═8d726e23-5a40-4a10-8405-7c9775a0c9d0
# ╠═2cdfdcbc-a6e7-444c-a9a0-47529ddb091a
# ╠═5f5785c4-b9be-4246-820d-3e3ccf841cc6
# ╠═2ce6cc8d-0582-4f08-9e2f-5317f2f5b2e0
# ╠═f5e2a904-5e2f-4dcd-8691-b7620d24a7f0
# ╠═35ce9514-521c-4985-a85e-f236f67417c9
# ╠═650e36ab-7ac8-4d0b-ad74-f5693c2c1586
# ╠═b2454e06-146c-4fc1-8928-70c7982922ec
# ╠═a8852b06-3f60-4bac-a610-8b0a1b8cf560
# ╠═67154f5a-e61e-4a5e-abf8-215bf8c8646e
# ╟─3f5069df-89e2-4ac7-96e4-f3c79936fdd9
# ╠═9c1284f7-f503-4c4d-be6a-e4e04043246a
# ╠═b36f7997-f29c-4a92-bb00-9e82ae7d2a05
# ╠═d4e1ac9e-84d7-4809-ac1c-c24a9263db84
# ╠═ca1d9705-b49e-43f8-bd17-ab4a7e28b578
# ╠═4e99725c-0895-4bc3-97ba-70eb675c7f47
# ╠═e8b8e4de-e548-447e-82c5-979124ce6cf3
# ╠═f87703c9-fa55-420a-b3cc-c6b947ef456a
# ╠═c1358b12-4f68-4485-8ee7-73ef672f3ed8
# ╠═fc699405-e1d6-48dd-aff3-a7ae7fb0114b
# ╠═9ed1bad6-9ac5-444f-8d02-974bdd955b09
# ╠═5d5bf742-89c2-45e6-81f6-3b22e90b7bd0
# ╠═60254d3a-97b0-4c7e-a6c3-ad9ea0d93c97
# ╠═2455abad-abf2-4b0a-bb3f-30d97de9a9ec
# ╠═4d5e8c8e-2225-43fd-a934-0c8dd9f21585
# ╠═af11d3d0-b93d-4916-baad-5c81cec45ab4
# ╠═28f3a13d-cee1-4415-be15-4cefc474207d
# ╠═7cd14da3-799d-4a53-a446-10a6842e5214
# ╠═2baddbea-169a-4558-9f47-1c43e1269f1f
# ╠═1d45213c-3921-4a67-a76c-f37db47a14ad
# ╠═3d423a59-b295-4026-8666-e08c38f01747
# ╠═0e771039-19b8-4cec-bb12-8ddcb4b6975e
# ╠═07e267b3-03e8-4b1f-af0a-6c5ea6c57684
# ╠═1d08e7de-026b-4b1c-ba29-d12ce9f2b4ef
# ╠═a6ff82d8-ac2f-4034-ad64-c3800d67ecfa
# ╠═a0cd088c-0412-4b1d-ba32-441f59e20c3c
# ╠═35864bca-4813-45ab-b579-ba9b40492189
# ╠═a835bb67-f4d0-4b42-8cda-de6b44936ca2
# ╠═5ac9e949-4b55-42f5-9dd4-1de517051112
# ╠═6a848cac-dc22-4419-8833-7279a4940834
# ╠═2041244d-1a95-4e9f-846b-0f64d5fca357
# ╠═4d9a6e4c-0b35-4a55-ac12-4c21e52603b7
# ╠═4a4d9cd9-b991-40c9-9e09-18ab10f0e7ba
# ╠═308b415e-8965-4d31-b374-078add023990
# ╠═4b001539-7208-4100-ae75-6c1232f540f8
# ╟─614fe5ba-f265-420d-b914-0ecadedfbacc
# ╠═750b1d30-79c4-4f5a-8df4-3b6fcffb7d49
# ╠═31c30c6c-d68a-48f7-a9f7-c6a6572367cb
# ╠═4e17762c-461b-43b5-b037-da144af61e27
# ╠═babe6b05-bb06-42b2-84aa-27b2c51d7226
# ╠═2042d4f1-b830-4d76-a38f-240506168618
# ╠═bab68f24-c8fb-4b8d-865f-5e9e760513ca
# ╠═9655e230-37d7-4504-8440-34078996a22d
# ╠═8455cb9b-e0f7-4832-a880-84a31d7b5d11
# ╟─75830e58-af79-4f94-b6a9-92df688355ce
# ╠═9452286f-e801-4db0-b884-f0c7b2373879
# ╠═36ce62de-bf89-4e9d-beb5-47da3cfebc3d
# ╠═dae9f84a-f7b9-4ec1-8399-d39410df2ab9
# ╠═fb473fb1-1e78-48c8-b20d-aa64f0d35ba2
# ╠═907cc295-5315-467f-8dc1-bec11b3b284a
# ╠═68ab6be9-2f77-4bc3-900a-9087b45a2a33
# ╠═86dd4b1d-51e9-4f41-b88e-fd51d639baa3
# ╠═207db8eb-8de1-4e29-8a95-679fcbc902f2
# ╟─5b150eed-0d08-4edc-a7a4-8e421bca66d0
# ╠═63eaadaa-9f2b-4de5-964e-7e1a8dfd0121
# ╠═2e6d1157-62db-4774-887d-4872a1b64d4d
# ╠═8e489020-256a-4014-9df3-99ef4bc35f6b
# ╠═71653a17-89dd-4cfd-8c5c-abd6c655c362
# ╠═cec30175-2096-42a6-965a-5cf5aec4fcaf
# ╠═5c758d90-a0b3-4e4e-a492-d96536ccbe6f
# ╠═0459ac5b-48bf-4645-8d91-a4c61c07578e
# ╠═4bd4af25-59ef-4008-b039-52ada489f81c
# ╠═639bd491-ccaa-4913-b4e9-c8f20977bce5
# ╠═799c1983-d4c4-4882-9024-29a492c7b2ba
# ╠═2858d507-0ce0-4a98-ad35-fa12cdb5fa89
# ╠═d1b784a5-ba76-46f4-8613-53503997ec05
# ╠═d41fc336-e9d1-4f00-8a68-fe2f26103a82
# ╠═0014b87c-e9e1-45a2-be18-ca587adcc1af
# ╠═560edf4e-fa8f-4c29-be6b-f51d17783f36
# ╠═0cc49e62-1191-43f6-976a-4b649b8348f6
# ╠═fe290c3a-7a5d-4d74-b5e3-495714598f8a
# ╠═3a8efacf-2936-4517-9fa3-538809d7c67f
# ╠═5e3e532c-10bc-4ba7-8a97-84c5b84e0b1b
# ╠═30c5b617-5d93-40fc-9b65-91f5dd1d2a1b
# ╠═f170dfc6-cdd2-4a13-9db7-3035f159d9b2
# ╠═c29f44ce-759e-4a52-bbf5-2cdc4136114e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
