### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ c14d5f69-8ed3-4ad7-b094-8de94394ad8e
using Revise

# ╔═╡ 3bbf5ec7-c772-4a1a-9eab-c0f26a99e94f
using PolaronMobility

# ╔═╡ f088c85f-4d12-4f4c-8e3e-318652b68a43
using Plots

# ╔═╡ 2b5f21a3-1130-47cd-9938-f49a60048e27
using DataFrames

# ╔═╡ b38c2c58-689b-48f9-b5b1-1499651e1928
using CSV

# ╔═╡ e9685260-9ece-4c3b-b110-096a7c39857b
using QuadGK

# ╔═╡ e495c578-ac0e-46c3-8f2f-aa0541cf3574
using Optim

# ╔═╡ 7f124d90-0430-11ec-129e-bf0ec7dbda06
# ESA Results

# ╔═╡ bf2bb3d0-6f5d-412a-8f9e-b10cb30560bd
begin
	ħ = 1.05457162825e-34
	eV = 1.602176487e-19
	me = 9.10938188e-31
	kB =  1.3806504e-23
	ϵ_0 = 8.854E-12
	amu = 1.660_539_066_60e-27
	Ha = 4.35974820e-18
	Bohr = 5.29177249e-11
	c = 2.99792458e8
end

# ╔═╡ 802b2d3b-d1d8-40c7-9f34-407780385457
α_range = 1:12

# ╔═╡ 30fe032c-f911-4f7e-aa24-1fb473ecf9e0
β_range = vcat([i for i in 0.1:0.1:8.0], [100.0])

# ╔═╡ 71355321-bb21-42d5-9966-60a9438b762c
Ω_range = vcat([1e-200], [i for i in 0.01:0.1:20.01])

# ╔═╡ a1161327-a8e8-44a0-8af7-b14d437ae14c
ω = 1

# ╔═╡ f40b9810-6a0e-4627-a290-80a446d101c8
m = 1

# ╔═╡ 95538cd0-3588-423a-ab1a-95915ee60828
begin
	# Equation 31: The <|X(t) - X(s)|^{-1}> * exp(-|t-w|) effective action.
	A_integrand(x, v, w) = (w^2 * x + (v^2 - w^2) / v * (abs(1 - exp(-v * x))))^(-0.5) * exp(-x)

	A(v, w, α) = π^(-0.5) * α * v * QuadGK.quadgk(x -> A_integrand(x, v, w), BigFloat(0.0), Inf)[1]

	# Equation 33: Lowest Free energy E = -B - A where B = -3/(4v)*(v-w)^2.
	free_energy(v, w, α) = (3 / (4 * v)) * (v - w)^2 - A(v, w, α)
end

# ╔═╡ 5eaec90c-f922-496e-a388-f40a49f38b26
begin
	# Equation 62d in Hellwarth.
	Y(x, v, β) = 1 / (1 - exp(-v * β)) * (1 + exp(-v * β) - exp(-v * x) - exp(v * (x -     β)))

	# Integrand of Equation 62c in Hellwarth.
	A_integrand(x, v, w, β) =  (exp(β - x) + exp(x)) / sqrt(1e-10 + abs(w^2 * x * (1 - x / β)     + Y(x, v, β) * (v^2 - w^2) / v))

	# Equation 62c in Hellwarth.
	A(v, w, α, β) = α * v / (sqrt(π) * (exp(BigFloat(β)) - 1)) * QuadGK.quadgk(x ->       A_integrand(x, v, w, β), BigFloat(0), BigFloat(β / 2))[1]

	# Equation 62b in Hellwarth. Equation 20 in Osaka.
	B(v, w, β) = 3 / β * (log(v / w) - 1 / 2 * log(2 * π * BigFloat(β)) - log(sinh(v *     BigFloat(β) / 2) / sinh(w * BigFloat(β) / 2)))

	# Equation 62e in Hellwarth. Equation 17 in Osaka.
	C(v, w, β) = 3 / 4 * (v^2 - w^2) / v * (coth(v * BigFloat(β) / 2) - 2 / (v * β))

	# Equation 62a in Hellwarth. In paragraph below Equation 22 in Osaka; has extra 1/	  β due to different definition of A, B & C.
	function free_energy(v, w, α, β)
		setprecision(BigFloat, 32)
		a = A(v, w, α, β)
		b = B(v, w, β)
		c = C(v, w, β)
		-(a + b + c)
	end
end

# ╔═╡ 5568b92e-5a20-4771-b440-4de309b6173c
function variation(α, β; v = 0.0, w = 0.0)

    initial = [v, w]
    
    # Limits of the optimisation.
    lower = [1, 1]
    upper = [200, 200]

    # Osaka Free Energy function to minimise.
    f(x) = free_energy(x[1], x[2], α, β)

    # Use Optim to optimise the free energy function w.r.t v and w.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff = :forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
    )

    # Get v and w values that minimise the free energy.
    v, w = Optim.minimizer(solution)
	
	println("v: $v, w: $w")

    # Return variational parameters that minimise the free energy.
    return v, w
end

# ╔═╡ 6fc1e642-77ef-49bb-895d-a4006b563f5e
function feynamn_mass(α, v, w, m, ω)
	j(τ) = 1 + (v * ω / w^2) * (1 - w^2 / v^2) * (1 - exp(-v * τ / ω)) / τ
	integrand(τ) = exp(-τ) * τ^(1/2) * j(τ)^(-3/2)
	return m * (1 + α * (v / w)^3 / (3 * π^(1/2)) * quadgk(τ -> integrand(τ), 0.0, Inf)[1])
end

# ╔═╡ ba42efa0-8f26-4ff1-a0ac-a56c626c7525
begin
	v = Array{Float64}(undef, length(β_range), length(α_range))
	w = Array{Float64}(undef, length(β_range), length(α_range))
	FE = Array{Float64}(undef, length(β_range), length(α_range))
	M = Array{Float64}(undef, length(β_range), length(α_range))
	for α in 1:length(α_range), β in 1:length(β_range)
		println("α: $(α_range[α]), β: $(β_range[β])")
		if β == 1
			v[β, α], w[β, α] = variation(α_range[α], β_range[β]; v = 4.0, w = 2.0)
		else
			v[β, α], w[β, α] = variation(α_range[α], β_range[β]; v = v[β-1, α], w = w[β-1, α])
		end
		FE[β, α] = free_energy(v[β, α], w[β, α], α_range[α], β_range[β])
		M[β, α] = feynamn_mass(α_range[α], v[β, α], w[β, α], m, ω)
		CSV.write("v_data.csv", DataFrames.DataFrame([[0.0, [i for i in α_range]...]'; [β_range v]], :auto))
		CSV.write("w_data.csv", DataFrames.DataFrame([[0.0, [i for i in α_range]...]'; [β_range w]], :auto))
		CSV.write("F_data.csv", DataFrames.DataFrame([[0.0, [i for i in α_range]...]'; [β_range FE]], :auto))
		CSV.write("M_data.csv", DataFrames.DataFrame([[0.0, [i for i in α_range]...]'; [β_range M]], :auto))
	end
end

# ╔═╡ c6a84269-a382-4fbd-a102-8a949fe531ec
v

# ╔═╡ 0e9a3b4a-4626-4efb-8925-f493f0ab5f34
w

# ╔═╡ 441cd179-c28f-4e23-ae23-e6a65d7edb46
FE

# ╔═╡ 3b93da85-80a2-466e-b28e-dcf22b4ff0fe
M

# ╔═╡ 2405eb4e-b81c-46fb-94a3-50456b4d7512
begin
	# Integrand of (31) in Feynman I (Feynman 1955, Physical Review, "Slow 	electrons...")
	fF(τ, v, w) = (abs(w^2 * τ + (v^2 - w^2) / v * (1 - exp(- v * τ))))^-0.5 * exp(-τ)
	
	# (31) in Feynman I
	AF(v,w,α) = π^(-0.5) * α * v * quadgk(τ -> fF(τ, v, w), 0, Inf)[1]
	
	# (33) in Feynman I
	F(v,w,α) = (3 / (4 * v)) * (v - w)^2 - AF(v, w, α)
end

# ╔═╡ ae352c06-69fe-41fa-ad4c-bafb23235012
function complex_conductivity(Ω, β, α, v, w)
	
	println("α = $α, Ω = $Ω")

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v)

    # FHIP1962, page 1009, eqn (35c).
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) * coth(β * v / 2) + x^2 / β - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36).
    S(x) = 2 * α / (3 * √π) * (exp(1im * x) + 2 * cos(x) / (exp(β) - 1)) / (D(x))^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(x) = (1 - exp(1im * Ω * x)) * imag(S(x)) / Ω
    1 / (-1im * Ω + 1im * QuadGK.quadgk(x -> integrand(x), 0.0, Inf)[1])
end

# ╔═╡ 6d015de5-3c44-400c-ba60-80d645a2789c
begin
	σ = Array{ComplexF64}(undef, length(Ω_range), length(β_range))
	for α in 1:length(α_range)
		for β in 1:length(β_range), Ω in 1:length(Ω_range)
			println("α: $(α_range[α]), β: $(β_range[β]), Ω: $(Ω_range[Ω])")
			σ[Ω, β] = complex_conductivity(Ω_range[Ω], β_range[β], α_range[α], v[β, α], w[β, α])
		end
		CSV.write("conductivity_data_$α.csv", DataFrames.DataFrame([[0.0, [i for i in β_range]...]'; [Ω_range σ]], :auto))
	end
end

# ╔═╡ 8f052fb6-65bd-4c5f-b931-af616548fb79
σ

# ╔═╡ 38f4dfbc-9a01-494c-8674-0d1b1cc324d3
contour_real = Plots.contourf(Ω_range, β_range, log.(abs.(real.(σ)))', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12)

# ╔═╡ 22c03bea-a0d5-4291-a568-f29c21d46443
begin
	plot(β_range[1:end-1], v[1:end-1])
	plot!(β_range[1:end-1], w[1:end-1])
end

# ╔═╡ 35d0a3c8-025c-4e63-a65d-0010baf5b61c
plot(β_range[1:end-1], M[1:end-1])

# ╔═╡ Cell order:
# ╠═7f124d90-0430-11ec-129e-bf0ec7dbda06
# ╠═c14d5f69-8ed3-4ad7-b094-8de94394ad8e
# ╠═3bbf5ec7-c772-4a1a-9eab-c0f26a99e94f
# ╠═f088c85f-4d12-4f4c-8e3e-318652b68a43
# ╠═2b5f21a3-1130-47cd-9938-f49a60048e27
# ╠═b38c2c58-689b-48f9-b5b1-1499651e1928
# ╠═e9685260-9ece-4c3b-b110-096a7c39857b
# ╠═e495c578-ac0e-46c3-8f2f-aa0541cf3574
# ╠═bf2bb3d0-6f5d-412a-8f9e-b10cb30560bd
# ╠═802b2d3b-d1d8-40c7-9f34-407780385457
# ╠═30fe032c-f911-4f7e-aa24-1fb473ecf9e0
# ╠═71355321-bb21-42d5-9966-60a9438b762c
# ╠═a1161327-a8e8-44a0-8af7-b14d437ae14c
# ╠═f40b9810-6a0e-4627-a290-80a446d101c8
# ╠═95538cd0-3588-423a-ab1a-95915ee60828
# ╠═5568b92e-5a20-4771-b440-4de309b6173c
# ╠═5eaec90c-f922-496e-a388-f40a49f38b26
# ╠═6fc1e642-77ef-49bb-895d-a4006b563f5e
# ╠═ba42efa0-8f26-4ff1-a0ac-a56c626c7525
# ╠═c6a84269-a382-4fbd-a102-8a949fe531ec
# ╠═0e9a3b4a-4626-4efb-8925-f493f0ab5f34
# ╠═441cd179-c28f-4e23-ae23-e6a65d7edb46
# ╠═3b93da85-80a2-466e-b28e-dcf22b4ff0fe
# ╠═2405eb4e-b81c-46fb-94a3-50456b4d7512
# ╠═ae352c06-69fe-41fa-ad4c-bafb23235012
# ╠═6d015de5-3c44-400c-ba60-80d645a2789c
# ╠═8f052fb6-65bd-4c5f-b931-af616548fb79
# ╠═38f4dfbc-9a01-494c-8674-0d1b1cc324d3
# ╠═22c03bea-a0d5-4291-a568-f29c21d46443
# ╠═35d0a3c8-025c-4e63-a65d-0010baf5b61c
