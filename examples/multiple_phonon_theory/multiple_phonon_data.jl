### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 1eb559b3-dbce-4caf-800f-96a2495dee1b
using Revise

# ╔═╡ 77d14c3f-bcaf-432c-9005-f5ec184d84b7
using QuadGK

# ╔═╡ ba2dadc8-d9a2-452e-9f94-406ccfc7404b
using Optim

# ╔═╡ c05cd489-68c8-418e-bbc2-a33bc42781f0
using Plots

# ╔═╡ 1d300e7d-5efb-4ad2-a826-c3f8a86083d8
using PolaronMobility

# ╔═╡ 6a99cd13-68cb-474c-96af-7d78ee4e7095
using DataFrames

# ╔═╡ 49663b33-1604-430b-a815-9eb985802987
using CSV

# ╔═╡ 314af942-8922-40d4-a92b-f75b715c276c
using LaTeXStrings

# ╔═╡ 6ff68ba5-df0c-42a0-8cf1-dedeced77e8d
using ColorSchemes

# ╔═╡ 29d62408-0629-4f39-909e-a6ba8696e955
begin
	# Physical constants
	const ħ = 1.05457162825e-34; # Reduced Planck's constant (kg m^2 s^{-1})
	const eV = 1.602176487e-19; # Electron Volt (kg m^2 s^{-2})
	const m_e = 9.10938188e-31; # Electron Mass (kg)
	const k_B = 1.3806504e-23; # Boltzmann's constant (kg m^2 K^{-1} s^2)
	const ϵ_0 = 8.854e-12 # Dielectric constant (C^2 N^{-1} m^{-2})
	const c = 2.99792458e8 # Speed of light (m s^{-1})
	const amu = 1.66053906660e-27 # Atomic Mass Unit (kg)
end

# ╔═╡ 35390bdf-e260-4566-8a97-7a0b5a3ef07a
MAPI= [
# 96.20813558773261 0.4996300522819191
# 93.13630357703363 1.7139631746083817
# 92.87834578121567 0.60108592692181
# 92.4847918585963 0.0058228799414729
# 92.26701437594754 0.100590086574602
# 89.43972834606603 0.006278895133832249
# 46.89209141511332 0.2460894564364346
# 46.420949316788 0.14174282581124137
# 44.0380222871706 0.1987196948553428
# 42.89702947649343 0.011159939465770681
# 42.67180170168193 0.02557751102757614
# 41.46971205834201 0.012555230726601503
# 37.08982543385215 0.00107488277468418
# 36.53555265689563 0.02126940080871224
# 30.20608114002676 0.009019481779712388
# 27.374810898415028 0.03994453721421388
# 26.363055017011728 0.05011922682554448
# 9.522966890022039 0.00075631870522737
4.016471586720514 0.08168931020200264
3.887605410774121 0.006311654262282101
3.5313112232401513 0.05353548710183397
2.755392921480459 0.021303020776321225
2.4380741812443247 0.23162784335484837
2.2490917637719408 0.2622203718355982
2.079632190634424 0.23382298607799906
2.0336707697261187 0.0623239656843172
1.5673011873879714 0.0367465760261409
1.0188379384951798 0.0126328938653956
1.0022960504442775 0.006817361620021601
0.9970130778462072 0.0103757951973341
0.9201781906386209 0.01095811116040592
0.800604081794174 0.0016830270365341532
0.5738689505255512 0.00646428491253749
# 0.022939578929507105 8.355742795827834e-06   # Acoustic modes!
# 0.04882611767873102 8.309858592685e-06
# 0.07575149723846182 2.778248540373041e-05
]

# ╔═╡ e0ce47a5-368a-4da7-83fc-9f921af0f747
function ϵ_ionic_mode(phonon_mode_freq, ir_activity, volume) # single ionic mode
    ω_j = 2π * phonon_mode_freq * 1e12 # angular phonon freq in Hz
    ϵ_mode = eV^2 * ir_activity / (3 * volume * ω_j^2 * amu) # single ionic mode
    return ϵ_mode / ϵ_0 # normalise with 1 / (4π ϵ_0)
end

# ╔═╡ 053d717a-2b63-40e7-9dc0-487e3892de98
function ϵ_total(freqs_and_ir_activity, volume) # total ionic contribution to dielectric
    phonon_freqs = freqs_and_ir_activity[:, 1] 
    ir_activity = freqs_and_ir_activity[:, 2]
    result = 0.0
    for (f, r) in zip(phonon_freqs, ir_activity)
        result += ϵ_ionic_mode(f, r, volume) # sum over all ionic contributions
    end
    return result
end

# ╔═╡ 86fe8a15-4d9e-48fa-8387-3cfaff7aa182
function frohlich_α_j(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff) # Frohlich alpha decomposed into phonon branch contributions
    Ry = eV^4 * m_e / (2 * ħ^2) # Rydberg energy
    ω = 2π * 1e12 * phonon_mode_freq # angular phonon freq (Hz)
    ϵ_static = ϵ_total + ϵ_optic # static dielectric. Calculate here instead of input so that ionic modes properly normalised.
    return (m_eff * Ry / (ħ * ω))^(1 / 2) * ϵ_ionic / (4π * ϵ_0) / (ϵ_optic * ϵ_static) # 1 / (4π ϵ_0) dielectric normalisation
end

# ╔═╡ 16e7d781-bc61-44b3-84a2-e8c8892bf773
function κ_i(i, v, w) # fictitious spring constant, multiple variational params
    κ = v[i]^2 - w[i]^2
    if length(v) > 1
        for j in 1:length(v)
            if j != i
                κ *= (v[j]^2 - w[i]^2) / (w[j]^2 - w[i]^2)
            end
        end
    end
    return κ
end

# ╔═╡ a78017d7-8835-4a8f-8dd3-f3a4a50eb37b
function h_i(i, v, w) # some vector relating to harmonic eigenmodes
    h = v[i]^2 - w[i]^2
    if length(v) > 1
        for j in 1:length(v)
            if j != i
                h *= (w[j]^2 - v[i]^2) / (v[j]^2 - v[i]^2)
            end
        end
    end
    return h
end

# ╔═╡ d7346684-6b68-4dfc-81da-12c435c4acf0
function C_ij(i, j, v, w) # generalised Feynman C variational parameter (inclusive of multiple v and w params)
    C = w[i] * κ_i(i, v, w) * h_i(j, v, w) / (4 * (v[j]^2 - w[i]^2))
    return C
end

# ╔═╡ 798a4f7a-a2a0-48de-82ad-f85b4199496a
function D_j(τ, β, v, w) # log of dynamic structure factor for polaron 
    D = τ * (1 - τ / β)
    for i in 1:length(v)
        if v[i] != w[i]
        D += (h_i(i, v, w) / v[i]^2) * (2 * sinh(v[i] * τ / 2) * sinh(v[i] * (β - τ) / 2) / (v[i] * sinh(v[i] * β / 2)) - τ * (1 - τ / β))
        end
    end
    return D
end

# ╔═╡ ec0744c0-2530-11ec-04c3-1b5d83413b78
function multi_free_energy(v_params, w_params, T, ϵ_optic, m_eff, volume, freqs_and_ir_activity)

    setprecision(BigFloat, 32) # Speed up. Stops potential overflows.

    # Extract phonon frequencies and ir activities.
    phonon_freqs = freqs_and_ir_activity[:, 1]
    ir_activity = freqs_and_ir_activity[:, 2]

    num_of_branches = length(phonon_freqs)
    
    # total dielectric contribution from all phonon branches (used as a normalisation)
    ϵ_tot = ϵ_total(freqs_and_ir_activity, volume)

    # Generalisation of B i.e. Equation 62c in Hellwarth.
    B_j_integrand(τ, β, v, w) = cosh(β / 2 - abs(τ)) / (sinh(β / 2) * sqrt(D_j(abs(τ), β, v, w)))
    B_j(β, α, v, w) = α / √π * quadgk(τ -> B_j_integrand(τ, β, v, w), 0.0, β / 2)[1]

    # Generalisation of C i.e. Equation 62e in Hellwarth.
    function C_j(β, v, w)
        s = 0.0
        for i in 1:length(v)
            for j in 1:length(v)
                s += C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2)  - 2 / (β * v[j]))
            end
        end
        3 * s / num_of_branches
    end

    # Generalisation of A i.e. Equation 62b in Hellwarth.
    function A_j(β, v, w)
        s = -log(2π * β) / 2
        for i in 1:length(v)
            if v[i] != w[i]
                s += log(v[i] / w[i]) - log(sinh(v[i] * β / 2) / sinh(w[i] * β / 2))
            end
        end
        3 / β * s / num_of_branches
    end
	
	F = 0.0
	for j in 1:num_of_branches
		ω_j = 2π * 1e12 * phonon_freqs[j] # angular phonon freq im 2π Hz
		β_j = BigFloat(ħ * ω_j / (k_B * T)) # reduced thermodynamic beta
		ϵ_ionic_j = ϵ_ionic_mode(phonon_freqs[j], ir_activity[j], volume) # ionic dielectric contribution for current phonon branch
		α_j = frohlich_α_j(ϵ_optic, ϵ_ionic_j, ϵ_tot, phonon_freqs[j], m_eff) # decomposed alpha for current phonon branch

		# F = -(A + B + C) in Hellwarth.
		F += -(B_j(β_j, α_j, v_params, w_params) + C_j(β_j, v_params, w_params) + A_j(β_j, v_params, w_params)) * ω_j # × ħω branch phonon energy
	end
		
    return F * ħ / eV * 1e3 # change to meV
end

# ╔═╡ f3b91e16-26d3-45c8-b28f-193afac15d7f
function multi_variation(T, ϵ_optic, m_eff, volume, freqs_and_ir_activity; initial_vw = false, N = 1) # N number of v and w params

    setprecision(BigFloat, 32) # Speed up. Stops potential overflows.

    if initial_vw isa Bool
		# Intial guess for v and w.
    	initial = sort(rand(2 * N)) .* 4.0 .+ 1.0 # initial guess around 4 and ≥ 1.
		
		# Limits of the optimisation.
		lower = repeat([0.1], 2 * N)
		upper = repeat([60.0], 2 * N)
	else
		# Intial guess for v and w.
		initial = sort(vcat(initial_vw...))
		
		# Limits of the optimisation.
		lower = repeat([0.1], 2 * N)
		upper = repeat([60.0], 2 * N)
	end
	
	println("Initial guess: ", initial)

	# Osaka Free Energy function to minimise.
	f(x) = multi_free_energy([x[2 * n] for n in 1:Int(N)], [x[2 * n - 1] for n in 1:Int(N)], T, ϵ_optic, m_eff, volume, freqs_and_ir_activity)

	# Use Optim to optimise the free energy function w.r.t v and w.
	solution = Optim.optimize(
		Optim.OnceDifferentiable(f, initial; autodiff = :forward),
		lower,
		upper,
		initial,
		Fminbox(BFGS()),
		Optim.Options(time_limit = 20.0),
	)

	# Get v and w params that minimised free energy.
	var_params = Optim.minimizer(solution)

	# Update matrices for v and w parameters.
	v_params = [var_params[2 * n] for n in 1:N]
	w_params = [var_params[2 * n - 1] for n in 1:N]

	# Show current v and w that minimise jth phonon branch.
	println("Variational parameters: ", var_params)

    # Return variational parameters that minimise the free energy.
    return v_params, w_params
end

# ╔═╡ a941ccc6-90cb-448a-91bb-f80233b630f8
function HellwarthAScheme(phonon_modes; T = 295)
	   
	phonon_mode_freqs = phonon_modes[:, 1]
	ir_activities = phonon_modes[:, 2]
	
	condition(f) = coth(π * f * 1e12 * ħ / (k_B * T)) / f - sum(ir_activities .* coth.(π .* phonon_mode_freqs .* 1e12 .* ħ ./ (k_B * T)) ./ phonon_mode_freqs) / sum(ir_activities)
	
	minimum_frequency = minimum(phonon_mode_freqs)
	maximum_frequency = maximum(phonon_mode_freqs)
	middle_frequency = (maximum_frequency + minimum_frequency) / 2
	print("\n")
			
	while (maximum_frequency - minimum_frequency) / 2 > 1e-6
		
		if sign(condition(middle_frequency)) == sign(condition(minimum_frequency))
			minimum_frequency = middle_frequency
			middle_frequency = (maximum_frequency + minimum_frequency) / 2
		else
			maximum_frequency = middle_frequency
			middle_frequency = (maximum_frequency + minimum_frequency) / 2
		end
	end
	
	return middle_frequency
end

# ╔═╡ 1372f247-4622-4388-b458-89aa6e695388
function HellwarthBScheme(LO)
    println("Hellwarth B Scheme... (athermal)")
    H58 = sum( LO[:,2] ./ LO[:,1].^2 )
    println("Hellwarth (58) summation: ",H58)

    H59 = sum( LO[:,2] ) # sum of total ir activity squarred
    println("Hellwarth (59) summation (total ir activity ^2): ", H59)
    println("Hellwarth (59) W_e (total ir activity ): ", sqrt(H59))

    omega = sqrt(H59 / H58)
    println("Hellwarth (61) Omega (freq): ",omega)

	return(omega)
end

# ╔═╡ c28f6d35-c3ed-4984-9c36-b282fb50001b
function ridgeline(x, y, z; shift = 1, jump = 10, ylabel = :none, xlabel = :none, palette = :twilight, linewidth = 1.5, color = :black, size = (595, 842), ymirror = false, fillalpha = 1, step = 1)
	if !ymirror
		p = plot(x, z[:, 1], fillrange = z[:, jump] .- shift, fillalpha = fillalpha, ylabel = ylabel, xlabel = xlabel, xticks = (0:maximum(x)/10:maximum(x), string.(x[1:Int(ceil((length(x)-1)/10)):end])), yticks = (z[1, 1:jump * step:end] .- shift .* range(0, stop = length(y[1:jump:end]) - 1, step = step), string.(y[1:jump * step:end])), legend = false, size = size, grid = false, ymirror = ymirror, palette = palette, tickfontsize = 12, labelfontsize = 12)
	else 
		p = plot(x, z[:, 1], fillrange = z[:, jump] .- shift, fillalpha = fillalpha, ylabel = ylabel, xlabel = xlabel, xticks = (0:maximum(x)/10:maximum(x), string.(x[1:Int((length(x))/10):end])), yticks = (z[end, 1:jump * step:end] .- shift .* range(0, stop = length(y[1:jump:end]) - 1, step = step), string.(y[1:jump * step:end])), legend = false, size = size, grid = false, ymirror = ymirror, palette = palette, tickfontsize = 12, labelfontsize = 12)
	end
	
	plot!(x, z[:, 1], color = color, linewidth = linewidth)
	plot!(x, z[:, jump] .- shift, color = color, linewidth = linewidth)

	for i in 1:Int(floor(length(y)/jump) - 1)
	plot!(x, z[:, i * jump] .- shift * i, fillrange = z[:, (i + 1) * jump] .- shift * (i + 1), fillalpha = fillalpha, palette = palette)
	plot!(x, z[:, i * jump] .- shift * i, color = color, linewidth = linewidth)
	plot!(x, z[:, (i + 1) * jump] .- shift * (i + 1), color = color, linewidth = linewidth)
	end
	return p	
end

# ╔═╡ 49480256-2dcc-4d10-a5a1-4202b243cfd2
T_range = 100:100:400

# ╔═╡ 4966ecf1-ef39-4da3-a8a5-9e47a6c6b0ec
phonon_mode_freqs = MAPI[:, 1]

# ╔═╡ 3321246b-76bd-4678-a5ae-332c3543a1ef
ir_activities = MAPI[:, 2]

# ╔═╡ b2944c8a-158e-461b-a83c-842925587179
A_scheme = [HellwarthAScheme(MAPI, T = i) for i in T_range]

# ╔═╡ fe4e5ef5-833a-4b61-9ea9-92fa52d7a21f
B_scheme = HellwarthBScheme(hcat(phonon_mode_freqs, ir_activities))

# ╔═╡ 8350c7a2-e9fd-4312-9561-44e2524734af
ϵ_ionic = [ϵ_ionic_mode(f, r, (6.29e-10)^3) for (f, r) in zip(phonon_mode_freqs, ir_activities)]

# ╔═╡ db2b4169-0737-44c1-9af5-aa2f2403f5db
ϵ_tot = sum(ϵ_ionic)

# ╔═╡ 2b6d00e1-317b-41cc-9489-14a1450012f4
α_j = [frohlich_α_j(4.5, ϵ_i, ϵ_tot, f, 0.12) for (ϵ_i, f) in zip(ϵ_ionic, phonon_mode_freqs)]

# ╔═╡ 3e1418c2-01a1-4651-8029-59f53264193d
begin
	scatter(phonon_mode_freqs, α_j, yaxis = :log, xticks = (phonon_mode_freqs, hcat(["$(round(i, digits = 2))" for i in phonon_mode_freqs]...)), label = L"\textbf{\alpha_j}", size = (500, 500), xlabel = "Phonon Frequencies (THz)", tickfontsize = 9, xrotation = 90)
	scatter!(phonon_mode_freqs, ir_activities, markershape = :diamond, label = "IR")
	plot!(phonon_mode_freqs, α_j ./ ir_activities, linewidth = 2, linestyle = :dash, label = "ratio")
end

# ╔═╡ d4f392ef-968a-439c-ba53-5707464446e7
begin
	multi_data = DataFrames.DataFrame(
			alpha = α_j,
			phonon_freqs = phonon_mode_freqs,
			ir_activities = ir_activities,
			ionic = ϵ_ionic
	    )
	CSV.write("multi_data.csv", multi_data)
end

# ╔═╡ b3c94ea5-d152-4705-9a30-f8d3e4da0874
α_eff = sum(α_j)

# ╔═╡ 3608af5e-7b44-4d50-9770-5e8afd4fca28
α_hellwarth_A = PolaronMobility.frohlichalpha.(4.5, 24.1, A_scheme .* 1e12, 0.12)

# ╔═╡ 1716f263-f8ad-40a1-a805-d8b7ea6f38df
α_hellwarth_B = PolaronMobility.frohlichalpha(4.5, 24.1, B_scheme * 1e12, 0.12)	

# ╔═╡ da629bec-b53c-4acf-a0ab-6057c3b973a9
β_j = [ħ * 2π * 1e12 * f / k_B / T for f in phonon_mode_freqs, T in T_range]

# ╔═╡ 0641f475-0f35-4f5c-abf0-85ec56f4e30c
begin
	multi_beta = DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [phonon_mode_freqs β_j]], :auto)
	CSV.write("multi_beta.csv", multi_beta)
end

# ╔═╡ c49a4dc3-2546-472f-9fde-b59ec3246237
β_j_avg = [sum(β_j[:, T]) / length(phonon_mode_freqs) for T in 1:length(T_range)]

# ╔═╡ 0e3cf449-6e2a-4353-853e-58d15693743d
β_hellwarth_A = [ħ * 2π / k_B / T_range[T] * A_scheme[T] * 1e12 for T in 1:length(T_range)]

# ╔═╡ 16040051-57d5-4c6c-aa13-8280ae964261
β_hellwarth_B = [ħ * 2π / k_B / T * B_scheme * 1e12 for T in T_range]

# ╔═╡ 32f32534-2a6b-4675-8508-890579058e07
begin
	var_params = []
	push!(var_params, multi_variation(T_range[1], 4.5, 0.12, (6.29e-10)^3, MAPI; N = 1))
	for i in 2:length(T_range)
		push!(var_params, multi_variation(T_range[i], 4.5, 0.12, (6.29e-10)^3, MAPI; initial_vw = var_params[i-1], N = 1))
	end
	var_params
end

# ╔═╡ de89f0aa-21fb-42d2-9a75-5b45938b1e65
var_A = PolaronMobility.feynmanvw.(α_hellwarth_A, β_hellwarth_A)

# ╔═╡ a74836cf-f9ab-41d6-9f01-5066f1afac34
var_B = PolaronMobility.feynmanvw.(α_hellwarth_B, β_hellwarth_B)

# ╔═╡ 2ea4fc88-5285-4c8e-baf9-f0b2f75b33e0
v_j = [var_params[T][1] for T in 1:length(T_range)]

# ╔═╡ e5b21223-f7a5-4560-b03d-f015cdc99c41
sum.(v_j) ./ 2

# ╔═╡ ff8beae7-50e2-46cb-a38d-eb65a8cedaa4
w_j = [var_params[T][2] for T in 1:length(T_range)]

# ╔═╡ 757c41b7-e32b-4a67-b02d-c3f25acc266f
v_A = [i[1] for i in var_A]

# ╔═╡ 368ee8d7-84a8-4711-a9f9-f08d3c03c0ec
w_A = [i[2] for i in var_A]

# ╔═╡ b032bda5-4920-4546-a80f-5491b4044f45
v_B = [i[1] for i in var_B]

# ╔═╡ 5dd88a96-75d2-409a-abf8-960c09c74188
w_B = [i[2] for i in var_B]

# ╔═╡ 83ec77a6-fcf8-4025-ac6d-cd7819312427
begin
	vw_temp = plot(T_range, v_A, label = "H-A v", legend = :bottomright, linewidth = 2, tickfontsize = 12, labelfontsize = 12, legendfontsize = 12, xlabel = "T(K)", ylabel = "v & w (THz)", size = (550, 550), minorgrid = true, linestyle = :dashdot, color = theme_palette(:default)[1])
	plot!(T_range, w_A, label = "H-A w", linewidth = 2, linestyle = :dashdot, color = theme_palette(:default)[2])
	plot!(T_range, v_B, label = "H-B v", linewidth = 2, linestyle = :dash, color = theme_palette(:default)[3])
	plot!(T_range, w_B, label = "H-B w", linewidth = 2, linestyle = :dash, color = theme_palette(:default)[4])
	plot!(T_range, sum.(v_j) ./ length(v_j[1]), label = "multi v", linewidth = 2, linestyle = :solid, color = theme_palette(:default)[1])
	plot!(T_range, sum.(w_j) ./ length(w_j[1]), label = "multi w", linewidth = 2, linestyle = :solid, color = theme_palette(:default)[2])
end

# ╔═╡ 3a52ff01-110b-4755-b43a-acc5b3b870aa
savefig(vw_temp, "vw_temp.pdf")

# ╔═╡ 8fafb718-d95a-4ca7-921f-4abe1ab846ac
F_j = [multi_free_energy(v_j[T], w_j[T], T_range[T], 4.5, 0.12, (6.29e-10)^3, MAPI) for T in 1:length(T_range)]

# ╔═╡ 552b01e3-fd72-489a-aa2b-ae00b6d37c05
Hellwarth_energy_A = F.(v_A, w_A, β_hellwarth_A, α_hellwarth_A) .* 1e3 .* ħ .* 2π .* A_scheme .* 1e12 ./ eV

# ╔═╡ a465eeb5-5b8f-4202-bd58-651867c9d423
Hellwarth_energy_B = F.(v_B, w_B, β_hellwarth_B, α_hellwarth_B) .* 1e3 .* ħ .* 2π .* B_scheme .* 1e12 ./ eV

# ╔═╡ 6b60ab62-6706-4c33-ad7a-efb11803eb6b
begin
	free_energy_temp = plot(T_range, Hellwarth_energy_A, linewidth = 2, linestyle = :dashdot, label = "H-A", tickfontsize = 12, labelfontsize = 12, legendfontsize = 12, xlabel = "T(K)", ylabel = "Free energy (meV)", minorgrid = true, size = (550, 550), legend = :topright)
	plot!(T_range, Hellwarth_energy_B, linewidth = 2, linestyle = :dash, label = "H-B")
	plot!(T_range, F_j, linewidth = 2, linestyle = :solid, label = "multi")
end

# ╔═╡ 7ae9767d-d18f-4188-8e59-dea17c9d5ac4
savefig(free_energy_temp, "free_energy_temp.pdf")

# ╔═╡ 523af3a6-7cc9-41cf-9acb-9a1d25086969
begin
	multi_vw = DataFrames.DataFrame(
			temperature = T_range,
			beta = β_j_avg,
			v = v_j,
			w = w_j,
			F = F_j
	    )
	CSV.write("multi_vw.csv", multi_vw)
end

# ╔═╡ 8db8a4fd-fc04-4d45-8d05-8f69eff611df
begin
	A_data = DataFrames.DataFrame(
		alpha = α_hellwarth_A,
		efffreq = A_scheme,
		temp = T_range,
		beta = β_hellwarth_A,
		v = v_A,
		w = w_A,
		F = Hellwarth_energy_A
    ) 
	CSV.write("A_data.csv", A_data)
end

# ╔═╡ 0e921a7a-a07d-4afc-a991-e06d9be96e3f
begin
	B_data = DataFrames.DataFrame(
			alpha = [α_hellwarth_B for i in 1:length(T_range)],
			efffreq = [B_scheme for i in 1:length(T_range)],
			temp = T_range,
			beta = β_hellwarth_B,
			v = v_B,
			w = w_B,
			F = Hellwarth_energy_B
		) 
	CSV.write("B_data.csv", B_data)
end

# ╔═╡ 2d4759ab-18d8-4368-a568-0da25474366c
function multi_conductivity(ν, β, α, v, w, ω, m_eff)
	z_integrand(t, β, ν) = (1 - exp(1im * 2π * ν * t)) * imag(cos(t - 1im * β / 2) / (  sinh(β / 2) * D_j(-1im * t, β, v, w)^(3/2)))
	z = 0.0
	for j in 1:length(ω)
		println("Photon frequency = $ν, Phonon mode frequency = $(ω[j] / 2π)")
		z += -1im * 2π * ν / length(ω) + 1im * 2 * α[j] * ω[j]^2 * quadgk(t -> z_integrand(t, β[j], ν / ω[j]), 0.0, Inf)[1] / (3 * √π * 2π * ν)
	end
	1 / z * eV * 100^2 / (m_eff * m_e * 1e12)
end

# ╔═╡ 2285006a-a115-414a-b788-b05fd6720bcd
ν_range = 0.01:0.01:3.51

# ╔═╡ 0fdcd049-8f06-4651-9864-a76647154680
begin
	σ_j = Array{ComplexF64}(undef, length(ν_range), length(T_range))
	for j in 1:length(T_range), i in 1:length(ν_range)
		println("ν: $(ν_range[i]) THz, T: $(T_range[j]) K")
		σ_j[i, j] = multi_conductivity(ν_range[i], β_j[:, j], α_j, v_j[1], w_j[1], phonon_mode_freqs .* 2π, 0.12)
		CSV.write("multi_conductivity_data.csv", DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [ν_range σ_j]], :auto))
	end
	σ_j
end

# ╔═╡ f09c4429-9415-49c1-b22d-a8867cbd4249
begin
	plot(ν_range, real.(σ_j), linewidth = 2, label = "Reσ", legend = :outerright, minorgrid = true)
	plot!(ν_range, imag.(σ_j), linewidth = 2, linestyle = :dash, label = "Imσ")
	plot!(ν_range, abs.(σ_j), linewidth = 2, linestyle = :dot, label = "|σ|")
end

# ╔═╡ 0178c326-eeca-4d7a-a111-58b589cba561
multi_contour_real = Plots.contourf(ν_range, T_range, real.(σ_j)', xlabel = "ν (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 600))

# ╔═╡ d2bbb4c9-98e5-4959-a36c-b3e48c37a03c
multi_contour_imag = Plots.contourf(ν_range, T_range, imag.(σ_j)', xlabel = "ν (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 600))

# ╔═╡ 0d78d4ef-6360-43b3-9952-3467410d8908
multi_contour_abs = Plots.contourf(ν_range, T_range, abs.(σ_j)', xlabel = "ν (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 600))

# ╔═╡ 8822d1c3-5bda-4e4e-9b0b-703c75d10816
begin
	savefig(multi_contour_real, "multi_contour_real.pdf")
	savefig(multi_contour_imag, "multi_contour_imag.pdf")
	savefig(multi_contour_abs, "multi_contour_abs.pdf")
end

# ╔═╡ 72d02e9b-d59e-498b-8947-3a86634b7eb4
multi_plot_temp_real = ridgeline(round.(ν_range, digits = 1), Int.(T_range), [i < 600 ? i : 600 for i in real.(σ_j)], jump = 5, shift = 10, xlabel = "ν (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ 0916b2de-ff08-4c60-9bb1-70efbea6f6f6
multi_plot_temp_imag = ridgeline(round.(ν_range, digits = 1), Int.(T_range), [i < 600 ? i : 600 for i in imag.(σ_j)], jump = 5, shift = 10, xlabel = "ν (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ dfd78449-6db6-4f6b-9ecf-00a202ad685d
multi_plot_temp_abs = ridgeline(round.(ν_range, digits = 1), Int.(T_range), [i < 600 ? i : 600 for i in abs.(σ_j)], jump = 5, shift = 10, xlabel = "ν (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ 7d101426-86d9-4925-b4f8-52c2d5415491
begin
	savefig(multi_plot_temp_real, "multi_plot_temp_real.pdf")
	savefig(multi_plot_temp_imag, "multi_plot_temp_imag.pdf")
	savefig(multi_plot_temp_abs, "multi_plot_temp_abs.pdf")
end

# ╔═╡ ad91591f-1ccf-4347-8156-1d0e2e090fdb
begin
	σ_hellwarth_A = Array{ComplexF64}(undef, length(ν_range), length(T_range))
	for j in 1:length(T_range)
		for i in 1:length(ν_range)
			println("ν: $(ν_range[i]) THz, T: $(T_range[j]) K")
			σ_hellwarth_A[i, j] = multi_conductivity(ν_range[i], β_hellwarth_A[j], α_hellwarth_A[j], v_A[1], w_A[1], A_scheme[j] * 2π, 0.12)
		end
		CSV.write("A_conductivity_data.csv", DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [ν_range σ_hellwarth_A]], :auto))
	end
	σ_hellwarth_A
end

# ╔═╡ 2c4a84f5-6839-4ea4-bc9c-4eb2f8e6281b
begin
	term = 40
	plot(ν_range, real.(σ_hellwarth_A)[:, term], linewidth = 2, label = "Reσ", legend = :outerright, minorgrid = true, ylims = (0, 500))
	plot!(ν_range, imag.(σ_hellwarth_A)[:, term], linewidth = 2, linestyle = :dash, label = "Imσ")
	plot!(ν_range, abs.(σ_hellwarth_A)[:, term], linewidth = 2, linestyle = :dot, label = "|σ|")
end

# ╔═╡ ee677bc3-2e9e-4e35-a66c-4c6052c4e32e
A_contour_real = Plots.contourf(ν_range, T_range, real.(σ_hellwarth_A)', xlabel = "ν (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 600))

# ╔═╡ 87154df5-6996-4b4b-9f63-45606205b945
A_contour_imag = Plots.contourf(ν_range, T_range, imag.(σ_hellwarth_A)', xlabel = "ν (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (-6, 600))

# ╔═╡ f4f84543-d7b8-43b2-93cb-0de8d7fa0a30
A_contour_abs = Plots.contourf(ν_range, T_range, abs.(σ_hellwarth_A)', xlabel = "ν (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true, scale = :log), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 600))

# ╔═╡ 005705cd-a857-431a-8956-8ba0035c7c05
begin
	savefig(A_contour_real, "A_contour_real.pdf")
	savefig(A_contour_imag, "A_contour_imag.pdf")
	savefig(A_contour_abs, "A_contour_abs.pdf")
end

# ╔═╡ 77a6559f-6c40-4673-9a8d-9e0c8fe3476f
A_plot_temp_real = ridgeline(round.(ν_range, digits = 1), Int.(T_range), [i < 600 ? i : 600 for i in abs.(σ_hellwarth_A)], jump = 5, shift = 10, xlabel = "ν (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ 498c610b-0646-4a84-ac40-ab8c6efc5278
A_plot_temp_imag = ridgeline(round.(ν_range, digits = 1), Int.(T_range), [i < 600 ? i : 600 for i in imag.(σ_hellwarth_A)], jump = 5, shift = 10, xlabel = "ν (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ bc47027d-21a0-4a02-914d-aba2381d2e01
A_plot_temp_abs = ridgeline(round.(ν_range, digits = 1), Int.(T_range), [i < 600 ? i : 600 for i in abs.(σ_hellwarth_A)], jump = 5, shift = 10, xlabel = "ν (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ 1e902322-d319-42c5-889c-c3048515674c
begin
	savefig(A_plot_temp_real, "A_plot_temp_real.pdf")
	savefig(A_plot_temp_imag, "A_plot_temp_imag.pdf")
	savefig(A_plot_temp_abs, "A_plot_temp_abs.pdf")
end

# ╔═╡ 905efcbf-b3f5-4da0-a45c-3204cd69c0ff
begin
	σ_hellwarth_B = Array{ComplexF64}(undef, length(ν_range), length(T_range))
	for j in 1:length(T_range)
		for i in 1:length(ν_range)
			println("ν: $(ν_range[i]) THz, T: $(T_range[j]) K")
			σ_hellwarth_B[i, j] = multi_conductivity(ν_range[i], β_hellwarth_B[j], α_hellwarth_B, v_B[1], w_B[1], B_scheme * 2π, 0.12)
		end
		CSV.write("B_conductivity_data.csv", DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [ν_range σ_hellwarth_B]], :auto))
	end
	σ_hellwarth_B
end

# ╔═╡ 38a23177-633a-4b2d-8083-5e8afd74796a
begin
	term_B = 2
	plot(ν_range, real.(σ_hellwarth_B)[:, term_B], linewidth = 2, label = "Reσ", legend = :outerright, minorgrid = true, ylims = (0, 600))
	plot!(ν_range, imag.(σ_hellwarth_B)[:, term_B], linewidth = 2, linestyle = :dash, label = "Imσ")
	plot!(ν_range, abs.(σ_hellwarth_B)[:, term_B], linewidth = 2, linestyle = :dot, label = "|σ|")
end

# ╔═╡ 64821290-9bbe-4a7b-af4c-79825f6d4015
B_contour_real = Plots.contourf(ν_range, T_range, real.(σ_hellwarth_B)', xlabel = "ν (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 600))

# ╔═╡ d292b963-3fd1-4d05-b72b-d71b66c5a919
B_contour_imag = Plots.contourf(ν_range, T_range, imag.(σ_hellwarth_B)', xlabel = "ν (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (-6, 600))

# ╔═╡ bfdc8a9b-9522-4f48-8b17-da79447994a8
B_contour_abs = Plots.contourf(ν_range, T_range, abs.(σ_hellwarth_B)', xlabel = "ν (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 600))

# ╔═╡ 273c8ea1-0f7c-42b4-8954-67671d157181
begin
	savefig(B_contour_real, "B_contour_real.pdf")
	savefig(B_contour_imag, "B_contour_imag.pdf")
	savefig(B_contour_abs, "B_contour_abs.pdf")
end

# ╔═╡ b96343f4-bb61-4f4b-9ec9-ad2f37a96cfa
B_plot_temp_real = ridgeline(round.(ν_range, digits = 1), Int.(T_range), [i < 600 ? i : 600 for i in real.(σ_hellwarth_A)], jump = 5, shift = 10, xlabel = "ν (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ b4ee14e9-a89f-4a29-b6ed-d777fe477b9a
B_plot_temp_imag = ridgeline(round.(ν_range, digits = 1), Int.(T_range), [i < 600 ? i : 600 for i in imag.(σ_hellwarth_B)], jump = 5, shift = 10, xlabel = "ν (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ f711553f-e567-437f-9c38-79e644bd18b9
B_plot_temp_abs = ridgeline(round.(ν_range, digits = 1), Int.(T_range), [i < 600 ? i : 600 for i in abs.(σ_hellwarth_B)], jump = 5, shift = 10, xlabel = "ν (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ 26e1a5da-6bdb-49b9-973b-768d901a2a60
begin
	savefig(B_plot_temp_real, "B_plot_temp_real.pdf")
	savefig(B_plot_temp_imag, "B_plot_temp_imag.pdf")
	savefig(B_plot_temp_abs, "B_plot_temp_abs.pdf")
end

# ╔═╡ Cell order:
# ╠═29d62408-0629-4f39-909e-a6ba8696e955
# ╠═1eb559b3-dbce-4caf-800f-96a2495dee1b
# ╠═77d14c3f-bcaf-432c-9005-f5ec184d84b7
# ╠═ba2dadc8-d9a2-452e-9f94-406ccfc7404b
# ╠═c05cd489-68c8-418e-bbc2-a33bc42781f0
# ╠═1d300e7d-5efb-4ad2-a826-c3f8a86083d8
# ╠═6a99cd13-68cb-474c-96af-7d78ee4e7095
# ╠═49663b33-1604-430b-a815-9eb985802987
# ╠═314af942-8922-40d4-a92b-f75b715c276c
# ╠═6ff68ba5-df0c-42a0-8cf1-dedeced77e8d
# ╠═35390bdf-e260-4566-8a97-7a0b5a3ef07a
# ╠═e0ce47a5-368a-4da7-83fc-9f921af0f747
# ╠═053d717a-2b63-40e7-9dc0-487e3892de98
# ╠═86fe8a15-4d9e-48fa-8387-3cfaff7aa182
# ╠═16e7d781-bc61-44b3-84a2-e8c8892bf773
# ╠═a78017d7-8835-4a8f-8dd3-f3a4a50eb37b
# ╠═d7346684-6b68-4dfc-81da-12c435c4acf0
# ╠═798a4f7a-a2a0-48de-82ad-f85b4199496a
# ╠═ec0744c0-2530-11ec-04c3-1b5d83413b78
# ╠═f3b91e16-26d3-45c8-b28f-193afac15d7f
# ╠═a941ccc6-90cb-448a-91bb-f80233b630f8
# ╠═1372f247-4622-4388-b458-89aa6e695388
# ╠═c28f6d35-c3ed-4984-9c36-b282fb50001b
# ╠═49480256-2dcc-4d10-a5a1-4202b243cfd2
# ╠═4966ecf1-ef39-4da3-a8a5-9e47a6c6b0ec
# ╠═3321246b-76bd-4678-a5ae-332c3543a1ef
# ╠═b2944c8a-158e-461b-a83c-842925587179
# ╠═fe4e5ef5-833a-4b61-9ea9-92fa52d7a21f
# ╠═8350c7a2-e9fd-4312-9561-44e2524734af
# ╠═db2b4169-0737-44c1-9af5-aa2f2403f5db
# ╠═2b6d00e1-317b-41cc-9489-14a1450012f4
# ╠═3e1418c2-01a1-4651-8029-59f53264193d
# ╠═d4f392ef-968a-439c-ba53-5707464446e7
# ╠═b3c94ea5-d152-4705-9a30-f8d3e4da0874
# ╠═3608af5e-7b44-4d50-9770-5e8afd4fca28
# ╠═1716f263-f8ad-40a1-a805-d8b7ea6f38df
# ╠═da629bec-b53c-4acf-a0ab-6057c3b973a9
# ╠═0641f475-0f35-4f5c-abf0-85ec56f4e30c
# ╠═c49a4dc3-2546-472f-9fde-b59ec3246237
# ╠═0e3cf449-6e2a-4353-853e-58d15693743d
# ╠═16040051-57d5-4c6c-aa13-8280ae964261
# ╠═32f32534-2a6b-4675-8508-890579058e07
# ╠═de89f0aa-21fb-42d2-9a75-5b45938b1e65
# ╠═a74836cf-f9ab-41d6-9f01-5066f1afac34
# ╠═2ea4fc88-5285-4c8e-baf9-f0b2f75b33e0
# ╠═e5b21223-f7a5-4560-b03d-f015cdc99c41
# ╠═ff8beae7-50e2-46cb-a38d-eb65a8cedaa4
# ╠═757c41b7-e32b-4a67-b02d-c3f25acc266f
# ╠═368ee8d7-84a8-4711-a9f9-f08d3c03c0ec
# ╠═b032bda5-4920-4546-a80f-5491b4044f45
# ╠═5dd88a96-75d2-409a-abf8-960c09c74188
# ╠═83ec77a6-fcf8-4025-ac6d-cd7819312427
# ╠═3a52ff01-110b-4755-b43a-acc5b3b870aa
# ╠═8fafb718-d95a-4ca7-921f-4abe1ab846ac
# ╠═552b01e3-fd72-489a-aa2b-ae00b6d37c05
# ╠═a465eeb5-5b8f-4202-bd58-651867c9d423
# ╠═6b60ab62-6706-4c33-ad7a-efb11803eb6b
# ╠═7ae9767d-d18f-4188-8e59-dea17c9d5ac4
# ╠═523af3a6-7cc9-41cf-9acb-9a1d25086969
# ╠═8db8a4fd-fc04-4d45-8d05-8f69eff611df
# ╠═0e921a7a-a07d-4afc-a991-e06d9be96e3f
# ╠═2d4759ab-18d8-4368-a568-0da25474366c
# ╠═2285006a-a115-414a-b788-b05fd6720bcd
# ╠═0fdcd049-8f06-4651-9864-a76647154680
# ╠═f09c4429-9415-49c1-b22d-a8867cbd4249
# ╠═0178c326-eeca-4d7a-a111-58b589cba561
# ╠═d2bbb4c9-98e5-4959-a36c-b3e48c37a03c
# ╠═0d78d4ef-6360-43b3-9952-3467410d8908
# ╠═8822d1c3-5bda-4e4e-9b0b-703c75d10816
# ╠═72d02e9b-d59e-498b-8947-3a86634b7eb4
# ╠═0916b2de-ff08-4c60-9bb1-70efbea6f6f6
# ╠═dfd78449-6db6-4f6b-9ecf-00a202ad685d
# ╠═7d101426-86d9-4925-b4f8-52c2d5415491
# ╠═ad91591f-1ccf-4347-8156-1d0e2e090fdb
# ╠═2c4a84f5-6839-4ea4-bc9c-4eb2f8e6281b
# ╠═ee677bc3-2e9e-4e35-a66c-4c6052c4e32e
# ╠═87154df5-6996-4b4b-9f63-45606205b945
# ╠═f4f84543-d7b8-43b2-93cb-0de8d7fa0a30
# ╠═005705cd-a857-431a-8956-8ba0035c7c05
# ╠═77a6559f-6c40-4673-9a8d-9e0c8fe3476f
# ╠═498c610b-0646-4a84-ac40-ab8c6efc5278
# ╠═bc47027d-21a0-4a02-914d-aba2381d2e01
# ╠═1e902322-d319-42c5-889c-c3048515674c
# ╠═905efcbf-b3f5-4da0-a45c-3204cd69c0ff
# ╠═38a23177-633a-4b2d-8083-5e8afd74796a
# ╠═64821290-9bbe-4a7b-af4c-79825f6d4015
# ╠═d292b963-3fd1-4d05-b72b-d71b66c5a919
# ╠═bfdc8a9b-9522-4f48-8b17-da79447994a8
# ╠═273c8ea1-0f7c-42b4-8954-67671d157181
# ╠═b96343f4-bb61-4f4b-9ec9-ad2f37a96cfa
# ╠═b4ee14e9-a89f-4a29-b6ed-d777fe477b9a
# ╠═f711553f-e567-437f-9c38-79e644bd18b9
# ╠═26e1a5da-6bdb-49b9-973b-768d901a2a60
