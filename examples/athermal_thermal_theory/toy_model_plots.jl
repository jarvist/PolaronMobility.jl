### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ ec4e1740-0b51-11ec-3fa2-d3c1d67f4572
using Revise

# ╔═╡ 9e67b5cd-f084-4bd4-9eb0-b94311184197
using CSV

# ╔═╡ b2a2fdb9-949a-498c-ac19-dab106582558
using DataFrames

# ╔═╡ 72553f5c-276d-4256-8354-04154017d346
using Plots

# ╔═╡ 67b93cb9-6335-448e-a09e-ead44db9d5c3
using LaTeXStrings

# ╔═╡ d20b171e-1bd3-44aa-b067-7bf1a3c0227e
Ω = real.(parse.(ComplexF64, (CSV.File("conductivity_data_6.csv") |> Tables.matrix)[3:end, 1]))

# ╔═╡ 8c909d72-7f28-49ed-9f6f-6b7f77e4e984
β = real.(parse.(ComplexF64, (CSV.File("conductivity_data_6.csv") |> Tables.matrix)[1, 2:end-1]))

# ╔═╡ 2980129c-4be6-4928-8bfe-b6b13eb0d50a
T = 1 ./ β

# ╔═╡ 29d4a2e3-692b-48b3-9ce0-fca83b55e8d3
v = real.((CSV.File("v_data.csv") |> Tables.matrix)[2:end-1, 2:end])

# ╔═╡ eab11453-2f25-493a-84e4-307231496653
v_0 = real.((CSV.File("v_data.csv") |> Tables.matrix)[end, 2:end])

# ╔═╡ 136a409d-9ea0-40ab-998a-f5bf954315ca
α = 1:12

# ╔═╡ a1bef2ba-0fb0-4a43-a09d-4e2e5ab6e069
w = real.((CSV.File("w_data.csv") |> Tables.matrix)[2:end-1, 2:end])

# ╔═╡ 671bed0e-082a-4b3a-9467-7be0ee2327df
w_0 = real.((CSV.File("w_data.csv") |> Tables.matrix)[end, 2:end])

# ╔═╡ a7431fd9-588a-41f1-b3d8-ef15837ccf75
F = real.((CSV.File("F_data.csv") |> Tables.matrix)[2:end-1, 2:end])

# ╔═╡ f7746913-7afa-4bce-b261-9532bd6f27cb
F_0 = real.((CSV.File("F_data.csv") |> Tables.matrix)[end, 2:end])

# ╔═╡ 780dd0f5-b922-4691-94b5-fa9881381c07
M = real.((CSV.File("M_data.csv") |> Tables.matrix)[2:end-1, 2:end])

# ╔═╡ d3f88ec2-6f32-4429-a34b-4db33fc55433
M_0 = real.((CSV.File("M_data.csv") |> Tables.matrix)[end, 2:end])

# ╔═╡ 39cff408-72c6-42ea-8172-d9d0fb9e94fc
σ = [parse.(ComplexF64, (CSV.File("conductivity_data_$i.csv") |> Tables.matrix)[3:end, 2:end-1]) for i in α]

# ╔═╡ 8e432e19-984e-4727-ab9c-25894e2585d8
σ_0 = hcat([parse.(ComplexF64, (CSV.File("conductivity_data_$i.csv") |> Tables.matrix)[3:end, end]) for i in α]...)

# ╔═╡ 83972098-4232-437e-9d2d-7620c17672f0
σ_0_dc = hcat([parse.(ComplexF64, (CSV.File("conductivity_data_$i.csv") |> Tables.matrix)[2, end]) for i in α]...)

# ╔═╡ d8762608-09eb-4ea0-9a84-5f46673a619d
σ_dc = hcat([parse.(ComplexF64, (CSV.File("conductivity_data_$i.csv") |> Tables.matrix)[2, 2:end-1]) for i in α]...)

# ╔═╡ 9ab20c3a-a0f5-4938-9339-720630fbd48b
dc_σ_real = Plots.contourf(α, reverse(T), reverse(log.(real.(σ_dc)), dims = 1), xlabel = "α", ylabel = "T / ω", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, xticks = (1:12, hcat(["$i" for i in 1:12]...)), yticks = ([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0], hcat([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]...)), yaxis = :log, right_margin = 4Plots.mm)

# ╔═╡ ed78a2eb-3313-4cae-b0d3-2118f5dff33f
savefig(dc_σ_real, "dc_mobility.pdf")

# ╔═╡ f524e982-6208-44bc-86df-ccfa851b0c92
athermal_theory = plot(α, [-F_0, M_0, v_0, w_0], yaxis = :log, legend = :topleft, linewidth = 3, linestyle = :auto, labels = [L"-\textrm{E_{gs}}" L"\textrm{m_p}" L"\textrm{v}" L"\textrm{w}" L"\textrm{\mu_0 (10^{6})}"], xlabel = "α", xticks = (1:12, ["$i" for i in 1:12]), ylabel = "Athermal Theory", tickfontsize = 12, labelfontsize = 12, legendfontsize = 12, size = (550, 550), yticks = ([1, 2, 4, 8, 16, 32, 64, 128, 256], ["1" "2" "4" "8" "16" "32" "64" "128" "256"]), ylims = (0.8, 256))


# ╔═╡ 9f09e259-014e-48b0-bd00-f7e2c9c897b3
plot(T, real.(σ_dc), ylims = (0.001, 10), linestyle = :dash, marker = :circle, xaxis = :log, yaxis = :log, legend = false)

# ╔═╡ ccda1394-45b4-4a9a-8c14-425cb0b81a87
log(v[end, 9])

# ╔═╡ fadad1b3-6da6-4724-be97-469d81641d8f
plot(α, real.(σ_dc)'[:, 1:10:end], linestyle = :dash, marker = :circle, legend = :outerright, yaxis = :log, label = hcat(["$(T[i])" for i in 1:10:length(T)]...))

# ╔═╡ 3918264a-f72b-4b4c-b67f-cd6fa4c6c846
savefig(athermal_theory, "athermal_theory.pdf")

# ╔═╡ c1d56910-a645-40ed-a5d8-0f7d5be19ae7
plot(Ω, hcat([real.(σ_0[:, i]) ./ maximum(real.(σ_0[:, i])) for i in 1:12]...), labels = hcat(["α = $i" for i in 1:12]...), linewidth = 2, legend = :outerright, xlabel = "Ω / ω", ylabel = "σ(Ω) [arb]", minorgrid= true)

# ╔═╡ 19065839-85a3-4768-830e-26f64bf4a278
hcat([real.(σ_0[:, i]) ./ maximum(real.(σ_0[:, i])) for i in 1:12]...)

# ╔═╡ 7c978142-5319-45f8-b854-efac81271fcc
function ridgeline(x, y, z; shift = 1, jump = 10, ylabel = :none, xlabel = :none, palette = :twilight, linewidth = 1.5, color = :black, size = (595, 842), ymirror = false, fillalpha = 1, step = 1)
	if !ymirror
		p = plot(x, z[:, 1], fillrange = z[:, 2 * jump] .- shift, fillalpha = fillalpha, ylabel = ylabel, xlabel = xlabel, xticks = (0:maximum(x)/10:maximum(x), string.(x[1:Int(floor((length(x)-1)/10)):end])), yticks = (z[1, 1:jump * step:end] .- shift .* range(0, stop = length(y[1:jump:end]) - 1, step = step), string.(y[1:jump * step:end])), legend = false, size = size, grid = false, ymirror = ymirror, palette = palette, tickfontsize = 12, labelfontsize = 12)
	else 
		p = plot(x, z[:, 1], fillrange = z[:, 2 * jump] .- 2 * shift, fillalpha = fillalpha, ylabel = ylabel, xlabel = xlabel, xticks = (0:maximum(x)/11:maximum(x), string.(reverse(reverse(x)[1:Int(floor(length(x)/11)):end]))), yticks = (z[end, 1:jump * step:end] .- shift .* range(0, stop = length(y[1:jump:end]) - 1, step = step), string.(y[1:jump * step:end])), legend = false, size = size, grid = false, ymirror = ymirror, palette = palette, tickfontsize = 12, labelfontsize = 12)
	end
	
	plot!(x, z[:, 1], color = color, linewidth = linewidth)
	plot!(x, z[:, 2 * jump] .- shift, color = color, linewidth = linewidth)

	for i in 2:Int(floor(length(y)/jump) - 1)
	plot!(x, z[:, i * jump] .- shift * (i - 1), fillrange = z[:, (i + 1) * jump] .- shift * i, fillalpha = fillalpha, palette = palette)
	plot!(x, z[:, i * jump] .- shift * (i - 1), color = color, linewidth = linewidth)
	plot!(x, z[:, (i + 1) * jump] .- shift * i, color = color, linewidth = linewidth)
	end
	return p	
end

# ╔═╡ b41bfff2-c899-491e-b382-e7e3b136edbd
ridgeline(log.(reverse(T)), α, reverse(log.(real.(σ_dc)), dims = 1), jump = 1, shift = 1, xlabel = "Ω / ω", ylabel = "T / ω", size = (550, 550), linewidth = 1, palette = :thermal, fillalpha = 0.5, step = 1)

# ╔═╡ 22a8f0c4-f15a-46eb-adbe-4ee2575233b0
zero_temp_α = ridgeline(round.(Ω, digits = 1), α, hcat([real.(σ_0[:, i]) ./ maximum(real.(σ_0[:, i])) for i in 1:12]...), jump = 1, shift = 0.2, xlabel = "Ω / ω", ylabel = "α", size = (550, 550), palette = cgrad(:tab20, 12, categorical = true), linewidth = 1.5, fillalpha = 0.35, step = 1)

# ╔═╡ a6e1cd52-3384-44a5-aa15-be0f13691eea
savefig(zero_temp_α, "zero_temp_α.pdf")

# ╔═╡ 0aa5786a-f7f4-43b1-becf-f35a848302db
plot(T, F, labels = hcat(["α = $i" for i in α]...), legend = :outerright, linewidth = 2, linestyle = :auto, xlabel = "T / ω", ylabel = "Free energy / ω", tickfontsize = 12, legendfontsize = 12, labelfontsize = 12, size = (700, 550), minorgrid = true)

# ╔═╡ c1e13990-899e-4c65-81ff-ca3d08a46851
plot(T, v .- w, labels = hcat(["α = $i" for i in α]...), legend = :outerright, linewidth = 2, linestyle = :auto, xlabel = "T / ω", ylabel = "Free energy / ω", tickfontsize = 12, legendfontsize = 12, labelfontsize = 12, size = (700, 550), minorgrid = true, xaxis = :log)

# ╔═╡ c49ef6c5-61b6-4a40-9362-2166c642638d
free_energy = Plots.contourf(α, reverse(T), reverse(F, dims = 1), xlabel = "α", ylabel = "T / ω", fill = cgrad(:thermal, rev = false, categorical = true, scale = :log), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, xticks = (1:12, hcat(["$i" for i in 1:12]...)), yticks = ([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0], hcat([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]...)), yaxis = :log, right_margin = 2Plots.mm)

# ╔═╡ d81b64a3-bf3d-443b-a712-783ae6f7eef3
savefig(free_energy, "free_energy.pdf")

# ╔═╡ 4020f1da-b186-48b9-b652-a480e20707a7
effective_mass = Plots.contourf(α, reverse(T), reverse(log.(M), dims = 1), xlabel = "α", ylabel = "T / ω", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, xticks = (1:12, hcat(["$i" for i in 1:12]...)), yticks = ([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0], hcat([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]...)), yaxis = :log)

# ╔═╡ ba9cc94a-26f3-45c8-b793-cc13ad62a28a
savefig(effective_mass, "effective_mass.pdf")

# ╔═╡ 2e68e769-79ff-4509-a9ec-3917237341bb
v_parameter = Plots.contourf(α, reverse(T), reverse(log.(v), dims = 1), xlabel = "α", ylabel = "T / ω", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, xticks = (1:12, hcat(["$i" for i in 1:12]...)), yticks = ([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0], hcat([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]...)), yaxis = :log)

# ╔═╡ 89b85619-ee78-46bf-8741-76d84741cc9f
savefig(v_parameter, "v_parameter.pdf")

# ╔═╡ db726be5-6572-4781-a407-1826a56dd1da
w_parameter = Plots.contourf(α, reverse(T), reverse(log.(w), dims = 1), xlabel = "α", ylabel = "T / ω", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, xticks = (1:12, hcat(["$i" for i in 1:12]...)), yticks = ([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0], hcat([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]...)), yaxis = :log)

# ╔═╡ a9cea3b1-55d8-4ac8-b140-120761e5daae
savefig(w_parameter, "w_parameter.pdf")

# ╔═╡ e4a79f98-c29a-4b9e-81fe-e9f55da2c71f
σ_contour_real = [Plots.contourf(Ω, reverse(T), reverse(log.(abs.(real.(σ[i])))', dims = 1), yaxis = :log, xlabel = "Ω / ω", ylabel = "T / ω", fill = cgrad(:viridis, rev = false, categorical = true, scale = :log), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, xticks = (0:4:20, hcat(["$i" for i in 0:4:20]...)), yticks = ([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0], hcat([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]...)), right_margin = 7Plots.mm) for i in α]

# ╔═╡ c39d711d-c97e-4c8b-9f8e-304ef2be54de
begin
	for i in α
		savefig(σ_contour_real[i], "conductivity_contour_plots/conductivity_contour_real_$i.pdf")
	end
end

# ╔═╡ f8f1085a-ce7f-4854-a273-f8f9db5f8a74
σ_contour_imag = [Plots.contourf(Ω, reverse(T), reverse(log.(abs.(imag.(σ[i])))', dims = 1), yaxis = :log, xlabel = "Ω / ω", ylabel = "T / ω", fill = cgrad(:viridis, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, xticks = (0:4:20, hcat(["$i" for i in 0:4:20]...)), yticks = ([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0], hcat([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]...)), right_margin = 7Plots.mm) for i in α]

# ╔═╡ c1d2a272-e4d4-4291-b424-f66cfde52135
begin
	for i in α
		savefig(σ_contour_imag[i], "conductivity_contour_plots/conductivity_contour_imag_$i.pdf")
	end
end

# ╔═╡ 67b1453c-e81a-4e7d-a7a5-a46ed60d4217
σ_contour_abs = [Plots.contourf(Ω, reverse(T), reverse(log.(abs.(abs.(σ[i])))', dims = 1), yaxis = :log, xlabel = "Ω / ω", ylabel = "T / ω", fill = cgrad(:viridis, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, xticks = (0:4:20, hcat(["$i" for i in 0:4:20]...)), yticks = ([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0], hcat([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]...)), right_margin = 7Plots.mm) for i in α]

# ╔═╡ 5e4f5ba2-21e7-427c-93d4-2ebf23d204d4
begin
	for i in α
		savefig(σ_contour_abs[i], "conductivity_contour_plots/conductivity_contour_abs_$i.pdf")
	end
end

# ╔═╡ 67248e4a-babc-4b69-8bf2-34d56d1c70f4
σ_plot_temp_real = [ridgeline(round.(Ω, digits = 1), round.(reverse(T), digits = 2), reverse(log.(abs.(real.(σ[i]))), dims = 2), jump = 2, shift = 0.8, xlabel = "Ω / ω", ylabel = "T / ω", size = (550, 550), linewidth = 1, palette = :thermal, fillalpha = 0.75, step = 8) for i in α]

# ╔═╡ 5559b8fd-c862-4fdd-9a8b-4c196d401b46
begin
	for i in α
		savefig(σ_plot_temp_real[i], "conductivity_contour_plots/conductivity_plot_temp_real_$i.pdf")
	end
end

# ╔═╡ 6c86267c-799c-4585-a5ad-a221e8990d16
σ_plot_freq_real = [ridgeline(round.(reverse(T), digits = 2), round.(Ω, digits = 1), reverse(log.(abs.(real.(σ[i]'))), dims = 1), jump = 4, shift = 0.8, ylabel = "Ω / ω", xlabel = "T / ω", size = (600, 600), linewidth = 1, palette = :twilight, fillalpha = 0.75, step = 8, ymirror = true) for i in α]

# ╔═╡ f872eee7-e2a8-4173-a149-3c68597d421a
begin
	for i in α
		savefig(σ_plot_freq_real[i], "conductivity_contour_plots/conductivity_plot_freq_real_$i.pdf")
	end
end

# ╔═╡ 4a31fd0c-3cc4-4680-91a0-cf52e6ca75e3
σ_plot_temp_imag = [ridgeline(round.(Ω, digits = 1), reverse(β[1:end]), reverse(log.(abs.(imag.(σ[i]))), dims = 2), jump = 2, shift = 0.8, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 1, palette = :thermal, fillalpha = 0.75, step = 9) for i in α]

# ╔═╡ 39bbd242-e6e9-419e-9bc4-a5f183c33e58
begin
	for i in α
		savefig(σ_plot_temp_imag[i], "conductivity_contour_plots/conductivity_plot_temp_imag_$i.pdf")
	end
end

# ╔═╡ c052d7ff-58ce-413c-8358-851180759fdb
σ_plot_freq_imag = [ridgeline(round.(reverse(T), digits = 2), round.(Ω, digits = 1), reverse(log.(abs.(imag.(σ[i]'))), dims = 1), jump = 4, shift = 0.8, ylabel = "Ω / ω", xlabel = "T / ω", size = (550, 550), linewidth = 1, palette = :twilight, fillalpha = 0.75, step = 10, ymirror = true) for i in α]

# ╔═╡ e5ebc4f2-be00-4d5b-bad1-a71c6b5cdd93
begin
	for i in α
		savefig(σ_plot_freq_imag[i], "conductivity_contour_plots/conductivity_plot_freq_imag_$i.pdf")
	end
end

# ╔═╡ f21050fa-9718-4c23-aff0-de5474d6c5f6
σ_plot_temp_abs = [ridgeline(round.(Ω, digits = 1), reverse(β[1:end]), reverse(log.(abs.(σ[i])), dims = 2), jump = 2, shift = 0.8, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 1, palette = :thermal, fillalpha = 0.75, step = 8) for i in α]

# ╔═╡ 71d9e25c-eec6-49bb-b43f-1927fb0618c8
begin
	for i in α
		savefig(σ_plot_temp_abs[i], "conductivity_contour_plots/conductivity_plot_temp_abs_$i.pdf")
	end
end

# ╔═╡ e3924eb5-8cfb-4db4-9d79-9f5dc7194a95
σ_plot_freq_abs = [ridgeline(round.(reverse(T), digits = 2), round.(Ω, digits = 1), reverse(log.(abs.(σ[i]')), dims = 1), jump = 4, shift = 0.8, ylabel = "Ω / ω", xlabel = "T / ω", size = (550, 550), linewidth = 1, palette = :twilight, fillalpha = 0.75, step = 8, ymirror = true) for i in α]

# ╔═╡ 87dfdce3-c94c-4f2a-9055-b699d37c1848
begin
	for i in α
		savefig(σ_plot_freq_abs[i], "conductivity_contour_plots/conductivity_plot_freq_abs_$i.pdf")
	end
end

# ╔═╡ Cell order:
# ╠═ec4e1740-0b51-11ec-3fa2-d3c1d67f4572
# ╠═9e67b5cd-f084-4bd4-9eb0-b94311184197
# ╠═b2a2fdb9-949a-498c-ac19-dab106582558
# ╠═72553f5c-276d-4256-8354-04154017d346
# ╠═67b93cb9-6335-448e-a09e-ead44db9d5c3
# ╠═d20b171e-1bd3-44aa-b067-7bf1a3c0227e
# ╠═8c909d72-7f28-49ed-9f6f-6b7f77e4e984
# ╠═2980129c-4be6-4928-8bfe-b6b13eb0d50a
# ╠═29d4a2e3-692b-48b3-9ce0-fca83b55e8d3
# ╠═eab11453-2f25-493a-84e4-307231496653
# ╠═136a409d-9ea0-40ab-998a-f5bf954315ca
# ╠═a1bef2ba-0fb0-4a43-a09d-4e2e5ab6e069
# ╠═671bed0e-082a-4b3a-9467-7be0ee2327df
# ╠═a7431fd9-588a-41f1-b3d8-ef15837ccf75
# ╠═f7746913-7afa-4bce-b261-9532bd6f27cb
# ╠═780dd0f5-b922-4691-94b5-fa9881381c07
# ╠═d3f88ec2-6f32-4429-a34b-4db33fc55433
# ╠═39cff408-72c6-42ea-8172-d9d0fb9e94fc
# ╠═8e432e19-984e-4727-ab9c-25894e2585d8
# ╠═83972098-4232-437e-9d2d-7620c17672f0
# ╠═d8762608-09eb-4ea0-9a84-5f46673a619d
# ╠═9ab20c3a-a0f5-4938-9339-720630fbd48b
# ╠═ed78a2eb-3313-4cae-b0d3-2118f5dff33f
# ╠═b41bfff2-c899-491e-b382-e7e3b136edbd
# ╠═f524e982-6208-44bc-86df-ccfa851b0c92
# ╠═9f09e259-014e-48b0-bd00-f7e2c9c897b3
# ╠═ccda1394-45b4-4a9a-8c14-425cb0b81a87
# ╠═fadad1b3-6da6-4724-be97-469d81641d8f
# ╠═3918264a-f72b-4b4c-b67f-cd6fa4c6c846
# ╠═c1d56910-a645-40ed-a5d8-0f7d5be19ae7
# ╠═22a8f0c4-f15a-46eb-adbe-4ee2575233b0
# ╠═a6e1cd52-3384-44a5-aa15-be0f13691eea
# ╠═19065839-85a3-4768-830e-26f64bf4a278
# ╠═7c978142-5319-45f8-b854-efac81271fcc
# ╠═0aa5786a-f7f4-43b1-becf-f35a848302db
# ╠═c1e13990-899e-4c65-81ff-ca3d08a46851
# ╠═c49ef6c5-61b6-4a40-9362-2166c642638d
# ╠═d81b64a3-bf3d-443b-a712-783ae6f7eef3
# ╠═4020f1da-b186-48b9-b652-a480e20707a7
# ╠═ba9cc94a-26f3-45c8-b793-cc13ad62a28a
# ╠═2e68e769-79ff-4509-a9ec-3917237341bb
# ╠═89b85619-ee78-46bf-8741-76d84741cc9f
# ╠═db726be5-6572-4781-a407-1826a56dd1da
# ╠═a9cea3b1-55d8-4ac8-b140-120761e5daae
# ╠═e4a79f98-c29a-4b9e-81fe-e9f55da2c71f
# ╠═c39d711d-c97e-4c8b-9f8e-304ef2be54de
# ╠═f8f1085a-ce7f-4854-a273-f8f9db5f8a74
# ╠═c1d2a272-e4d4-4291-b424-f66cfde52135
# ╠═67b1453c-e81a-4e7d-a7a5-a46ed60d4217
# ╠═5e4f5ba2-21e7-427c-93d4-2ebf23d204d4
# ╠═67248e4a-babc-4b69-8bf2-34d56d1c70f4
# ╠═5559b8fd-c862-4fdd-9a8b-4c196d401b46
# ╠═6c86267c-799c-4585-a5ad-a221e8990d16
# ╠═f872eee7-e2a8-4173-a149-3c68597d421a
# ╠═4a31fd0c-3cc4-4680-91a0-cf52e6ca75e3
# ╠═39bbd242-e6e9-419e-9bc4-a5f183c33e58
# ╠═c052d7ff-58ce-413c-8358-851180759fdb
# ╠═e5ebc4f2-be00-4d5b-bad1-a71c6b5cdd93
# ╠═f21050fa-9718-4c23-aff0-de5474d6c5f6
# ╠═71d9e25c-eec6-49bb-b43f-1927fb0618c8
# ╠═e3924eb5-8cfb-4db4-9d79-9f5dc7194a95
# ╠═87dfdce3-c94c-4f2a-9055-b699d37c1848
