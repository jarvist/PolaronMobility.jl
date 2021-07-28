### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ ed9ab030-e61d-11eb-32f3-5be220868ae7
using Revise

# ╔═╡ d51f978f-a2c8-410a-8836-f57fc362eed1
using CSV

# ╔═╡ ffa81468-a1cd-44f3-aa5a-c970e56b7f27
using DataFrames

# ╔═╡ 4c98c769-0aed-4b8a-bf92-fb06268100a9
using Plots

# ╔═╡ 88fe20b4-303a-494e-95cc-e8183b2f474e
Ω = real.(parse.(ComplexF64, (CSV.File("C:/Users/neutr/OneDrive - Imperial College London/PhD/Code/PolaronMobility.jl/examples/MAPI_data/multi_conductivity_data.csv") |> Tables.matrix)[2:end, 1]))

# ╔═╡ a76c9858-0fe7-45c2-9b5c-a8cf7f800c08
T = real.(parse.(ComplexF64, (CSV.File("C:/Users/neutr/OneDrive - Imperial College London/PhD/Code/PolaronMobility.jl/examples/MAPI_data/multi_conductivity_data.csv") |> Tables.matrix)[1, 2:end]))

# ╔═╡ 9d039885-65fa-4f43-ac15-ba927a68be60
multi_conductivity_data = parse.(ComplexF64, (CSV.File("C:/Users/neutr/OneDrive - Imperial College London/PhD/Code/PolaronMobility.jl/examples/MAPI_data/multi_conductivity_data.csv") |> Tables.matrix)[2:end, 2:end])

# ╔═╡ 88c443f2-f6ea-464e-9190-aa158a0f7f50
A_conductivity_data = parse.(ComplexF64, (CSV.File("C:/Users/neutr/OneDrive - Imperial College London/PhD/Code/PolaronMobility.jl/examples/MAPI_data/A_conductivity_data.csv") |> Tables.matrix)[2:end, 2:end])

# ╔═╡ db784d36-b859-4022-b585-0730bf081a65
B_conductivity_data = parse.(ComplexF64, (CSV.File("C:/Users/neutr/OneDrive - Imperial College London/PhD/Code/PolaronMobility.jl/examples/MAPI_data/B_conductivity_data.csv") |> Tables.matrix)[2:end, 2:end])

# ╔═╡ 329f1be8-967d-4c5b-8536-f6f1008f2a9c
B_contour = Plots.contourf(Ω, T, log.(real.(B_conductivity_data).* 2π * 2.25)', xlabel = "Ω (THz)", ylabel = "T (K)", clims = (-1.1, 8.3))

# ╔═╡ 26ec12fd-5c2d-458c-8282-b411ef5184a8
savefig(B_contour, "B_conductivity_contour.png")

# ╔═╡ 1e89541e-0635-47d9-96f6-99dceef27d00
multi_contour = Plots.contourf(Ω, T, log.(real.(multi_conductivity_data).* 2π)', xlabel = "Ω (THz)", ylabel = "T (K)", clims = (-1.1, 8.3))

# ╔═╡ 512e223b-249d-48ff-9eec-4282f4658c76
savefig(multi_contour, "multi_conductivity_contour.png")

# ╔═╡ b80aedaf-b264-47b3-a679-1c73ce1e8a0f
log.(real.(multi_conductivity_data))'

# ╔═╡ c3281d7d-3d8e-4283-ad14-0fbe30759990
multi_plo = Plots.plot(Ω, real.(multi_conductivity_data)[:, 1:25:end] .* 2π, legend = true, label = hcat(["T = $i K" for i in T[1:25:end]]...), yaxis = :log)

# ╔═╡ 4bab65b0-a561-4355-8075-8cd3d8ba09f1
Plots.plot(Ω, real.(B_conductivity_data)[:, 1:25:end] .* 2π * 2.25, legend = true, label = hcat(["T = $i K" for i in T[1:25:end]]...), yaxis = :log, xlabel = "Ω (THz)", ylabel = "Hellwarth Mobility (Vs/kg)")

# ╔═╡ da3417bb-2570-41df-ac54-bdf1a46299d9
Plots.plot(T, real.(B_conductivity_data)[1:25:end, :]' .* 2.5 .* 2π, legend = true, label = hcat(["Ω = $i THz" for i in Ω[1:25:end]]...), xlabel = "T (K)", ylabel = "Hellwarth Mobility (cm^2/Vs)", yaxis = :log)

# ╔═╡ 67223f4b-4bc5-4f73-bb8c-80ebd9a16f39
Plots.plot(T, real.(multi_conductivity_data)[1:25:end, :]' .* 2π, legend = true, label = hcat(["Ω = $i THz" for i in Ω[1:25:end]]...), xlabel = "T (K)", ylabel = "Multi Mobility (cm^2/Vs)", yaxis = :log)

# ╔═╡ 72cf81af-42ac-4a4f-a2a2-6221548b047b
begin
	multi_plot_temp = Plots.plot(Ω, log.(real.(multi_conductivity_data[:, 1])), legend = false, color = "black", ylabel = "T (K)", xlabel = "Ω (THz)")
	for i in 1:3:length(T)
		plot!(Ω, log.(real.(multi_conductivity_data[:, i])) .- i/10, color = "black")
	end
	multi_plot_temp
end

# ╔═╡ d224e64d-2205-4159-a990-8bdce9d9c9be
savefig(multi_plot_temp, "multi_conductivity_temp_zeppelin.png")

# ╔═╡ 96eb31a1-d7bb-4d16-9b14-5c48845f5b54
begin
	multi_plot_freq = Plots.plot(T, log.(real.(multi_conductivity_data[1, :])), legend = false, color = "black", xlabel = "T (K)", ylabel = "Ω (THz)")
	for i in 1:1:length(Ω)
		plot!(T, log.(real.(multi_conductivity_data[i, :])) .- i/5, color = "black")
	end
	multi_plot_freq
end

# ╔═╡ 7d4c7313-3ba2-406c-8115-5937f66b2e0e
savefig(multi_plot_freq, "multi_conductivity_freq_zeppelin.png")

# ╔═╡ 4d1875ff-6f21-43e4-8c59-066450f84707
begin
	B_plot_temp = Plots.plot(Ω, log.(real.(B_conductivity_data[:, 1])), legend = false, color = "black", ylabel = "T (K)", xlabel = "Ω (THz)")
	for i in 1:3:length(T)
		plot!(Ω, log.(real.(B_conductivity_data[:, i])) .- i/10, color = "black")
	end
	B_plot_temp
end

# ╔═╡ becba837-6712-429e-b136-bdacd34ad5fc
savefig(B_plot_temp, "B_conductivity_temp_zeppelin.png")

# ╔═╡ a1dab1bc-672e-49a6-bb37-92b0d60ca896
begin
	B_plot_freq = Plots.plot(T, log.(real.(B_conductivity_data[1, :])), legend = false, color = "black", xlabel = "T (K)", ylabel = "Ω (THz)")
	for i in 1:1:length(Ω)
		plot!(T, log.(real.(B_conductivity_data[i, :])) .- i/5, color = "black")
	end
	B_plot_freq
end

# ╔═╡ 7ecb4d64-1b56-4dd1-ab08-3399b846a704
savefig(B_plot_freq, "B_conductivity_freq_zeppelin.png")

# ╔═╡ Cell order:
# ╠═ed9ab030-e61d-11eb-32f3-5be220868ae7
# ╠═d51f978f-a2c8-410a-8836-f57fc362eed1
# ╠═ffa81468-a1cd-44f3-aa5a-c970e56b7f27
# ╠═4c98c769-0aed-4b8a-bf92-fb06268100a9
# ╠═88fe20b4-303a-494e-95cc-e8183b2f474e
# ╠═a76c9858-0fe7-45c2-9b5c-a8cf7f800c08
# ╠═9d039885-65fa-4f43-ac15-ba927a68be60
# ╠═88c443f2-f6ea-464e-9190-aa158a0f7f50
# ╠═db784d36-b859-4022-b585-0730bf081a65
# ╠═329f1be8-967d-4c5b-8536-f6f1008f2a9c
# ╠═26ec12fd-5c2d-458c-8282-b411ef5184a8
# ╠═1e89541e-0635-47d9-96f6-99dceef27d00
# ╠═512e223b-249d-48ff-9eec-4282f4658c76
# ╠═b80aedaf-b264-47b3-a679-1c73ce1e8a0f
# ╠═c3281d7d-3d8e-4283-ad14-0fbe30759990
# ╠═4bab65b0-a561-4355-8075-8cd3d8ba09f1
# ╠═da3417bb-2570-41df-ac54-bdf1a46299d9
# ╠═67223f4b-4bc5-4f73-bb8c-80ebd9a16f39
# ╠═72cf81af-42ac-4a4f-a2a2-6221548b047b
# ╠═d224e64d-2205-4159-a990-8bdce9d9c9be
# ╠═96eb31a1-d7bb-4d16-9b14-5c48845f5b54
# ╠═7d4c7313-3ba2-406c-8115-5937f66b2e0e
# ╠═4d1875ff-6f21-43e4-8c59-066450f84707
# ╠═becba837-6712-429e-b136-bdacd34ad5fc
# ╠═a1dab1bc-672e-49a6-bb37-92b0d60ca896
# ╠═7ecb4d64-1b56-4dd1-ab08-3399b846a704
