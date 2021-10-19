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

# ╔═╡ c819bc8f-cc85-4edb-8289-5884f7fbb8f7
using LaTeXStrings

# ╔═╡ e2d229a0-eaae-481c-9ee1-94d641fc3014
using ColorSchemes

# ╔═╡ 88fe20b4-303a-494e-95cc-e8183b2f474e
Ω = real.(parse.(ComplexF64, (CSV.File("multi_conductivity_data.csv") |> Tables.matrix)[2:end, 1]))

# ╔═╡ a76c9858-0fe7-45c2-9b5c-a8cf7f800c08
T = real.(parse.(ComplexF64, (CSV.File("multi_conductivity_data.csv") |> Tables.matrix)[1, 2:end]))

# ╔═╡ 6d08ec24-51d5-4105-bddf-443bc8496456
A_freq = (CSV.File("B_data.csv") |> Tables.matrix)[:, 2]

# ╔═╡ 5fa3101f-d24b-41f1-9a77-f5b7782ec4a3
B_freq = (CSV.File("B_data.csv") |> Tables.matrix)[:, 2]

# ╔═╡ 78156817-9c43-44c0-90e9-55e1100538a7
multi_freqs = real.((CSV.File("multi_data.csv") |> Tables.matrix)[1:end, 2])

# ╔═╡ 4bf7a882-1306-415a-8760-6191e5284306
multi_α = real.((CSV.File("multi_data.csv") |> Tables.matrix)[1:end, 1])

# ╔═╡ 8b91053e-28cc-47ce-8551-72e3e3c8ee93
A_α = real.((CSV.File("A_data.csv") |> Tables.matrix)[1:end, 1])

# ╔═╡ fd24b1c2-a6c2-4328-854b-ada51d40c56a
B_α = real.((CSV.File("B_data.csv") |> Tables.matrix)[1:end, 1])

# ╔═╡ db273a6b-bd14-4935-92dd-1271def823ba
ir_activity = real.((CSV.File("multi_data.csv") |> Tables.matrix)[1:end, end])

# ╔═╡ 830ce1e4-80f4-4dfa-a214-d7d0c84cfbbb
begin
	scatter(multi_freqs, multi_α, yaxis = :log, xticks = (multi_freqs, hcat(["$(round(i, digits = 2))" for i in multi_freqs]...)), label = L"\textbf{\alpha_j}", size = (500, 500), xlabel = "Phonon Frequencies (THz)", tickfontsize = 9, xrotation = 90)
	scatter!(multi_freqs, ir_activity, markershape = :diamond, label = "IR")
	plot!(multi_freqs, multi_α ./ ir_activity, linewidth = 2, linestyle = :dash, label = "ratio")
end

# ╔═╡ 4b5dc310-b4f6-4bd3-88e5-f082ea8ee47e
F_multi = real.((CSV.File("multi_vw.csv") |> Tables.matrix)[1:end, 5])

# ╔═╡ f9e47e6f-fdc5-42b6-9ea2-46c801f8cf6b
F_multi_avg = [F_multi[i] for i in 1:length(T)]

# ╔═╡ e2a502b1-cca9-40c2-a107-2b577f0662dc
F_A = real.((CSV.File("A_data.csv") |> Tables.matrix)[1:end, end])

# ╔═╡ 92f361d6-2052-45bd-a2a1-0d19d04eec1b
F_B = real.((CSV.File("B_data.csv") |> Tables.matrix)[1:end, end])

# ╔═╡ b1473f00-5aa9-4ace-a0f9-9e5cb9747af5
multi_energy_temp = plot(T, F_multi, label = hcat(["$(round(i, digits = 3)) THz" for i in multi_freqs]...), legend = :outerright, linewidth = 2, size = (650, 500), minorgrid = true, ylabel = "Free energy (meV)", xlabel = "T (K)", tickfontsize = 9, legendfontsize = 9)

# ╔═╡ 4c4d1b61-6710-478c-a8f5-067c219be342
savefig(multi_energy_temp, "multi_energy_temp.pdf")

# ╔═╡ 9493c1f8-d2b4-42b3-9af5-2d69ee7499bf
begin
	free_energy_temp = plot(T, F_A, linewidth = 2, linestyle = :dashdot, label = "H-A", tickfontsize = 12, labelfontsize = 12, legendfontsize = 12, xlabel = "T(K)", ylabel = "Free energy (meV)", minorgrid = true, size = (550, 550), legend = :topright)
	plot!(T, F_B, linewidth = 2, linestyle = :dash, label = "H-B")
	plot!(T, F_multi_avg, linewidth = 2, linestyle = :solid, label = "multi")
end

# ╔═╡ d4a13a12-e1ca-4ded-9c49-3c33b1cb69a0
savefig(free_energy_temp, "free_energy_temp.pdf")

# ╔═╡ 40fc4b41-c835-41a3-b468-765d4bbd0042
begin
	v_A = real.((CSV.File("A_data.csv") |> Tables.matrix)[1:end, end-2])
	w_A = real.((CSV.File("A_data.csv") |> Tables.matrix)[1:end, end-1])
	v_B = real.((CSV.File("B_data.csv") |> Tables.matrix)[1:end, end-2])
	w_B = real.((CSV.File("B_data.csv") |> Tables.matrix)[1:end, end-1])
end

# ╔═╡ aa1b7720-91ee-48e5-a97e-f8e3a499ad2f
v_multi = real.((CSV.File("multi_vw.csv") |> Tables.matrix)[1:end, 3])

# ╔═╡ 10b31119-9d2f-445d-9353-9e6788a464d3
w_multi = real.((CSV.File("multi_vw.csv") |> Tables.matrix)[1:end, 4])

# ╔═╡ eef00978-0153-4031-b748-d9e1e5e2f5c5
multi_v_temp = plot(T, v_multi .* (multi_freqs)', label = hcat(["$(round(i, digits = 3)) THz" for i in multi_freqs]...), legend = :outerright, linewidth = 1.5, size = (650, 500), minorgrid = true, ylabel = "v (THz)", xlabel = "T (K)", tickfontsize = 9, legendfontsize = 9)

# ╔═╡ 45e2101a-010d-404d-aaf6-c68d28536677
multi_w_temp = plot(T, w_multi .* multi_freqs', label = hcat(["$(round(i, digits = 3)) THz" for i in multi_freqs]...), legend = :outerright, linewidth = 1.5, size = (650, 500), minorgrid = true, ylabel = "w (THz)", xlabel = "T (K)", tickfontsize = 9, legendfontsize = 9)

# ╔═╡ 4da848c3-23dd-4640-bab3-a8b4b4784639
plot(T, (v_multi .- w_multi) .* multi_freqs', label = hcat(["$(round(i, digits = 3)) THz" for i in multi_freqs]...), legend = :outerright, linewidth = 1.5, size = (650, 500), minorgrid = true, ylabel = "v - w (THz)", xlabel = "T (K)", tickfontsize = 9, legendfontsize = 9)

# ╔═╡ 8f23b6a9-9cb8-47b5-bec1-b7f4f2fe050b
begin
	savefig(multi_v_temp, "multi_v_temp.pdf")
	savefig(multi_w_temp, "multi_w_temp.pdf")
end

# ╔═╡ e7631e4b-dc2c-4b0c-b7a9-7a92dbe0fd1e
begin
	vw_temp = plot(T, v_A, label = "H-A v", legend = :bottomright, linewidth = 2, tickfontsize = 12, labelfontsize = 12, legendfontsize = 12, xlabel = "T(K)", ylabel = "v & w (THz)", size = (550, 550), minorgrid = true, linestyle = :dashdot, color = theme_palette(:default)[1])
	plot!(T, w_A, label = "H-A w", linewidth = 2, linestyle = :dashdot, color = theme_palette(:default)[2])
	plot!(T, v_B, label = "H-B v", linewidth = 2, linestyle = :dash, color = theme_palette(:default)[3])
	plot!(T, w_B, label = "H-B w", linewidth = 2, linestyle = :dash, color = theme_palette(:default)[4])
	plot!(T, v_multi, label = "multi v", linewidth = 2, linestyle = :solid, color = theme_palette(:default)[1])
	plot!(T, w_multi, label = "multi w", linewidth = 2, linestyle = :solid, color = theme_palette(:default)[2])
end

# ╔═╡ 32c47e89-a988-4586-840b-0ac952e1d916
savefig(vw_temp, "vw_temp.pdf")

# ╔═╡ dd0f2762-4cbb-4fc4-b79e-a70b1006905e
begin
	plot(T, v_A, label = "H-A v", legend = :bottomright, linewidth = 2, tickfontsize = 12, labelfontsize = 12, legendfontsize = 12, xlabel = "T(K)", ylabel = "v & w (THz)", size = (550, 550), minorgrid = true, linestyle = :dashdot, color = theme_palette(:default)[1])
	plot!(T, w_A, label = "H-A w", linewidth = 2, linestyle = :dashdot, color = theme_palette(:default)[2])
	plot!(T, v_B, label = "H-B v", linewidth = 2, linestyle = :dash, color = theme_palette(:default)[3])
	plot!(T, w_B, label = "H-B w", linewidth = 2, linestyle = :dash, color = theme_palette(:default)[4])
end

# ╔═╡ 9d039885-65fa-4f43-ac15-ba927a68be60
multi_conductivity_data = parse.(ComplexF64, (CSV.File("multi_conductivity_data.csv") |> Tables.matrix)[2:end, 2:end])

# ╔═╡ 88c443f2-f6ea-464e-9190-aa158a0f7f50
A_conductivity_data = parse.(ComplexF64, (CSV.File("A_conductivity_data.csv") |> Tables.matrix)[2:end, 2:end])

# ╔═╡ db784d36-b859-4022-b585-0730bf081a65
B_conductivity_data = parse.(ComplexF64, (CSV.File("B_conductivity_data.csv") |> Tables.matrix)[2:end, 2:end])

# ╔═╡ 656db2a0-91d6-447b-95a8-2a2981979645
A_contour_real = Plots.contourf(Ω, T, real.(A_conductivity_data)', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 500))

# ╔═╡ a586ecee-fd54-44e7-8dea-fa1f2a6829df
A_contour_imag = Plots.contourf(Ω, T, imag.(A_conductivity_data)', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), clims = (-6, 500), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm)

# ╔═╡ 2aeb1ede-b786-494f-9f15-0a8382943c29
A_contour_abs = Plots.contourf(Ω, T, abs.(A_conductivity_data)', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true),  clims = (0, 500), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm)

# ╔═╡ 88c2fd6b-a8cc-436a-b9e9-5412c19b42bd
begin
	savefig(A_contour_real, "A_contour_real.pdf")
	savefig(A_contour_imag, "A_contour_imag.pdf")
	savefig(A_contour_abs, "A_contour_abs.pdf")
end

# ╔═╡ 329f1be8-967d-4c5b-8536-f6f1008f2a9c
B_contour_real = Plots.contourf(Ω, T, real.(B_conductivity_data)', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 500))

# ╔═╡ ba30cda1-9a2e-4427-9cde-3a5e3e7e4dcf
B_contour_imag = Plots.contourf(Ω, T, imag.(B_conductivity_data)', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (-6, 500))

# ╔═╡ e30b9ee1-0eac-44ad-bb6e-4c243672468e
B_contour_abs = Plots.contourf(Ω, T, abs.(B_conductivity_data)', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 500))

# ╔═╡ 26ec12fd-5c2d-458c-8282-b411ef5184a8
begin
	savefig(B_contour_real, "B_contour_real.pdf")
	savefig(B_contour_imag, "B_contour_imag.pdf")
	savefig(B_contour_abs, "B_contour_abs.pdf")
end

# ╔═╡ 1e89541e-0635-47d9-96f6-99dceef27d00
multi_contour_real = Plots.contourf(Ω, T, real.(multi_conductivity_data)', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 500))

# ╔═╡ 6b617a87-f84e-4818-be31-5806d2e8894d
multi_contour_imag = Plots.contourf(Ω, T, imag.(multi_conductivity_data)', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (-6, 500))

# ╔═╡ bc190237-ddff-4abc-bf29-0e23d10fad0d
multi_contour_abs = Plots.contourf(Ω, T, abs.(multi_conductivity_data)', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true),  linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 12, labelfontsize = 12, right_margin = 7Plots.mm, clims = (0, 500))

# ╔═╡ 512e223b-249d-48ff-9eec-4282f4658c76
begin
	savefig(multi_contour_real, "multi_contour_real.pdf")
	savefig(multi_contour_imag, "multi_contour_imag.pdf")
	savefig(multi_contour_abs, "multi_contour_abs.pdf")
end

# ╔═╡ 0ec007f8-9b46-4d0c-a8d2-c12ff85b0f5a
AB_contour_real = Plots.contourf(Ω, T, log.(abs.(real.(A_conductivity_data .- B_conductivity_data)))', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true, scale = :log), clims = (-16, 5), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 10, right_margin = 7Plots.mm)

# ╔═╡ 7a68a5ca-2c8f-40bc-ac7b-8a996f48b41b
AB_contour_imag = Plots.contourf(Ω, T, log.(abs.(imag.(A_conductivity_data .- B_conductivity_data)))', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true, scale = :log), clims = (-16, 5), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 10, right_margin = 7Plots.mm)

# ╔═╡ 3a14ec9b-c726-40da-a4ed-ad22ee89e8d3
AB_contour_abs = Plots.contourf(Ω, T, log.(abs.(A_conductivity_data .- B_conductivity_data))', xlabel = "Ω (THz)", ylabel = "T (K)", fill = cgrad(:thermal, rev = false, categorical = true, scale = :log), clims = (-16, 5), linewidth = 1.3, color = :grey10, size = (550, 500), tickfontsize = 10, right_margin = 7Plots.mm)

# ╔═╡ 6071f2a4-ec70-4a12-b40a-84062ae2bb6a
begin
	savefig(AB_contour_real, "AB_contour_real.pdf")
	savefig(AB_contour_imag, "AB_contour_imag.pdf")
	savefig(AB_contour_abs, "AB_contour_abs.pdf")
end

# ╔═╡ b80aedaf-b264-47b3-a679-1c73ce1e8a0f
real.(multi_conductivity_data)[:, 291]

# ╔═╡ c30386e2-e275-4e8c-bc46-e043ce3484b2
real.(A_conductivity_data)[:, 291]

# ╔═╡ 585d49c7-0b54-4884-922c-84e9ae387837
real.(B_conductivity_data)[:, 291]

# ╔═╡ c3281d7d-3d8e-4283-ad14-0fbe30759990
begin
	multi_plot = Plots.plot(Ω[1:250], real.(multi_conductivity_data[1:250, 1:1:30]), legend = false, label = hcat(["T = $i K" for i in T[1:1:30]]...), minorgrid = true, linewidth = 1.5, size = (550, 500), ylabel = L"\mathrm{Multi\ Mobility\ } (cm^2/Vs)", xlabel = L"\nu\ (THz)", tickfontsize = 9, legendfontsize = 9, ylims = (0, 1500), palette = palette(:thermal, 30))
end

# ╔═╡ aec2c420-5d16-4130-914b-0ba94f4cd9ca
red = hcat([(real.(multi_conductivity_data[i, 1]) .- real.(multi_conductivity_data[i-1, 10])) for i in 2:350]...)

# ╔═╡ 32773157-09b8-4feb-ae44-7f41e244b101
imd = hcat([(imag.(multi_conductivity_data[i, 10]) .- imag.(multi_conductivity_data[i-1, 10])) for i in 2:350]...)

# ╔═╡ 41b5f3a6-1fb7-4cab-b712-57816e8209c5
begin
	plot(Ω[1:250], red[1:250], size = (550, 500), minorgrid = true)
	plot!(Ω[1:250], -imd[1:250])
end

# ╔═╡ a64cfd0a-67f6-4941-8ba4-2e26b3361040
begin
	multi_plot2 = Plots.plot(Ω[1:250], imag.(multi_conductivity_data[1:250, 1:1:30]), legend = false, label = hcat(["T = $i K" for i in T[1:1:30]]...), minorgrid = true, linewidth = 1.5, size = (550, 500), ylabel = L"\mathrm{Multi\ Mobility\ } (cm^2/Vs)", xlabel = L"\nu\ (THz)", tickfontsize = 9, legendfontsize = 9, ylims = (0, 1500), palette = palette(:thermal, 30))
end

# ╔═╡ 835feade-6750-48ef-baab-be9dce52500e
Plots.plot(Ω ./ 2π, abs.(real.(A_conductivity_data))[:, 400], legend = true, label = hcat(["T = $i K" for i in T[1]]...), yaxis = :log, xlabel = L"\nu\ (THz)", ylabel = L"\mathrm{Hellwarth\ A\ Mobility\ } (cm^2/Vs)", linewidth = 1.5, minorgrid = true, ylims = (5, 10^2.8), size = (500, 500), tickfontsize = 9, legendfontsize = 9)

# ╔═╡ 4bab65b0-a561-4355-8075-8cd3d8ba09f1
Plots.plot(Ω, abs.(real.(B_conductivity_data))[:, 1:60:end], legend = true, label = hcat(["T = $i K" for i in T[1:60:end]]...), yaxis = :log, xlabel = L"\Omega\ (THz)", ylabel = L"\mathrm{Hellwarth\ B\ Mobility\ } (cm^2/Vs)", linewidth = 1.5, minorgrid = true, ylims = (5, 10^2.8), size = (500, 500), tickfontsize = 9, legendfontsize = 9)

# ╔═╡ f4f970e4-6113-4782-b6dc-9072254c5b4e
Plots.plot(T, real.(multi_conductivity_data)[1:25:end, :]', legend = :topright, label = hcat(["Ω = $i THz" for i in Ω[1:25:end]]...), xlabel = L"\textrm{T\ (K)}", ylabel = L"\textbf{Multi\ Mobility\ (cm^2/Vs)}", yaxis = :log, minorgrid = true, linewidth = 1.5, ylims = (5, 10^4.5), size = (500, 500), tickfontsize = 9, legendfontsize = 9, xaxis = :log, xlims = (10, 400))

# ╔═╡ 86ff3a78-852f-424a-af6d-32e886726d0a
plot!(T, T.^(-0.2) * 18^2, linestyle = :dash, linewidth = 3, label  = "T^(-0.2)")

# ╔═╡ efca960d-a7d6-42df-8196-7e00ae0ac024
Plots.plot(T, real.(A_conductivity_data)[1:25:end, :]', legend = :topright, label = hcat(["Ω = $i THz" for i in Ω[1:25:end]]...), xlabel = L"\textrm{T\ (K)}", ylabel = L"\textrm{Hellwarth\ A\ Mobility\ (cm^2/Vs)}", yaxis = :log, minorgrid = true, linewidth = 1.5, ylims = (5, 10^5), size = (500, 500), tickfontsize = 9, legendfontsize = 9)

# ╔═╡ da3417bb-2570-41df-ac54-bdf1a46299d9
Plots.plot(T, real.(B_conductivity_data)[1:25:end, :]', legend = true, label = hcat(["Ω = $i THz" for i in Ω[1:25:end]]...), xlabel = L"\textrm{T\ (K)}", ylabel = L"\textrm{Hellwarth\ B\ Mobility\ (cm^2/Vs)}", yaxis = :log, minorgrid = true, linewidth = 1.5, ylims = (5, 10^5), size = (500, 500), tickfontsize = 9, legendfontsize = 9, xaxis = :log, xlims = (10, 400))

# ╔═╡ 5c685b4e-8702-4b74-90e8-25ffd53e5582
plot!(T, T.^(-0.5) * 50^2, linestyle = :dash, linewidth = 3, label  = "T^(-0.2)", legend = :outerright, size = (700, 500))

# ╔═╡ 9bd4dd63-c784-4a1e-86fa-935bd8c25ad2
function ridgeline(x, y, z; shift = 1, jump = 10, ylabel = :none, xlabel = :none, palette = :twilight, linewidth = 1.5, color = :black, size = (595, 842), ymirror = false, fillalpha = 1, step = 1)
	if !ymirror
		p = plot(x, z[:, 1], fillrange = z[:, jump] .- shift, fillalpha = fillalpha, ylabel = ylabel, xlabel = xlabel, xticks = (0:maximum(x)/10:maximum(x), string.(x[1:Int((length(x)-1)/10):end])), yticks = (z[1, 1:jump * step:end] .- shift .* range(0, stop = length(y[1:jump:end]) - 1, step = step), string.(y[1:jump * step:end])), legend = false, size = size, grid = false, ymirror = ymirror, palette = palette, tickfontsize = 12, labelfontsize = 12)
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

# ╔═╡ 7ae0d8a1-d2ad-4599-bce1-f71865fa7497
multi_plot_temp_real = ridgeline(round.(Ω, digits = 1), Int.(T), log.(real.(multi_conductivity_data)), jump = 10, shift = 0.6, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ 23c44da4-8468-4987-8d5e-1069ffb41d13
multi_plot_temp_imag = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(imag.(multi_conductivity_data))), jump = 10, shift = 0.6, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 7)

# ╔═╡ 077a2c40-9aae-4045-89da-99f1db15fae7
multi_plot_temp_abs = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(multi_conductivity_data)), jump = 10, shift = 0.6, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ fbf4e5d9-c20b-4953-8477-c0144e34e9e9
begin
	savefig(multi_plot_temp_real, "multi_plot_temp_real.pdf")
	savefig(multi_plot_temp_imag, "multi_plot_temp_imag.pdf")
	savefig(multi_plot_temp_abs, "multi_plot_temp_abs.pdf")
end

# ╔═╡ bfd794a2-abf9-499e-92fe-aa02f262171d
A_plot_temp_real = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(real.(A_conductivity_data))), jump = 10, shift = 0.6, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ cbdee07e-77e5-4113-b099-4feb30d5d74a
A_plot_temp_imag = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(imag.(A_conductivity_data))), jump = 10, shift = 0.6, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 7)

# ╔═╡ 2b4876d2-6d1e-4774-89c1-3ff8be380e11
A_plot_temp_abs = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(A_conductivity_data)), jump = 10, shift = 0.6, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ 79ddc218-4710-466a-a5ba-976431cd7165
begin
	savefig(A_plot_temp_real, "A_plot_temp_real.pdf")
	savefig(A_plot_temp_imag, "A_plot_temp_imag.pdf")
	savefig(A_plot_temp_abs, "A_plot_temp_abs.pdf")
end

# ╔═╡ 9b8183fb-3a65-463b-bae5-62a8aa9f06ac
B_plot_temp_real = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(real.(B_conductivity_data))), jump = 10, shift = 0.6, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ ed8f7ce5-c4eb-44da-a747-56ac7c152760
B_plot_temp_imag = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(imag.(B_conductivity_data))), jump = 10, shift = 0.6, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 7)

# ╔═╡ 80cde1da-7d6a-4fea-a2c2-89b2a0a6cd0d
B_plot_temp_abs = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(B_conductivity_data)), jump = 10, shift = 0.6, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ 4cfb19f7-6852-4224-9b09-54a844624cde
begin
	savefig(B_plot_temp_real, "B_plot_temp_real.pdf")
	savefig(B_plot_temp_imag, "B_plot_temp_imag.pdf")
	savefig(B_plot_temp_abs, "B_plot_temp_abs.pdf")
end

# ╔═╡ ccfa612b-9cf0-4e4c-a69c-e3b002e5fea5
AB_plot_temp_real = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(real.(A_conductivity_data .- B_conductivity_data))), jump = 5, shift = 0.3, xlabel = "Ω (THz)", ylabel = "T (K)", size = (550, 550), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 5)

# ╔═╡ 39eeac1c-5db8-4d6c-8552-d308b5199d5e
AB_plot_temp_imag = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(imag.(A_conductivity_data .- B_conductivity_data))), jump = 5, shift = 0.3, xlabel = "Ω (THz)", ylabel = "T (K)", size = (500, 500), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 4)

# ╔═╡ a8ab6e8a-2d89-4383-9ab6-e49990712489
AB_plot_temp_abs = ridgeline(round.(Ω, digits = 1), Int.(T), log.(abs.(A_conductivity_data .- B_conductivity_data)), jump = 5, shift = 0.3, xlabel = "Ω (THz)", ylabel = "T (K)", size = (500, 500), linewidth = 0.5, palette = :thermal, fillalpha = 0.75, step = 4)

# ╔═╡ ba68e7aa-963f-4cdb-a141-b4ee9ff697db
begin
	savefig(AB_plot_temp_real, "AB_plot_temp_real.pdf")
	savefig(AB_plot_temp_imag, "AB_plot_temp_imag.pdf")
	savefig(AB_plot_temp_abs, "AB_plot_temp_abs.pdf")
end

# ╔═╡ 9e0e8eea-f5e9-4d91-80cd-a1e7847933e6
multi_plot_freq_real = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(real.(multi_conductivity_data)))', jump = 6, shift = 0.4, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (550, 550), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 7)

# ╔═╡ 3e185071-6fda-4dc6-8e84-4900a106b798
multi_plot_freq_imag = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(imag.(multi_conductivity_data)))', jump = 6, shift = 0.4, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (550, 550), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 10)

# ╔═╡ 05516d71-0f45-4e97-97b0-8a1819aed270
multi_plot_freq_abs = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(multi_conductivity_data))', jump = 6, shift = 0.4, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (500, 550), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 7)

# ╔═╡ c3f82311-03fa-4cb0-9dca-4cf495b5e8ca
begin
	savefig(multi_plot_freq_real, "multi_plot_freq_real.pdf")
	savefig(multi_plot_freq_imag, "multi_plot_freq_imag.pdf")
	savefig(multi_plot_freq_abs, "multi_plot_freq_abs.pdf")
end

# ╔═╡ c3bdb1f6-889b-4c84-9857-88f619d6e617
A_plot_freq_real = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(real.(A_conductivity_data)))', jump = 6, shift = 0.6, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (550, 550), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 7)

# ╔═╡ 428b7174-9d3f-4762-8941-6feabcfbb4e0
A_plot_freq_imag = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(imag.(A_conductivity_data)))', jump = 6, shift = 0.6, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (550, 550), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 10)

# ╔═╡ 5900719d-e804-4a7d-a142-38fddf56a085
A_plot_freq_abs = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(A_conductivity_data))', jump = 6, shift = 0.6, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (550, 550), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 7)

# ╔═╡ ff9a0b04-b007-468e-acba-677a7f25a340
begin
	savefig(A_plot_freq_real, "A_plot_freq_real.pdf")
	savefig(A_plot_freq_imag, "A_plot_freq_imag.pdf")
	savefig(A_plot_freq_abs, "A_plot_freq_abs.pdf")
end

# ╔═╡ bab7f306-94d2-4360-8c0d-08b8a619046d
B_plot_freq_real = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(real.(B_conductivity_data)))', jump = 6, shift = 0.6, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (550, 550), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 7)

# ╔═╡ 0879342a-c552-4bc6-8dce-21f020abeb25
B_plot_freq_imag = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(imag.(B_conductivity_data)))', jump = 6, shift = 0.6, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (550, 550), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 10)

# ╔═╡ 8acd98e2-23ef-478b-97c2-533e0e2bef98
B_plot_freq_abs = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(B_conductivity_data))', jump = 6, shift = 0.6, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (550, 550), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 7)

# ╔═╡ 34411f71-63bc-4c91-9cd4-486e69f5876a
begin
	savefig(B_plot_freq_real, "B_plot_freq_real.pdf")
	savefig(B_plot_freq_imag, "B_plot_freq_imag.pdf")
	savefig(B_plot_freq_abs, "B_plot_freq_abs.pdf")
end

# ╔═╡ aff41df3-144f-452e-b059-c306000ad10d
AB_plot_freq_real = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(real.(A_conductivity_data .- B_conductivity_data)))', jump = 3, shift = 1, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (500, 500), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 5)

# ╔═╡ 9288099f-4e2a-43d6-ba2c-1debb17c4901
AB_plot_freq_imag = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(imag.(A_conductivity_data .- B_conductivity_data)))', jump = 3, shift = 1, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (500, 500), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 5)

# ╔═╡ a844da51-68ad-442f-9517-bbcde197c38e
AB_plot_freq_abs = ridgeline(Int.(T), round.(Ω, digits = 1), log.(abs.(A_conductivity_data .- B_conductivity_data))', jump = 3, shift = 1, ylabel = "Ω (THz)", xlabel = "T (K)", ymirror = true, size = (500, 500), linewidth = 0.75, palette = :twilight, fillalpha = 0.75, step = 5)

# ╔═╡ f47f71d3-870a-47b8-9ab2-426b8ceacf72
begin
	savefig(AB_plot_freq_real, "AB_plot_freq_real.pdf")
	savefig(AB_plot_freq_imag, "AB_plot_freq_imag.pdf")
	savefig(AB_plot_freq_abs, "AB_plot_freq_abs.pdf")
end

# ╔═╡ 8b7c9fd6-2289-47c6-ad73-2e60e867979e
plot(multi_freqs, F_multi[:, 1:1:20], markershape = :circle, legend = :bottomright, label = hcat(["$(T[i]) K" for i in 1:1:20]...), xticks = (multi_freqs, hcat(["$(round(i, digits = 2))" for i in multi_freqs]...)), xrotation = 90, tickfontsize = 9, ylabel = "Free energy (meV)", xlabel = "Phonon frequencies (THz)", linestyle = :dot)

# ╔═╡ 509b518e-07d9-4b5d-8938-b4a9ec188d8f
ridgeline(round.(multi_freqs, digits = 1), T[4:400], F_multi[:, 4:400], ymirror = true)

# ╔═╡ 18edac61-b8e1-4119-b49d-56d2b9ba2251
begin
	plot(multi_freqs, v_multi[:, 1:50:400] .* multi_freqs, markershape = :circle, legend = :outerright, label = hcat(["$(T[i]) K" for i in 1:50:400]...), xticks = (multi_freqs, hcat(["$(round(i, digits = 2))" for i in multi_freqs]...)), ylims = (0, 65), xrotation = 90, tickfontsize = 9, ylabel = "v & w (THz)", xlabel = "Phonon frequencies (THz)", linestyle = :dot)
	plot!(multi_freqs, w_multi[:, 1:50:400] .* multi_freqs, markershape = :diamond, legend = :outerright, label = :none, xticks = (multi_freqs, hcat(["$(round(i, digits = 2))" for i in multi_freqs]...)), ylims = (0, 65), xrotation = 90, tickfontsize = 9, linestyle = :dash, color = hcat([theme_palette(:default)[i] for i in 1:8]...))
end

# ╔═╡ 6dae57f3-cfd4-4598-bc1b-ea25f8d39565
plot(multi_freqs, (v_multi[:, 1:50:400] .- w_multi[:, 1:50:400]) .* multi_freqs ./ ir_activity, markershape = :circle, legend = :outerright, label = hcat(["$(T[i]) K" for i in 1:50:400]...), xticks = (multi_freqs, hcat(["$(round(i, digits = 2))" for i in multi_freqs]...)), xrotation = 90, tickfontsize = 9, ylabel = "(v - w) / IR", xlabel = "Phonon frequencies (THz)", linestyle = :dot)

# ╔═╡ Cell order:
# ╠═ed9ab030-e61d-11eb-32f3-5be220868ae7
# ╠═d51f978f-a2c8-410a-8836-f57fc362eed1
# ╠═ffa81468-a1cd-44f3-aa5a-c970e56b7f27
# ╠═4c98c769-0aed-4b8a-bf92-fb06268100a9
# ╠═c819bc8f-cc85-4edb-8289-5884f7fbb8f7
# ╠═88fe20b4-303a-494e-95cc-e8183b2f474e
# ╠═a76c9858-0fe7-45c2-9b5c-a8cf7f800c08
# ╠═6d08ec24-51d5-4105-bddf-443bc8496456
# ╠═5fa3101f-d24b-41f1-9a77-f5b7782ec4a3
# ╠═78156817-9c43-44c0-90e9-55e1100538a7
# ╠═4bf7a882-1306-415a-8760-6191e5284306
# ╠═8b91053e-28cc-47ce-8551-72e3e3c8ee93
# ╠═fd24b1c2-a6c2-4328-854b-ada51d40c56a
# ╠═db273a6b-bd14-4935-92dd-1271def823ba
# ╠═830ce1e4-80f4-4dfa-a214-d7d0c84cfbbb
# ╠═4b5dc310-b4f6-4bd3-88e5-f082ea8ee47e
# ╠═f9e47e6f-fdc5-42b6-9ea2-46c801f8cf6b
# ╠═e2a502b1-cca9-40c2-a107-2b577f0662dc
# ╠═92f361d6-2052-45bd-a2a1-0d19d04eec1b
# ╠═b1473f00-5aa9-4ace-a0f9-9e5cb9747af5
# ╠═4c4d1b61-6710-478c-a8f5-067c219be342
# ╠═9493c1f8-d2b4-42b3-9af5-2d69ee7499bf
# ╠═d4a13a12-e1ca-4ded-9c49-3c33b1cb69a0
# ╠═40fc4b41-c835-41a3-b468-765d4bbd0042
# ╠═aa1b7720-91ee-48e5-a97e-f8e3a499ad2f
# ╠═10b31119-9d2f-445d-9353-9e6788a464d3
# ╠═eef00978-0153-4031-b748-d9e1e5e2f5c5
# ╠═45e2101a-010d-404d-aaf6-c68d28536677
# ╠═4da848c3-23dd-4640-bab3-a8b4b4784639
# ╠═8f23b6a9-9cb8-47b5-bec1-b7f4f2fe050b
# ╠═e7631e4b-dc2c-4b0c-b7a9-7a92dbe0fd1e
# ╠═32c47e89-a988-4586-840b-0ac952e1d916
# ╠═dd0f2762-4cbb-4fc4-b79e-a70b1006905e
# ╠═9d039885-65fa-4f43-ac15-ba927a68be60
# ╠═88c443f2-f6ea-464e-9190-aa158a0f7f50
# ╠═db784d36-b859-4022-b585-0730bf081a65
# ╠═656db2a0-91d6-447b-95a8-2a2981979645
# ╠═a586ecee-fd54-44e7-8dea-fa1f2a6829df
# ╠═2aeb1ede-b786-494f-9f15-0a8382943c29
# ╠═88c2fd6b-a8cc-436a-b9e9-5412c19b42bd
# ╠═329f1be8-967d-4c5b-8536-f6f1008f2a9c
# ╠═ba30cda1-9a2e-4427-9cde-3a5e3e7e4dcf
# ╠═e30b9ee1-0eac-44ad-bb6e-4c243672468e
# ╠═26ec12fd-5c2d-458c-8282-b411ef5184a8
# ╠═1e89541e-0635-47d9-96f6-99dceef27d00
# ╠═6b617a87-f84e-4818-be31-5806d2e8894d
# ╠═bc190237-ddff-4abc-bf29-0e23d10fad0d
# ╠═512e223b-249d-48ff-9eec-4282f4658c76
# ╠═0ec007f8-9b46-4d0c-a8d2-c12ff85b0f5a
# ╠═7a68a5ca-2c8f-40bc-ac7b-8a996f48b41b
# ╠═3a14ec9b-c726-40da-a4ed-ad22ee89e8d3
# ╠═6071f2a4-ec70-4a12-b40a-84062ae2bb6a
# ╠═b80aedaf-b264-47b3-a679-1c73ce1e8a0f
# ╠═c30386e2-e275-4e8c-bc46-e043ce3484b2
# ╠═585d49c7-0b54-4884-922c-84e9ae387837
# ╠═c3281d7d-3d8e-4283-ad14-0fbe30759990
# ╠═aec2c420-5d16-4130-914b-0ba94f4cd9ca
# ╠═32773157-09b8-4feb-ae44-7f41e244b101
# ╠═41b5f3a6-1fb7-4cab-b712-57816e8209c5
# ╠═a64cfd0a-67f6-4941-8ba4-2e26b3361040
# ╠═835feade-6750-48ef-baab-be9dce52500e
# ╠═4bab65b0-a561-4355-8075-8cd3d8ba09f1
# ╠═f4f970e4-6113-4782-b6dc-9072254c5b4e
# ╠═86ff3a78-852f-424a-af6d-32e886726d0a
# ╠═efca960d-a7d6-42df-8196-7e00ae0ac024
# ╠═da3417bb-2570-41df-ac54-bdf1a46299d9
# ╠═5c685b4e-8702-4b74-90e8-25ffd53e5582
# ╠═e2d229a0-eaae-481c-9ee1-94d641fc3014
# ╠═9bd4dd63-c784-4a1e-86fa-935bd8c25ad2
# ╠═7ae0d8a1-d2ad-4599-bce1-f71865fa7497
# ╠═23c44da4-8468-4987-8d5e-1069ffb41d13
# ╠═077a2c40-9aae-4045-89da-99f1db15fae7
# ╠═fbf4e5d9-c20b-4953-8477-c0144e34e9e9
# ╠═bfd794a2-abf9-499e-92fe-aa02f262171d
# ╠═cbdee07e-77e5-4113-b099-4feb30d5d74a
# ╠═2b4876d2-6d1e-4774-89c1-3ff8be380e11
# ╠═79ddc218-4710-466a-a5ba-976431cd7165
# ╠═9b8183fb-3a65-463b-bae5-62a8aa9f06ac
# ╠═ed8f7ce5-c4eb-44da-a747-56ac7c152760
# ╠═80cde1da-7d6a-4fea-a2c2-89b2a0a6cd0d
# ╠═4cfb19f7-6852-4224-9b09-54a844624cde
# ╠═ccfa612b-9cf0-4e4c-a69c-e3b002e5fea5
# ╠═39eeac1c-5db8-4d6c-8552-d308b5199d5e
# ╠═a8ab6e8a-2d89-4383-9ab6-e49990712489
# ╠═ba68e7aa-963f-4cdb-a141-b4ee9ff697db
# ╠═9e0e8eea-f5e9-4d91-80cd-a1e7847933e6
# ╠═3e185071-6fda-4dc6-8e84-4900a106b798
# ╠═05516d71-0f45-4e97-97b0-8a1819aed270
# ╠═c3f82311-03fa-4cb0-9dca-4cf495b5e8ca
# ╠═c3bdb1f6-889b-4c84-9857-88f619d6e617
# ╠═428b7174-9d3f-4762-8941-6feabcfbb4e0
# ╠═5900719d-e804-4a7d-a142-38fddf56a085
# ╠═ff9a0b04-b007-468e-acba-677a7f25a340
# ╠═bab7f306-94d2-4360-8c0d-08b8a619046d
# ╠═0879342a-c552-4bc6-8dce-21f020abeb25
# ╠═8acd98e2-23ef-478b-97c2-533e0e2bef98
# ╠═34411f71-63bc-4c91-9cd4-486e69f5876a
# ╠═aff41df3-144f-452e-b059-c306000ad10d
# ╠═9288099f-4e2a-43d6-ba2c-1debb17c4901
# ╠═a844da51-68ad-442f-9517-bbcde197c38e
# ╠═f47f71d3-870a-47b8-9ab2-426b8ceacf72
# ╠═8b7c9fd6-2289-47c6-ad73-2e60e867979e
# ╠═509b518e-07d9-4b5d-8938-b4a9ec188d8f
# ╠═18edac61-b8e1-4119-b49d-56d2b9ba2251
# ╠═6dae57f3-cfd4-4598-bc1b-ea25f8d39565
