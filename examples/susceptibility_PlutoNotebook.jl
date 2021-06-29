### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 8336bf6e-0fa8-11eb-1f48-37cfe47a6d0f
using Pkg

# ╔═╡ 8ae04a5a-0fa8-11eb-2821-8dc46ee73269
Pkg.activate("../")

# ╔═╡ 66923a50-0fa8-11eb-1547-cf8e289818b0
using Revise

# ╔═╡ 3f9bdb02-7028-449b-93fc-01c4d3632db4
using PolaronMobility

# ╔═╡ 2024f235-aca8-427b-8579-2c869a0ba901
using QuadGK

# ╔═╡ f0529460-0fa8-11eb-1bec-1feb30645d42
using Plots

# ╔═╡ ffff1d5b-e9db-4da4-b01c-8c8075023997
using Gnuplot

# ╔═╡ d263806b-0851-4854-ba06-df1071b3d38b
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s

# ╔═╡ c6606f73-3c31-4df8-8585-81698966008a
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2

# ╔═╡ 8b357378-449b-45cb-aa8e-fa2a5bf5feed
T=10

# ╔═╡ c3fd3779-7c5c-4370-8378-66c4b0d963af
MAPIe=polaronmobility(T, 4.5, 24.1, 2.25E12, 0.12)

# ╔═╡ 386fd0a4-0faa-11eb-14e3-a3cdf3f3cb5d
#α=5
α=MAPIe.α[1]

# ╔═╡ a86c04ae-0faa-11eb-3565-7bd3108906b2
#v,w=feynmanvw(α) # Variational solution, athermal (original) action

# ╔═╡ 4306f2f4-0faa-11eb-1d32-dd8ade1d930a
f=2.25E12

# ╔═╡ b9a25625-38bd-43d7-a4b7-5a0da5a28772
ω=2π*f

# ╔═╡ 947530ce-e248-47af-ae49-34d698312d12
βred=ħ*ω/(kB*T)

# ╔═╡ b3e10abe-0faa-11eb-308e-91c4c0886b8b
v,w=feynmanvw(α, βred) # Variational solution, thermal (Osaka) action

# ╔═╡ 4b26ec8c-0faa-11eb-2af1-e3f5e39fb358
meff=0.12

# ╔═╡ 34a7c84e-0faa-11eb-08b6-83f4ab86b689
nurange=0.0:0.1:3.5 #I know this says nu, but due to the rescaling below, this is actually THz

# ╔═╡ ad10a50c-0fa8-11eb-30a9-1b1bb4fc26b7
oldX=PolaronMobility.ImX(nurange,v,w,βred,α,ω,meff*PolaronMobility.MassElectron)

# ╔═╡ 9ca217bf-985f-48fc-8a35-913c18fd79fc
plot(oldX.nu, oldX.ImX, label="ImX")

# ╔═╡ e45514ef-e24f-46ed-9844-b525e1042a80
βred1=ħ*2π*1E12/(kB*T)

# ╔═╡ 6ffe0fba-4146-43d4-90bd-dc325da74482
feynmanvw(α,βred1)

# ╔═╡ 2f09ce0a-8243-4e3f-9a3a-48ddf0a5b12a
βred2=ħ*2π*2.25E12/(kB*T)

# ╔═╡ 172d0553-d260-491a-9b94-79b865e3ee9a
feynmanvw(α,βred2)

# ╔═╡ 28040f04-9d08-449d-bc6c-079fae287657
X=[PolaronMobility.χ(nu, βred, α, v,w) for nu in nurange]

# ╔═╡ d721c0a3-234d-4824-9619-8a8159be9220
begin
	plot(nurange.*f/1E12, imag.(X), label="ImX")
	plot!(nurange.*f/1E12, real.(X), label="ReX")
end

# ╔═╡ 73ac4039-4675-4dce-8eab-b6635481406e
begin
	# re-scaling to frequency directly within the selection of nurange
	#  i.e. now nurange is in THz
	X1=[PolaronMobility.χ(nu, βred1, α,feynmanvw(α,βred1)...)/nu for nu in nurange./(1.0)]
	X2=[PolaronMobility.χ(nu, βred2, α,feynmanvw(α,βred2)...)/nu for nu in nurange./(2.25)]
end

# ╔═╡ 399bb77e-13a9-41dc-b432-a9bc605af25a
begin
	plot(nurange, imag.(X1), label="ImX 1")
	plot!(nurange, real.(X1), label="ReX 1")
	plot!(nurange, imag.(X2), label="ImX 2")
	plot!(nurange, real.(X2), label="ReX 2")
end

# ╔═╡ d33d2c9d-ac56-47ac-a0ad-b769c2fcf245
begin
	plot(nurange, imag.(X1+1.5*X2), label="ImX")
	plot!(nurange, real.(X1+1.5*X2), label="ReX")
	xlabel!("Frequeny (THz)")
	ylabel!("Χ(ν)/ν")
end

# ╔═╡ 3474df75-bd91-4803-acf3-bd6f697a6163
begin
	@gp nurange imag.(0.4*X1+0.6*X2) "w l t 'Im'"
	@gp :- nurange real.(0.4*X1+0.6*X2) "w l t 'Re'"
	@gp :- "set xlabel 'Frequency (THz)' "
	@gp :- "set ylabel 'Χ(ν)/ν' "
end

# ╔═╡ 04a89c80-7ecc-49fe-90f1-e1e081824d87
Gnuplot.save(term="pngcairo size 1024,768", output="~/tmp_plot.png")

# ╔═╡ 7e796fb8-1149-44fa-9076-9d7c2641be03
#Notable Mode f= 2.4380741812443247e12    α_partial=0.33509445526579384
#Notable Mode f= 2.249091763771941e12    α_partial=0.46413279655805495
#Notable Mode f= 2.0796321906344238e12    α_partial=0.5034016572517709
#Notable Mode f= 2.0336707697261187e12    α_partial=0.14188851068346853
#Notable Mode f= 1.5673011873879714e12    α_partial=0.16044621957165453
#Notable Mode f= 1.0188379384951798e12    α_partial=0.16189515169321733
#Notable Mode f= 9.970130778462072e11    α_partial=0.14036635422200763
#Notable Mode f= 9.20178190638621e11    α_partial=0.18115470952505333
#Notable Mode f= 5.738689505255511e11    α_partial=0.3479224374216951
#Sum alpha: 2.658367465100708
0.33509445526579384+0.46413279655805495+0.5034016572517709+0.14188851068346853

# ╔═╡ cafcff56-6e56-4b92-8ef0-87c5df309c7a
0.16044621957165453+0.16189515169321733+0.14036635422200763+0.18115470952505333+0.3479224374216951

# ╔═╡ 8530958f-7390-470a-9b61-ce88046bc084
1.5/2.5

# ╔═╡ 41398cba-bb88-4d4a-a817-1de8f551cc9c
function σ(χ)
	ϵr=χ+1
	n=√ϵr # assuming relative permeability of MAPI is 1...
	σ=im*(1-n)
	σ
end

# ╔═╡ c21ad7eb-9d9b-4202-b4ad-b1206dc94b3c
σ(1+0im)

# ╔═╡ 775c2be5-5aa4-41e2-9956-8e4131bff73b
σ(0+1im)

# ╔═╡ dab2fc71-d8d8-4e6d-9c43-d85688ca7558
σ1=σ.([PolaronMobility.χ(nu, βred1, α,feynmanvw(α,βred1)...)/nu for nu in nurange./(1.0)])


# ╔═╡ 9b8b2293-2395-4285-86fe-1b7fb7895297
σ2=σ.([PolaronMobility.χ(nu, βred2, α,feynmanvw(α,βred2)...)/nu for nu in nurange./(2.25)])

# ╔═╡ 53303534-308e-448f-9051-2bc3d859dc14
Gnuplot.options.mime

# ╔═╡ cdba35db-b258-433c-9387-42b505e2efb4
begin
	@gp nurange imag.(0.4*σ1+0.6*σ2) "w lp t 'Im'"
	@gp :- nurange real.(0.4*σ1+0.6*σ2) "w lp t 'Re'"
	@gp :- "set xlabel 'Frequency (THz)' "
	@gp :- "set ylabel 'σ' "
	
	Gnuplot.save("~/complex_conductivity_MAPI.gp")
	Gnuplot.save(term="pngcairo size 1024,768 fontscale 2 pointscale 2 linewidth 2", output="~/complex_conductivity_MAPI.png")

end

# ╔═╡ e0c76116-8c8a-4053-86c7-2f7c48f808bf
begin
	@gp nurange imag.(σ1) "w lp t '1THz phonons Im'"
	@gp :- nurange real.(σ1) "w lp t '1THz phonons Re'"
	@gp :- nurange imag.(σ2) "w lp t '2.25THz phonons Im'"
	@gp :- nurange real.(σ2) "w lp t '2.25THz phonons Re'"
	@gp :- nurange imag.(0.4*σ1+0.6*σ2) "w lp t 'Combined Im'"
	@gp :- nurange real.(0.4*σ1+0.6*σ2) "w lp t 'Combined Re'"

	@gp :- "set xlabel 'Frequency (THz)' "
	@gp :- "set ylabel 'σ' "
end

# ╔═╡ 1880d08f-f808-4690-8ad1-af1616cda0bf
# Multiple branch free energy optimisation

# ╔═╡ e39b7ce9-5396-4065-ba6a-eac83e626794
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
# 0.07575149723846182 2.78248540373041e-05
]

# ╔═╡ 3f2969ce-53c4-44ce-b442-c6c1e4c8e63c
v_params, w_params = PolaronMobility.multi_variation(0.01, 4.5, 0.12, (6.29e-10)^3, MAPI; N = 1)

# ╔═╡ 9f8e375a-e45d-430e-89f0-e62573f2f45b
F_j = [PolaronMobility.multi_free_energy(v_params[j, :], w_params[j, :], 10, 4.5, 0.12, (6.29e-10)^3, MAPI, j) for j in 1:length(MAPI[:, 1])] 

# ╔═╡ b0fc3a52-b4da-4ac4-a36e-e6229cd1f20e
F_total = sum(F_j)

# ╔═╡ d36a215e-9411-4332-82cf-caa87a56c7f4
# multiple branch susceptibility and conductivity

# ╔═╡ c2e007f4-dac0-4be5-ab55-10186fb34d61
phonon_mode_freqs = MAPI[:, 1]

# ╔═╡ ecd96809-9808-46fa-9405-884aeb8c7e6e
ir_activities = MAPI[:, 2]

# ╔═╡ f0d30402-ee5f-4cf9-b4bd-af35a102f71d
ϵ_ionic = [PolaronMobility.ϵ_ionic_mode(f, r, (6.29e-10)^3) for (f, r) in zip(phonon_mode_freqs, ir_activities)]

# ╔═╡ e2faa878-0ad1-42f7-99e8-53699c1c3468
ϵ_total = sum(ϵ_ionic)

# ╔═╡ 2d161b89-5594-4601-bb8e-c1748ade6269
α_j = [PolaronMobility.frohlich_α_j(4.5, ϵ_i, ϵ_total, f, 0.12) for (ϵ_i, f) in zip(ϵ_ionic, phonon_mode_freqs)]

# ╔═╡ 928e7f3b-6eab-4081-94bf-a4fa3a953f35
α_eff = sum(α_j)

# ╔═╡ 1dd3fcda-6818-4601-9fa4-7b75a07ee8d8
βr = [ħ * 2π * 1e12 * f / kB / 40 for f in phonon_mode_freqs]

# ╔═╡ 9b410a25-1c09-4bcb-b06d-63c928ab5bbc
Ω_range = 0.01:0.01:2.5

# ╔═╡ d20fd376-b5a6-48ce-9380-fccf420bdddd
Plots.PyPlotBackend()

# ╔═╡ c69faa18-759e-465c-b720-d8c6d429780d
χ = [PolaronMobility.multi_impedence(Ω, βr, α_j, v_params, w_params, phonon_mode_freqs, 0.12) for Ω in Ω_range]

# ╔═╡ e052454b-f885-4d76-8eb1-baedeff24e1e
begin
	χ_plot = plot(Ω_range, (imag.(χ)), xlabel = "Frequency (THz)", ylabel = "z", label = "Imz", minorgrid = true, markershape = :circle)
	plot!(χ_plot, Ω_range, (real.(χ)), label = "Rez", markershape = :circle)
end

# ╔═╡ dbe241ba-5b5c-45d0-9d63-3ccefc82766f
σ_new = PolaronMobility.multi_conductivity(Ω_range, χ)

# ╔═╡ ec793073-968b-447b-a445-9182f1cfdf0f
begin
	σ_plot = plot(Ω_range, (real.(σ_new)), xlabel = "Frequency (THz)", ylabel = "σ", label = "Reσ", legend = :bottomright, minorgrid = true, markershape = :circle)
	plot!(σ_plot, Ω_range, (imag.(σ_new)), label = "Imσ", markershape = :circle)
end

# ╔═╡ 36f09479-cf29-4a0f-9f97-81169abbcd9b
σ_abs_plot = plot(Ω_range, (sqrt.(abs.(σ_new))), xlabel = "Frequency (THz)", ylabel = "abs(σ)", legend = false, minorgrid = true, markershape = :circle)

# ╔═╡ Cell order:
# ╠═d263806b-0851-4854-ba06-df1071b3d38b
# ╠═c6606f73-3c31-4df8-8585-81698966008a
# ╠═66923a50-0fa8-11eb-1547-cf8e289818b0
# ╠═8336bf6e-0fa8-11eb-1f48-37cfe47a6d0f
# ╠═8ae04a5a-0fa8-11eb-2821-8dc46ee73269
# ╠═3f9bdb02-7028-449b-93fc-01c4d3632db4
# ╠═2024f235-aca8-427b-8579-2c869a0ba901
# ╠═8b357378-449b-45cb-aa8e-fa2a5bf5feed
# ╠═c3fd3779-7c5c-4370-8378-66c4b0d963af
# ╠═386fd0a4-0faa-11eb-14e3-a3cdf3f3cb5d
# ╠═a86c04ae-0faa-11eb-3565-7bd3108906b2
# ╠═947530ce-e248-47af-ae49-34d698312d12
# ╠═b3e10abe-0faa-11eb-308e-91c4c0886b8b
# ╠═4306f2f4-0faa-11eb-1d32-dd8ade1d930a
# ╠═b9a25625-38bd-43d7-a4b7-5a0da5a28772
# ╠═4b26ec8c-0faa-11eb-2af1-e3f5e39fb358
# ╠═34a7c84e-0faa-11eb-08b6-83f4ab86b689
# ╠═ad10a50c-0fa8-11eb-30a9-1b1bb4fc26b7
# ╠═f0529460-0fa8-11eb-1bec-1feb30645d42
# ╠═9ca217bf-985f-48fc-8a35-913c18fd79fc
# ╠═e45514ef-e24f-46ed-9844-b525e1042a80
# ╠═6ffe0fba-4146-43d4-90bd-dc325da74482
# ╠═2f09ce0a-8243-4e3f-9a3a-48ddf0a5b12a
# ╠═172d0553-d260-491a-9b94-79b865e3ee9a
# ╠═28040f04-9d08-449d-bc6c-079fae287657
# ╠═d721c0a3-234d-4824-9619-8a8159be9220
# ╠═73ac4039-4675-4dce-8eab-b6635481406e
# ╠═399bb77e-13a9-41dc-b432-a9bc605af25a
# ╠═d33d2c9d-ac56-47ac-a0ad-b769c2fcf245
# ╠═ffff1d5b-e9db-4da4-b01c-8c8075023997
# ╠═3474df75-bd91-4803-acf3-bd6f697a6163
# ╠═04a89c80-7ecc-49fe-90f1-e1e081824d87
# ╠═7e796fb8-1149-44fa-9076-9d7c2641be03
# ╠═cafcff56-6e56-4b92-8ef0-87c5df309c7a
# ╠═8530958f-7390-470a-9b61-ce88046bc084
# ╠═41398cba-bb88-4d4a-a817-1de8f551cc9c
# ╠═c21ad7eb-9d9b-4202-b4ad-b1206dc94b3c
# ╠═775c2be5-5aa4-41e2-9956-8e4131bff73b
# ╠═dab2fc71-d8d8-4e6d-9c43-d85688ca7558
# ╠═9b8b2293-2395-4285-86fe-1b7fb7895297
# ╠═53303534-308e-448f-9051-2bc3d859dc14
# ╠═cdba35db-b258-433c-9387-42b505e2efb4
# ╠═e0c76116-8c8a-4053-86c7-2f7c48f808bf
# ╠═1880d08f-f808-4690-8ad1-af1616cda0bf
# ╠═e39b7ce9-5396-4065-ba6a-eac83e626794
# ╠═3f2969ce-53c4-44ce-b442-c6c1e4c8e63c
# ╠═9f8e375a-e45d-430e-89f0-e62573f2f45b
# ╠═b0fc3a52-b4da-4ac4-a36e-e6229cd1f20e
# ╠═d36a215e-9411-4332-82cf-caa87a56c7f4
# ╠═c2e007f4-dac0-4be5-ab55-10186fb34d61
# ╠═ecd96809-9808-46fa-9405-884aeb8c7e6e
# ╠═f0d30402-ee5f-4cf9-b4bd-af35a102f71d
# ╠═e2faa878-0ad1-42f7-99e8-53699c1c3468
# ╠═2d161b89-5594-4601-bb8e-c1748ade6269
# ╠═928e7f3b-6eab-4081-94bf-a4fa3a953f35
# ╠═1dd3fcda-6818-4601-9fa4-7b75a07ee8d8
# ╠═9b410a25-1c09-4bcb-b06d-63c928ab5bbc
# ╠═d20fd376-b5a6-48ce-9380-fccf420bdddd
# ╠═c69faa18-759e-465c-b720-d8c6d429780d
# ╠═e052454b-f885-4d76-8eb1-baedeff24e1e
# ╠═dbe241ba-5b5c-45d0-9d63-3ccefc82766f
# ╠═ec793073-968b-447b-a445-9182f1cfdf0f
# ╠═36f09479-cf29-4a0f-9f97-81169abbcd9b
