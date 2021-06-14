### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 8336bf6e-0fa8-11eb-1f48-37cfe47a6d0f
using Pkg

# ╔═╡ 8ae04a5a-0fa8-11eb-2821-8dc46ee73269
Pkg.activate("../")

# ╔═╡ 66923a50-0fa8-11eb-1547-cf8e289818b0
using Revise

# ╔═╡ 93125894-0fa8-11eb-0ff6-211feff97889
using PolaronMobility

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
	
	Gnuplot.save("/home/jarvist/complex_conductivity_MAPI.gp")
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


# ╔═╡ ccab0b3e-c9f7-4081-b037-290eaf44292f


# ╔═╡ Cell order:
# ╠═d263806b-0851-4854-ba06-df1071b3d38b
# ╠═c6606f73-3c31-4df8-8585-81698966008a
# ╠═66923a50-0fa8-11eb-1547-cf8e289818b0
# ╠═8336bf6e-0fa8-11eb-1f48-37cfe47a6d0f
# ╠═8ae04a5a-0fa8-11eb-2821-8dc46ee73269
# ╠═93125894-0fa8-11eb-0ff6-211feff97889
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
# ╠═ccab0b3e-c9f7-4081-b037-290eaf44292f
