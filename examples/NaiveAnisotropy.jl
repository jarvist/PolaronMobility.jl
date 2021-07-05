### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ a6049db0-dd73-11eb-2757-05d672d68b4e
using Revise

# ╔═╡ 728b1147-6022-4787-9e75-0990b797a7a3
using PolaronMobility

# ╔═╡ 51cc807b-cf2c-459c-a1b7-85c00f73b145
begin
	ħ = 1.05457162825e-34
	eV = 1.602176487e-19
	me = 9.10938188e-31
	kB =  1.3806504e-23
	ϵ_0 = 8.854E-12
	amu = 1.660_539_066_60e-27
	Ha = 4.35974820e-18
	Bohr = 5.29177249e-11
end

# ╔═╡ b8654fc8-0b81-4ee2-866d-bc88679028e0
data = [
	"Material" "ϵ∞" "ϵ0" "m^*_⟂" "m^*_z" "ω_LO" "ZPR" "∂^2Σ/∂k_2⟂" "∂^2Σ/∂k_2z" "m^pol_⟂" "m^pol_z" "m^*_⟂/m^pol_⟂" "m^*_z/m^pol_z" "a_p" "a_p⟂" "a_pz";
	"AlAs"	9.49 	11.51 	0.24276 	0.89650 	47.3	-8.8 	-0.14802 	-0.02315 	0.25181 0.91550 	0.96407 	0.97925 	105.67 	123.26 	77.68;
	"AlP"	8.12 	10.32 	0.25190		0.80934		59.9 	-14.0 	-0.17791	-0.03398 	0.26372 0.83222		0.95519		0.97250		88.57	101.64	67.26;		
	"AlSb"	12.02 	13.35 	0.22156		1.14183		39.8	-3.6	-0.08089	-0.00782	0.22561	1.15212		0.98208		0.99107		157.37 	190.81 	107.04;		
	"BAs"	9.81 	9.89 	0.21642		1.09408		84.4	-0.5	-0.00557	-0.00055	0.21668	1.09474		0.99879		0.99939		1665.87	2015.38	1138.17;	
	"BN"	4.52 	6.69	0.29866		0.89524		161.0	-67.9	-0.26823	-0.05659	0.32467	0.94301		0.91989		0.94934		30.38	34.58	23.45;
	"CdS"	6.21 	10.24 	0.11773		0.11773		34.4	-14.9	-0.61225	-0.61225	0.12688	0.12688		0.92792		0.92792		503.93	0.0		0.0;		
	"CdSe"	7.83 	11.78 	0.05115		0.05115		23.6	-5.5	-0.75785	-0.75785	0.05322	0.05322		0.96124		0.96123		1716.44	0.0		0.0;	
	"CdTe"	8.89 	12.37 	0.05178		0.05178		19.1	-3.7	-0.61860	-0.61860	0.05349	0.05349		0.96797		0.96797		2294.84	0.0		0.0;		
	"GaAs"	15.31 	17.55 	0.00912		0.00912		33.5	-0.5	-0.29326	-0.29326	0.00914	0.00914		0.99733		0.99733		49463.62	0.0		0.0;	
	"GaN"	6.13 	11.00 	0.14355		0.14355		86.0	-29.6	-0.39962	-0.39965	0.15229	0.15228		0.94263		0.94263		362.68	362.68	362.69;		
	"GaP"	10.50 	12.53 	0.22958		1.06158		48.6	-7.4	-0.13113	-0.01483	0.23670	1.07856		0.96990		0.98426		95.73	114.62	66.78;		
	"SiC"	6.97 	10.30 	0.22814		0.67709		117.0	-32.7	-0.23217	-0.04969	0.24090	0.70067		0.94703		0.96636		62.48	71.04	48.33;		
	"ZnS"	5.97 	9.40 	0.16715		0.16715		40.6	-18.6	-0.45614	-0.45614	0.18094	0.18094		0.92376		0.92376		368.04	0.0		0.0;		
	"ZnSe"	7.35 	10.73 	0.08932		0.08932		29.3	-8.1	-0.51498	-0.51498	0.09362	0.09362		0.95400		0.95400		982.25	0.0		0.0;		
	"ZnTe"	9.05 	 11.99 	0.07644		0.07644		24.1	-4.3	-0.38801	-0.38801	0.07877	0.07877		0.97034		0.97034		1815.51	0.0		0.0;		
	"BaO"	4.21 	92.43 	0.38040		1.19717		47.3	-132.0	-1.40050	-0.27555	0.81412	1.78651		0.46725		0.67012		7.01	8.02	5.35;		
	"CaO"	3.77 	16.67 	0.44286		1.42415		66.8	-153.9	-0.99507	-0.18981	0.79178	1.95175		0.55932		0.72968		6.42	7.37	4.88;		
	"Li2O"	2.9 	7.8 	0.43735		0.84989		86.3	-171.7	-0.82404	-0.32275	0.68378	1.17113		0.63960		0.72570		13.49	14.59	11.53;		
	"MgO"	3.23 	11.14 	0.33954		0.33954		84.5	-137.3	-0.79790	-0.79790	0.46571	0.46571		0.72908		0.72908		50.37	0.0		0.0;		
	"SrO"	3.77 	20.91 	0.40701		1.22525		55.4	-141.0	-1.18792	-0.24910	0.78801	1.76349		0.51650		0.69479		7.31	8.33	5.64;		
]

# ╔═╡ 3d1f6989-12f1-44ee-8dcb-8a06cc83f6b1
ϵ_optic = data[2:end, 2]

# ╔═╡ 4bb07926-4ed8-4f36-9d4e-561f3693a801
ϵ_static = data[2:end, 3]

# ╔═╡ 3f7c4d03-3f7a-4d7a-af1a-fed238ea4ce7
phonon_freqs = data[2:end, 6] ./ 1e3 .* eV ./ ħ ./ 2π # Hz

# ╔═╡ 7af80fc4-e2a9-433c-bc98-f2497045bf3e
m_perp = data[2:end, 4] ./ Ha ./ Bohr^2 .* ħ^2 ./ me # me

# ╔═╡ 5cfd7351-ca1b-4757-8652-fc7147693101
m_z = data[2:end, 5] ./ Ha ./ Bohr^2 .* ħ^2 ./ me # me

# ╔═╡ 74598315-7ac3-4aa8-9e3c-f38bfe1e7205
α_perp = [frohlichalpha(ϵ∞, ϵ0, f, m) for (ϵ∞, ϵ0, f, m) in zip(ϵ_optic, ϵ_static, phonon_freqs, m_perp)]

# ╔═╡ 1b524ae3-60a4-47d6-8535-10e77abd2571
α_z = [frohlichalpha(ϵ∞, ϵ0, f, m) for (ϵ∞, ϵ0, f, m) in zip(ϵ_optic, ϵ_static, phonon_freqs, m_z)]

# ╔═╡ be0fd60e-633a-40e3-b07e-8bf98cd243f6
# Athermal Theory

# ╔═╡ 53103b1e-c993-480c-bcbc-004d6b56ac3a
begin
	var_perp_athermal = feynmanvw.(α_perp)
	v_perp_athermal = [i[1] for i in var_perp_athermal]
	w_perp_athermal = [i[2] for i in var_perp_athermal]
end

# ╔═╡ a932ff04-a4d5-4cab-baf5-79c6a763a132
v_perp_athermal

# ╔═╡ 0e89abc0-9332-42df-aec9-77f94f07c241
w_perp_athermal

# ╔═╡ add797ac-99c9-4ea8-8fda-bff4524fb4e0
begin
	var_z_athermal = feynmanvw.(α_z)
	v_z_athermal = [i[1] for i in var_z_athermal]
	w_z_athermal = [i[2] for i in var_z_athermal]
end

# ╔═╡ 64382d9f-2575-4521-9fb2-94bdef83833a
v_z_athermal

# ╔═╡ b5391e54-1b0f-4e73-9489-79c803f66acf
w_z_athermal

# ╔═╡ be65dc13-6c3b-4ecb-83d5-ed2f079b18c6
F_perp_athermal = F.(v_perp_athermal, w_perp_athermal, α_perp) .* 1e3 .* ħ .* 2π .* phonon_freqs ./ eV

# ╔═╡ b2e5b2d2-65bd-45e7-8f3f-ef82f343582b
F_z_athermal = F.(v_z_athermal, w_z_athermal, α_z) .* 1e3 .* ħ .* 2π .* phonon_freqs ./ eV

# ╔═╡ f8313eab-8792-4854-9288-f6f803e68e03
F_avg_athermal = (F_perp_athermal .* 2 .+ F_z_athermal) ./ 3

# ╔═╡ 72e71b1b-ca34-4f64-a561-a4c3e5fcb758
m_pol_perp_athermal = m_perp .* (1 .+ (v_perp_athermal.^2 .- w_perp_athermal.^2) ./ w_perp_athermal.^2) .* Ha .* Bohr^2 ./ ħ^2 .* me

# ╔═╡ 0bcdc70b-e564-4874-94ea-c10b7e638522
m_pol_z_athermal = m_z .* (1 .+ (v_z_athermal.^2 .- w_z_athermal.^2) ./ w_z_athermal.^2) .* Ha .* Bohr^2 ./ ħ^2 .* me

# ╔═╡ 73b1fe6d-771e-4eed-be9d-cc8eee02a7e1
m_perp_athermal_ratio = m_perp ./ m_pol_perp_athermal

# ╔═╡ 9a5479da-addc-4d74-a546-36105778798c
m_z_athermal_ratio = m_z ./ m_pol_z_athermal

# ╔═╡ fab8e16b-473f-4231-b165-87d027c40fb6
size_perp_athermal = sqrt.(3 ./ me ./ m_perp ./ (v_perp_athermal.^2 - w_perp_athermal.^2) .* v_perp_athermal ./ 2π ./ phonon_freqs .* ħ) ./ Bohr ./ 2

# ╔═╡ 41b02e23-0430-4e5d-bba8-b8de54258911
size_z_athermal = sqrt.(3 ./ me ./ m_z ./ (v_z_athermal.^2 - w_z_athermal.^2) .* v_z_athermal ./ 2π ./ phonon_freqs .* ħ) ./ Bohr ./2

# ╔═╡ a965eff9-e19e-46aa-bd20-f735429cd273
size_avg_athermal = (size_perp_athermal.^2 .* size_z_athermal).^(1/3)

# ╔═╡ 411691a5-53c7-40e5-8115-6d0ed0be07cb
# Thermal Theory

# ╔═╡ 66d7a406-7061-4714-9197-c65db6ef75cf
T = 300 # K

# ╔═╡ 50c5c581-2500-4746-8554-a51598c7aeea
β_red = ħ * 2π / kB / T .* phonon_freqs

# ╔═╡ 8fd12478-4320-4ffd-b686-45ba77b13e87
begin
	var_perp_thermal = feynmanvw.(α_perp, β_red)
	v_perp_thermal = [i[1] for i in var_perp_thermal]
	w_perp_thermal = [i[2] for i in var_perp_thermal]
end

# ╔═╡ 21ccb4da-ca25-4e50-9a3b-9eea67ff2699
begin
	var_z_thermal = feynmanvw.(α_z, β_red)
	v_z_thermal = [i[1] for i in var_z_thermal]
	w_z_thermal = [i[2] for i in var_z_thermal]
end

# ╔═╡ e9aa5e25-9467-427e-bd4b-7ad7d4906c22
F_perp_thermal = F.(v_perp_thermal, w_perp_thermal, β_red, α_perp) .* 1e3 .* ħ .* 2π .* phonon_freqs ./ eV

# ╔═╡ aa91544c-f776-4e2f-8fcf-e8c821bf5ac8
F_z_thermal = F.(v_z_thermal, w_z_thermal, β_red, α_z) .* 1e3 .* ħ .* 2π .* phonon_freqs ./ eV

# ╔═╡ 76d164e5-0e87-4aea-a40b-5e329371e848
F_avg_thermal = (F_perp_thermal .* 2 .+ F_z_thermal) / 3

# ╔═╡ ece97292-8325-4bf6-a18b-16f0a5be4621
m_pol_perp_thermal = m_perp ./ (1 .+ (v_perp_thermal.^2 .- w_perp_thermal.^2) ./ w_perp_thermal.^2 .* Ha .* Bohr^2 ./ ħ^2 .* me)

# ╔═╡ 67d60382-30d6-4b7f-9774-76a4d05e2098
m_pol_z_thermal = m_z ./ (1 .+ (v_z_thermal.^2 .- w_z_thermal.^2) ./ w_z_thermal.^2 .* Ha .* Bohr^2 ./ ħ^2 .* me)

# ╔═╡ 16d1552c-a7fe-4775-ac98-f97add75d1e4
m_perp_thermal_ratio = m_perp ./ m_pol_perp_thermal

# ╔═╡ f84c0d7f-13f7-49a0-9943-4f5c56527368
m_z_thermal_ratio = m_z ./ m_pol_z_thermal

# ╔═╡ 1b0d8352-7384-4355-ad13-2aae6aa6cf87
size_perp_thermal = sqrt.(3 ./ me ./ m_perp ./ (v_perp_thermal.^2 - w_perp_thermal.^2) .* v_perp_thermal ./ 2π ./ phonon_freqs .* ħ) ./ Bohr ./ 2

# ╔═╡ c274dd67-7e0c-417a-a972-26886fd32573
size_z_thermal = sqrt.(3 ./ me ./ m_z ./ (v_z_thermal.^2 - w_z_thermal.^2) .* v_z_thermal ./ 2π ./ phonon_freqs .* ħ) ./ Bohr ./ 2

# ╔═╡ 7e3aac26-62bf-410c-981f-654d918624bd
size_avg_thermal = (size_perp_thermal.^2 .* size_z_thermal).^(1/3)

# ╔═╡ Cell order:
# ╠═a6049db0-dd73-11eb-2757-05d672d68b4e
# ╠═728b1147-6022-4787-9e75-0990b797a7a3
# ╠═51cc807b-cf2c-459c-a1b7-85c00f73b145
# ╠═b8654fc8-0b81-4ee2-866d-bc88679028e0
# ╠═3d1f6989-12f1-44ee-8dcb-8a06cc83f6b1
# ╠═4bb07926-4ed8-4f36-9d4e-561f3693a801
# ╠═3f7c4d03-3f7a-4d7a-af1a-fed238ea4ce7
# ╠═7af80fc4-e2a9-433c-bc98-f2497045bf3e
# ╠═5cfd7351-ca1b-4757-8652-fc7147693101
# ╠═74598315-7ac3-4aa8-9e3c-f38bfe1e7205
# ╠═1b524ae3-60a4-47d6-8535-10e77abd2571
# ╠═be0fd60e-633a-40e3-b07e-8bf98cd243f6
# ╠═53103b1e-c993-480c-bcbc-004d6b56ac3a
# ╠═a932ff04-a4d5-4cab-baf5-79c6a763a132
# ╠═0e89abc0-9332-42df-aec9-77f94f07c241
# ╠═add797ac-99c9-4ea8-8fda-bff4524fb4e0
# ╠═64382d9f-2575-4521-9fb2-94bdef83833a
# ╠═b5391e54-1b0f-4e73-9489-79c803f66acf
# ╠═be65dc13-6c3b-4ecb-83d5-ed2f079b18c6
# ╠═b2e5b2d2-65bd-45e7-8f3f-ef82f343582b
# ╠═f8313eab-8792-4854-9288-f6f803e68e03
# ╠═72e71b1b-ca34-4f64-a561-a4c3e5fcb758
# ╠═0bcdc70b-e564-4874-94ea-c10b7e638522
# ╠═73b1fe6d-771e-4eed-be9d-cc8eee02a7e1
# ╠═9a5479da-addc-4d74-a546-36105778798c
# ╠═fab8e16b-473f-4231-b165-87d027c40fb6
# ╠═41b02e23-0430-4e5d-bba8-b8de54258911
# ╠═a965eff9-e19e-46aa-bd20-f735429cd273
# ╠═411691a5-53c7-40e5-8115-6d0ed0be07cb
# ╠═66d7a406-7061-4714-9197-c65db6ef75cf
# ╠═50c5c581-2500-4746-8554-a51598c7aeea
# ╠═8fd12478-4320-4ffd-b686-45ba77b13e87
# ╠═21ccb4da-ca25-4e50-9a3b-9eea67ff2699
# ╠═e9aa5e25-9467-427e-bd4b-7ad7d4906c22
# ╠═aa91544c-f776-4e2f-8fcf-e8c821bf5ac8
# ╠═76d164e5-0e87-4aea-a40b-5e329371e848
# ╠═ece97292-8325-4bf6-a18b-16f0a5be4621
# ╠═67d60382-30d6-4b7f-9774-76a4d05e2098
# ╠═16d1552c-a7fe-4775-ac98-f97add75d1e4
# ╠═f84c0d7f-13f7-49a0-9943-4f5c56527368
# ╠═1b0d8352-7384-4355-ad13-2aae6aa6cf87
# ╠═c274dd67-7e0c-417a-a972-26886fd32573
# ╠═7e3aac26-62bf-410c-981f-654d918624bd
