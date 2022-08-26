### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ d09dec61-ab12-4424-bfaf-30f155c90ef9
using Revise

# ╔═╡ 89fcb3a0-db48-11eb-30ac-b501d7fcd55f
using PolaronMobility

# ╔═╡ 4a7ed3cb-88e4-4e86-90c1-285849994ab9
using Plots

# ╔═╡ 8c54eedc-7c8f-4c5d-879e-62743e302801
using DataFrames

# ╔═╡ 3edfd42b-4b8a-4080-b617-740ef3fbf76b
using CSV

# ╔═╡ 27af9453-8044-4b1f-85fd-9e2c5dd104eb
using QuadGK

# ╔═╡ 97aa8c27-f2fe-418c-a662-906dbc126fb6
begin
    ħ = 1.05457162825e-34
    eV = 1.602176487e-19
    me = 9.10938188e-31
    kB = 1.3806504e-23
    ϵ_0 = 8.854E-12
    amu = 1.660_539_066_60e-27
    Ha = 4.35974820e-18
    Bohr = 5.29177249e-11
end

# ╔═╡ d08e66f9-636b-4377-b653-389ee7b78a24
T_range = 100:100:400

# ╔═╡ 50c1a877-807e-498e-87f9-0e909080a452
MAPI = [
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

# ╔═╡ e349bf0a-554f-46ae-8e88-8a7b80d0c4cd
# Multiple Phonon Theory

# ╔═╡ 4cd925a7-048a-4696-b2b9-6a6b394e0528
B_scheme = PolaronMobility.HellwarthBScheme(MAPI)

# ╔═╡ b04ccf74-853d-4ab2-886f-5b6c21cc9a51
phonon_mode_freqs = MAPI[:, 1]

# ╔═╡ 4e349b51-0590-4f98-b777-f13db8050b6d
ir_activities = MAPI[:, 2]

# ╔═╡ 1c6fff5f-839a-43e1-acfe-040df2ed4822
ϵ_ionic = [PolaronMobility.ϵ_ionic_mode(f, r, (6.29e-10)^3) for (f, r) in zip(phonon_mode_freqs, ir_activities)]

# ╔═╡ 2f524fa1-c08b-4c11-b05a-e2d016ce6ec0
ϵ_total = sum(ϵ_ionic)

# ╔═╡ 8510da63-37d7-45d0-9012-8dd17745a9e8
α_j = [PolaronMobility.frohlich_α_j(4.5, ϵ_i, ϵ_total, f, 0.12) for (ϵ_i, f) in zip(ϵ_ionic, phonon_mode_freqs)]

# ╔═╡ d9fca06f-ee9c-4790-a059-98667dcd8428
begin
    multi_data = DataFrames.DataFrame(
        alpha=α_j,
        phonon_freqs=phonon_mode_freqs,
        ir_activities=ir_activities,
        ionic=ϵ_ionic
    )
    CSV.write("multi_data.csv", multi_data)
end

# ╔═╡ 04e16581-ca4f-4456-bea2-2e8bf44a08f8
α_eff = sum(α_j)

# ╔═╡ ba626ebc-d7da-4eca-8149-907d724a2ef5
var_params = [PolaronMobility.multi_variation(i, 4.5, 0.12, (6.29e-10)^3, MAPI; N=1) for i in T_range]

# ╔═╡ af26a6c8-9c97-4130-bd54-c9fe40e58879
[var_params[i][1] for i in 1:length(T_range)]

# ╔═╡ dd41fbf9-b4b9-46b9-abbf-96eddf59b0ad
v_j = [var_params[i][1][j] for j in 1:length(MAPI[:, 1]), i in 1:length(T_range)]

# ╔═╡ 2975a6b8-bcec-4e4e-be94-5a3897c650db
begin
    multi_v = DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [phonon_mode_freqs v_j]], :auto)
    CSV.write("multi_v.csv", multi_v)
end

# ╔═╡ 748d6916-6125-4f9a-adbf-cee3b25dc3ba
w_j = [var_params[i][2][j] for j in 1:length(MAPI[:, 1]), i in 1:length(T_range)]

# ╔═╡ a304f008-b3fa-4a92-a8c1-9ce2f7b4b3df
begin
    multi_w = DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [phonon_mode_freqs w_j]], :auto)
    CSV.write("multi_w.csv", multi_w)
end

# ╔═╡ 5652b186-0e8d-469c-bfea-bf665f75c06a
F_j = [PolaronMobility.multi_free_energy(var_params[i][1][j, :], var_params[i][2][j, :], T_range[i], 4.5, 0.12, (6.29e-10)^3, MAPI, j) for j in 1:length(MAPI[:, 1]), i in 1:length(T_range)]

# ╔═╡ 2e0eb8f9-3b6b-4464-b715-258c96e22422
begin
    multi_energy = DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [phonon_mode_freqs F_j]], :auto)
    CSV.write("multi_energy.csv", multi_energy)
end

# ╔═╡ 143b1c45-ce4b-4f3c-9674-ff181c64fe90
F_total = [sum(F_j[:, i]) for i in 1:length(T_range)]

# ╔═╡ 3900999b-fdfa-48ae-afd4-bca5d03008e3
Ω_range = 0.01:0.1:20.01

# ╔═╡ 30d37afa-7b70-479f-8dcc-860af723a96e
β_red = [ħ / kB / T * f * 2π * 1e12 for f in MAPI[:, 1], T in T_range]

# ╔═╡ 9c140cb6-286b-4e25-a3dc-829b1eea9dfb
begin
    multi_beta = DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [phonon_mode_freqs β_red]], :auto)
    CSV.write("multi_beta.csv", multi_beta)
end

# ╔═╡ 65a76ecd-1dcd-4d5c-88be-1b57d31369da
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

# ╔═╡ 6a793986-844a-4a82-84b0-b71a48aa8242
function D_j(τ, β, v, w) # log of dynamic structure factor for polaron 
    D = τ * (1 - τ / β)
    for i in 1:length(v)
        if v[i] != w[i]
            D += (h_i(i, v, w) / v[i]^2) * (2 * sinh(v[i] * τ / 2) * sinh(v[i] * (β - τ) / 2) / (v[i] * sinh(v[i] * β / 2)) - τ * (1 - τ / β))
        end
    end
    return D
end

# ╔═╡ 146bcc43-15fc-45de-bb9b-f16a68d84e41
function multi_memory_function(Ω, β, α, v, w, f, m_eff)

    m_e = 9.10938188e-31 # Electron Mass (kg)
    eV = 1.602176487e-19 # Electron Volt (kg m^2 s^{-2})

    # FHIP1962, page 1009, eqn (36).
    S(t, β_j, v_j, w_j) = cos(t - 1im * β_j / 2) / sinh(β_j / 2) / D_j(-1im * t, β_j, v_j, w_j)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a). Scale Frequency Ω by phonon branch frequency f_j.
    integrand(t, β_j, v_j, w_j, Ω) = (1 - exp(-1im * Ω * t)) * imag(S(t, β_j, v_j, w_j))

    impedence = 0.0
    for j in 1:length(f) # sum over phonon branches
        impedence += 1im * (Ω - 2 * α[j] * (f[j])^2 * quadgk(t -> integrand(t, β[j], v[j], w[j], Ω * 2π / f[j]), 0.0, Inf)[1] / (3 * √π * Ω))
    end
    return impedence / eV * m_eff * m_e * 1e12 / 100^2
end

# ╔═╡ 3c7564b0-b9fe-49e8-9ad7-eaa80e9841a7
begin
    Z = Array{ComplexF64}(undef, length(Ω_range), length(T_range))
    for j in 1:length(T_range), i in 1:length(Ω_range)
        println("Ω: $(Ω_range[i]), T: $(T_range[j])")
        Z[i, j] = multi_memory_function(Ω_range[i], β_red[:, j], α_j, var_params[j][1], var_params[j][2], [i for i in phonon_mode_freqs], 0.12)
        CSV.write("multi_conductivity_data.csv", DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [Ω_range 1 ./ Z]], :auto))
    end
end

# ╔═╡ 0c424ddb-aecd-4fc6-9f26-54bb1fc9f099
Z

# ╔═╡ fff9cb8d-00aa-426c-bb66-4b7b21871c36
Plots.plot(Ω_range, real.(1 ./ Z), minorgrid=true, label=hcat(["T = $i K" for i in T_range]...), yaxis=:log)

# ╔═╡ 8ace5afe-29da-4bd1-9599-b51c03b67b6e
Plots.plot(Ω_range, imag.(1 ./ Z), minorgrid=true, label=hcat(["T = $i K" for i in T_range]...))

# ╔═╡ fab3ad29-61f6-45dd-b63f-4ea40c960afb
Plots.plot(Ω_range, abs.(1 ./ Z), minorgrid=true, label=hcat(["T = $i K" for i in T_range]...))

# ╔═╡ 3fb0f337-c6dd-4f56-a694-00d1c2f403bd
# Hellwarth Theory

# ╔═╡ 53117290-63b5-469c-93a9-bab8a40b2658
# B Scheme

# ╔═╡ a6b7c896-04d8-4704-8173-e8895ff84cbd
α_hellwarth_B = PolaronMobility.frohlichalpha(4.5, 24.1, B_scheme * 1e12, 0.12)

# ╔═╡ 6d322516-d247-4bd0-bcd5-9b8371b47750
β_red_hellwarth_B = [ħ * 2π / kB / T * B_scheme * 1e12 for T in T_range]

# ╔═╡ ae756cc9-f4ca-4073-abff-c56ecd1fdade
var_B = PolaronMobility.feynmanvw.(α_hellwarth_B, β_red_hellwarth_B)

# ╔═╡ 4c90b122-1bb9-460e-ac68-608ca8074a50
v_B = [i[1] for i in var_B]

# ╔═╡ 83068836-36d9-42aa-9867-7555d89e0d42
w_B = [i[2] for i in var_B]

# ╔═╡ 8082688f-ae0f-4aad-9bb9-71852be0fe00
Hellwarth_energy_B = F.(v_B, w_B, β_red_hellwarth_B, α_hellwarth_B) .* 1e3 .* ħ .* 2π .* B_scheme .* 1e12 ./ eV

# ╔═╡ 466354b5-7d93-4e2b-ab6d-f7a01112fbfb
begin
    B_data = DataFrames.DataFrame(
        alpha=[α_hellwarth_B for i in 1:length(T_range)],
        efffreq=[B_scheme for i in 1:length(T_range)],
        temp=T_range,
        beta=β_red_hellwarth_B,
        v=v_B,
        w=w_B,
        F=Hellwarth_energy_B
    )
    CSV.write("B_data.csv", B_data)
end

# ╔═╡ b1cf4fe2-8dda-4882-b39f-f23adb2738d1
begin
    ZB = Array{ComplexF64}(undef, length(Ω_range), length(T_range))
    for j in 1:length(T_range), i in 1:length(Ω_range)
        println("Ω: $(Ω_range[i]), T: $(T_range[j])")
        ZB[i, j] = multi_memory_function(Ω_range[i], β_red_hellwarth_B[j], [α_hellwarth_B], v_B[j], w_B[j], [B_scheme / 2π], 0.12) * 2π
        CSV.write("B_conductivity_data.csv", DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [Ω_range 1 ./ ZB]], :auto))
    end
end

# ╔═╡ a15e4d30-9ed2-492a-bfa4-0608180738c1
ZB

# ╔═╡ 2015d16c-e9f8-4072-834b-a317223b854f
Plots.plot(Ω_range, real.(1 ./ ZB), minorgrid=true, label=hcat(["T = $i K" for i in T_range]...), yaxis=:log)

# ╔═╡ 313a5d22-9a70-4c20-acb2-402c1683fbf8
Plots.plot(Ω_range, imag.(1 ./ ZB), minorgrid=true, label=hcat(["T = $i K" for i in T_range]...))

# ╔═╡ b8061ceb-1610-481b-986b-c5b82b9678db
Plots.plot(Ω_range, abs.(1 ./ ZB), minorgrid=true, label=hcat(["T = $i K" for i in T_range]...))

# ╔═╡ c7cdee43-9ac6-490b-9466-367a1f136766
# A Scheme

# ╔═╡ bb33a5c2-555b-4078-b278-90e5ae864fed
A_scheme = reverse([PolaronMobility.HellwarthAScheme(MAPI, T=i) for i in T_range])

# ╔═╡ 44604f2b-e5fd-438d-87df-0f10f874d00d
α_hellwarth_A = PolaronMobility.frohlichalpha.(4.5, 24.1, A_scheme .* 1e12, 0.12)

# ╔═╡ e858791b-e869-490f-9af3-2826f480d722
β_red_hellwarth_A = [ħ * 2π / kB / T_range[i] * A_scheme[i] * 1e12 for i in 1:length(T_range)]

# ╔═╡ 1f47a9e8-c2de-4963-ad58-aca71dbbb2ac
var_A = PolaronMobility.feynmanvw.(α_hellwarth_A, β_red_hellwarth_A)

# ╔═╡ 105665d5-ba53-4ddd-86df-dc40d5fefa8f
v_A = [i[1] for i in var_A]

# ╔═╡ 4468bf47-7668-4309-8f67-59fee2a1cf43
w_A = [i[2] for i in var_A]

# ╔═╡ 4bae3659-666b-4ae7-924a-b20d98a43148
Hellwarth_energy_A = F.(v_A, w_A, β_red_hellwarth_A, α_hellwarth_A) .* 1e3 .* ħ .* 2π .* A_scheme .* 1e12 ./ eV

# ╔═╡ 6e03cbc9-c6cd-440b-8652-fa1be4d7934a
begin
    A_data = DataFrames.DataFrame(
        alpha=α_hellwarth_A,
        efffreq=A_scheme,
        temp=T_range,
        beta=β_red_hellwarth_A,
        v=v_A,
        w=w_A,
        F=Hellwarth_energy_A
    )
    CSV.write("A_data.csv", A_data)
end

# ╔═╡ 53ac6ee6-f5d1-4923-9352-9d1cd96817fc
begin
    ZA = Array{ComplexF64}(undef, length(Ω_range), length(T_range))
    for j in 1:length(T_range), i in 1:length(Ω_range)
        println("Ω: $(Ω_range[i]), T: $(T_range[j])")
        ZA[i, j] = multi_memory_function(Ω_range[i], β_red_hellwarth_A[j], α_hellwarth_A[j], v_A[j], w_A[j], A_scheme[j] / 2π, 0.12) * 2π
        CSV.write("A_conductivity_data.csv", DataFrames.DataFrame([[0.0, [i for i in T_range]...]'; [Ω_range 1 ./ ZA]], :auto))
    end
end

# ╔═╡ 9ec88954-f8a6-4ae0-a800-3723de537d10
ZA

# ╔═╡ e129c967-858e-46cb-b6ce-29c238f45277
Plots.plot(Ω_range, real.(1 ./ ZA), label=hcat(["T = $i K" for i in T_range]...), minorgrid=true)

# ╔═╡ c8cdaa0e-05cb-4b9a-99dc-6932312ec56e
Plots.plot(Ω_range, imag.(1 ./ ZA), label=hcat(["T = $i K" for i in T_range]...), minorgrid=true)

# ╔═╡ 0f1066e0-d13a-4e1d-b9c0-c51509f25eb6
Plots.plot(Ω_range, abs.(1 ./ ZA), label=hcat(["T = $i K" for i in T_range]...), minorgrid=true)

# ╔═╡ Cell order:
# ╠═d09dec61-ab12-4424-bfaf-30f155c90ef9
# ╠═89fcb3a0-db48-11eb-30ac-b501d7fcd55f
# ╠═4a7ed3cb-88e4-4e86-90c1-285849994ab9
# ╠═8c54eedc-7c8f-4c5d-879e-62743e302801
# ╠═3edfd42b-4b8a-4080-b617-740ef3fbf76b
# ╠═97aa8c27-f2fe-418c-a662-906dbc126fb6
# ╠═d08e66f9-636b-4377-b653-389ee7b78a24
# ╠═50c1a877-807e-498e-87f9-0e909080a452
# ╠═e349bf0a-554f-46ae-8e88-8a7b80d0c4cd
# ╠═4cd925a7-048a-4696-b2b9-6a6b394e0528
# ╠═b04ccf74-853d-4ab2-886f-5b6c21cc9a51
# ╠═4e349b51-0590-4f98-b777-f13db8050b6d
# ╠═1c6fff5f-839a-43e1-acfe-040df2ed4822
# ╠═2f524fa1-c08b-4c11-b05a-e2d016ce6ec0
# ╠═8510da63-37d7-45d0-9012-8dd17745a9e8
# ╠═d9fca06f-ee9c-4790-a059-98667dcd8428
# ╠═04e16581-ca4f-4456-bea2-2e8bf44a08f8
# ╠═ba626ebc-d7da-4eca-8149-907d724a2ef5
# ╠═af26a6c8-9c97-4130-bd54-c9fe40e58879
# ╠═dd41fbf9-b4b9-46b9-abbf-96eddf59b0ad
# ╠═2975a6b8-bcec-4e4e-be94-5a3897c650db
# ╠═748d6916-6125-4f9a-adbf-cee3b25dc3ba
# ╠═a304f008-b3fa-4a92-a8c1-9ce2f7b4b3df
# ╠═5652b186-0e8d-469c-bfea-bf665f75c06a
# ╠═2e0eb8f9-3b6b-4464-b715-258c96e22422
# ╠═143b1c45-ce4b-4f3c-9674-ff181c64fe90
# ╠═3900999b-fdfa-48ae-afd4-bca5d03008e3
# ╠═30d37afa-7b70-479f-8dcc-860af723a96e
# ╠═9c140cb6-286b-4e25-a3dc-829b1eea9dfb
# ╠═27af9453-8044-4b1f-85fd-9e2c5dd104eb
# ╠═65a76ecd-1dcd-4d5c-88be-1b57d31369da
# ╠═6a793986-844a-4a82-84b0-b71a48aa8242
# ╠═146bcc43-15fc-45de-bb9b-f16a68d84e41
# ╠═3c7564b0-b9fe-49e8-9ad7-eaa80e9841a7
# ╠═0c424ddb-aecd-4fc6-9f26-54bb1fc9f099
# ╠═fff9cb8d-00aa-426c-bb66-4b7b21871c36
# ╠═8ace5afe-29da-4bd1-9599-b51c03b67b6e
# ╠═fab3ad29-61f6-45dd-b63f-4ea40c960afb
# ╠═3fb0f337-c6dd-4f56-a694-00d1c2f403bd
# ╠═53117290-63b5-469c-93a9-bab8a40b2658
# ╠═a6b7c896-04d8-4704-8173-e8895ff84cbd
# ╠═6d322516-d247-4bd0-bcd5-9b8371b47750
# ╠═ae756cc9-f4ca-4073-abff-c56ecd1fdade
# ╠═4c90b122-1bb9-460e-ac68-608ca8074a50
# ╠═83068836-36d9-42aa-9867-7555d89e0d42
# ╠═8082688f-ae0f-4aad-9bb9-71852be0fe00
# ╠═466354b5-7d93-4e2b-ab6d-f7a01112fbfb
# ╠═b1cf4fe2-8dda-4882-b39f-f23adb2738d1
# ╠═a15e4d30-9ed2-492a-bfa4-0608180738c1
# ╠═2015d16c-e9f8-4072-834b-a317223b854f
# ╠═313a5d22-9a70-4c20-acb2-402c1683fbf8
# ╠═b8061ceb-1610-481b-986b-c5b82b9678db
# ╠═c7cdee43-9ac6-490b-9466-367a1f136766
# ╠═bb33a5c2-555b-4078-b278-90e5ae864fed
# ╠═44604f2b-e5fd-438d-87df-0f10f874d00d
# ╠═e858791b-e869-490f-9af3-2826f480d722
# ╠═1f47a9e8-c2de-4963-ad58-aca71dbbb2ac
# ╠═105665d5-ba53-4ddd-86df-dc40d5fefa8f
# ╠═4468bf47-7668-4309-8f67-59fee2a1cf43
# ╠═4bae3659-666b-4ae7-924a-b20d98a43148
# ╠═6e03cbc9-c6cd-440b-8652-fa1be4d7934a
# ╠═53ac6ee6-f5d1-4923-9352-9d1cd96817fc
# ╠═9ec88954-f8a6-4ae0-a800-3723de537d10
# ╠═e129c967-858e-46cb-b6ce-29c238f45277
# ╠═c8cdaa0e-05cb-4b9a-99dc-6932312ec56e
# ╠═0f1066e0-d13a-4e1d-b9c0-c51509f25eb6
