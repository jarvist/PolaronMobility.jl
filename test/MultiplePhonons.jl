# Check values of the new multimode polaron codes

@testset verbose = true "MultiplePhonons" begin

    # Physical constants
    ħ = 1.05457162825e-34          # kg m2 / s
    eV = 1.602176487e-19           # kg m2 / s2
    me = 9.10938188e-31            # kg
    kB = 1.3806504e-23            # kg m2 / K s2
    ε_0 = 8.854E-12                 # Units: C2N−1m−2, permittivity of free space
    c = 3e8

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
        # 0.022939578929507105 8.355742795827834e-05   # Acoustic modes!
        # 0.04882611767873102 8.309858592685e-06
        # 0.07575149723846182 2.778248540373041e-05
    ]

    volume = (6.29e-10)^3
    ϵ_optic = 4.5
    m_eff = 0.12
    phonon_freq = MAPI[:, 1]
    ir_activity = MAPI[:, 2]
    ϵ_static = 24.1

    # Ionic dielectric and decomposed alpha

    ω = 2π .* phonon_freq
    ϵ_ionic = [ϵ_ionic_mode(i, j, volume) for (i, j) in zip(phonon_freq, ir_activity)]

    α = [frohlichalpha(ϵ_optic, i, sum(ϵ_ionic), j, m_eff) for (i, j) in zip(ϵ_ionic, phonon_freq)]

    @testset "Ionic dielectric and decomposed alphas" begin

        @test ϵ_ionic ≈ [0.2999680470756664, 0.0247387647244569, 0.2543132184018061, 0.16621617310133838, 2.3083204422506296, 3.0707979601813267, 3.2026782087486407, 0.892674135624958, 0.8861579096771846, 0.7209278375756829, 0.40199759805819046, 0.6183279038315278, 0.7666391525823296, 0.1555444994147378, 1.1627710200840813] rtol = 1e-3

        @test α ≈ [0.03401013445306177, 0.002850969846883158, 0.03075081562607006, 0.02275292381052607, 0.33591418423943553, 0.46526818717696034, 0.5046331098089347, 0.14223560721522646, 0.16083871312929882, 0.1622911897190622, 0.09123913334086006, 0.14070972715961402, 0.18159786148348117, 0.039500396977008016, 0.3487735470060312] rtol = 1e-3

        @test sum(α) ≈ 2.663366500992453 rtol = 1e-3

        println('\n', "Ionic dielectric and alphas:")
        println("-------------------------------")
        println("  ω (THz)    ϵ_ionic     αs")
        println("-------------------------------")
        display(hcat(ω, ϵ_ionic, α))
        println("\n-------------------------------")
    end

    # Variations

    athermal_energy(v, w) = frohlich_energy(v, w, α, ω)
    v_0, w_0, E_0, A_0, B_0, C_0 = vw_variation(athermal_energy, 4, 3) # Athermal

    β = ħ / (kB * 300) * 1e12
    thermal_energy(v, w) = frohlich_energy(v, w, α, ω, β)
    v, w, E = vw_variation(thermal_energy, v_0, w_0) # Thermal

    @testset "Multiple mode variations" begin

        @test v_0[1] ≈ 3.292283619446986 rtol = 1e-3
        @test w_0[1] ≈ 2.679188425097246 rtol = 1e-3

        @test v[1] ≈ 35.19211042393129 rtol = 1e-3
        @test w[1] ≈ 32.454157668863225 rtol = 1e-3

        println("\nVariational parameters:")
        println("Athermal v = $(v_0[1]) | athermal w = $(w_0[1])")
        println("300K v = $(v[1]) | 300K w = $(w[1])")
    end

    # Energies

    E_0 *= 1000 * ħ / eV * 1e12 # Enthalpy

    E *=  1000 * ħ / eV * 1e12 # Free energy

    @testset "Multiple mode energies" begin

        @test E_0 ≈ -19.51150570496287 rtol = 1e-3

        @test E ≈ -42.79764110613318 rtol = 1e-3

        println("\nPolaron energies")
        println("0K energy: $E_0 meV")
        println("300K energy: $E meV")
    end

    # Mobility

    μ = frohlich_mobility(v, w, α, ω, β) * eV / (1e12 * me * m_eff) * 100^2

    @testset "Multiple mode mobility" begin

        @test μ ≈ 160.39270615550004 rtol = 1e-3

        println("\nMobility at 300K: $μ")
    end

    # Impedence 

    Z_dc = frohlich_complex_impedence(0, v, w, α, ω, β) # DC limit 

    Z = frohlich_complex_impedence(50, v, w, α, ω, β) # AC limit

    @testset "Multiple mode impedence" begin

        @test Z_dc ≈ 91.38457766169273 + 0.0im rtol = 1e-3

        @test Z ≈ 76.18909903384355 - 23.457492887508455im rtol = 1e-3

        println("\nComplex Impedence:")
        println("DC limit: $Z_dc")
        println("At 50 THz: $Z")
    end

    # Conductivity

    σ_dc = frohlich_complex_conductivity(0, v, w, α, ω, β) # DC limit

    σ = frohlich_complex_conductivity(50, v, w, α, ω, β) # AC limit

    @testset "Multiple mode conductivity" begin

        @test σ_dc ≈ 0.0109427654598571 - 0.0im rtol = 1e-3

        @test σ ≈ 0.011988781430646566 + 0.003691167879729921im rtol = 1e-3

        println("\nComplex Conductivity:")
        println("DC limit: $σ_dc")
        println("At 50 THz: $σ\n")
    end

    # Single Mode Tests

    Hellwarth_B_freq = HellwarthBScheme(hcat(phonon_freq, ir_activity))

    @testset "Hellwarth effective frequency" begin

        @test Hellwarth_B_freq ≈ 2.25 rtol = 1e-3
    end

    @testset "Single effective mode MAPI" begin

        MAPIs = material(ϵ_optic, ϵ_static, m_eff, Hellwarth_B_freq)

        singlemode_polaron = frohlichpolaron(MAPIs, [0, 300], 3, verbose = true)

        display(singlemode_polaron)

        @test singlemode_polaron.α ≈ 2.393156008589176 rtol = 1e-3
        @test singlemode_polaron.v ≈ [3.3086321909253815, 19.847527941504232] rtol = 1e-3
        @test singlemode_polaron.w ≈ [2.6634018089397014, 16.948211346396874] rtol = 1e-3
        @test singlemode_polaron.F ≈ [-5.571343967852482, -8.585426348039439] rtol = 1e-3
        @test singlemode_polaron.χ ≈ [1.8185958394302872 + 2.8215456392868816im, -2.6431955087396575 + 16.602938172029397im] rtol = 1e-3
        @test singlemode_polaron.z ≈ [0.3385854767144258 - 0.5782315007316344im, 1.9923525806435276 - 0.0428165389512411im] rtol = 1e-3
        @test singlemode_polaron.σ ≈ [0.7541017043762221 + 1.287844252674551im, 0.5016874943626112 + 0.010781486345548893im] rtol = 1e-3
        @test singlemode_polaron.μ ≈ [Inf, 0.4873858565327472] rtol = 1e-3
    end

    # Multimode Tests

    @testset "Multiple mode (15) MAPI" begin

        MAPIm = material(ϵ_optic, ϵ_static, m_eff, phonon_freq, ir_activity, volume)

        multimode_polaron = frohlichpolaron(MAPIm, [0, 300], 3, verbose = true)

        display(multimode_polaron)

        @test multimode_polaron.α ≈ [0.03401013640864639, 0.0028509700108140822, 0.030750817394243672, 0.02275292511882046, 0.33591420355451984, 0.46526821392990697, 0.5046331388253664, 0.14223561539378177, 0.16083872237753374, 0.16229119905081466, 0.0912391385871153, 0.14070973525043112, 0.18159787192536828, 0.03950039924828304, 0.3487735670605295] rtol = 1e-3
        @test sum(multimode_polaron.α) ≈ 2.663366654136176 rtol = 1e-3
        @test multimode_polaron.v ≈ [3.292268504574272, 35.19207894065109] rtol = 1e-3
        @test multimode_polaron.w ≈ [2.679193078382525, 32.45420279136673] rtol = 1e-3
        @test multimode_polaron.F ≈ [-4.719036156508013, -10.357755243151143] rtol = 1e-3
        @test multimode_polaron.χ ≈ [1.0521976989565067 + 2.089861861859462im, -2.1193321531426945 + 14.077489512387219im] rtol = 1e-3
        @test multimode_polaron.z ≈ [0.2507834234231354 - 0.4862637238747808im, 1.6892987414864662 - 0.10568014162287666im] rtol = 1e-3
        @test multimode_polaron.σ ≈ [0.8377746271072964 + 1.624427182564005im, 0.5896539523503038 + 0.036887917845742565im] rtol = 1e-3
        @test multimode_polaron.μ ≈ [Inf, 0.5729621483218189] rtol = 1e-3
    end

    print('\n')

end

