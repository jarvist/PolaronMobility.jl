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

    ω = phonon_freq .* 2π
    ϵ_ionic = [ϵ_ionic_mode(i, j, volume) for (i, j) in zip(phonon_freq, ir_activity)]

    α = [frohlichalpha(ϵ_optic, i, sum(ϵ_ionic), j, m_eff) for (i, j) in zip(ϵ_ionic, phonon_freq)]

    Hellwarth_B_freq = HellwarthBScheme(hcat(phonon_freq, ir_activity))

    @testset "Ionic dielectric and decomposed alphas" begin

        @test ϵ_ionic ≈ [0.299961772906185, 0.024738247285335784, 0.2543078991544757, 0.16621269650291148, 2.3082721611293127, 3.07073373098414, 3.2026112211275524, 0.8926554643402222, 0.8861393746866522, 0.7209127585581645, 0.40198918982576065, 0.6183149708071556, 0.7666231174611848, 0.1555412460263828, 1.1627466994188664] rtol = 1e-3

        @test α ≈ [0.03400925262932455, 0.002850895926184945, 0.030750018310811124, 0.0227523338667158, 0.33590547456785985, 0.4652561235806758, 0.5046200255485525, 0.14223191929288398, 0.1608345428607191, 0.16228698179028636, 0.09123676766854395, 0.1407060788006921, 0.1815931529662271, 0.03949937280029587, 0.34876450391350566] rtol = 1e-3

        @test sum(α) ≈ 2.6632974445232787 rtol = 1e-3

        println('\n', "Ionic dielectric and alphas:")
        println("-------------------------------")
        println("  ω (THz)    ϵ_ionic     αs")
        println("-------------------------------")
        display(hcat(ω, ϵ_ionic, α))
        println("\n-------------------------------")
    end

    # Variations

    v_0, w_0, E_0, A_0, B_0, C_0 = feynmanvw(α, ω) # Athermal

    β = ħ / (kB * 300) * 1e12
    v, w, E, A, B, C = feynmanvw(α, ω, β) # Thermal

    @testset "Multiple mode variations" begin

        @test v_0[1] / Hellwarth_B_freq ≈ 16.283245452806543 rtol = 1e-3
        @test w_0[1] / Hellwarth_B_freq ≈ 12.93471218329606 rtol = 1e-3

        @test v[1] / Hellwarth_B_freq ≈ 123.90465526352604 rtol = 1e-3
        @test w[1] / Hellwarth_B_freq ≈ 106.7711782241155 rtol = 1e-3

        println("\nVariational parameters:")
        println("Athermal v = $(v_0[1]) | athermal w = $(w_0[1])")
        println("300K v = $(v[1]) | 300K w = $(w[1])")
    end

    # Energies

    E_0 *= 1000 * ħ / eV * 1e12 # Enthalpy

    E *=  1000 * ħ / eV * 1e12 # Free energy

    @testset "Multiple mode energies" begin

        @test E_0 ≈ -19.54914965269833 rtol = 1e-3

        @test E ≈ -63.315841394023664 rtol = 1e-3

        println("\nPolaron energies")
        println("0K energy: $E_0 meV")
        println("300K energy: $E meV")
    end

    # Mobility

    μ = frohlich_mobility(v, w, α, ω, β) * eV / (1e12 * me * m_eff) * 100^2

    @testset "Multiple mode mobility" begin

        @test μ ≈ 145.35917679067404  rtol = 1e-3

        println("\nMobility at 300K: $μ")
    end

    # Impedence 

    Z_dc = frohlich_complex_impedence(0, v, w, α, ω, β) # DC limit 

    Z = frohlich_complex_impedence(50, v, w, α, ω, β) # AC limit

    @testset "Multiple mode impedence" begin

        @test Z_dc ≈ 100.83185205569625 + 0.0im  rtol = 1e-3

        @test Z ≈ 82.51429805705051 - 18.759832224634735im rtol = 1e-3

        println("\nComplex Impedence:")
        println("DC limit: $Z_dc")
        println("At 50 THz: $Z")
    end

    # Conductivity

    σ_dc = frohlich_complex_conductivity(0, v, w, α, ω, β) # DC limit

    σ = frohlich_complex_conductivity(50, v, w, α, ω, β) # AC limit

    @testset "Multiple mode conductivity" begin

        @test σ_dc ≈ 0.009917501063529335 - 0.0im rtol = 1e-3

        @test σ ≈ 0.011523473106500403 + 0.0026198904579370222im rtol = 1e-3

        println("\nComplex Conductivity:")
        println("DC limit: $σ_dc")
        println("At 50 THz: $σ\n")
    end

    # Single Mode Tests

    @testset "Hellwarth effective frequency" begin

        @test Hellwarth_B_freq ≈ 2.25 rtol = 1e-3
    end

    @testset "Single effective mode MAPI" begin

        MAPIs = material(ϵ_optic, ϵ_static, m_eff, Hellwarth_B_freq)

        singlemode_polaron = frohlichpolaron(MAPIs, [0, 300], 3, verbose = true)

        display(singlemode_polaron)

        @test singlemode_polaron.α ≈ 2.393156008589176 rtol = 1e-3
        @test singlemode_polaron.v ≈ [7.449619526820557, 44.68801621691182] rtol = 1e-3
        @test singlemode_polaron.w ≈ [5.996837775164936, 38.160005552975356] rtol = 1e-3
        @test singlemode_polaron.F ≈ [-5.5702304614638845, -16.23756639407282] rtol = 1e-3
        @test singlemode_polaron.χ ≈ [1.8185958394302872 + 2.8215456392868816im, -2.6431955087396575 + 16.602938172029397im] rtol = 1e-3
        @test singlemode_polaron.z ≈ [2.821550935533812 - 4.818576894764461im, 16.60369379393479 - 0.356819708523505im] rtol = 1e-3
        @test singlemode_polaron.σ ≈ [0.09049281752137871 + 0.15454145950697004im, 0.06019975971109561 + 0.0012937157827582215im] rtol = 1e-3
        @test singlemode_polaron.μ ≈ [Inf, 0.058486537513714805] rtol = 1e-3
    end

    # Multimode Tests

    @testset "Multiple mode (15) MAPI" begin

        MAPIm = material(ϵ_optic, ϵ_static, m_eff, phonon_freq, ir_activity, volume)

        multimode_polaron = frohlichpolaron(MAPIm, [0, 300], 3, verbose = true)

        display(multimode_polaron)

        @test multimode_polaron.α ≈ [0.03401013640864639, 0.0028509700108140822, 0.030750817394243672, 0.02275292511882046, 0.33591420355451984, 0.46526821392990697, 0.5046331388253664, 0.14223561539378177, 0.16083872237753374, 0.16229119905081466, 0.0912391385871153, 0.14070973525043112, 0.18159787192536828, 0.03950039924828304, 0.3487735670605295] rtol = 1e-3
        @test sum(multimode_polaron.α) ≈ 2.663366654136176 rtol = 1e-3
        @test multimode_polaron.v ≈ [5.835086623282841, 44.40143374273715] rtol = 1e-3
        @test multimode_polaron.w ≈ [4.635142386961883, 38.26163944180995] rtol = 1e-3
        @test multimode_polaron.F ≈ [-4.728298068562315, -15.318577804140983]  rtol = 1e-3
        @test multimode_polaron.χ ≈ [1.0083819418968094 + 2.262844921017003im, -2.5135009156247414 + 15.484273072268337im] rtol = 1e-3
        @test multimode_polaron.z ≈ [2.262844921017003 - 4.008381941896809im, 15.484273072268337 - 0.4864990843752586im] rtol = 1e-3
        @test multimode_polaron.σ ≈ [0.10680047179649949 + 0.18918533857934283im, 0.0645179673929453 + 0.002027084637162277im] rtol = 1e-3
        @test multimode_polaron.μ ≈ [Inf, 0.062313430372364004] rtol = 1e-3
    end

    print('\n')

end

