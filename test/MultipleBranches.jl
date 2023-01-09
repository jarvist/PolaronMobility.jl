@testset "MultipleBranches" begin

    # ((freq THz)) ((IR Activity / e^2 amu^-1))
    # These data from MAPbI3-Cubic_PeakTable.csv
    # https://github.com/WMD-group/Phonons/tree/master/2015_MAPbI3/SimulatedSpectra
    # Data published in Brivio2015 (PRB)
    # https://doi.org/10.1103/PhysRevB.92.144308
    MAPI = [
        96.20813558773261 0.4996300522819191
        93.13630357703363 1.7139631746083817
        92.87834578121567 0.60108592692181
        92.4847918585963 0.0058228799414729
        92.26701437594754 0.100590086574602
        89.43972834606603 0.006278895133832249
        46.89209141511332 0.2460894564364346
        46.420949316788 0.14174282581124137
        44.0380222871706 0.1987196948553428
        42.89702947649343 0.011159939465770681
        42.67180170168193 0.02557751102757614
        41.46971205834201 0.012555230726601503
        37.08982543385215 0.00107488277468418
        36.53555265689563 0.02126940080871224
        30.20608114002676 0.009019481779712388
        27.374810898415028 0.03994453721421388
        26.363055017011728 0.05011922682554448
        9.522966890022039 0.00075631870522737
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
        #0.022939578929507105 8.355742795827834e-05   # Acoustic modes!
        #0.04882611767873102 8.309858592685e-06
        #0.07575149723846182 2.778248540373041e-05
    ]

    vol = (6.29E-10)^3
    ϵ_o = 4.5
    meff = 0.12

    ϵ_i = IRtoDielectric(MAPI, vol)
    #ϵ_ired=ϵ_i/ϵ_0
    #ϵ_s=ϵ_ired + ϵ_o

    ϵ_s = sum(ϵ_i) + ϵ_o # total (static) dielectric = sum of ionic, and optical

    println("Sum of ionic dielectric: $(ϵ_s)")

    IRtoalpha(MAPI, volume=vol, ϵ_o=ϵ_o, ϵ_s=ϵ_s, meff=meff)


    println()
    splat = DielectricFromIRmode.(eachrow(MAPI), volume=vol)
    println(splat)
    println("Sum of dieletric: ", sum(splat))

    f_dielectric = hcat(MAPI[:, 1], splat)
    println("f_dielectric: ", f_dielectric)

    alphas = frohlichPartial.(eachrow(f_dielectric), ϵ_o=ϵ_o, ϵ_s=ϵ_o + sum(splat), meff=meff)
    println("alphas: ", alphas)
    println("sum alphas: ", sum(alphas))

    println("Feynman v,w for alphas: ", feynmanvw.(3.1, 3.0, alphas))

    mobilityproblem = hcat(alphas, feynmanvw.(3.1, 3.0, alphas), MAPI[:, 1])
    println("Specify mobility problem: ", mobilityproblem)

    inverse_μ = Hellwarth1999mobilityRHS.(eachrow(mobilityproblem), meff, 300)

    μ = sum(inverse_μ)^-1
    @printf("\n\tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs", μ, μ * 100^2)
    println()

    ħ = 1.05457162825e-34
    kB = 1.3806504e-23

    MAPI = [
        5 20.0
        0.5 0.67
    ] # fictious two-frequency MAPI-a-like

    splat = DielectricFromIRmode.(eachrow(MAPI), volume=vol)
    println(splat)
    println("Sum of dieletric: ", sum(splat))

    f_dielectric = hcat(MAPI[:, 1], splat)

    for T in 10:1:500
        alphas = frohlichPartial.(eachrow(f_dielectric), ϵ_o=ϵ_o, ϵ_s=ϵ_o + sum(splat), meff=meff)
        βred = ħ * MAPI[:, 1] * 1E12 * 2π / (kB * T)

        mobilityproblem = hcat(alphas, feynmanvw.(3.1, 3.0, alphas, βred), MAPI[:, 1])
        inverse_μ = Hellwarth1999mobilityRHS.(eachrow(mobilityproblem), meff, T)
        μ = sum(inverse_μ)^-1

        @printf("\nT= %d μ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs", T, μ, μ * 100^2)
    end

    println()
    # A bit cyclical, but this is extrapolating back to if MAPI had a single LO
    # mode, at the point where the Hellwarth1999 averaging approximation puts it.
    MAPI_singlemode = [
    2.25 1.67502750447212449732]

    # By comparison to 2017 PRB, this should be 24.1-4.5 = 19.6
    #  Nb: however, that dielectric is from my 2014 NanoLetters, whereas the
    #  2015PRB above is a different calculation, with different convergence.
    ϵ_i = IRtoDielectric(MAPI_singlemode, vol)
    ϵ_s = sum(ϵ_i) + ϵ_o
    println("Sum of ionic dielectric: $(ϵ_s)")
    αmode_MAPIe = IRtoalpha(MAPI_singlemode, volume=vol, ϵ_o=ϵ_o, ϵ_s=ϵ_s, meff=meff)

    α_MAPIe = frohlichalpha(4.5, 24.1, 2.25, meff)

    @test αmode_MAPIe ≈ α_MAPIe atol = 0.01


    f = IRtoalpha(MAPI, volume=vol, ϵ_o=ϵ_o, ϵ_s=ϵ_s, meff=meff)

end # @testset

