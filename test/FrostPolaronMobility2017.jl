# Cross-check against data published in Frost 2017 PRB. 

@testset "FrostPolaronMobility2017" begin

    ϵ_optic = 4.5
    electron_m_eff = 0.12
    hole_m_eff = 0.15
    phonon_freq = 2.25
    ϵ_static = 24.1

    temperature = 300

    MAPIe = material(ϵ_optic, ϵ_static, electron_m_eff, phonon_freq)
    MAPIh = material(ϵ_optic, ϵ_static, hole_m_eff, phonon_freq)

    MAPIe_polaron = frohlichpolaron(MAPIe, temperature)
    MAPIh_polaron = frohlichpolaron(MAPIh, temperature)

    addunits!(MAPIe_polaron)
    addunits!(MAPIh_polaron)

    # Hellwarth mobility
    @test ustrip(MAPIe_polaron.μ |> u"cm^2/V/s") ≈ 136.42 rtol=0.02

    # Test variational parameters
    @test ustrip(MAPIe_polaron.v) / phonon_freq ≈ 19.86 rtol=0.02
    @test ustrip(MAPIe_polaron.w) / phonon_freq ≈ 16.96 rtol=0.02

    # Same for the MAPI holes @ 300 K
    @test ustrip(MAPIh_polaron.μ |> u"cm^2/V/s") ≈ 94.15 rtol=0.02
    @test ustrip(MAPIh_polaron.v) / phonon_freq ≈ 20.09  rtol=0.02
    @test ustrip(MAPIh_polaron.w) / phonon_freq ≈ 16.81  rtol=0.02

end

