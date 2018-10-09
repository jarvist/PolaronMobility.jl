# Cross-check against data published in Frost 2017 PRB. 

@testset "FrostPolaronMobility2017" begin

T=300

#Ts,Kμs, Hμs, FHIPμs, ks, Ms, As, Bs, Cs, Fs, Taus
#effectivemass=0.12 # the bare-electron band effective-mass.
# --> 0.12 for electrons and 0.15 for holes, in MAPI. See 2014 PRB.
# MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=
MAPIe=polaronmobility(T, 4.5, 24.1, 2.25E12, 0.12)
MAPIh=polaronmobility(T, 4.5, 24.1, 2.25E12, 0.15)

# Hellwarth mobility
@test MAPIe.Hμ[1] ≈ 136.42 rtol=0.02
# Test variational parameters
@test MAPIe.v[1] ≈ 19.86 rtol=0.02
@test MAPIe.w[1] ≈ 16.96 rtol=0.02

# Same for the MAPI holes @ 300 K
@test MAPIh.Hμ[1] ≈ 94.15 rtol=0.02
@test MAPIh.v[1] ≈ 20.09  rtol=0.02
@test MAPIh.w[1] ≈ 16.81  rtol=0.02

end

