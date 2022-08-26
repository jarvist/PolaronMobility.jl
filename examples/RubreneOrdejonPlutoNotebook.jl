### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 3df3da8e-2bc9-44f5-bbe0-2cabdf3a2a81
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate("../")

    using PolaronMobility
    using QuadGK
end


# ╔═╡ 69e7d9cd-7900-474d-aeba-06b38e62d756
using Gnuplot

# ╔═╡ 4df7f468-544d-11ec-06a3-0b757c39e918
# Rubrene Ordejon
# Rubrene electron-phonon matrix elements

# Ordejón, P., Boskovic, D., Panhans, M., Ortmann, F., 2017. Ab initio study of
# electron-phonon coupling in rubrene. Phys. Rev. B 96, 035202.
# https://doi.org/10.1103/PhysRevB.96.035202

# ╔═╡ 78d3ef38-fd4a-4957-bc9b-4523cb641df4
begin
    CM1inTHz = 33.356
    hbar = ħ = 1.0545718E-34
    q = 1.6E-19
    const Boltzmann = const kB = 1.3806504e-23
    const me = MassElectron = 9.10938188e-31
end

# ╔═╡ c4d0bfad-d837-4814-aa60-d71bd41ce28b
# Table II
# Units: cm^-1 meV | convert to THz meV
EffectiveH = [1208.9 106.8] # Holstein

# ╔═╡ 43f24329-1a03-4012-9b31-cf72d273c3d6
cm1tomeV(cm1) = 1000 * 2π * ħ * cm1 / CM1inTHz * 1E12 / q

# ╔═╡ 139b3aed-b4f8-4705-8eff-d1e0d6ea9d12
cm1tomeV(EffectiveH[1]) # I got 149.88 meV by hand 

# ╔═╡ e1573de4-450a-4ed8-a9db-003fedcc9cc3
g0 = EffectiveH[2] / cm1tomeV(EffectiveH[1]) # I got 0.71 by hand

# ╔═╡ cc9153d2-86d6-4272-a266-053165197d99
EffectiveP = [117.9 21.9]

# ╔═╡ d303b91b-567c-4705-b51e-5f32c91b3d80
gi = EffectiveP[2] / cm1tomeV(EffectiveP[1])

# ╔═╡ f192deb9-5de9-40e1-959b-a518aee7a6a3
α = gi + g0 # I got 2.21 by hand; this is essentially the Frohlich alpha

# ╔═╡ b1c36a20-9c5b-4fc2-9907-8061ed78d88f
v, w = PolaronMobility.feynmanvw(α)

# ╔═╡ 7f6aedfa-700b-47b5-8030-14d1448c4737
# We've just solved the Feynman polaron theory!

# ╔═╡ 9829e516-9240-48b2-b4f6-9eff06bf5018
EffectiveFreq = (g0 * EffectiveH[1] + gi * EffectiveP[1]) / α

# ╔═╡ edacbf0c-959f-4a2e-b406-afe7fdfa98bc
ω = 1E12 * EffectiveFreq / CM1inTHz

# ╔═╡ 6d6d286a-2f36-4762-803b-59299d7373a1
EffectiveMass(; a=1E-9, J=0.030q) = ħ^2 / (2J * a^2)
# Tight-binding model, then cos(2θ) expansion around extremum

# ╔═╡ 88f68d4c-4782-47da-ab8f-6f7ea728f23c


# ╔═╡ 79c95049-89a3-42cc-962a-e5266d8902a1
# Rubrene crystals have an orthorhombic unit cell with lattice parameters a = 1.44 nm, b = 0.718 nm, and c = 2.69 nm. 
# https://www.lehigh.edu/~inlo/rubrene.html#:~:text=Rubrene%20crystals%20have%20an%20orthorhombic,the%20a%20and%20c%20directions.
EffectiveMass(a=0.718E-9, J=0.134 * q) / MassElectron

# ╔═╡ 58828fe1-12c4-499d-837d-ae45049758e7
mb = 0.88 * MassElectron

# ╔═╡ ac1f76d8-1cdd-4de3-b1e1-049ba33374d2
βred = (q / 1000 * cm1tomeV(EffectiveFreq)) / (kB * 300)

# ╔═╡ 691ace17-eb19-444a-9807-887156a813aa


# ╔═╡ 5f62b499-3363-4d0f-8b70-6dabe5e579ed
cm1tomeV(EffectiveH[1])

# ╔═╡ 755da976-d29f-4c82-9645-7acef0a140c5
EffectiveH[1]

# ╔═╡ 78695321-c088-4e99-987c-fca646aa7af3
# Hellwarth1999 - directly do contour integration in Feynman1962, for
# finite temperature DC mobility
# Hellwarth1999 Eqn (2) and (1) - These are going back to the general
# (pre low-T limit) formulas in Feynman1962.  to evaluate these, you
# need to do the explicit contour integration to get the polaron
# self-energy
R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)

# ╔═╡ 76794859-df83-4c35-bb51-3f34774798e0
#b=R*βred/sinh(b*βred*v/2) # This self-references b! What on Earth?
# OK! I now understand that there is a typo in Hellwarth1999 and
# Biaggio1997. They've introduced a spurious b on the R.H.S. compared to
# the original, Feynman1962:
b = R * βred / sinh(βred * v / 2) # Feynman1962 version; page 1010, Eqn (47b)

# ╔═╡ 10adfb5e-d9ed-4f5b-bc1e-98b8a452e886
a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))

# ╔═╡ 1e2a33b0-624c-41a5-8c94-5f1f51fbe835
k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)

# ╔═╡ 4a36c74a-d44d-4a92-9d24-35ec22d8ae57
K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

# ╔═╡ 6f49c66c-82b8-428c-8420-4646a6394a9e
#Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
RHS = α / (3 * sqrt(π)) * βred^(5 / 2) / sinh(βred / 2) * (v^3 / w^3) * K

# ╔═╡ 7255cc2c-4680-4cde-93da-1ba2154e2d22
μ = RHS^-1 * (q) / (ω * mb) #, "m^2/Vs"

# ╔═╡ e6fe559e-46f1-4336-9513-7306c32f2515
μ * 100^2, "cm^2/Vs"

# ╔═╡ ed7cb9bb-4b83-4ed6-a8fd-9936eb02643b
function FHIPmob(T; v=v, w=w, ω=ω)
    #βred = (q/1000 * cm1tomeV(EffectiveP[1]) ) / (kB*T)
    βred = ħ * ω / (kB * T)
    R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)
    b = R * βred / sinh(βred * v / 2) # Feynman1962 version; page 1010, Eqn (47b)
    a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))
    k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)
    K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)
    #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
    RHS = α / (3 * sqrt(π)) * βred^(5 / 2) / sinh(βred / 2) * (v^3 / w^3) * K
    μ = RHS^-1 * (q) / (ω * mb) #, "m^2/Vs"
    μ = μ * 100^2 # convert to cm^2/Vs, from SI
end

# ╔═╡ 084ac351-17e6-4916-b0b7-12aa593984b2
Trange = 90:10:400

# ╔═╡ 322962f3-c1f5-47d1-882b-f4efed0c2fa0
mobs = [FHIPmob(T) for T = Trange]

# ╔═╡ 5e996c83-88dc-4569-a435-4d00684de2af
FHIPmob(300)

# ╔═╡ 19b69d34-367e-4794-abb1-3f76e1356561
begin
    @gp Trange mobs " w lp title 'Rubrene - b axis'"
    @gp :- "set yrange [0:]"
    @gp :- "set xlabel 'Temperature (K)'"
    @gp :- "set ylabel 'Mobility (cm^2/Vs)'"
end

# ╔═╡ bbe9a49b-6495-4c13-91d0-3a4dbd0a2800
begin
    @gp Trange mobs " w lp title 'Rubrene - b axis'"
    #@gp :- "set yrange [0:]"
    @gp :- "set xlabel 'Temperature (K)'"
    @gp :- "set ylabel 'Mobility (cm^2/Vs)'"
    @gp :- "set logscale y"
end

# ╔═╡ 444d7454-48a9-44dd-87a3-f4cf1a495072
@printf("\n\tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs", μ, μ * 100^2)

# ╔═╡ b0f72692-d76c-44f1-91aa-2e9b7b882d1b
# Jan's COF
Jperylenecoff = 0.204

# ╔═╡ 0c043a99-18bf-4b75-95fa-ce8aea8f8bef
relevantmode = 200 #cm^-1

# ╔═╡ 06ccb204-694d-4c64-9ecd-742c009fb972
perylenecofmass = EffectiveMass(a=1.0E-9, J=Jperylenecoff * q)

# ╔═╡ 2ba19614-f2ef-4dca-9f61-b95908b7c81c
perylenecofmass / MassElectron

# ╔═╡ 27ee4adf-b7d9-41e0-b0a0-f414422b7b73
perylenecofg0 = 128 / cm1tomeV(relevantmode)

# ╔═╡ ed6d3a1a-3fcd-4ef1-a57f-a394895d4b29
pv, pw = feynmanvw(perylenecofg0)

# ╔═╡ fbee75ec-96f0-41b6-8c14-a724872476fb
FHIPmob(300; v=pv, w=pw, ω=relevantmode / CM1inTHz * 1E12 * 2π)

# ╔═╡ 7d7d3ff6-30d4-477a-8f03-e4c08ac8a608
md"""
# Stuff below here for later
"""

# ╔═╡ ba64c0a2-8cd4-421a-b494-a70f1f92ca0f
EffectiveHGirlando = [1277 99]

# ╔═╡ d30b10f4-2c62-4596-aa73-e4f2796b0541
EffectivePGirlando = [77 20]

# ╔═╡ 04b0a526-b808-4c61-a565-6c60e12be2b7
# Table IV
# ωp      ωpg0  ωpgi
# cm^-1    meV   meV
PerModeOrdejon =
    [57.8 -1.7 0.85
        59.6 1.4 -0.83
        89.0 1.6 -4.8
        107.3 -0.14 2.8
        139.1 -2.3 -3.7
        639.1 -7.5 1.0
        1011.2 -3.6 -0.04
        1344.7 19.8 0.04
        1593.3 -42.0 -0.12
    ] .* [1 / CM1inTHz 1 1]

# ╔═╡ b38b707f-6b97-424a-a22a-cf2fdebcfb31
PerModeGirlando =
    [37.4 -0.9 3.4
        66.6 1.6 -6.6
        86.7 -0.6 -9.3
        106.3 0 -4.4
        125.1 1.4 -4.7
        631.2 -10.8 1.3
        1002.3 24.6 0
        1348.6 49.9 0
        1593.8 -45.6 1.6
    ] .* [1 / CM1inTHz 1 1]

# ╔═╡ Cell order:
# ╠═4df7f468-544d-11ec-06a3-0b757c39e918
# ╠═3df3da8e-2bc9-44f5-bbe0-2cabdf3a2a81
# ╠═78d3ef38-fd4a-4957-bc9b-4523cb641df4
# ╠═c4d0bfad-d837-4814-aa60-d71bd41ce28b
# ╠═43f24329-1a03-4012-9b31-cf72d273c3d6
# ╠═139b3aed-b4f8-4705-8eff-d1e0d6ea9d12
# ╠═e1573de4-450a-4ed8-a9db-003fedcc9cc3
# ╠═cc9153d2-86d6-4272-a266-053165197d99
# ╠═d303b91b-567c-4705-b51e-5f32c91b3d80
# ╠═f192deb9-5de9-40e1-959b-a518aee7a6a3
# ╠═b1c36a20-9c5b-4fc2-9907-8061ed78d88f
# ╠═7f6aedfa-700b-47b5-8030-14d1448c4737
# ╠═9829e516-9240-48b2-b4f6-9eff06bf5018
# ╠═edacbf0c-959f-4a2e-b406-afe7fdfa98bc
# ╠═6d6d286a-2f36-4762-803b-59299d7373a1
# ╠═88f68d4c-4782-47da-ab8f-6f7ea728f23c
# ╠═79c95049-89a3-42cc-962a-e5266d8902a1
# ╠═58828fe1-12c4-499d-837d-ae45049758e7
# ╠═ac1f76d8-1cdd-4de3-b1e1-049ba33374d2
# ╠═691ace17-eb19-444a-9807-887156a813aa
# ╠═5f62b499-3363-4d0f-8b70-6dabe5e579ed
# ╠═755da976-d29f-4c82-9645-7acef0a140c5
# ╠═78695321-c088-4e99-987c-fca646aa7af3
# ╠═76794859-df83-4c35-bb51-3f34774798e0
# ╠═10adfb5e-d9ed-4f5b-bc1e-98b8a452e886
# ╠═1e2a33b0-624c-41a5-8c94-5f1f51fbe835
# ╠═4a36c74a-d44d-4a92-9d24-35ec22d8ae57
# ╠═6f49c66c-82b8-428c-8420-4646a6394a9e
# ╠═7255cc2c-4680-4cde-93da-1ba2154e2d22
# ╠═e6fe559e-46f1-4336-9513-7306c32f2515
# ╠═ed7cb9bb-4b83-4ed6-a8fd-9936eb02643b
# ╠═084ac351-17e6-4916-b0b7-12aa593984b2
# ╠═322962f3-c1f5-47d1-882b-f4efed0c2fa0
# ╠═5e996c83-88dc-4569-a435-4d00684de2af
# ╠═69e7d9cd-7900-474d-aeba-06b38e62d756
# ╠═19b69d34-367e-4794-abb1-3f76e1356561
# ╠═bbe9a49b-6495-4c13-91d0-3a4dbd0a2800
# ╠═444d7454-48a9-44dd-87a3-f4cf1a495072
# ╠═b0f72692-d76c-44f1-91aa-2e9b7b882d1b
# ╠═0c043a99-18bf-4b75-95fa-ce8aea8f8bef
# ╠═06ccb204-694d-4c64-9ecd-742c009fb972
# ╠═2ba19614-f2ef-4dca-9f61-b95908b7c81c
# ╠═27ee4adf-b7d9-41e0-b0a0-f414422b7b73
# ╠═ed6d3a1a-3fcd-4ef1-a57f-a394895d4b29
# ╠═fbee75ec-96f0-41b6-8c14-a724872476fb
# ╠═7d7d3ff6-30d4-477a-8f03-e4c08ac8a608
# ╠═ba64c0a2-8cd4-421a-b494-a70f1f92ca0f
# ╠═d30b10f4-2c62-4596-aa73-e4f2796b0541
# ╠═04b0a526-b808-4c61-a565-6c60e12be2b7
# ╠═b38b707f-6b97-424a-a22a-cf2fdebcfb31
