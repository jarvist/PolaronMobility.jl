# FiniteTemperatures.jl
# Checks finite temperature algorithms for ImX and ReX with MAPI data.
include("../src/PolaronMobility.jl")
plotly()
Plots.PlotlyBackend()
struct susceptibility
    nu
    ImX
    μ
end
Susceptibility()=susceptibility([],[],[])

# Physical constants
const T = 1
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2
const ϵ_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space
const amu = 1.660_539_066_60e-27 # kg

MAPIe = polaronmobility(T, 4.5, 24.1, 2.25E12, 0.12)
#MAPIh = polaronmobility(T, 4.5, 24.1, 2.25E12, 0.15)

s = Susceptibility()

# Yup, this is a bit horrid. Agreed?
v = MAPIe.v[1]
w = MAPIe.w[1]
βred = MAPIe.βred[1]
α = MAPIe.α[1]
ω = MAPIe.ω[1]
mb = MAPIe.mb[1]
@show(βred, v, w, α)
Ω_range = 1:0.5:50
println("Integrating Imχ for Ω = $Ω_range range...")
ImagX = [PolaronMobility.ℑχ(Ω, βred, α, v, w) for Ω in Ω_range]
μ = [x^-1 * (q) / (ω * mb) for x in ImagX]
# s = ImX(nu, v, w, βred,α,ω,mb)
append!(s.nu, Ω_range)
append!(s.ImX, ImagX)
append!(s.μ, μ)

println("Loading Plots for plotting...")
using Plots

p1 = plot(s.nu, s.ImX, label="ImX",
         markersize=3,marker=:downtriangle, xlab="nu (units Omega)", ylab="ImX")
yaxis!(:log10)
display(p1)
# savefig("MAPIe-ImX.png")
p2 = plot(s.nu, s.μ, label="mu",
         markersize=3,marker=:uptriangle, xlab="nu (units Omega)",ylab="Mob")
yaxis!(:log10)
display(p2)

# savefig("MAPIe-mu.png")

println("That's me!")
