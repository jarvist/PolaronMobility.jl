push!(LOAD_PATH,"../src/") # load module from local directory
using PolaronMobility

const T=300

#Ts,Kμs, Hμs, FHIPμs, ks, Ms, As, Bs, Cs, Fs, Taus
#effectivemass=0.12 # the bare-electron band effective-mass.
# --> 0.12 for electrons and 0.15 for holes, in MAPI. See 2014 PRB.
# MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=
println("Calculating polaron for T=$T for MAPIe params...")
MAPIe=polaronmobility(T, 4.5, 24.1, 2.25E12, 0.12)
#MAPIh=polaronmobility(T, 4.5, 24.1, 2.25E12, 0.15)

#function ImX(nurange,v,w,βred,α,ω,mb)

# Yup, this is a bit horrid.
v=MAPIe.v[1]
w=MAPIe.w[1]
βred=MAPIe.βred[1]
α=MAPIe.α[1]
ω=MAPIe.ω[1]
mb=MAPIe.mb[1]

nu=0:0.1:10
println("Integrating ImX for nu=$nu range...")
s=ImX(nu, v, w, βred,α,ω,mb)

println("Loading Plots for plotting...")
using Plots

plot( s.nu,s.ImX,label="ImX",
         markersize=3,marker=:downtriangle, xlab="nu (units Omega)",ylab="ImX")
savefig("MAPIe-ImX.png")
plot( s.nu,s.μ,label="mu",
         markersize=3,marker=:uptriangle, xlab="nu (units Omega)",ylab="Mob")
savefig("MAPIe-mu.png")

println("That's me!")
