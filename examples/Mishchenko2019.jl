# Mishchenko2019.jl - reproducing Mishchenko et al. 2019 PRL, with the Feynman
# variational polaron technique

println("Loading PolaronMobility...")
using PolaronMobility
println("Loading Plots for plotting...")
using Plots
gr(size=(500,375))


#function ImX(nurange,v,w,βred,α,ω,mb)

function savefigs(name)
    println("Saving PNG")
    savefig("$(name).png")
    println("Saving PDF")
    savefig("$(name).pdf")
# high res conversion of vector PDF -> suitable, anti-aliased, PNG for slides
    println("Convert PDF to high res PNG")
    run(`convert -density 300 -resize 1600x $(name).pdf $(name).pdf.png`)
end

# Specify problem 
βred=5 # temperature, in thermodynamic Beta units
mb=0.12*9.1093837015E-31
ω=1

# FHIP1962 reproduction
# Fig 2 FHIP1962
α=5
#in [6] #1:1:10
v,w = feynmanvw(α)
println("α: $(α) v: $(v) w: $(w)")

nu=3.5:0.2:22 
println("Integrating ImX for nu=$nu range...")
s=ImX(nu, v, w, βred,α,ω,mb)
plot( s.nu,s.ImX,label="ImX",
         markersize=3,marker=:downtriangle, xlab="nu (units Omega)",ylab="ImX")
savefigs("FHIP1962_Fig2")

# Fig 3 FHIP1962
α=7
#in [6] #1:1:10
v,w = feynmanvw(α)
println("α: $(α) v: $(v) w: $(w)")

nu=3.5:0.2:28
println("Integrating ImX for nu=$nu range...")
s=ImX(nu, v, w, βred,α,ω,mb)
plot( s.nu,s.ImX,label="ImX",
         markersize=3,marker=:downtriangle, xlab="nu (units Omega)",ylab="ImX")



# Athernmal Feynman variational technique, to compare best to the DiagMC
# results, though nb. they have some kind of T-dep in the diagrams

# Specify problem
#mb=1
βred=8 # temperature, in thermodynamic Beta units
ω=1

nurange=0.0:0.05:20

for α in [6] #1:1:10
   # athermal action
    v,w = feynmanvw(α)
    plot(xlab="nu (units Omega)",ylab="Mob", ylims=(0,1.0E12))

    for T in [8,4,2,1,0.5,0.25,0.125]
        βred=1/T
        @time s=ImX(nurange, v, w, βred,α,ω,mb)
        plot!( s.nu, s.μ, label="T=$(T)")
    end
    
    savefigs("Mishchenko-Fig4-AthermalAction_alpha_$(α)")

    plot(xlab="nu (units Omega)",ylab="Mob", ylims=(0,3.0E12))
    
    for T in [8,4,2,1,0.5,0.25,0.125]
        βred=1/T
        # Osaka finite temperature action
        v,w=feynmanvw(α, βred, verbose=true) # temperature dependent Action

        @time s=ImX(nurange, v, w, βred,α,ω,mb)
        plot!( s.nu, s.μ, label="T=$(T) v=$(round(v,sigdigits=2)) w=$(round(w,sigdigits=2))")
    end

    savefigs("Mishchenko-Fig4-OsakaFiniteTemperatureAction_alpha_$(α)")
end

println("That's me!")
