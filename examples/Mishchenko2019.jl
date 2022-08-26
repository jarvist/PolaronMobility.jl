# Mishchenko2019.jl - reproducing Mishchenko et al. 2019 PRL, with the Feynman
# variational polaron technique

println("Loading PolaronMobility...")
using PolaronMobility
println("Loading Plots for plotting...")
using Plots
gr(size=(500, 375))

using Printf
import QuadGK.quadgk # one ring to integrate them all...

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
βred = 5 # temperature, in thermodynamic Beta units
mb = 0.12 * 9.1093837015E-31
ω = 1
const eV = const q = const ElectronVolt = 1.602176487e-19;      # kg m2 / s2
const Boltzmann = const kB = 1.3806504e-23;                  # kg m2 / K s2
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s

function FHIP1962_Fig2()
    # FHIP1962 reproduction
    # Fig 2 FHIP1962
    α = 5
    #in [6] #1:1:10
    v, w = feynmanvw(α)
    println("α: $(α) v: $(v) w: $(w)")

    nu = 3.5:0.2:22
    println("Integrating ImX for nu=$nu range...")
    s = ImX(nu, v, w, βred, α, ω, mb)
    plot(s.nu, s.ImX, label="ImX",
        markersize=3, marker=:downtriangle, xlab="nu (units Omega)", ylab="ImX")
    savefigs("FHIP1962_Fig2")
end

function FHIP1962_Fig3()
    # Fig 3 FHIP1962
    α = 7
    #in [6] #1:1:10
    v, w = feynmanvw(α)
    println("α: $(α) v: $(v) w: $(w)")

    nu = 3.5:0.2:28
    println("Integrating ImX for nu=$nu range...")
    s = ImX(nu, v, w, βred, α, ω, mb)
    plot(s.nu, s.ImX, label="ImX",
        markersize=3, marker=:downtriangle, xlab="nu (units Omega)", ylab="ImX")
end


# Athernmal Feynman variational technique, to compare best to the DiagMC
# results, though nb. they have some kind of T-dep in the diagrams

function Mishchenko_FigureFour()
    # Specify problem
    #mb=1
    βred = 8 # temperature, in thermodynamic Beta units
    ω = 1

    nurange = 0.0:0.05:20
    Trange = [8, 4, 2, 1, 0.5, 0.25, 0.125]

    #using Distributed
    #addprocs(6)
    #@distributed 

    for α in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20] #1:1:10
        # athermal action
        v, w = feynmanvw(α)
        plot(xlab="nu (units Omega) alpha=$(α)", ylab="Mob", ylims=(0, 1.0E12))

        for T in Trange
            βred = 1 / T
            @time s = ImX(nurange, v, w, βred, α, ω, mb)
            plot!(s.nu, s.μ, label="T=$(T)")
        end

        savefigs("Mishchenko-Fig4-AthermalAction_alpha_$(α)")

        plot(xlab="nu (units Omega) alpha=$(α)", ylab="Mob", ylims=(0, 3.0E12))

        for T in Trange
            βred = 1 / T
            # Osaka finite temperature action
            v, w = feynmanvw(α, βred, verbose=true) # temperature dependent Action

            @time s = ImX(nurange, v, w, βred, α, ω, mb)
            plot!(s.nu, s.μ, label="T=$(T) v=$(round(v,sigdigits=2)) w=$(round(w,sigdigits=2))")
        end

        savefigs("Mishchenko-Fig4-OsakaFiniteTemperatureAction_alpha_$(α)")
    end
end
#wait()

function HellwarthMobility(v, w, βred, α, ω)
    # Hellwarth1999 - directly do contour integration in Feynman1962, for
    # finite temperature DC mobility
    # Hellwarth1999 Eqn (2) and (1) - These are going back to the general
    # (pre low-T limit) formulas in Feynman1962.  to evaluate these, you
    # need to do the explicit contour integration to get the polaron
    # self-energy
    R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)

    b = R * βred / sinh(βred * v / 2) # Feynman1962 version; page 1010, Eqn (47b)

    a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))
    k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)
    K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

    #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
    RHS = α / (3 * sqrt(π)) * βred^(5 / 2) / sinh(βred / 2) * (v^3 / w^3) * K
    #μ=RHS^-1 * (q)/(ω*mb)
    #@printf("\n\tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)
    RHS^-1
end


function Mischenko_MobilityComparison(α; Trange=0.1:0.01:25)
    println("Mishchenko Feynman athermal Hellwarth explicit integration mobility with α= $(α)")
    βred = 8
    ω = 1 #(50*kB)/ħ # in artifical units of T
    fout = open("MishchenkoMobility_$(α).dat", "w")
    @printf(fout, "# T v w mob(SI)")

    for T in Trange
        β = ω / T #  ħ*ω/(kB*T)

        # Athermal variational
        v, w = feynmanvw(α)
        μ = HellwarthMobility(v, w, β, α, ω)
        @printf(fout, "%f    %f %f %f ", T, v, w, μ)

        # Thermal theory
        v, w = feynmanvw(α, β)
        μ = HellwarthMobility(v, w, β, α, ω)
        @printf(fout, "%f %f %f ", v, w, μ)
        @printf(fout, "\n")
    end
    close(fout)
end

for α in [2.4, 4, 6, 8, 10, 12, 14] # Mishchenko2019 Fig 2 and 3
    Mischenko_MobilityComparison(α)
end

println("That's me!")

