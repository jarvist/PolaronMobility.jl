# CheckAlpha.jl
#   - check units in Frohlich / Feynman alpha, by testing against lit. values 

push!(LOAD_PATH,"../src/") # load module from local directory
using PolaronMobility 

" Copy and pasted out of a Jupyter notebook; this calculates 'alpha' parameters
for various materials, as a comparison to the literature used when figuring out
the oft-quoted units. "
function checkalpha()
    println(" Alpha-parameter, Cross check 'frohlichalpha()' fn vs. literature values.\n")
    println("NaCl Frohlich paper α=",frohlichalpha(2.3, 5.6, (4.9E13/(2*pi)), 1.0))
    println(" \t should be ~about 5 (Feynman1955)")
    println("CdTe  α=",frohlichalpha(7.1,   10.4,  5.08E12, 0.095))
    println("\t Stone 0.39 / Devreese 0.29 ")
    println("GaAs  α=",frohlichalpha(10.89, 12.9,  8.46E12, 0.063))
    println("\t Devreese 0.068 ")

    println()
    println("Values which were once upon a time of interest, back in 2014 before we had good phonons. Still interesting for sensitivity analysis / scaling.")
    println("Guess at PCBM: 4.0, 6.0 ; α=",frohlichalpha(4.0,6.0, 1E12, 50))
    println("MAPI:")
    println("MAPI  4.5, 24.1, 9THz ; α=",frohlichalpha(4.5,   24.1,  9.0E12,    0.12))
    println("MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=",frohlichalpha(4.5,   24.1,  2.25E12,    0.12))
    println("MAPI  6.0, 25.7, 9THz ; α=",frohlichalpha(6.0,   25.7,  9.0E12,    0.12))
    println("MAPI  6.0, 36.0, 9THz ; α=",frohlichalpha(6.0,   36,  9.0E12,    0.12))
    println("MAPI  6.0, 36.0, 1THz ; α=",frohlichalpha(6.0,   36,  1.0E12,    0.12))
end
checkalpha()

