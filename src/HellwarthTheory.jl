# HellwarthTheory.jl

#####
"""
	HellwarthBScheme(LO)

    Multiple phonon mode reduction to a single effective frequency.
	Hellwarth et al. 1999 PRB, 'B scheme'; the athermal method.
    Averaging procedure is constructed by considering the average effect of the action of multiple branches.

	Follows Eqn (58) in this paper, assuming typo on LHS, should actually be W_e.
"""
function HellwarthBScheme(LO)
    println("Hellwarth B Scheme... (athermal)")
    H58 = sum( (LO[:,2].^2)./ LO[:,1].^2 )
    println("Hellwarth (58) summation: ",H58)

    H59 = sum( LO[:,2].^2 ) # sum of total ir activity squarred
    println("Hellwarth (59) summation (total ir activity ^2): ", H59)
    println("Hellwarth (59) W_e (total ir activity ): ", sqrt(H59))

    omega = sqrt(H59 / H58)
    println("Hellwarth (61) Omega (freq): ",omega)

	return(omega)
end

# More complex scheme, involving thermodynamic Beta
# Hellwarth(50), RHS
"""
	HellwarthAScheme(LO,T=295)

    Multiple phonon mode reduction to a single effective frequency.
	Temperature dependent, defaults to T=295 K.

    Follows Hellwarth et al. 1999 PRB 'A' scheme, Eqn 50 RHS.

	UNTESTED AND UNCERTAIN CODE.
"""
function HellwarthAScheme(LO; T=295)
    println("Hellwarth A scheme...T=$T K")
    omega=LO[:,1] # No unit conversion
    #omega=LO[:,1].*2*pi*1E12 # THz --> Hz
    #omega=LO[:,1].*2*pi*0.02998*1E12 # cm^-1 --> Hz
    β=(omega.*ħ)/(kB*T) #assuming units SI 
    H50 = sum( ((LO[:,2].^2).*coth.(β))./omega )
    println("Hellwarth (50) summation: ",H50)

    H51= sum( LO[:,2].^2 ) # sum of total ir activity squarred
    println("Hellwarth (51) summation (total ir activity ^2): ", H51)
    println("Hellwarth (51) W_e (total ir activity ): ", sqrt(H51))

    # OK; so this is deriving Omega / coth(Beta/2)
    omegacoth=H51/H50
    println("omegacoth: ",omegacoth)

    solnfreq=0.0 #required for Julia 0.5 so it realises variable still required for return.
    # Very primitive manner to decouple Omega from both sides of the eqn.
	# Should really rewrite as a bisection (at least!)
    maxfreq=maximum(omega)
    for freq in 0.1:maxfreq/2500:maxfreq
        pseudo_omega=omegacoth/coth(freq * ħ/(2*kB*T))
        if freq>pseudo_omega
            println("freq: $freq pseudo-omega: $pseudo_omega")
            println("freq: ",freq/(2*pi*0.02998*1E12), " pseudo-omega: ",pseudo_omega/(2*pi*0.02998*1E12))
            solnfreq=freq
            break
        end
    end
	return solnfreq
end

