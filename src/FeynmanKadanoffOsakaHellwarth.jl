# FeynmanKadanoffOsakaHellwarth.jl
# Codes by Jarvist Moore Frost, 2017
# Calculate Polaron Mobility - by a Osaka/Hellwarth variational solution to the Feynman model
# If you run this, it should construct the model, solve for varying temperature, then produce plots as .pngs in the local directory.
# These codes were developed with Julia 0.5.0, and requires the Optim and Plots packages.

module FeynmanKadanoffOsakaHellwarth
export feynmanalpha, polaronmobility
export HellwarthBScheme, HellwarthAScheme

##### load in library routines... #####
# one-dimensional numerical integration in Julia using adaptive Gauss-Kronrod quadrature
import QuadGK.quadgk

# Plot figures with Plots, which defaults to Pyplot backend
using Plots
#gr()
#default(size=(800,600)) # For the .png file output
# Using the powerful Julia Optim package to optimise the variational parameters
using Optim

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space

# via Feynman 1955
#   http://dx.doi.org/10.1103/PhysRev.97.660
"""
feynmanalpha(ε_Inf,ε_S,freq,m_eff)

    Calculates the Frohlich alpha parameter, for a given dielectric constant,
    frequency (f) of phonon in Hertz, and effective mass (in units of the
    bare electron mass).
"""
function feynmanalpha(ε_Inf,ε_S,freq,m_eff)
    ω=freq*2*pi #frequency to angular velocity
    # Note to self - you've introduced the 4*pi factor into the dielectric constant; 
    # This gives numeric agreement with literature values, but I'm uncertain of the justification.
    # Such a factor also seems to be necessary for calculation of Exciton binding energies (see below).
    # Maybe scientists in the 50s assumed that the Epsilon subsumed the 4*pi factor?
    α=0.5* (1/ε_Inf - 1/ε_S)/(4*pi*ε_0) * (q^2/(hbar*ω)) * (2*me*m_eff*ω/hbar)^0.5
    return (α)
end

#####
# Hellwarth multiple phonon branch reduction

# Most simple scheme
# Hellwarth (58), assuming further typo on LHS, actually should be W_e
"""
	HellwarthBScheme(LO)

	Multiple phonon mode reduction to a single effective frequency.
	Hellwarth et al. 1999 PRB, 'B scheme'.
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
	UNTESTED AND UNCERTAIN CODE.

	Follows Hellwarth et al. 1999 PRB, Eqn 50 RHS.
"""
function HellwarthAScheme(LO,T=295)
    println("Hellwarth A scheme...T=$T K")
    β=LO[:,1].*2*pi*1E12*ħ/(kB*T) #assuming units of THz
    H50 = sum( ((LO[:,2].^2).*coth.(β))./LO[:,1] )
    println("Hellwarth (50) summation: ",H50)

    H51= sum( LO[:,2].^2 ) # sum of total ir activity squarred
    println("Hellwarth (51) summation (total ir activity ^2): ", H51)
    println("Hellwarth (51) W_e (total ir activity ): ", sqrt(H51))

    # OK; so this is deriving Omega / coth(Beta/2)
    omegacoth=H51/H50
    println("omegacoth: ",omegacoth)

    freq=0.0 #required for Julia 0.5 so it realises variable still required for return.
    # Very primitive manner to decouple Omega from both sides of the eqn.
	# Should really rewrite as a bisection (at least!)
    for freq in 0.1:0.1:20
        pseudo_omega=omegacoth*coth(freq * 2*pi*1E12*ħ/(2*kB*T))
        if freq>pseudo_omega
            println("freq: $freq pseudo-omega: $pseudo_omega")
            break
        end
    end
	return(freq)
end



#####
# Set up equations for the polaron free energy, which we will variationally improve upon

# Hellwarth 1999 PRB - Part IV; T-dep of the Feynman variation parameter
# Originally a Friday afternoon of hacking to try and implement the T-dep electron-phonon coupling from the above PRB
# Which was unusually successful! And more or less reproduced Table III

# Define Osaka's free-energies (Hellwarth1999 version) as Julia functions
# Equation numbers follow above Hellwarth 1999 PRB
# 62b
A(v,w,β)=3/β*( log(v/w) - 1/2*log(2*π*β) - log(sinh(v*β/2)/sinh(w*β/2)))
# 62d
Y(x,v,β)=1/(1-exp(-v*β))*(1+exp(-v*β)-exp(-v*x)-exp(v*(x-β)))
# 62c integrand
f(x,v,w,β)=(exp(β-x)+exp(x))/(w^2*x*(1-x/β)+Y(x,v,β)*(v^2-w^2)/v)^(1/2)
# 62c
B(v,w,β,α) = α*v/(sqrt(π)*(exp(β)-1)) * quadgk(x->f(x,v,w,β),0,β/2)[1]
#62e
C(v,w,β)=3/4*(v^2-w^2)/v * (coth(v*β/2)-2/(v*β))

F(v,w,β,α)=-(A(v,w,β)+B(v,w,β,α)+C(v,w,β)) #(62a)

# Can now evaluate, e.g.
# F(v,w,β,α)=F(7.2,6.5,1.0,1.0)
# BUT - this is just the objective function! (The temperature-dependent free-energy.) Not the optimised parameters.
# Also there's a scary numeric integration (quadgk) buried within...

#####
# OK, this was all in the global scope, but has now been put within a function so it can be called for varying parameters
function polaronmobility(fileprefix,Trange, ε_Inf, ε_S,  freq,    effectivemass; figures::Bool=true)
    @printf("Calculating polaron mobility for %s ...\n",fileprefix)

    # Internally we have 'mb' for the 'band mass' in SI units, of the effecitve-mass of the electron
    mb=effectivemass*MassElectron 
    ω = (2*pi)*freq # angular-frequency
        
    α=feynmanalpha(ε_Inf, ε_S,  freq,    effectivemass)
    #α=2.395939683378253 # Hard coded; from MAPI params, 4.5, 24.1, 2.25THz, 0.12me

    @printf("Polaron mobility input parameters: ε_Inf=%f ε_S=%f freq=%g α=%f \n",ε_Inf, ε_S, freq, α )
    @printf("Derived params in SI: ω =%g mb=%g \n",ω ,mb)

    # Dump to log file
    f=open("$fileprefix.dat","w")
    @printf(f,"#Model parameters: ε_Inf=%f ε_S=%f freq=%g mb=%f \n",ε_Inf, ε_S, freq, effectivemass)
    @printf(f,"#Params in SI: ω =%g mb=%g \n",ω ,mb)
    @printf(f,"#Alpha parameter: α = %f  \n",α )
    close(f)


    # Initial v,w to use
    initial=[7.1,6.5]
    # Main use of these bounds is stopping v or w going negative, at which you get a NaN error as you are evaluating log(-ve Real)
    lower=[0.1,0.1]
    upper=[100.0,100.0]

    # Empty arrays for storing data 
    # Surely some better way of doing this ф_ф 
    Ts=[]
    Kμs=[]
    Hμs=[]
    FHIPμs=[]
    ks=[]
    Ms=[]
    As=[]
    Bs=[]
    Cs=[]
    Fs=[]
    Taus=[]
    vs=[]
    ws=[]
    βreds=[]
    rfsis=[]

    # We define βred as the subsuming the energy of the phonon; i.e. kbT c.f. ħω
    for T in Trange #10:10:4000
        β=1/(kB*T)
        βred=ħ*ω*β
        append!(βreds,βred)
        @printf("T: %f β: %.2g βred: %.2g ħω  = %.2g meV\t",T,β,βred, 1000.0*ħ*ω  / q)
        myf(x) = F(x[1],x[2],βred,α) # Wraps the function so just the two variational params are exposed
        res=optimize(DifferentiableFunction(myf), initial, lower, upper, Fminbox(); 
                optimizer=BFGS, optimizer_o=(Optim.Options(autodiff=true)))
        minimum=Optim.minimizer(res)
        #show(Optim.converged(res)) # All came out as 'true'
        
        v=minimum[1]
        w=minimum[2]
            
        @printf("\n VariationalParams v= %.2f = %.2g Hz \t w= %.2f = %.2g Hz\t",v,v*ω /(2*pi),w,w*ω /(2*pi))
        
        # From 1962 Feynman, definition of v and w in terms of the coupled Mass and spring-constant
        # See Page 1007, just after equation (18)
        # Units of M appear to be 'electron masses'
        # Unsure of units for k, spring coupling constant
        k=(v^2-w^2)
        M=(v^2-w^2)/w^2
        append!(ks,k)
        append!(Ms,M)

        @printf(" M=%f k=%f\t",M,k)
        
		# (46) in Feynman1955
		meSmallAlpha(α )=α /6 + 0.025*α ^2
		# (47) In Feynman1955
		meLargeAlpha(α )=16*α ^4 / (81*π ^4)
		#meLargeAlpha(α )=202*(α /10)^4
		println("\n Feynman1955(46,47): meSmallAlpha(α)=",meSmallAlpha(α)," meLargeAlpha(α)=",meLargeAlpha(α))

        @printf("\n Feynman1962: Large alpha, v/w = %.2f  =~approx~= alpha^2 = %.2f ",v/w,α ^2)

        # Schultz1959 - rather nicely he actually specifies everything down into units!
        # just before (2.4) in Schultz1959
        mu=((v^2-w^2)/v^2)
        # (2.4)
        rf=sqrt(3/(2*mu*v))
        # (2.4) SI scaling inferred from units in (2.5a) and Table II
        rfsi=rf*sqrt(2*me*ω )
        @printf("\n Schultz1959(2.4): rf= %g (int units) = %g m [SI]",rf,rfsi )
        append!(rfsis,rfsi)

		rf=(3/(0.44*α ))^0.5
        rfsi=rf*sqrt(2*me*ω )
        @printf("\n Schultz1959(2.5a), Feynman alpha->0 expansion: rf= %g (int units) = %g m [SI]",rf,rfsi )
		rf=3*(pi/2)^0.5 * α 
        rfsi=rf*sqrt(2*me*ω )
        @printf("\n Schultz1959(2.5a), Feynman alpha>-Inf expansion: rf= %g (int units) = %g m [SI]",rf,rfsi )
        
        # Schultz1959 - Between (5.7) and (5.8) - resonance of Feynman SHM system
        phononfreq=sqrt(k/M)
        @printf("\n Schultz1959 (5.7-5.8) fixed-e: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            phononfreq,phononfreq*ω /(2*pi), phononfreq*hbar*ω *1000/q)
 
        phononfreq=sqrt(k/mu) # reduced mass
        @printf("\n Schultz1959: (5.7-5.8) reducd mass: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            phononfreq,phononfreq*ω /(2*pi), phononfreq*hbar*ω *1000/q)
        
        @printf("\n Schultz1959: electronfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            sqrt(k/1),sqrt(k/1)*ω /(2*pi), sqrt(k/1)*hbar*ω *1000/q)
        @printf("\n Schultz1959: combinedfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            sqrt(k/(1+M)),sqrt(k/(1+M))*ω /(2*pi), sqrt(k/(1+M))*hbar*ω *1000/q)

        # Devreese1972: 10.1103/PhysRevB.5.2367
        # p.2371, RHS.
        @printf("\n Devreese1972: (Large Alpha) Franck-Condon frequency = %.2f", 4*α ^2/(9*pi))

        # F(v,w,β,α)=-(A(v,w,β)+B(v,w,β,α)+C(v,w,β)) #(62a) - Hellwarth 1999
        @printf("\n Polaron Free Energy: A= %f B= %f C= %f F= %f",A(v,w,βred),B(v,w,βred,α),C(v,w,βred),F(v,w,βred,α))
        @printf("\t = %f meV",1000.0 * F(v,w,βred,α) * ħ*ω  / q)
        append!(As,A(v,w,βred))
        append!(Bs,B(v,w,βred,α))
        append!(Cs,C(v,w,βred))
        append!(Fs,F(v,w,βred,α))
        

        # FHIP 
        #    - low-T mobility, final result of Feynman1962
        # [1.60] in Devreese2016 page 36; 6th Edition of Frohlich polaron notes (ArXiv)
        # I believe here β is in SI (expanded) units
        μ=(w/v)^3 * (3*q)/(4*mb*ħ*ω^2*α*β) * exp(ħ*ω*β) * exp((v^2-w^2)/(w^2*v))
        @printf("\n\tμ(FHIP)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)
        append!(Ts,T)
        append!(FHIPμs,μ*100^2)
        

        # Kadanoff
        #     - low-T mobility, constructed around Boltzmann equation. 
        #     - Adds factor of 3/(2*beta) c.f. FHIP, correcting phonon emission behaviour
        # [1.61] in Devreese2016 - Kadanoff's Boltzmann eqn derived mob
        μ=(w/v)^3 * (q)/(2*mb*ω*α) * exp(ħ*ω*β) * exp((v^2-w^2)/(w^2*v))
        @printf("\n\tμ(Kadanoff)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)    

        append!(Kμs,μ*100^2)

        ######
        # OK, now deep-diving into Kadanoff1963 itself to extract Boltzmann equation components
        # Particularly right-hand-side of page 1367
        #
        # From 1963 Kadanoff, (9), define eqm. number of phonons (just from T and phonon omega)
        Nbar=(exp(βred) - 1)^-1
        @printf("\n\t Eqm. Phonon. pop. Nbar: %f ",Nbar)
        @printf("\n\texp(Bred): %f exp(-Bred): %f exp(Bred)-1: %f",exp(βred),exp(-βred),exp(βred)-1)
        Nbar=exp(-βred) 
        #Note - this is only way to get Kadanoff1963 to be self-consistent with
        #FHIP, and later statements (Devreese) of the Kadanoff mobility. 
        #It suggests that Kadanoff used the wrong identy for Nbar in 23(b) for the
        #Gamma0 function, and should have used a version -1 to account for phonon
        #statistics

        #myv=sqrt(k*(1+M)/M) # cross-check maths between different papers
        #@printf("\nv: %f myv: %f\n",v,myv)

        # Between 23 and 24 in Kadanoff 1963, for small momenta skip intergration --> Gamma0
        Gamma0=2*α * Nbar * (M+1)^(1/2) * exp(-M/v)
        Gamma0*=ω  #* (ω *hbar) # Kadanoff 1963 uses hbar=omega=mb=1 units
            # Factor of omega to get it as a rate relative to phonon frequency
            # Factor of omega*hbar to get it as a rate per energy window
        μ=q/( mb*(M+1) * Gamma0 ) #(25) Kadanoff 1963, with SI effective mass
        @printf("\n\tμ(Kadanoff [Eqn. 25]) = %f m^2/Vs \t = %.2f cm^2/Vs",μ,μ*100^2)
        @printf("\n\tGamma0 = %g rad/s = %g /s Tau=1/Gamma0 = %g = %f ps",
            Gamma0, Gamma0/(2*pi), 2*pi/Gamma0, 2*pi*1E12/Gamma0)
        append!(Taus, 2*pi*1E12/Gamma0)

        # Hellwarth1999 Eqn (2) and (1) - These are going back to the general (pre low-T limit) formulas in Feynman1962.
        # to evaluate these, you need to do the explicit contour integration to get the polaron self-energy
        R=(v^2-w^2)/(w^2*v) # inline, page 300 just after Eqn (2)
        #b=R*βred/sinh(b*βred*v/2) # This self-references b! What on Earth?
        # OK! I now understand that there is a typo in Hellwarth1999 and
        # Biaggio1997. They've introduced a spurious b on the R.H.S. compared to
        # the original, Feynman1962...
        b=R*βred/sinh(βred*v/2) # Feynman1962 version; page 1010, Eqn (47b)
        a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2))
        k(u,a,b,v) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u) # integrand in (2)
        K=quadgk(u->k(u,a,b,v),0,Inf)[1] # numerical quadrature integration of (2)
        
        #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS= α/(3*sqrt(π)) * βred^(5/2) / sinh(βred/2) * (v^3/w^3) * K
        μ=RHS^-1 * (q)/(ω*mb)
        @printf("\n\tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)    
        append!(Hμs,μ*100^2)
        
        #Hellwarth1999/Biaggio1997, b=0 version... 'Setting b=0 makes less than 0.1% error'
        R=(v^2-w^2)/(w^2*v) # inline, page 300 just after Eqn (2)
        b=0 
        a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2))
        k(u,a,b,v) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u) # integrand in (2)
        K=quadgk(u->k(u,a,b,v),0,Inf)[1] # numerical quadrature integration of (2)

        #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS= α/(3*sqrt(π)) * βred^(5/2) / sinh(βred/2) * (v^3/w^3) * K
        μ=RHS^-1 * (q)/(ω*mb)
        @printf("\n\tμ(Hellwarth1999,b=0)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)
        @printf("\n\tError due to b=0; %f",(100^2*μ-Hμs[length(Hμs)])/(100^2*μ))
        #append!(Hμs,μ*100^2)
        
        
        @printf("\n\n")
        
        # Recycle previous variation results (v,w) as next guess
        initial=[v,w] # Caution! Might cause weird sticking in local minima
        append!(vs,v)
        append!(ws,w)
    end

    println("Saving data to $fileprefix.dat ...")
    f=open("$fileprefix.dat","a")
    @printf(f,"# %s \n# Ts, βreds, Kμs, Hμs, FHIPμs, vs, ws, ks, Ms, As, Bs, Cs, Fs, Taus, rfsis\n",fileprefix)
    @printf(f,"#  1    2     3    4     5     6   7    8  9   10  11  12  13  14 15\n") # columns for GNUPLOT etc.
    for i in 1:length(Ts)
        @printf(f,"%d %03f %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
        Ts[i], βreds[i], Kμs[i], Hμs[i], FHIPμs[i], 
        vs[i], ws[i],
        ks[i], Ms[i], As[i], Bs[i], Cs[i], Fs[i], 
        Taus[i], rfsis[i])
    end
    close(f)

    if figures # only if asked to plot... 
    println("OK - everything calculated and stored. Now plotting..")

    #####
    ## Mass vs. Temperature plot
    plot(Ts,Ms,label="Phonon effective-mass",markersize=3,marker=:rect,xlab="Temperature (K)",ylab="Phonon effective-mass",ylim=(0,1.2))

    savefig("$fileprefix-mass.png")
    savefig("$fileprefix-mass.eps")

    #####
    ## Relaxationtime vs. Temperature plot
    plot(Ts,Taus,label="Kadanoff relaxation time (ps)",markersize=3,marker=:rect,xlab="Temperature (K)",ylab="Relaxation time (ps)",ylim=(0,1.2))

    savefig("$fileprefix-tau.png")
    savefig("$fileprefix-tau.eps")

    ## Mass + relaxation time vs. Temperature plot
    plot(Ts,Ms,label="Phonon effective-mass (m\$_b\$)",markersize=3,marker=:rect,
        xlab="Temperature (K)",ylab="Effective-mass / relaxation time",ylim=(0,1.2))
    plot!(Ts,Taus,label="Kadanoff relaxation time (ps)",markersize=3,marker=:diamond,
        xlab="Temperature (K)",ylab="Relaxation time (ps)",ylim=(0,1.2))

    savefig("$fileprefix-mass-tau.png")
    savefig("$fileprefix-mass-tau.eps")

    #####
    ## Spring Constants vs. Temperature plot
    plot(Ts,ks,label="Polaron spring-constant",markersize=3, marker=:uptriangle, xlab="Temperature (K)",ylab="Spring-constant",)

    savefig("$fileprefix-spring.png")
    savefig("$fileprefix-spring.eps")

    #####
    ## Variation Energy vs. Temperature plots
    plot( Ts,As,label="A",markersize=3,marker=:downtriangle, xlab="Temperature (K)",ylab="Polaron free-energy")
    plot!(Ts,Bs,label="B",markersize=3,marker=:diamond)
    plot!(Ts,Cs,label="C",markersize=3,marker=:uptriangle)
    plot!(Ts,Fs,label="F",markersize=3,marker=:rect)
    #plot!(Ts,Fs,label="F=-(A+B+C)",markersize=3,marker=:rect)

    savefig("$fileprefix-variational.png")
    savefig("$fileprefix-variational.eps")

    #####
    ## Polaron radius vs. Temperature
    plot(Ts,rfsis.*10^10,label="Schultz Feynman radius",xlab="Temperature (K)",ylab="Polaron Radius (Angstrom)")
    savefig("$fileprefix-radius.png")
    savefig("$fileprefix-radius.eps")

    #####
    ## Calculated mobility comparison plot
    plot(Ts,Kμs,label="Kadanoff",markersize=3,marker=:rect,xlab="Temperature (K)",ylab="Mobility (cm\$^2\$/Vs)",ylims=(0,1000))
    plot!(Ts,FHIPμs,label="FHIP",markersize=3,marker=:diamond)
    plot!(Ts,Hμs,label="Hellwarth1999",markersize=3,marker=:uptriangle)

    savefig("$fileprefix-mobility-calculated.png")
    savefig("$fileprefix-mobility-calculated.eps")
    end

    return(Ts,Kμs, Hμs, FHIPμs, ks, Ms, As, Bs, Cs, Fs, Taus,rfsis)
end

end

