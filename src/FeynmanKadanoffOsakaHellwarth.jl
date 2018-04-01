# FeynmanKadanoffOsakaHellwarth.jl
# Codes by Jarvist Moore Frost, 2017
# Calculate Polaron Mobility - by a Osaka/Hellwarth variational solution to the Feynman model

# These codes were developed with Julia 0.5.0 - Julia 0.6.2, and require the Optim and Plots packages.

export Polaron # Type to hold the data
export feynmanalpha, feynmanvw, F, polaronmobility, savepolaron, plotpolaron
export HellwarthBScheme, HellwarthAScheme
export ImX

##### load in library routines... #####
# one-dimensional numerical integration in Julia using adaptive Gauss-Kronrod quadrature
import QuadGK.quadgk

# Using the powerful Julia Optim package to optimise the variational parameters
using Optim

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space

"""
feynmanalpha(ε_Inf,ε_S,freq,m_eff)

    Calculates the Frohlich alpha parameter, for a given dielectric constant,
    frequency (f) of phonon in Hertz, and effective mass (in units of the
    bare electron mass).

    See Feynman 1955:
    http://dx.doi.org/10.1103/PhysRev.97.660

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

# Hellwarth et al. 1999 PRB - Part IV; T-dep of the Feynman variation parameter

# Originally a Friday afternoon of hacking to try and implement the T-dep electron-phonon coupling from the above PRB
# Which was unusually successful! And more or less reproduced Table III

# Define Osaka's free-energies (Hellwarth1999 version) as Julia functions
# Equation numbers follow above Hellwarth et al. 1999 PRB
# 62b
A(v,w,β)=3/β*( log(v/w) - 1/2*log(2*π*β) - log(sinh(v*β/2)/sinh(w*β/2)))
# 62d
Y(x,v,β)=1/(1-exp(-v*β))*(1+exp(-v*β)-exp(-v*x)-exp(v*(x-β)))
# 62c integrand
f(x,v,w,β)=(exp(β-x)+exp(x))/(w^2*x*(1-x/β)+Y(x,v,β)*(v^2-w^2)/v)^(1/2)
# 62c
B(v,w,β,α) = α*v/(sqrt(π)*(exp(β)-1)) * quadgk(x->f(x,v,w,β),0,β/2)[1]
# 62e
C(v,w,β)=3/4*(v^2-w^2)/v * (coth(v*β/2)-2/(v*β))
# 62a
F(v,w,β,α)=-(A(v,w,β)+B(v,w,β,α)+C(v,w,β)) 

# Can now evaluate, e.g.
# F(v,w,β,α)=F(7.2,6.5,1.0,1.0)
# BUT - this is just the objective function! (The temperature-dependent free-energy.) Not the optimised parameters.
# Also there's a scary numeric integration (quadgk) buried within...

# In Julia we have 'Multiple dispatch', so let's just construct the Feynman (athermal) 
# energy as the same signature, but without the thermodynamic beta

# Integrand of (31) in Feynman I (Feynman 1955, Physical Review, "Slow electrons...")
fF(τ,v,w)=(w^2 * τ + (v^2-w^2)/v*(1-exp(-v*τ)))^-0.5 * exp(-τ)
# (31) in Feynman I
AF(v,w,α)=π^(-0.5) * α*v * quadgk(τ->fF(τ,v,w),0,Inf)[1]
# (33) in Feynman I
F(v,w,α)=(3/(4*v))*(v-w)^2-AF(v,w,α)

# Let's wrap this in a simple function
"""
    feynmanvw(α)

    Calculate v and w variational polaron parameters (Feynman original athermal action), 
    for the supplied α Frohlich coupling.
	Returns v,w. 
"""
function feynmanvw(α)
    initial=[7.0,6.0]
    # Main use of these bounds is stopping v or w going negative, at which you get a NaN error as you are evaluating log(-ve Real)
    lower=[0.1,0.1]
    upper=[1000.0,1000.0]

    myf(x) = F(x[1],x[2],α) # Wraps the function so just the two variational params are exposed

    res=optimize(OnceDifferentiable(myf, initial; autodiff = :forward), initial, lower, upper, Fminbox(); optimizer=BFGS)
    # allow_f_increases=true,  - increases stability of algorith, but makes it more likely to crash as it steps outside Fminbox(), le sigh.
    
    v,w=Optim.minimizer(res)
    
    return v,w
end

# Structure to store data of polaron solution + other parameters, for each temperature
struct Polaron
    T
    # Mobilities 
    Kμ; Hμ; FHIPμ
    # Spring constant and renormalised (phonon-drag) mass
    k; M
    # Osaka free energy components (A,B,C) and total (F). See Hellwarth et al. 1999 PRB Part IV
    A; B; C; F
    # Relaxation time from Kadanoff Boltzmann transport equation
    Tau
    # Raw variational parameters
    v; w
    # Reduced thermodynamic beta
    βred
    # Feynman polaron radius (Schultz), in SI units. Then also the small-alpha asymptotic approx
    rfsi; rfsmallalpha
    # Setup of simulation. These parameters are sent to the function. 
    # Alpha = Frohlich alpha
    α
    # Ban effective mass
    mb
    # Effective dielectric frequency
    ω
end
Polaron()=Polaron([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]) # structure initialisation 

"""
    polaronmobility(Trange, 
                    ε_Inf, ε_S, freq, effectivemass; 
                    verbose::Bool=false)

    Solves the Feynman polaron problem variationally with finite temperature
    Osaka energies.  From the resulting v, and w parameters, calculates polaron
    structure (wave function size, etc.).  Uses FHIP, Kadanoff (Boltzmann
    relaxation time) and Hellwarth direct contour integration to predict
    a temperature-dependent mobility for the material system. 
    Input is a temperature range (e.g. 10:50:1000), 
    reduced dielectric constants (e.g. 5, 20),
    characteristic dielectric phonon frequency (e.g. 2.25E12) - units Hertz
    bare-band effective-mass (e.g. 012) - units electron mass.

    Returns a structure of type Polaron, containing arrays of useful
    information.  Also prints a lot of information to the standard out - which
    may be more useful if you're just inquiring as to a particular data point,
    rather than plotting a temperature-dependent parameter. 

    As an example, to calculate the electron polaron in MAPI at 300 K:
    polaronmobility(300, 4.5, 24.1, 2.25E12, 0.12)
"""
function polaronmobility(Trange, ε_Inf, ε_S, freq, effectivemass; verbose::Bool=false)
    println("\n\nPolaron mobility for system ε_Inf=$ε_Inf, ε_S=$ε_S, freq=$freq, 
                 effectivemass=$effectivemass; with Trange $Trange ...")

    # Internally we we use 'mb' for the 'band mass' in SI units, of the effecitve-mass of the electron
    mb=effectivemass*MassElectron 
    ω = (2*pi)*freq # angular-frequency, of the phonon mode 
        
    α=feynmanalpha(ε_Inf, ε_S,  freq,    effectivemass)
    #α=2.395939683378253 # Hard coded; from MAPI params, 4.5, 24.1, 2.25THz, 0.12me

    @printf("Polaron mobility input parameters: ε_Inf=%f ε_S=%f freq=%g α=%f \n",ε_Inf, ε_S, freq, α )
    @printf("Derived params in SI: ω =%g mb=%g \n",ω ,mb)

    # Initial v,w to use
    initial=[7.1,6.5]
    # Main use of these bounds is stopping v or w going negative, at which you get a NaN error as you are evaluating log(-ve Real)
    lower=[0.1,0.1]
    upper=[100.0,100.0]

    # Empty struct for storing data 
    # A slightly better way of doing this ф_ф ...
    p=Polaron()
    # populate data structure with (athermal) parameters supplied...
    append!(p.α,α) # appending so as not to mess with type immutability
    append!(p.mb,mb)
    append!(p.ω,ω)

    # We define βred as the subsuming the energy of the phonon; i.e. kbT c.f. ħω
    for T in Trange 
        β=1/(kB*T)
        βred=ħ*ω*β
        append!(p.βred,βred)
        @printf("T: %f β: %.3g βred: %.3g ħω  = %.3g meV\t",T,β,βred, 1E3* ħ*ω  / q)
        myf(x) = F(x[1],x[2],βred,α) # Wraps the function so just the two variational params are exposed
        # Now updated to use Optim > 0.7.8 call signature (Julia >0.6 only) 
        res=optimize(OnceDifferentiable(myf, initial; autodiff = :forward), initial, lower, upper, Fminbox(); 
            optimizer=BFGS)
        # allow_f_increases=true,  
        #  - increases stability of algorithm, particular when calling it for
        #  a single, extreme value of temperature, but makes it more likely to
        #  crash as it steps outside Fminbox(), le sigh.

        minimum=Optim.minimizer(res)
        print("\tConverged? : ",Optim.converged(res) ) # All came out as 'true'
        if verbose
            println()
            show(res)
        end
        
        v=minimum[1] # unpack these from minimised results
        w=minimum[2]
            
        @printf("\n VariationalParams v= %.2f  w= %.2f",v,w)

        # From 1962 Feynman, definition of v and w in terms of the coupled Mass and spring-constant
        # See Page 1007, just after equation (18)
        # Units of M appear to be 'electron masses'
        # Unsure of units for k, spring coupling constant
        k=(v^2-w^2)
        M=(v^2-w^2)/w^2
        append!(p.k,k)
        append!(p.M,M)

        @printf("\t||\t M=%f k=%f\t",M,k)
       
        @printf("\n Polaron frequency (SI) v=  %.2g Hz \t w=  %.2g Hz\t",v*ω /(2*pi),w*ω /(2*pi))

		# (46) in Feynman1955
		meSmallAlpha(α )=α /6 + 0.025*α ^2
		# (47) In Feynman1955
		meLargeAlpha(α )=16*α ^4 / (81*π ^4)
		#meLargeAlpha(α )=202*(α /10)^4
		@printf("\n Feynman1955(46,47): meSmallAlpha(α)= %.3f meLargeAlpha(α)= %.3f",meSmallAlpha(α),meLargeAlpha(α))
        @printf("\n Feynman1962: Approximate ~ Large alpha limit, v/w = %.2f  =~approx~= alpha^2 = %.2f ",v/w,α ^2)

        # POLARON SIZE
        @printf("\n POLARON SIZE (rf), following Schultz1959. (s.d. of Gaussian polaron ψ )")
        # Schultz1959 - rather nicely he actually specifies everything down into units!
        # just before (2.4) in Schultz1959
        mu=((v^2-w^2)/v^2)
        # (2.4)
        rf=sqrt(3/(2*mu*v))
        # (2.4) SI scaling inferred from units in (2.5a) and Table II
        rfsi=rf*sqrt(2*me*ω )
        @printf("\n\t Schultz1959(2.4): rf= %g (int units) = %g m [SI]",rf,rfsi )
        append!(p.rfsi,rfsi)

        scale=sqrt(2*mb*ω) # Note we're using mb; 
        #band effective-mass in SI units (i.e. meff*melectron)
		
        rfa=(3/(0.44*α ))^0.5 # As given in Schultz1959(2.5a), but that 0.44 is actually 4/9
        @printf("\n\t Schultz1959(2.5a) with 0.44: Feynman α→0 expansion: rfa= %g (int units) = %g m [SI]",rfa,scale*rfa )
        rfa=(3/((4/9)*α ))^0.5 # Rederived from Feynman1955, 8-8-2017; Yellow 2017.B Notebook pp.33-34	
        @printf("\n\t Schultz1959(2.5a) with 4/9 re-derivation: Feynman α→0 expansion: rfa= %g (int units) = %g m [SI]",rfa,scale*rfa )
        append!(p.rfsmallalpha,scale*rfa)

        rfb=3*(pi/2)^0.5 * α 
        @printf("\n\t Schultz1959(2.5b): Feynman α→∞ expansion: rf= %g (int units) = %g m [SI]",rfb,scale*rfb )
        
        # Schultz1959 - Between (5.7) and (5.8) - resonance of Feynman SHM system
        phononfreq=sqrt(k/M)
        @printf("\n\t Schultz1959(5.7-5.8): fixed-e: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            phononfreq,phononfreq*ω /(2*pi), phononfreq*hbar*ω *1000/q)
 
        phononfreq=sqrt(k/mu) # reduced mass
        @printf("\n\t Schultz1959(5.7-5.8): reducd mass: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            phononfreq,phononfreq*ω /(2*pi), phononfreq*hbar*ω *1000/q)
        
        @printf("\n\t Schultz1959: electronfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            sqrt(k/1),sqrt(k/1)*ω /(2*pi), sqrt(k/1)*hbar*ω *1000/q)
        @printf("\n\t Schultz1959: combinedfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            sqrt(k/(1+M)),sqrt(k/(1+M))*ω /(2*pi), sqrt(k/(1+M))*hbar*ω *1000/q)

        # Devreese1972: 10.1103/PhysRevB.5.2367
        # p.2371, RHS.
        @printf("\n Devreese1972: (Large Alpha) Franck-Condon frequency = %.2f", 4*α ^2/(9*pi))

        # F(v,w,β,α)=-(A(v,w,β)+B(v,w,β,α)+C(v,w,β)) #(62a) - Hellwarth 1999
        @printf("\n Polaron Free Energy: A= %f B= %f C= %f F= %f",A(v,w,βred),B(v,w,βred,α),C(v,w,βred),F(v,w,βred,α))
        @printf("\t = %f meV",1000.0 * F(v,w,βred,α) * ħ*ω  / q)
        append!(p.A,A(v,w,βred))
        append!(p.B,B(v,w,βred,α))
        append!(p.C,C(v,w,βred))
        append!(p.F,F(v,w,βred,α))
        

        # FHIP 
        #    - low-T mobility, final result of Feynman1962
        # [1.60] in Devreese2016 page 36; 6th Edition of Frohlich polaron notes (ArXiv)
        # I believe here β is in SI (expanded) units
        @printf("\nPolaron Mobility theories:")
        μ=(w/v)^3 * (3*q)/(4*mb*ħ*ω^2*α*β) * exp(ħ*ω*β) * exp((v^2-w^2)/(w^2*v))
        @printf("\n\tμ(FHIP)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)
        append!(p.T,T)
        append!(p.FHIPμ,μ*100^2)
        

        # Kadanoff
        #     - low-T mobility, constructed around Boltzmann equation. 
        #     - Adds factor of 3/(2*beta) c.f. FHIP, correcting phonon emission behaviour
        # [1.61] in Devreese2016 - Kadanoff's Boltzmann eqn derived mob
        μ=(w/v)^3 * (q)/(2*mb*ω*α) * exp(ħ*ω*β) * exp((v^2-w^2)/(w^2*v))
        @printf("\n\tμ(Kadanoff,via Devreese2016)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)    

        append!(p.Kμ,μ*100^2)

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
        #It suggests that Kadanoff used the wrong identy for Nbar in 23(b) for
        #the Gamma0 function, and should have used a version with the -1 to
        #account for Bose / phonon statistics

        #myv=sqrt(k*(1+M)/M) # cross-check maths between different papers
        #@printf("\nv: %f myv: %f\n",v,myv)

        # Between 23 and 24 in Kadanoff 1963, for small momenta skip intergration --> Gamma0
        Gamma0=2*α * Nbar * (M+1)^(1/2) * exp(-M/v)
        Gamma0*=ω  #* (ω *hbar) # Kadanoff 1963 uses hbar=omega=mb=1 units
            # Factor of omega to get it as a rate relative to phonon frequency
            # Factor of omega*hbar to get it as a rate per energy window
        μ=q/( mb*(M+1) * Gamma0 ) #(25) Kadanoff 1963, with SI effective mass
        @printf("\n\tμ(Kadanoff1963 [Eqn. 25]) = %f m^2/Vs \t = %.2f cm^2/Vs",μ,μ*100^2)
        @printf("\n\tGamma0 = %g rad/s = %g /s \n\tTau=1/Gamma0 = %g = %f ps",
            Gamma0, Gamma0/(2*pi), 2*pi/Gamma0, 2*pi*1E12/Gamma0)
        Eloss=hbar*ω * Gamma0/(2*pi) # Simply Energy * Rate
        @printf("\n\tEnergy Loss = %g J/s = %g meV/ps",Eloss,Eloss * 1E3 / (q*1E12) ) 
        append!(p.Tau, 2*pi*1E12/Gamma0) # Boosted into ps ?


        # Hellwarth1999 - directly do contour integration in Feynman1962, for finite temperature DC mobility

        # Hellwarth1999 Eqn (2) and (1) - These are going back to the general (pre low-T limit) formulas in Feynman1962.
        # to evaluate these, you need to do the explicit contour integration to get the polaron self-energy
        R=(v^2-w^2)/(w^2*v) # inline, page 300 just after Eqn (2)
        
        #b=R*βred/sinh(b*βred*v/2) # This self-references b! What on Earth?
        # OK! I now understand that there is a typo in Hellwarth1999 and
        # Biaggio1997. They've introduced a spurious b on the R.H.S. compared to
        # the original, Feynman1962:
        b=R*βred/sinh(βred*v/2) # Feynman1962 version; page 1010, Eqn (47b)

        a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2))
        k(u,a,b,v) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u) # integrand in (2)
        K=quadgk(u->k(u,a,b,v),0,Inf)[1] # numerical quadrature integration of (2)
        
        #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS= α/(3*sqrt(π)) * βred^(5/2) / sinh(βred/2) * (v^3/w^3) * K
        μ=RHS^-1 * (q)/(ω*mb)
        @printf("\n\tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)    
        append!(p.Hμ,μ*100^2)
        
        # Hellwarth1999/Biaggio1997, b=0 version... 'Setting b=0 makes less than 0.1% error'
        # So let's test this
        R=(v^2-w^2)/(w^2*v) # inline, page 300 just after Eqn (2)
        b=0 
        a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2))
        k(u,a,b,v) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u) # integrand in (2)
        K=quadgk(u->k(u,a,b,v),0,Inf)[1] # numerical quadrature integration of (2)

        #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS= α/(3*sqrt(π)) * βred^(5/2) / sinh(βred/2) * (v^3/w^3) * K
        μ=RHS^-1 * (q)/(ω*mb)
        @printf("\n\tμ(Hellwarth1999,b=0)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)
        @printf("\n\tError due to b=0; %f",(100^2*μ-p.Hμ[length(p.Hμ)])/(100^2*μ))
        #append!(Hμs,μ*100^2)
        
        #ImX(v,w,βred,α,ω,mb)

        # Recycle previous variation results (v,w) as next guess
        initial=[v,w] # Caution! Might cause weird sticking in local minima
        append!(p.v,v)
        append!(p.w,w)
    end

    @printf("\n") # blank line at end of spiel.
    return(p)
end

"""
    savepolaron(fileprefix, p::Polaron)

Saves data from polaron 'p' into file "fileprefix". 
This is a simple space-delimited text file, with each entry a separate temperature, for plotting with Gnuplot or similar. 

Structure of file is written to the header:
# Ts, βreds, Kμs, Hμs, FHIPμs, vs, ws, ks, Ms, As, Bs, Cs, Fs, Taus, rfsis
# 1    2     3    4     5      6   7   8   9  10  11  12  13    14     15
"""
function savepolaron(fileprefix, p::Polaron) 
    println("Saving data to $fileprefix.dat ...")
    f=open("$fileprefix.dat","w")
    
    @printf(f,"# %s\n",fileprefix) # put name / material at header
    @printf(f,"# Params in SI: ω =%g mb=%g \n",p.ω[1] ,p.mb[1])
    @printf(f,"# Alpha parameter: α = %f  \n",p.α[1] )

    @printf(f,"# Ts, βreds, Kμs, Hμs, FHIPμs, vs, ws, ks, Ms, As, Bs, Cs, Fs, Taus, rfsis\n")
    @printf(f,"#  1    2     3    4     5      6   7   8   9  10  11  12  13    14     15\n") # columns for GNUPLOT etc.

    for i in 1:length(p.T)
        @printf(f,"%d %03f %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
        p.T[i], p.βred[i], p.Kμ[i], p.Hμ[i], p.FHIPμ[i], 
        p.v[i], p.w[i],
        p.k[i], p.M[i], p.A[i], p.B[i], p.C[i], p.F[i], 
        p.Tau[i], p.rfsi[i])
    end
    close(f)
end

