# FeynmanKadanoffOsakaHellwarth.jl
# Codes by Jarvist Moore Frost, 2017
# Calculate Polaron Mobility - by a Osaka/Hellwarth variational solution to the Feynman model
# If you run this, it should construct the model, solve for varying temperature, then produce plots as .pngs in the local directory.
# These codes were developed with Julia 0.5.0, and requires the Optim and Plots packages.

module FeynmanKadanoffOsakaHellwarth
export feynmanalpha, polaronmobility, savepolaron, plotpolaron
export HellwarthBScheme, HellwarthAScheme

##### load in library routines... #####
# one-dimensional numerical integration in Julia using adaptive Gauss-Kronrod quadrature
import QuadGK.quadgk

# Plot figures with Plots, which defaults to Pyplot backend
using Plots
#pyplot()
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

# Struct to store data
struct Polaron
    T
    Kμ
    Hμ
    FHIPμ
    k
    M
    A
    B
    C
    F
    Tau
    v
    w
    βred
    rfsi
end
Polaron()=Polaron([],[],[],[],[],[],[],[],[],[],[],[],[],[],[])


#####
# OK, this was all in the global scope, but has now been put within a function so it can be called for varying parameters
function polaronmobility(fileprefix,Trange, ε_Inf, ε_S,  freq,    effectivemass; figures::Bool=true, verbose::Bool=false)
    println("Polaron mobility called for $fileprefix with Trange $Trange ...")
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

    # Empty struct for storing data 
    # A slightly better way of doing this ф_ф ...
    p=Polaron()

    # We define βred as the subsuming the energy of the phonon; i.e. kbT c.f. ħω
    for T in Trange 
        β=1/(kB*T)
        βred=ħ*ω*β
        append!(p.βred,βred)
        @printf("T: %f β: %.2g βred: %.2g ħω  = %.2g meV\t",T,β,βred, 1000.0*ħ*ω  / q)
        myf(x) = F(x[1],x[2],βred,α) # Wraps the function so just the two variational params are exposed
        # OK; this was as working on Julia 0.5; before the great Optim update
        # For Julia 0.6; pin the Optim package to this old interface: Pkg.pin("Optim",v"0.7.8")
        res=optimize(DifferentiableFunction(myf), initial, lower, upper, Fminbox(); 
            optimizer=BFGS, 
            optimizer_o=(Optim.Options(allow_f_increases=true, autodiff=true)))


        minimum=Optim.minimizer(res)
        print("\tConverged? : ",Optim.converged(res) ) # All came out as 'true'
        if verbose
            println()
            show(res)
        end
        
        v=minimum[1]
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
        @printf("\n Schultz1959(2.4): rf= %g (int units) = %g m [SI]",rf,rfsi )
        append!(p.rfsi,rfsi)

        scale=sqrt(2*mb*ω) # Note we're using mb; 
        #band effective-mass in SI units (i.e. meff*melectron)
		
        rfa=(3/(0.44*α ))^0.5 # As given in Schultz1959(2.5a), but that 0.44 is actually 4/9
        @printf("\n Schultz1959(2.5a) with 0.44: Feynman α→0 expansion: rfa= %g (int units) = %g m [SI]",rfa,scale*rfa )
        rfa=(3/((4/9)*α ))^0.5 # Rederived from Feynman1955, 8-8-2017; Yellow 2017.B Notebook pp.33-34	
        @printf("\n Schultz1959(2.5a) with 4/9 re-derivation: Feynman α→0 expansion: rfa= %g (int units) = %g m [SI]",rfa,scale*rfa )
        rfb=3*(pi/2)^0.5 * α 
        @printf("\n Schultz1959(2.5b): Feynman α→∞ expansion: rf= %g (int units) = %g m [SI]",rfb,scale*rfb )
        
        # Schultz1959 - Between (5.7) and (5.8) - resonance of Feynman SHM system
        phononfreq=sqrt(k/M)
        @printf("\n Schultz1959(5.7-5.8): fixed-e: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            phononfreq,phononfreq*ω /(2*pi), phononfreq*hbar*ω *1000/q)
 
        phononfreq=sqrt(k/mu) # reduced mass
        @printf("\n Schultz1959(5.7-5.8): reducd mass: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
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
        append!(p.A,A(v,w,βred))
        append!(p.B,B(v,w,βred,α))
        append!(p.C,C(v,w,βred))
        append!(p.F,F(v,w,βred,α))
        

        # FHIP 
        #    - low-T mobility, final result of Feynman1962
        # [1.60] in Devreese2016 page 36; 6th Edition of Frohlich polaron notes (ArXiv)
        # I believe here β is in SI (expanded) units
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
        @printf("\n\tGamma0 = %g rad/s = %g /s Tau=1/Gamma0 = %g = %f ps",
            Gamma0, Gamma0/(2*pi), 2*pi/Gamma0, 2*pi*1E12/Gamma0)
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

    return(p)
end

"""
function ImX(v,w,βred,α,ω,mb)
 
Impedance in (47a) from Feynman1962, directly solving freq dep without taking 
Hellwarth1999 limit of v->0
"""
function ImX(v,w,βred,α,ω,mb)
        @printf("Ain't no Imaginary Impedance like a Feynman 1962 ImX...")
    # Feynman, I love you - but using Nu, v; Omega, w in the same paper + formulas, for similar objects?!
    for nu=0:0.5:20
        R=(v^2-w^2)/(w^2*v) # inline, page 300 just after Eqn (2)
        b=R*βred/sinh(βred*v/2) # Feynman1962 version; page 1010, Eqn (47b)
        a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2))
        k(u,a,b,v,nu) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u)*cos(nu*u) # integrand with cos(vu) term, as (47a) 
        @printf("Numerical integration of Feynman1962(42a): nu=%.2f ",nu)
        @time n=quadgk(u->k(u,a,b,v,nu),0,Inf,maxevals=10^9,reltol=0.1) # numerical quadrature integration of (2)
        K=n[1]
        err=n[2]
        @printf(" quadgk: K=%g err=%g\n",K,err)

        # Full 47a constructed here
        ImX= 2*α/(3*sqrt(π)) * βred^(3/2) * (sinh(βred*nu/2))/sinh(βred/2) * (v^3/w^3) * K

        μ=ImX^-1 * (q)/(ω*mb)

        @printf(" %.3f %g %g\n",nu,ImX,μ)
    end

    @printf("\n\n")
end

function savepolaron(fileprefix, p::Polaron) 
    println("Saving data to $fileprefix.dat ...")
    f=open("$fileprefix.dat","a")
    @printf(f,"# %s \n# Ts, βreds, Kμs, Hμs, FHIPμs, vs, ws, ks, Ms, As, Bs, Cs, Fs, Taus, rfsis\n",fileprefix)
    @printf(f,"#  1    2     3    4     5     6   7    8  9   10  11  12  13  14 15\n") # columns for GNUPLOT etc.
    for i in 1:length(p.T)
        @printf(f,"%d %03f %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
        p.T[i], p.βred[i], p.Kμ[i], p.Hμ[i], p.FHIPμ[i], 
        p.v[i], p.w[i],
        p.k[i], p.M[i], p.A[i], p.B[i], p.C[i], p.F[i], 
        p.Tau[i], p.rfsi[i])
    end
    close(f)
end

function plotpolaron(fileprefix, p::Polaron)  
    println("Plotting polaron to $fileprefix...")

    #####
    ## Mass vs. Temperature plot
    plot(p.T,p.M,label="Phonon effective-mass",markersize=3,marker=:rect,xlab="Temperature (K)",ylab="Phonon effective-mass",ylim=(0,1.2))

    savefig("$fileprefix-mass.png")
    savefig("$fileprefix-mass.eps")

    #####
    ## Relaxationtime vs. Temperature plot
    plot(p.T,p.Tau,label="Kadanoff relaxation time (ps)",markersize=3,marker=:rect,xlab="Temperature (K)",ylab="Relaxation time (ps)",ylim=(0,1.2))

    savefig("$fileprefix-tau.png")
    savefig("$fileprefix-tau.eps")

    ## Mass + relaxation time vs. Temperature plot
    plot(p.T,p.M,label="Phonon effective-mass (m\$_b\$)",markersize=3,marker=:rect,
        xlab="Temperature (K)",ylab="Effective-mass / relaxation time",ylim=(0,1.2))
    plot!(p.T,p.Tau,label="Kadanoff relaxation time (ps)",markersize=3,marker=:diamond,
        xlab="Temperature (K)",ylab="Relaxation time (ps)",ylim=(0,1.2))

    savefig("$fileprefix-mass-tau.png")
    savefig("$fileprefix-mass-tau.eps")

    #####
    ## Spring Constants vs. Temperature plot
    plot(p.T,p.k,label="Polaron spring-constant",markersize=3, marker=:uptriangle, xlab="Temperature (K)",ylab="Spring-constant",)

    savefig("$fileprefix-spring.png")
    savefig("$fileprefix-spring.eps")

    #####
    ## Variation Energy vs. Temperature plots
    plot( p.T,p.A,label="A",markersize=3,marker=:downtriangle, xlab="Temperature (K)",ylab="Polaron free-energy")
    plot!(p.T,p.B,label="B",markersize=3,marker=:diamond)
    plot!(p.T,p.C,label="C",markersize=3,marker=:uptriangle)
    plot!(p.T,p.F,label="F",markersize=3,marker=:rect)
    #plot!(Ts,Fs,label="F=-(A+B+C)",markersize=3,marker=:rect)

    savefig("$fileprefix-variational.png")
    savefig("$fileprefix-variational.eps")

    #####
    ## Polaron radius vs. Temperature
    plot(p.T,p.rfsi.*10^10,label="Schultz Feynman radius",xlab="Temperature (K)",ylab="Polaron Radius (Angstrom)")
    savefig("$fileprefix-radius.png")
    savefig("$fileprefix-radius.eps")

    #####
    ## Calculated mobility comparison plot
    plot(p.T,p.Kμ,label="Kadanoff",markersize=3,marker=:rect,xlab="Temperature (K)",ylab="Mobility (cm\$^2\$/Vs)",ylims=(0,1000))
    plot!(p.T,p.FHIPμ,label="FHIP",markersize=3,marker=:diamond)
    plot!(p.T,p.Hμ,label="Hellwarth1999",markersize=3,marker=:uptriangle)

    savefig("$fileprefix-mobility-calculated.png")
    savefig("$fileprefix-mobility-calculated.eps")
end

end

