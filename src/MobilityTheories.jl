# MobilityTheories.jl

"""
    polaronmobility(Trange, ε_Inf, ε_S, freq, effectivemass; verbose::Bool=false)

Solves the Feynman polaron problem variationally with finite temperature Osaka energies. From the resulting v, and w parameters, calculates polaron structure (wave function size, etc.).  Uses FHIP, Kadanoff (Boltzmann relaxation time) and Hellwarth direct contour integration to predict a temperature-dependent mobility for the material system.

# Arguments  
- `Trange::range`: temperature range.
- `ε_Inf`: reduced optical dielectric constant.
- `ε_S`: reduced static dielectric constant.
- `freq`: characteristic dielectric phonon frequency (THz).
- `effectivemass`: bare-band effective mass (mₑ).

Returns a structure of type `Polaron`, containing arrays of useful
information.  Also prints a lot of information to the standard out - which
may be more useful if you're just inquiring as to a particular data point,
rather than plotting a temperature-dependent parameter.

As an example, to calculate the electron polaron in MAPI at 300 K:
# Examples
```jldoctest
polaronmobility(300, 4.5, 24.1, 2.25, 0.12)
```
"""
function polaronmobility(Trange, ε_Inf, ε_S, freq, effectivemass; verbose::Bool=false)
    println("\n\nPolaron mobility for system ε_Inf=$ε_Inf, ε_S=$ε_S, freq=$freq,
                 effectivemass=$effectivemass; with Trange $Trange ...")

    # Internally we we use 'mb' for the 'band mass' in SI units, of the effecitve-mass of the electron
    mb = effectivemass * MassElectron
    ω = (2 * pi) * freq * 1e12 # angular-frequency, of the phonon mode

    α = frohlichalpha(ε_Inf, ε_S, freq, effectivemass)
    #α=2.395939683378253 # Hard coded; from MAPI params, 4.5, 24.1, 2.25THz, 0.12me

    v = 7.1 # starting guess for v,w variational parameters
    w = 6.5

    @printf("Polaron mobility input parameters: ε_Inf=%f ε_S=%f freq=%g α=%f \n", ε_Inf, ε_S, freq, α)
    @printf("Derived params in SI: ω =%g mb=%g \n", ω, mb)

    # Empty struct for storing data
    # A slightly better way of doing this ф_ф ...
    p = OldPolaron()
    # populate data structure with (athermal) parameters supplied...
    append!(p.α, α) # appending so as not to mess with type immutability
    append!(p.mb, mb)
    append!(p.ω, ω)

    # We define βred as the subsuming the energy of the phonon; i.e. kbT c.f. ħω
    for T in Trange
        β = 1 / (kB * T)
        βred = ħ * ω * β
        append!(p.βred, βred)
        @printf("T: %f β: %.3g βred: %.3g ħω  = %.3g meV\t", T, β, βred, 1E3 * ħ * ω / q)

        if T == 0
            v, w = feynmanvw(v, w, α, 1.0)
        else
            v, w = feynmanvw(v, w, α, 1.0, βred)
        end

        @printf("\n Polaron Parameters:  v= %.4f  w= %.4f", v, w)

        # From 1962 Feynman, definition of v and w in terms of the coupled Mass and spring-constant
        # See Page 1007, just after equation (18)
        # Units of M appear to be 'electron masses'
        # Unsure of units for k, spring coupling constant
        k = (v^2 - w^2)
        M = (v^2 - w^2) / w^2
        append!(p.k, k)
        append!(p.M, M)

        @printf("  ||   M=%f  k=%f\t", M, k)
        @printf("\n Bare-band effective mass: %f Polaron effective mass: %f Polaron mass enhancement: %f%%", effectivemass, effectivemass * (1 + M), M * 100)

        @printf("\n Polaron frequency (SI) v= %.2g Hz  w= %.2g Hz",
            v * ω / (2 * pi), w * ω / (2 * pi))

        # (46) in Feynman1955
        meSmallAlpha(α) = α / 6 + 0.025 * α^2
        # (47) In Feynman1955
        meLargeAlpha(α) = 16 * α^4 / (81 * π^4)
        #meLargeAlpha(α )=202*(α /10)^4
        if (verbose) # asymptotic solutions - not that interesting when you have the actual ones!
            @printf("\n Feynman1955(46,47): meSmallAlpha(α)= %.3f meLargeAlpha(α)= %.3f",
                meSmallAlpha(α), meLargeAlpha(α))
            @printf("\n Feynman1962: Approximate ~ Large alpha limit, v/w = %.2f  =~approx~= alpha^2 = %.2f ",
                v / w, α^2)
        end

        # POLARON SIZE
        @printf("\n Polaron size (rf), following Schultz1959. (s.d. of Gaussian polaron ψ )")
        # Schultz1959 - rather nicely he actually specifies everything down into units!
        # just before (2.4) in Schultz1959
        mu = ((v^2 - w^2) / v^2)
        # (2.4)
        rf = sqrt(3 / (2 * mu * v))
        # (2.4) SI scaling inferred from units in (2.5a) and Table II
        rfsi = rf * sqrt(2 * me * ω)
        @printf("\n\t Schultz1959(2.4): rf= %g (int units) = %g m [SI]", rf, rfsi)
        append!(p.rfsi, rfsi)

        if (verbose)
            scale = sqrt(2 * mb * ω) # Note we're using mb;
            #band effective-mass in SI units (i.e. meff*melectron)

            rfa = (3 / (0.44 * α))^0.5 # As given in Schultz1959(2.5a), but that 0.44 is actually 4/9
            @printf("\n\t Schultz1959(2.5a) with 0.44: Feynman α→0 expansion: rfa= %g (int units) = %g m [SI]", rfa, scale * rfa)
            rfa = (3 / ((4 / 9) * α))^0.5 # Rederived from Feynman1955, 8-8-2017; Yellow 2017.B Notebook pp.33-34
            @printf("\n\t Schultz1959(2.5a) with 4/9 re-derivation: Feynman α→0 expansion: rfa= %g (int units) = %g m [SI]", rfa, scale * rfa)
            append!(p.rfsmallalpha, scale * rfa)

            rfb = 3 * (pi / 2)^0.5 * α
            @printf("\n\t Schultz1959(2.5b): Feynman α→∞ expansion: rf= %g (int units) = %g m [SI]", rfb, scale * rfb)

            # Schultz1959 - Between (5.7) and (5.8) - resonance of Feynman SHM system
            phononfreq = sqrt(k / M)
            @printf("\n\t Schultz1959(5.7-5.8): fixed-e: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
                phononfreq, phononfreq * ω / (2 * pi), phononfreq * hbar * ω * 1000 / q)

            phononfreq = sqrt(k / mu) # reduced mass
            @printf("\n\t Schultz1959(5.7-5.8): reducd mass: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
                phononfreq, phononfreq * ω / (2 * pi), phononfreq * hbar * ω * 1000 / q)

            @printf("\n\t Schultz1959: electronfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
                sqrt(k / 1), sqrt(k / 1) * ω / (2 * pi), sqrt(k / 1) * hbar * ω * 1000 / q)
            @printf("\n\t Schultz1959: combinedfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
                sqrt(k / (1 + M)), sqrt(k / (1 + M)) * ω / (2 * pi), sqrt(k / (1 + M)) * hbar * ω * 1000 / q)

            # Devreese1972: 10.1103/PhysRevB.5.2367
            # p.2371, RHS.
            @printf("\n Devreese1972: (Large Alpha) Franck-Condon frequency = %.2f", 4 * α^2 / (9 * pi))
        end

        # F(v,w,β,α)=-(A(v,w,β)+B(v,w,β,α)+C(v,w,β)) #(62a) - Hellwarth 1999
        if (T > 0)
            @printf("\n Polaron Free Energy: A= %f B= %f C= %f F= %f", A(v, w, βred), B(v, w, βred, α), C(v, w, βred), F(v, w, βred, α)[1])
            @printf("\t = %f meV", 1000.0 * F(v, w, βred, α)[1] * ħ * ω / q)
            append!(p.A, A(v, w, βred))
            append!(p.B, B(v, w, βred, α))
            append!(p.C, C(v, w, βred))
            append!(p.F, F(v, w, βred, α))
        else # Athermal case; Enthalpy
            @printf("\n Polaron Enthalpy: F= %f = %f meV \n", F(v, w, α, 1.0), 1_000 * F(v, w, α, 1.0) * ħ * ω / q)

            return # return early, as if T=0, all mobility theories = infinite / fall over
        end

        # FHIP
        #    - low-T mobility, final result of Feynman1962
        # [1.60] in Devreese2016 page 36; 6th Edition of Frohlich polaron notes (ArXiv)
        # I believe here β is in SI (expanded) units
        @printf("\n Polaron Mobility theories:")
        μ = (w / v)^3 * (3 * q) / (4 * mb * ħ * ω^2 * α * β) * exp(ħ * ω * β) * exp((v^2 - w^2) / (w^2 * v))
        @printf("\n\tμ(FHIP)= %f m^2/Vs \t= %.2f cm^2/Vs", μ, μ * 100^2)
        append!(p.T, T)
        append!(p.FHIPμ, μ * 100^2)


        # Kadanoff
        #     - low-T mobility, constructed around Boltzmann equation.
        #     - Adds factor of 3/(2*beta) c.f. FHIP, correcting phonon emission behaviour
        # [1.61] in Devreese2016 - Kadanoff's Boltzmann eqn derived mob
        μ = (w / v)^3 * (q) / (2 * mb * ω * α) * exp(ħ * ω * β) * exp((v^2 - w^2) / (w^2 * v))
        @printf("\n\tμ(Kadanoff,via Devreese2016)= %f m^2/Vs \t= %.2f cm^2/Vs", μ, μ * 100^2)

        append!(p.Kμ, μ * 100^2)

        ######
        # OK, now deep-diving into Kadanoff1963 itself to extract Boltzmann equation components
        # Particularly right-hand-side of page 1367
        #
        # From 1963 Kadanoff, (9), define eqm. number of phonons (just from T and phonon omega)
        Nbar = (exp(βred) - 1)^-1
        @printf("\n\t\tEqm. Phonon. pop. Nbar: %f ", Nbar)
        if (verbose)
            @printf("\n\texp(Bred): %f exp(-Bred): %f exp(Bred)-1: %f", exp(βred), exp(-βred), exp(βred) - 1)
        end
        Nbar = exp(-βred)
        #Note - this is only way to get Kadanoff1963 to be self-consistent with
        #FHIP, and later statements (Devreese) of the Kadanoff mobility.
        #It suggests that Kadanoff used the wrong identy for Nbar in 23(b) for
        #the Gamma0 function, and should have used a version with the -1 to
        #account for Bose / phonon statistics

        #myv=sqrt(k*(1+M)/M) # cross-check maths between different papers
        #@printf("\nv: %f myv: %f\n",v,myv)

        # Between 23 and 24 in Kadanoff 1963, for small momenta skip intergration --> Gamma0
        Gamma0 = 2 * α * Nbar * (M + 1)^(1 / 2) * exp(-M / v)
        Gamma0 *= ω  #* (ω *hbar) # Kadanoff 1963 uses hbar=omega=mb=1 units
        # Factor of omega to get it as a rate relative to phonon frequency
        # Factor of omega*hbar to get it as a rate per energy window
        μ = q / (mb * (M + 1) * Gamma0) #(25) Kadanoff 1963, with SI effective mass

        if (verbose) # these are cross-checks
            @printf("\n\tμ(Kadanoff1963 [Eqn. 25]) = %f m^2/Vs \t = %.2f cm^2/Vs", μ, μ * 100^2)
            @printf("\n\t\t Eqm. Phonon. pop. Nbar: %f ", Nbar)
        end

        @printf("\n\t\tGamma0 = %g rad/s = %g /s ",
            Gamma0, Gamma0 / (2 * pi))
        @printf(" \n\t\tTau=1/Gamma0 = %g s = %f ps",
            2 * pi / Gamma0, 2 * pi * 1E12 / Gamma0)
        Eloss = hbar * ω * Gamma0 / (2 * pi) # Simply Energy * Rate
        @printf("\n\t\tEnergy Loss = %g J/s = %g meV/ps", Eloss, Eloss * 1E3 / (q * 1E12))
        append!(p.Tau, 2 * pi * 1E12 / Gamma0) # Boosted into ps ?


        # Hellwarth1999 - directly do contour integration in Feynman1962, for
        # finite temperature DC mobility
        # Hellwarth1999 Eqn (2) and (1) - These are going back to the general
        # (pre low-T limit) formulas in Feynman1962.  to evaluate these, you
        # need to do the explicit contour integration to get the polaron
        # self-energy
        R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)

        #b=R*βred/sinh(b*βred*v/2) # This self-references b! What on Earth?
        # OK! I now understand that there is a typo in Hellwarth1999 and
        # Biaggio1997. They've introduced a spurious b on the R.H.S. compared to
        # the original, Feynman1962:
        b = R * βred / sinh(βred * v / 2) # Feynman1962 version; page 1010, Eqn (47b)

        a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))
        k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)
        K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

        #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS = α / (3 * sqrt(π)) * βred^(5 / 2) / sinh(βred / 2) * (v^3 / w^3) * K
        μ = RHS^-1 * (q) / (ω * mb)
        @printf("\n\tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs", μ, μ * 100^2)
        append!(p.Hμ, μ * 100^2)

        if (verbose)
            # Hellwarth1999/Biaggio1997, b=0 version... 'Setting b=0 makes less than 0.1% error'
            # So let's test this
            R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)
            b = 0
            a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))
            #k(u,a,b,v) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u) # integrand in (2)
            K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

            #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
            RHS = α / (3 * sqrt(π)) * βred^(5 / 2) / sinh(βred / 2) * (v^3 / w^3) * K
            μ = RHS^-1 * (q) / (ω * mb)
            @printf("\n\tμ(Hellwarth1999,b=0)= %f m^2/Vs \t= %.2f cm^2/Vs", μ, μ * 100^2)
            @printf("\n\tError due to b=0; %f", (100^2 * μ - p.Hμ[length(p.Hμ)]) / (100^2 * μ))
            #append!(Hμs,μ*100^2)
        end
        @printf("\n") # blank line at end of spiel.

        # Recycle previous variation results (v,w) as next guess
        initial = [v, w] # Caution! Might cause weird sticking in local minima

        append!(p.v, v)
        append!(p.w, w)

    end

    return (p)
end

"""
polaron is the outer constructor for the Polaron type. This function evaluates model data for the polaron, i.e. unitless and not material specific. 
"""
function polaron(αrange, Trange, Ωrange, ω, v_guesses, w_guesses; verbose=false)

    # v_guesses and w_guesses are initial values for v and w (including many v and w parameters).
    # These guesses are generally not needed unless instabilities are found in the minimisation and better initial values improve stability.

    # Get the length of any arrays etc.
    num_α = size(αrange, 1)
    num_T = length(Trange)
    num_Ω = length(Ωrange)
    num_ω = length(ω)

    println(αrange)

    # For multiple variational modes, ensure that the number of v and w parameters is the same.
    @assert length(v_guesses) == length(w_guesses) "v and w guesses must be the same length."
    num_vw = length(v_guesses)

    # Instantiate all the polaron data that will go into the Polaron type.
    p = Dict(
        "α"     => αrange,                                          # alphas
        "αeff"  => sum(αrange, dims=2),                             # alphas sums
        "T"     => Trange,                                          # temperatures
        "ω"     => ω,                                               # phonon frequencies
        "β"     => Matrix{Float64}(undef, num_T, num_ω),            # betas
        "Ω"     => Ωrange,                                          # photon frequencies
        "v0"    => Matrix{Float64}(undef, num_α, num_vw),           # v ground state params
        "w0"    => Matrix{Float64}(undef, num_α, num_vw),           # w ground state params
        "F0"    => Vector{Float64}(undef, num_α),                   # ground state energies
        "A0"    => Vector{Float64}(undef, num_α),                   # A ground state parameter
        "B0"    => Vector{Float64}(undef, num_α),                   # B ground state parameter
        "C0"    => Vector{Float64}(undef, num_α),                   # C ground state parameter
        "v"     => Array{Float64,3}(undef, num_T, num_α, num_vw),   # v params
        "w"     => Array{Float64,3}(undef, num_T, num_α, num_vw),   # w params
        "F"     => Matrix{Float64}(undef, num_T, num_α),            # energies
        "A"     => Matrix{Float64}(undef, num_T, num_α),            # A parameter
        "B"     => Matrix{Float64}(undef, num_T, num_α),            # B parameter
        "C"     => Matrix{Float64}(undef, num_T, num_α),            # C parameter
        "κ"     => Array{Float64,3}(undef, num_T, num_α, num_vw),   # spring constants
        "M"     => Array{Float64,3}(undef, num_T, num_α, num_vw),   # fictitious masses
        "R"     => Array{Float64,3}(undef, num_T, num_α, num_vw),   # polaron radii
        "z"     => Array{ComplexF64,3}(undef, num_Ω, num_T, num_α), # complex impedences
        "σ"     => Array{ComplexF64,3}(undef, num_Ω, num_T, num_α), # complex conductivities
        "μ"     => Matrix{Float64}(undef, num_T, num_α)             # mobilities
    )

    # Print the phonon frequencies. 
    if verbose
        # Instantiate IO. Compact and limit reduce printed info.
        io = IOContext(stdout, :compact => true, :limit => true) 

        # Hide cursor whilst printing live data.
        print(io, "\e[?25l") 

        println(io, "\e[K-----------------------------------------------------------------------")
        println(io, "\e[K                         Polaron Information:                          ")
        println(io, "\e[K-----------------------------------------------------------------------")

        # Different formatting for single vs multipelf frequencies (limits an array to head and tail to limit prints).
        if num_ω == 1
            println(io, "\e[KPhonon frequencies         | ω = ", ω, " ω₀")
        else
            println(io, "\e[KPhonon frequencies         | ω = ", join(round.(first(ω, 2), digits=1), ", ")..., " ... ", join(round.(last(ω, 2), digits=1), ", ")..., " ω₀")
        end
        process = 1     # Counter for total iterations.
        αprocess = 1    # Counter for αrange iterations. 
    end

    for j in axes(αrange, 1)    # αrange loop.

        # Reduce αrange to either a Vector of values per phonon mode, or a single scalar. 
        α = reduce_array(αrange[j, :])

        # Print α coupling parameters.
        if verbose
            if num_ω == 1
                println(io, "\e[KFröhlich coupling          | α = ", join(α, ", ")...)
            else
                println(io, "\e[KFröhlich coupling          | α = ", join(round.(first(α, 2), digits=3), ", ")..., " ... ", join(round.(last(α, 2), digits=3), ", ")...)
            end
        end

        # Extract the ground-state, athermal polaron properties (energy (enthalpy) and variational parameters v and w).
        # w is also the frequency of oscillation of the SHM trial system composed of the bare particle and fictitous mass.
        # A, B, C are components of the total energy: A is the bare electron energy, B the electron-phonon interaction energy, C is the energy of the harmonic trial system.
        v_gs, w_gs, F_gs, A_gs, B_gs, C_gs = feynmanvw(v_guesses, w_guesses, α, ω)

        # Update the guesses to keep them close-ish to the true solutions during loops over alphas.
        v_guesses, w_guesses = v_gs, w_gs

        # Store the athermal data.
        p["v0"][j, :] .= v_gs
        p["w0"][j, :] .= w_gs
        p["F0"][j] = F_gs
        p["A0"][j] = A_gs
        p["B0"][j] = B_gs
        p["C0"][j] = C_gs

        # Print ground-state data.
        if verbose
            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[K              Ground State Information: [$(αprocess[]) / $(num_α) ($(round(αprocess[] / (num_α) * 100, digits=1)) %)]")
            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[KGS variational parameter   | v₀ = ", v_gs, " ω₀")
            println(io, "\e[KGS variational parameter   | w₀ = ", w_gs, " ω₀")
            println(io, "\e[KGS Energy                  | E₀ = ", F_gs, " ħω₀")
            println(io, "\e[KGS Electron energy         | A₀ = ", A_gs, " ħω₀")
            println(io, "\e[KGS Interaction energy      | B₀ = ", B_gs, " ħω₀")
            println(io, "\e[KGS Trial energy            | C₀ = ", C_gs, " ħω₀")
            Tprocess = 1    # Counter for Trange.
            αprocess += 1   # Increment αrange iteration.
        end

        # Calculate and store fictitious spring constants. See after Eqn. (18), pg. 1007 of Feynman1962. Thermal
        κ_gs = v_gs .^ 2 .- w_gs .^ 2

        # Calculate and store fictitious masses. Athermal.
        M_gs = κ_gs ./ w_gs .^ 2

        # Small coupling (α → 0) polaron mass approximation. Eqn. (46) in Feynman1955.
        M_small = α ./ 6 .+ 0.025 .* α .^2

        # Large coupling (α → ∞) polaron mass approximation. Eqn. (47) In Feynman1955.
        M_large = 16 .* α .^4 ./ (81 * π^4)

        # Approximate large coupling asymptotic limit for the polaron mass. Feynman1962. Athermal
        M_asymp_gs = v_gs ./ w_gs

        # Reduced mass of particle and fictitious mass system. Before eqn. (2.4) in Schultz1959. Athermal.
        M_reduced_gs = v_gs .^2 .- w_gs .^2 / v_gs .^ 2

        # Calculate and store polaron radii. Approximates the polaron wavefunction as a Gaussian and relates the size to the standard deviation. Eqn. (2.4) in Schultz1959. Athermal.
        R_gs = sqrt.(3 .* v_gs ./ (v_gs .^ 2 .- w_gs .^ 2) .^ 2)

        # Small coupling (α → 0) polaron radii approximiation. Eqn. (2.5a) in Schultz1959.
        R_small = sqrt.(3 ./ (4 ./ 9 .* α)) 

        # Large coupling (α → ∞) polaron radii approximiation. Eqn. (2.5b) in Schultz1959.
        R_large = 3 * √(π / 2) .* α

        # Franck-Condon (FC) frequency in large coupling (α → ∞) limit. RHS of pg. 2371 in Devreese1972.
        Ω_FC = 4 / 9π .* α .^2

        for i in eachindex(Trange)  # Temperatures loop.
            T = Trange[i]

            # Print temperature.
            if verbose
                println(io, "\e[K-----------------------------------------------------------------------")
                println(io, "\e[K         Finite Temperature Information: [$(Tprocess[]) / $(num_T) ($(round(Tprocess[] / (num_T) * 100, digits=1)) %)]")
                println(io, "\e[K-----------------------------------------------------------------------")
                println(io, "\e[KTemperatures               | T = ", T, " K")
            end

            # Calculate the reduced (unitless) thermodynamic betas for each phonon mode.
            β = ω ./ T
            p["β"][i, :] .= β

            # Print thermodynamic betas.
            if verbose
                if num_ω == 1
                    println(io, "\e[KReduced thermodynamic      | β = ", β)
                else
                    println(io, "\e[KReduced thermodynamic      | β = ", join(round.(first(β, 2), digits=3), ", ")..., " ... ", join(round.(last(β, 2), digits=3), ", ")...)
                end
            end

            # Calculate thermal polaron properties (energy (Gibbs free energy) and variational parameters v and w).
            # w is also the frequency of oscillation of the SHM trial system composed of the bare particle and fictitous mass.
            # A, B, C are components of the total energy: A is the bare electron energy, B the electron-phonon interaction energy, C is the energy of the harmonic trial system.
            v, w, F, A, B, C = feynmanvw(v_guesses, w_guesses, α, ω, β)

            # Update the guesses to keep them close-ish to the true solutions during loops over temperatures.
            v_guesses, w_guesses = v, w

            # Store thermal data.
            p["v"][i, j, :] .= v
            p["w"][i, j, :] .= w
            p["F"][i, j] = F
            p["A"][i, j] = A
            p["B"][i, j] = B
            p["C"][i, j] = C

            # Print thermal data.
            if verbose
                println(io, "\e[KVariational parameter      | v = ", v, " ω₀")
                println(io, "\e[KVariational parameter      | w = ", w, " ω₀")
                println(io, "\e[KFree energy                | F = ", F, " ħω₀")
                println(io, "\e[KElectron energy            | A = ", A, " ħω₀")
                println(io, "\e[KInteraction energy         | B = ", B, " ħω₀")
                println(io, "\e[KTrial energy               | C = ", C, " ħω₀")
            end

            # Calculate and store fictitious spring constants. See after Eqn. (18), pg. 1007 of Feynman1962. Thermal
            κ = v .^ 2 .- w .^ 2
            p["κ"][i, j, :] .= κ

            # Print spring constants.
            if verbose
                println(io, "\e[K-----------------------------------------------------------------------")
                println(io, "\e[K                      Trial System Information:                        ")
                println(io, "\e[K-----------------------------------------------------------------------")
                println(io, "\e[KFictitious spring constant | κ = ", κ, " m₀/ω₀²")
            end

            # Calculate and store fictitious masses. Thermal.
            M = κ ./ w .^ 2
            p["M"][i, j, :] .= M

            # Print masses.
            if verbose
                println(io, "\e[KFictitious mass            | M = ", M, " m₀.")
            end

            # Approximate large coupling asymptotic limit for the polaron mass. Feynman1962. Thermal
            M_asymp = v ./ w

            # Reduced mass of particle and fictitious mass system. Before eqn. (2.4) in Schultz1959. Athermal.
            M_reduced = v .^2 .- w .^2 / v .^ 2

            # Calculate and store polaron radii.
            R = sqrt.(3 .* v ./ (v .^ 2 .- w .^ 2) .^ 2)
            p["R"][i, j, :] .= R

            # Print radii.
            if verbose
                println(io, "\e[KPolaron radius             | R = ", R, " √(ħ/2m₀ω₀)")
            end

            # Calculate and store the DC mobiliies.
            μ = polaron_mobility(v, w, α, ω, β)
            p["μ"][i, j] = μ

            # Print DC mobilities.
            if verbose
                println(io, "\e[KMobility                   | μ = ", μ, " q/m₀ω₀")
                Ωprocess = 1    # Counter for Ωrange.
                Tprocess += 1   # Increment Trange iterator.
            end

            for k in eachindex(Ωrange)  # E-field frequencies loop. 
                Ω = Ωrange[k]

                # Print E-field frequency.
                if verbose
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K           Linear Reponse Information: [$(Ωprocess[]) / $(num_Ω) ($(round(Ωprocess[] / (num_Ω) * 100, digits=1)) %)]")
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[KElectric field frequency   | Ω = ", Ω, " ω₀")
                end

                # Calculate and store polaron memory functions (akin to self energy).
                χ = polaron_memory_function(v, w, α, ω, β, Ω)

                # Print memory function.
                if verbose
                    println(io, "\e[KMemory function            | χ = ", χ .|> y -> ComplexF64.(y), " ω₀m₀V₀/q²")
                end

                # Calculate and store polaron complex impedances.
                z = -im * Ω + im * χ
                p["z"][k, i, j] = z

                # Print complex impedances.
                if verbose
                    println(io, "\e[KComplex impedance          | z = ", z .|> y -> ComplexF64.(y), " ω₀m₀V₀/q²")
                end

                # Calculate and store polaron complex conductivities.
                σ = 1 / z
                p["σ"][k, i, j] = σ

                # Print complex conductivities and show total algorithm progress.
                if verbose
                    println(io, "\e[KComplex conductivity       | σ = ", σ .|> y -> ComplexF64.(y), " q²/ω₀m₀V₀")
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K[Total Progress: $(process[]) / $(num_α * num_T * num_Ω) ($(round(process[] / (num_α * num_T * num_Ω) * 100, digits=1)) %)]")
                    print(io, "\e[9F")
                    Ωprocess += 1   # Increment Ωrange iterator.
                    process += 1    # Increment total iterator.
                end
            end
            if verbose print(io, "\e[18F")end   # Move up 18 lines and erase.
        end 
        if verbose print(io, "\e[10F") end      # Move up 10 lines and erase.
    end
    if verbose print(io, "\e[5F\e[?25h") end    # Move up 5 lines and erase. Show cursor.

    # Extract polaron data as Vector ready for Polaron type.
    polaron_data = [
        p["α"],     # alphas
        p["αeff"],  # alphas sums
        p["T"],     # temperatures
        p["ω"],     # phonon frequencies
        p["β"],     # betas
        p["Ω"],     # photon frequencies
        p["v0"],    # v params ground state
        p["w0"],    # w params ground state
        p["F0"],    # energies ground state
        p["A0"],    # A parameter, bare electron energy ground state
        p["B0"],    # B parameter, interaction energy ground state
        p["C0"],    # C parameter, trial system free energy ground state
        p["v"],     # v params
        p["w"],     # w params
        p["F"],     # energies
        p["A"],     # A parameter, bare electron energy
        p["B"],     # B parameter, interaction energy
        p["C"],     # C parameter, trial system free energy
        p["κ"],     # spring constants
        p["M"],     # fictitious masses
        p["R"],     # polaron radii
        p["z"],     # complex impedences
        p["σ"],     # complex conductivities
        p["μ"]      # mobilities
    ]

    # Return Polaron type containing all generated data over any coupling strengths, temperatures and frequencies.
    return Polaron(polaron_data...)
end

"""
Set default guesses for variational paramaters. Reduces user inputs.
"""
polaron(αrange, Trange, Ωrange, ω; verbose=false) = polaron(αrange, Trange, Ωrange, ω, 3.11, 2.87; verbose=verbose)

"""
Single alpha parameter. polaron() expects alpha parameters to be in a Vector.
"""
polaron(α::Real, Trange, Ωrange, ω, v_guesses, w_guesses; verbose=false) = polaron([α], Trange, Ωrange, ω, v_guesses, w_guesses; verbose=verbose)

"""
DC (zero frequency) limited version. Can just exclude a frequency input.
"""
polaron(αrange, Trange, ω; verbose=false) = polaron(αrange, Trange, 0, ω; verbose=verbose)

"""
Reduce to one single, unitless phonon mode ω=1 (polaron units).
"""
polaron(αrange, Trange; verbose=false) = polaron(αrange, Trange, 1; verbose=verbose)

"""
Only alpha parameter(s) just calculates ground-state and room temperature DC polarons.
"""
polaron(αrange; verbose=false) = polaron(αrange, 298, 1; verbose=verbose)

"""
Material specific constructors that use material specific parameters to parameterise the polaron.
Material data is inputted through the `Material` type.
Returns all data in either SI units or other common, suitable units otherwise.
"""
function polaron(material::Material, Trange, Ωrange, v_guesses, w_guesses; verbose=false)

    # Show material data.
    if verbose
        display(material)
    end
    
    # Unitful parameters. These are our measuring sticks!
    ω₀ = 2π * 1e12
    m₀ = me * 0.12
    r₀ = sqrt(ħ / (m₀ * ω₀))
    T₀ = ħ * ω₀ / kB

    # Extract material data from Material type.
    m_eff = material.mb
    volume = material.volume
    phonon_freqs = material.freqs

    # Generate polaron data from the arbitrary model constructor.
    p = polaron(material.α', Trange ./ T₀, Ωrange, phonon_freqs, v_guesses, w_guesses, verbose=verbose)

    # Add SI or otherise suitable units to the polaron data.
    F0_unit = p.F0 .* ħ * ω₀ / q * 1e3                                               # meV
    A0_unit = p.A0 .* ħ * ω₀ / q * 1e3                                               # meV
    B0_unit = p.B0 .* ħ * ω₀ / q * 1e3                                               # meV
    C0_unit = p.C0 .* ħ * ω₀ / q * 1e3                                               # meV
    F_unit = p.F .* ħ * ω₀ / q * 1e3                                                 # meV
    A_unit = p.A .* ħ * ω₀ / q * 1e3                                                 # meV
    B_unit = p.B .* ħ * ω₀ / q * 1e3                                                 # meV
    C_unit = p.C .* ħ * ω₀ / q * 1e3                                                 # meV
    M_units = p.M .* m_eff                                                           # electron mass units
    R_unit = [R * sqrt(ħ / 2 / m_eff / me / ω / 1e12) * 1e10 for R in p.R, ω in p.ω] # Angstroms
    z_unit = p.z .* (m₀ * r₀^2 * ω₀ / q^2)                                           # Ohms
    σ_unit = p.σ ./ (m₀ * r₀^2 * ω₀ / q^2)                                           # Siemens
    μ_unit = p.μ ./ m₀ / ω₀ * q * 1e4                                                # cm^2/Vs

    # Return material-specific, unitful Polaron type.
    return Polaron(p.α, p.αeff, p.T, phonon_freqs, p.β, p.Ω, p.v0, p.w0, F0_unit, A0_unit, B0_unit, C0_unit, p.v, p.w, F_unit, A_unit, B_unit, C_unit, p.κ, M_units, R_unit, z_unit, σ_unit, μ_unit)
end

"""
Set default guesses for variational paramaters. Reduces user inputs.
"""
polaron(material::Material, Trange, Ωrange; verbose=false) = polaron(material, Trange, Ωrange, 3.11, 2.87; verbose=verbose)

"""
DC (zero frequency) limit for material specific code. Can just exclude a frequency input.
"""
polaron(material::Material, Trange; verbose=false) = polaron(material, Trange, 0; verbose=verbose)

"""
Only material data calculates room-temperature DC polarons.
"""
polaron(material::Material; verbose=false) = polaron(material, 298; verbose=verbose)

"""
    save_polaron(p::Polaron, prefix)

Saves data from 'polaron' into file "prefix".
This is a .jdl file for storing the polaron data whilst preserving types. Allows for saving multidimensional arrays that sometimes arise in the polaron data.
Each parameter in the NewPolaron type is saved as a dictionary entry. E.g. NewPolaron.α is saved under JLD.load("prefix.jld")["alpha"].
"""
function save_polaron(polaron::Polaron, prefix)

    println("Saving polaron data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "alpha", polaron.α,
        "temperature", polaron.T,
        "beta", polaron.β,
        "phonon freq", polaron.ω,
        "v", polaron.v,
        "w", polaron.w,
        "spring", polaron.κ,
        "mass", polaron.M,
        "energy", polaron.F,
        "efield freq", polaron.Ω,
        "impedance", polaron.Z,
        "conductivity", polaron.σ,
        "mobility", polaron.μ
    )

    println("... Polaron data saved.")
end

"""
    load_polaron(p::NewPolaron, prefix)

Loads data from file "polaron_file_path" into a NewPolaron type.
"""
function load_polaron(polaron_file_path)

    println("Loading polaron data from $polaron_file_path ...")

    data = JLD.load("$polaron_file_path")

    polaron = NewPolaron(
        data["alpha"],
        data["temperature"],
        data["beta"],
        data["phonon freq"],
        data["v"],
        data["w"],
        data["spring"],
        data["mass"],
        data["energy"],
        data["efield freq"],
        data["impedance"],
        data["conductivity"],
        data["mobility"]
    )

    println("... Polaron loaded.")

    return polaron
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
    f = open("$fileprefix.dat", "w")

    @printf(f, "# %s\n", fileprefix) # put name / material at header
    @printf(f, "# Params in SI: ω =%g mb=%g \n", p.ω[1], p.mb[1])
    @printf(f, "# Alpha parameter: α = %f  \n", p.α[1])

    @printf(f, "# Ts, βreds, Kμs, Hμs, FHIPμs, vs, ws, ks, Ms, As, Bs, Cs, Fs, Taus, rfsis\n")
    @printf(f, "#  1    2     3    4     5      6   7   8   9  10  11  12  13    14     15\n") # columns for GNUPLOT etc.

    for i in 1:length(p.T)
        @printf(f, "%d %03f %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
            p.T[i], p.βred[i], p.Kμ[i], p.Hμ[i], p.FHIPμ[i],
            p.v[i], p.w[i],
            p.k[i], p.M[i], p.A[i], p.B[i], p.C[i], p.F[i],
            p.Tau[i], p.rfsi[i])
    end
    close(f)
end

"""
    Hellwarth1999mobilityRHS((α, (v, w) ,f), effectivemass, T)

Calculates the DC mobility using Hellwarth et al. 1999 Eqn. (2).

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
function Hellwarth1999mobilityRHS((α, (v, w), f), effectivemass, T)
    mb = effectivemass * MassElectron
    ω = f * 1e12 * 2π
    βred = ħ * ω / (kB * T)

    R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)
    b = R * βred / sinh(βred * v / 2) # Feynman1962 version; page 1010, Eqn (47b)
    a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))
    k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)
    K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

    # Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
    RHS = α / (3 * sqrt(π)) * βred^(5 / 2) / sinh(βred / 2) * (v^3 / w^3) * K
    μ = RHS^(-1) * q / (ω * mb)

    return 1 / μ
end

function Hellwarth_mobility(β, α, v, w; ω=ω)
    R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)
    b = R * β / sinh(β * v / 2) # Feynman1962 version; page 1010, Eqn (47b)
    a = sqrt((β / 2)^2 + R * β * coth(β * v / 2))
    k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)
    K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

    # Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
    RHS = α / (3 * sqrt(π)) * β^(5 / 2) / sinh(β / 2) * (v^3 / w^3) * K
    μ = RHS^(-1)
    return μ
end
