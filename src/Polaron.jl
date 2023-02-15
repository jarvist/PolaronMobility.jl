# Polaron.jl

"""
    Polaron(x...)

Type for storing the polaron information, `x...`.
"""
struct Polaron
    α       # Fröhlich coupling.
    αeff    # Effective Fröhlich coupling summed for multiple modes.
    T       # Temperature.
    ω       # Phonon mode frequencies.
    β       # Reduced thermodynamic beta ħω₀/kBT.
    Ω       # Electric field frequency.
    v0      # Athermal variational parameter v.
    w0      # Athermal variational parameter w.
    F0      # Polaron athermal energy (enthalpy).
    A0      # Bare electron enthalpy.
    B0      # ⟨S⟩ₜ interaction enthalpy.
    C0      # ⟨Sₜ⟩ₜ enthalpy of trial system.
    Fs      # Small alpha (α→0) approximate energy.
    Fl      # Large alpha (α→∞) approximate energy.
    κ0      # Fictitious spring constant.
    M0      # Athermal fictitious mass.
    Ms      # Small alpha (α→0) approximate fictitious mass.
    Ml      # Large alpha (α→∞) approximate fictitious mass.
    M0a     # Athermal asymptotic approximate fictitious mass (v0/w0).
    M0r     # Athermal reduced mass of particle + fictitious particle system.
    R0      # Athermal polaron radius (s. d. of a Gaussian wavefunction).
    Rs      # Small alpha (α→0) approximate polaron radius.
    Rl      # Large alpha (α→∞) approximate polaron radius.
    ΩFC     # Approximate Franck-Condon peak frequency.
    v       # Thermal variational parameter v.
    w       # Thermal variational parameter w.
    F       # Polaron free energy.
    A       # Bare electron free energy.
    B       # ⟨S⟩ₜ interaction energy.
    C       # ⟨Sₜ⟩ₜ free energy of trial system.
    κ       # Fictitious spring constant.
    M       # Thermal fictitious mass.
    Ma      # Thermal asymptotic approximate fictitious mass (v/w).
    Mr      # Thermal reduced mass of particle + fictitious particle system.
    R       # Schultz polaron radius.
    μ       # DC mobility.
    μFHIP   # FHIP low-temperature DC mobiltiy approximation.
    μD      # Kadanoff DC mobility, low-temperature approximation using Devreese2016.
    μK      # Kadanoff DC mobility, low-temperature approximation using Kadanoff1963.
    τ       # Relaxation time from Kadanoff Boltzmann transport equation.
    μH      # Hellwarth DC mobility.
    μH0     # Hellwarth DC mobility with b=0.
    χ       # Memory function or polaron self energy.
    z       # Complex impedence.
    σ       # Complex conductivity.

    function Polaron(x...)
        new(reduce_array.(x)...)
    end
end

"""
    polaron(αrange, Trange, Ωrange, ω, v_guesses, w_guesses; verbose=false)

Outer constructor for the Polaron type. This function evaluates model data for the polaron, i.e. unitless and not material specific. 

# Examples
```jldoctest
julia> polaron(6, 300, 3, 1.0, 3.6, 2.8)
```
"""
function polaron(αrange, Trange, Ωrange, ω, v_guesses, w_guesses; verbose=false)

    # v_guesses and w_guesses are initial values for v and w (including many v and w parameters).
    # These guesses are generally not needed unless instabilities are found in the minimisation and better initial values improve stability.

    # Get the length of any arrays etc.
    num_α = size(αrange, 1)
    num_T = length(Trange)
    num_Ω = length(Ωrange)
    num_ω = length(ω)

    # For multiple variational modes, ensure that the number of v and w parameters is the same.
    @assert length(v_guesses) == length(w_guesses) "v and w guesses must be the same length."
    num_vw = length(v_guesses)

    # Instantiate all the polaron data that will go into the Polaron type.
    p = Dict(
        "α"     => αrange, # Alphas.
        "αeff"  => sum(αrange, dims=2), # Alphas sums.
        "T"     => Trange, # Temperatures.
        "ω"     => ω, # Phonon frequencies.
        "β"     => Matrix{Any}(undef, num_T, num_ω), # Betas.
        "Ω"     => Ωrange, # Photon frequencies.
        "v0"    => Matrix{Float64}(undef, num_α, num_vw), # Athermal v params.
        "w0"    => Matrix{Float64}(undef, num_α, num_vw), # Athermal w params.
        "F0"    => Vector{Unitful.Energy}(undef, num_α), # Athermal energies.
        "A0"    => Vector{Unitful.Energy}(undef, num_α), # Athermal A parameter.
        "B0"    => Vector{Unitful.Energy}(undef, num_α), # Athermal B parameter.
        "C0"    => Vector{Unitful.Energy}(undef, num_α), # Athermal C parameter.
        "Fs"    => Vector{Unitful.Energy}(undef, num_α), # Small alpha (α→0) approximate energy.
        "Fl"    => Vector{Unitful.Energy}(undef, num_α), # Large alpha (α→∞) approximate energy.
        "κ0"    => Matrix{Any}(undef, num_α, num_vw), # Fictitious spring constant.
        "M0"    => Matrix{Unitful.Mass}(undef, num_α, num_vw), # Athermal fictitious mass.
        "Ms"    => Vector{Unitful.Mass}(undef, num_α), # Small alpha (α→0) approximate fictitious mass.
        "Ml"    => Vector{Unitful.Mass}(undef, num_α), # Large alpha (α→∞) approximate fictitious mass.
        "M0a"   => Matrix{Unitful.Mass}(undef, num_α, num_vw), # Athermal asymptotic approximate fictitious mass (v0/w0).
        "M0r"   => Matrix{Unitful.Mass}(undef, num_α, num_vw), # Athermal reduced mass of particle + fictitious particle system.
        "R0"    => Matrix{Unitful.Length}(undef, num_α, num_vw), # Athermal polaron radius (s. d. of a Gaussian wavefunction).
        "Rs"    => Vector{Unitful.Length}(undef, num_α), # Small alpha (α→0) approximate polaron radius.
        "Rl"    => Vector{Unitful.Length}(undef, num_α), # Large alpha (α→∞) approximate polaron radius.
        "ΩFC"   => Vector{Unitful.Frequency}(undef, num_α), # Large alpha (α → ∞)approximate Franck-Condon peak frequency.
        "v"     => Array{Float64,3}(undef, num_T, num_α, num_vw), # v params.
        "w"     => Array{Float64,3}(undef, num_T, num_α, num_vw), # w params.
        "F"     => Matrix{Unitful.Energy}(undef, num_T, num_α), # Energies.
        "A"     => Matrix{Unitful.Energy}(undef, num_T, num_α), # A parameter.
        "B"     => Matrix{Unitful.Energy}(undef, num_T, num_α), # B parameter.
        "C"     => Matrix{Unitful.Energy}(undef, num_T, num_α), # C parameter.
        "κ"     => Array{Any,3}(undef, num_T, num_α, num_vw), # Spring constants.
        "M"     => Array{Unitful.Mass,3}(undef, num_T, num_α, num_vw), # Fictitious masses.
        "Ma"    => Array{Unitful.Mass,3}(undef, num_T, num_α, num_vw),  # Thermal asymptotic approximate fictitious mass (v/w).
        "Mr"    => Array{Unitful.Mass,3}(undef, num_T, num_α, num_vw),  # Thermal reduced mass of particle + fictitious particle system.
        "R"     => Array{Unitful.Length,3}(undef, num_T, num_α, num_vw), # Polaron radii.
        "μ"     => Matrix{Any}(undef, num_T, num_α), # Mobilities.
        "μFHIP" => Matrix{Any}(undef, num_T, num_α), # FHIP low-temperature DC mobiltiy approximation.
        "μD"    => Matrix{Any}(undef, num_T, num_α), # Kadanoff DC mobility, low-temperature approximation using Devreese2016.
        "μK"    => Matrix{Any}(undef, num_T, num_α), # Kadanoff DC mobility, low-temperature approximation.
        "τ"     => Matrix{Any}(undef, num_T, num_α), # relaxation time from Kadanoff Boltzmann transport equation.
        "μH"    => Matrix{Any}(undef, num_T, num_α), # Hellwarth DC mobility.
        "μH0"   => Matrix{Any}(undef, num_T, num_α), # Hellwarth DC mobility with b=0.
        "χ"     => Array{Any,3}(undef, num_Ω, num_T, num_α), # 'Memory function' or self energy.
        "z"     => Array{Any,3}(undef, num_Ω, num_T, num_α), # Complex impedences.
        "σ"     => Array{Any,3}(undef, num_Ω, num_T, num_α), # Complex conductivities.
    )

    # Print the phonon frequencies. 
    if verbose
        # Instantiate IO. Compact and limit reduce printed info.
        io = IOContext(stdout, :compact => true, :limit => true) 
        process = 1     # Counter for total iterations.
        αprocess = 1    # Counter for αrange iterations. 
    end

    for j in axes(αrange, 1)    # αrange loop.

        # Reduce αrange to either a Vector of values per phonon mode, or a single scalar. 
        α = reduce_array(αrange[j, :])

        # Effective alpha parameter is the sum over all contributions per phonon mode.
        αeff = sum(α)

        if verbose
            # Hide cursor whilst printing live data.
            print(io, "\e[?25l") 

            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[K               Polaron Information: [$(αprocess[]) / $(num_α) ($(round(αprocess[] / (num_α) * 100, digits=1)) %)]")
            println(io, "\e[K-----------------------------------------------------------------------")

            # Different formatting for single vs multipelf frequencies (limits an array to head and tail to limit prints).
            if num_ω == 1
                println(io, "\e[KPhonon frequencies             | ω = ", ω)
            else
                println(io, "\e[KPhonon frequencies             | ω = ", join(round.(first(ω, 2), digits=1), ", ")..., " ... ", join(round.(last(ω, 2), digits=1), ", ")...)
            end
        end

        # Print α coupling parameters.
        if verbose
            if num_ω == 1
                println(io, "\e[KFröhlich coupling              | αeff = ", αeff, " | α = ", join(α, ", ")...)
            else
                println(io, "\e[KFröhlich coupling              | αeff = ", αeff, " | α = ", join(round.(first(α, 2), digits=3), ", ")..., " ... ", join(round.(last(α, 2), digits=3), ", ")...)
            end
        end

        # Small alpha (α → 0) approximate energy.
        F_small = (-αeff - αeff^2 / 81) .* E0_pu .|> u"meV"
        p["Fs"][j] = F_small

        # Print small alpha energy.
        if verbose
            println(io, "\e[KSmall α→0 energy               | Fₛ = ", F_small)
        end

        # Large alpha (α → ∞) approximate energy.
        F_large = (-αeff^2 / 3π - 3 * log(2) - 3 / 4) .* E0_pu .|> u"meV"
        p["Fl"][j] = F_large

        # Print large alpha energy.
        if verbose
            println(io, "\e[KLarge α→∞ energy               | Fₗ = ", F_large)
        end

        # Small coupling (α → 0) polaron mass approximation. Eqn. (46) in Feynman1955.
        M_small = (αeff ./ 6 .+ 0.025 .* αeff .^2) .* m0_pu
        p["Ms"][j] = M_small

        # Print small alpha fictitious mass.
        if verbose
            println(io, "\e[KSmall α→0 fictitious mass      | Mₛ = ", M_small)
        end

        # Large coupling (α → ∞) polaron mass approximation. Eqn. (47) In Feynman1955.
        M_large = 16 .* αeff .^4 ./ (81 * π^4) .* m0_pu
        p["Ml"][j] = M_large

        # Print large alpha fictitious mass.
        if verbose
            println(io, "\e[KLarge α→∞ fictitious mass      | Mₗ = ", M_large)
        end

        # Small coupling (α → 0) polaron radii approximiation. Eqn. (2.5a) in Schultz1959.
        R_small = sqrt.(3 ./ (4 ./ 9 .* αeff))  .* a0_pu
        p["Rs"][j] = R_small

        # Print small alpha polaron radius.
        if verbose
            println(io, "\e[KSmall α→0 polaron radius       | Rₛ = ", R_small)
        end

        # Large coupling (α → ∞) polaron radii approximiation. Eqn. (2.5b) in Schultz1959.
        R_large = 3 * √(π / 2) .* αeff * a0_pu
        p["Rl"][j] = R_large

        # Print large alpha polaron radius.
        if verbose
            println(io, "\e[KLarge α→∞ polaron radius       | Rₗ = ", R_large)
        end

        # Franck-Condon (FC) frequency in large coupling (α → ∞) limit. RHS of pg. 2371 in Devreese1972.
        Ω_FC = 4 / 9π .* αeff .^2 * ω0_pu
        p["ΩFC"][j] = Ω_FC

        # Print large alpha Franck-Condon peak frequency.
        if verbose
            println(io, "\e[KLarge α→∞ FC peak freq.        | ΩFC = ", Ω_FC)
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
        p["F0"][j] = F_gs .* E0_pu .|> u"meV"
        p["A0"][j] = A_gs .* E0_pu .|> u"meV"
        p["B0"][j] = B_gs .* E0_pu .|> u"meV"
        p["C0"][j] = C_gs .* E0_pu .|> u"meV"

        # Print athermal variational parameter and energy data.
        if verbose
            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[K                      Zero Temperature Information:                    ")
            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[KVariational parameter          | v₀ = ", v_gs)
            println(io, "\e[KVariational parameter          | w₀ = ", w_gs)
            println(io, "\e[KEnergy                         | E₀ = ", F_gs .* E0_pu .|> u"meV")
            println(io, "\e[KElectron energy                | A₀ = ", A_gs .* E0_pu .|> u"meV")
            println(io, "\e[KInteraction energy             | B₀ = ", B_gs .* E0_pu .|> u"meV")
            println(io, "\e[KTrial energy                   | C₀ = ", C_gs .* E0_pu .|> u"meV")
            Tprocess = 1    # Counter for Trange.
            αprocess += 1   # Increment αrange iteration.
        end

        # Calculate and store fictitious spring constants. See after Eqn. (18), pg. 1007 of Feynman1962. Thermal
        κ_gs = (v_gs .^ 2 .- w_gs .^ 2) .* m0_pu * ω0_pu^2
        p["κ0"][j, :] .= κ_gs

        # Print athermal fictitious spring constant.
        if verbose
            println(io, "\e[KFictitious spring constant     | κ₀ = ", κ_gs)
        end

        # Calculate and store fictitious masses. Athermal.
        M_gs = κ_gs ./ w_gs .^ 2 ./ ω0_pu^2
        p["M0"][j, :] .= M_gs

        # Print athermal fictitious mass.
        if verbose
            println(io, "\e[KFictitious mass                | M₀ = ", M_gs)
        end

        # Approximate large coupling asymptotic limit for the polaron mass. Feynman1962. Athermal
        M_asymp_gs = v_gs ./ w_gs .* m0_pu
        p["M0a"][j, :] .= M_asymp_gs

        # Print athermal asymptotic fictitious mass.
        if verbose
            println(io, "\e[KFictitious mass (asymptotic)   | M₀ₐ = ", M_asymp_gs)
        end

        # Reduced mass of particle and fictitious mass system. Before eqn. (2.4) in Schultz1959. Athermal.
        M_reduced_gs = (v_gs .^2 .- w_gs .^2) / v_gs .^ 2 .* m0_pu
        p["M0r"][j, :] .= M_reduced_gs

        # Print athermal reduced mass.
        if verbose
            println(io, "\e[KReduced mass                   | μ₀ᵣ = ", M_reduced_gs)
        end

        # Calculate and store polaron radii. Approximates the polaron wavefunction as a Gaussian and relates the size to the standard deviation. Eqn. (2.4) in Schultz1959. Athermal.
        R_gs = sqrt.(3 .* v_gs ./ (v_gs .^ 2 .- w_gs .^ 2) .^ 2) .* a0_pu
        p["R0"][j, :] .= R_gs

        # Print athermal polaron radius.
        if verbose
            println(io, "\e[KPolaron radius                 | R₀ = ", R_gs)
        end

        for i in eachindex(Trange)  # Temperatures loop.
            T = Trange[i] * u"K"

            # Print temperature.
            if verbose
                println(io, "\e[K-----------------------------------------------------------------------")
                println(io, "\e[K         Finite Temperature Information: [$(Tprocess[]) / $(num_T) ($(round(Tprocess[] / (num_T) * 100, digits=1)) %)]")
                println(io, "\e[K-----------------------------------------------------------------------")
                println(io, "\e[KTemperatures                   | T = ", T)
            end

            # Calculate the reduced (unitless) thermodynamic betas for each phonon mode.
            β = 1 ./ T .* ħ_pu ./ k_pu
            p["β"][i, :] .= β

            # Print thermodynamic betas.
            if verbose
                println(io, "\e[KReduced thermodynamic          | β = ", β)
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
            p["F"][i, j] = F .* E0_pu .|> u"meV"
            p["A"][i, j] = A .* E0_pu .|> u"meV"
            p["B"][i, j] = B .* E0_pu .|> u"meV"
            p["C"][i, j] = C .* E0_pu .|> u"meV"

            # Print thermal data.
            if verbose
                println(io, "\e[KVariational parameter          | v = ", v)
                println(io, "\e[KVariational parameter          | w = ", w)
                println(io, "\e[KFree energy                    | F = ", F .* E0_pu .|> u"meV")
                println(io, "\e[KElectron energy                | A = ", A .* E0_pu .|> u"meV")
                println(io, "\e[KInteraction energy             | B = ", B .* E0_pu .|> u"meV")
                println(io, "\e[KTrial energy                   | C = ", C .* E0_pu .|> u"meV")
            end

            # Calculate and store fictitious spring constants. See after Eqn. (18), pg. 1007 of Feynman1962. Thermal
            κ = (v .^ 2 .- w .^ 2) .* m0_pu * ω0_pu^2
            p["κ"][i, j, :] .= κ

            # Print spring constants.
            if verbose
                println(io, "\e[KFictitious spring constant     | κ = ", κ)
            end

            # Calculate and store fictitious masses. Thermal.
            M = κ ./ w .^ 2 ./ ω0_pu^2
            p["M"][i, j, :] .= M

            # Print masses.
            if verbose
                println(io, "\e[KFictitious mass                | M = ", M)
            end

            # Approximate large coupling asymptotic limit for the polaron mass. Feynman1962. Thermal
            M_asymp = v ./ w .* m0_pu
            p["Ma"][i, j, :] .= M_asymp

            # Print asymptotic masses.
            if verbose
                println(io, "\e[KFictitious mass (asymptotic)   | Mₐ = ", M_asymp)
            end

            # Reduced mass of particle and fictitious mass system. Before eqn. (2.4) in Schultz1959. Athermal.
            M_reduced = (v .^2 .- w .^2) / v .^ 2 .* m0_pu
            p["Mr"][i, j, :] .= M_reduced

            # Print redcued masses.
            if verbose
                println(io, "\e[KReduced mass                   | μᵣ = ", M_reduced)
            end

            # Calculate and store polaron radii.
            R = sqrt.(3 .* v ./ (v .^ 2 .- w .^ 2) .^ 2) .* a0_pu
            p["R"][i, j, :] .= R

            # Print polaron radius.
            if verbose
                println(io, "\e[KPolaron radius                 | R = ", R)
            end

            if verbose
                println("\e[K-----------------------------------------------------------------------")
                println("\e[K                      DC Mobility Information:                         ")
                println("\e[K-----------------------------------------------------------------------")
            end

            # Calculate and store the DC mobiliies.
            μ = polaron_mobility(v, w, α, ω, β) / m0_pu * e_pu |> u"cm^2/V/s"
            p["μ"][i, j] = μ 

            # Print DC mobilities.
            if verbose
                println(io, "\e[KFinite temperature mobility    | μ = ", μ)
            end

            # FHIP low-temperature mobility, final result of Feynman1962.
            # Eqn. (1.60) in Devreese2016 page 36; 6th Edition of Frohlich polaron notes (ArXiv).
            μ_FHIP = FHIP_mobility_lowT(v, w, α, ω, β) / m0_pu * e_pu |> u"cm^2/V/s"
            p["μFHIP"][i, j] = μ_FHIP

            # Print low-temperature FHIP mobility.
            if verbose
                println(io, "\e[KFHIP low-temp. mobility        | μFHIP = ", μ_FHIP)
            end

            # Kadanoff low-temperaure mobility, constructed around Boltzmann equation.
            # Adds factor of 3 / (2 * β) c.f. FHIP, correcting phonon emission behaviour.
            # Provide also the Kadanoff mobility that is consistent with the
            # FHIP, and later statements (Devreese) of the Kadanoff mobility.
            # It suggests that Kadanoff used the wrong identy for Nbar in Eqn. (23b) for
            # the Γ₀ function, and should have used a version with the -1 to
            # account for Bose / phonon statistics!
            μ_Kadanoff_Devreese, μ_Kadanoff, rel_time = Kadanoff_mobility_lowT(v, w, α, ω, β) ./ m0_pu .* e_pu .|> u"cm^2/V/s"
            p["μD"][i, j] = μ_Kadanoff_Devreese
            p["μK"][i, j] = μ_Kadanoff
            p["τ"][i, j] = rel_time .* m0_pu ./ e_pu .|> u"ns"

            # Print low-temperature Kadanoff mobility (both original and Devreese-corrected versions).
            if verbose
                println(io, "\e[KDevreese low-temp. mobility    | μDev = ", μ_Kadanoff_Devreese)
                println(io, "\e[KKadanoff low-temp. mobility    | μKad = ", μ_Kadanoff)
            end

            #C alculates the DC mobility using Hellwarth et al. 1999 Eqn. (2).
            # Directly performs contour integration in Feynman1962, for finite temperature DC mobility.
            # Eqns. (2) and (1) are going back to the general (pre low-T limit) formulas in Feynman1962.  
            # To evaluate these, you need to do the explicit contour integration to get the polaron self-energy.
            # See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
            μ_Hellwarth, μ_Hellwarth_b0 = Hellwarth_mobility(v, w, α, ω, β) ./ m0_pu .* e_pu ./ ω0_pu^2 .|> u"cm^2/V/s"
            p["μH"][i, j] = μ_Hellwarth
            p["μH0"][i, j] = μ_Hellwarth_b0

            # Print Hellwarth mobility (both with b and b=0) and Kadanoff relaxation time.
            if verbose
                println(io, "\e[KHellwarth mobility             | μHel = ", μ_Hellwarth)
                println(io, "\e[KHellwarth mobility (b=0)       | μHel₀ = ", μ_Hellwarth_b0)
                println(io, "\e[KKadanoff relaxation time       | τ = ", rel_time .* m0_pu ./ e_pu .|> u"ns")
                Ωprocess = 1    # Counter for Ωrange.
                Tprocess += 1   # Increment Trange iterator.
            end

            for k in eachindex(Ωrange)  # E-field frequencies loop. 
                Ω = Ωrange[k] .* ω0_pu

                # Print E-field frequency.
                if verbose
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K         Frequency Reponse Information: [$(Ωprocess[]) / $(num_Ω) ($(round(Ωprocess[] / (num_Ω) * 100, digits=1)) %)]")
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[KElectric field frequency       | Ω = ", Ω)
                end

                # Calculate and store polaron memory functions (akin to self energy).
                χ = polaron_memory_function(v, w, α, ω, β, Ω) .|> u"THz2π" 
                p["χ"][k, i, j] = χ

                # Print memory function.
                if verbose
                    println(io, "\e[KMemory function                | χ = ", χ)
                end

                # Calculate and store polaron complex impedances.

                z = (-im * Ω + im * χ) ./ ω0_pu .* u"Ω" .|> u"Ω"
                p["z"][k, i, j] = z 

                # Print complex impedances.
                if verbose
                    println(io, "\e[KComplex impedance              | z = ", z)
                end

                # Calculate and store polaron complex conductivities.
                σ = 1 / z .|> u"mS"
                p["σ"][k, i, j] = σ 

                # Print complex conductivities and show total algorithm progress.
                if verbose
                    println(io, "\e[KComplex conductivity           | σ = ", σ)
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K[Total Progress: $(process[]) / $(num_α * num_T * num_Ω) ($(round(process[] / (num_α * num_T * num_Ω) * 100, digits=1)) %)]")
                    print(io, "\e[9F")
                    Ωprocess += 1   # Increment Ωrange iterator.
                    process += 1    # Increment total iterator.
                end
            end
            if verbose print(io, "\e[26F")end   # Move up 26 lines and erase.
        end 
        if verbose print(io, "\e[26F") end   # Move up 26 lines and erase. 
    end
    if verbose print("\e[?25h") end # Show cursor again.

    # Extract polaron data as Vector ready for Polaron type.
    polaron_data = [
        p["α"], # Alphas.
        p["αeff"], # Alphas sums.
        p["T"], # Temperatures.
        p["ω"], # Phonon frequencies.
        p["β"], # Betas.
        p["Ω"], # Photon frequencies.
        p["v0"], # Athermal v params.
        p["w0"], # Athermal w params.
        p["F0"], # Athermal energies.
        p["A0"], # Athermal A parameter.
        p["B0"], # Athermal B parameter.
        p["C0"], # Athermal C parameter.
        p["Fs"], # Small alpha (α→0) approximate energy.
        p["Fl"], # Large alpha (α→∞) approximate energy.
        p["κ0"], # Fictitious spring constant.
        p["M0"], # Athermal fictitious mass.
        p["Ms"], # Small alpha (α→0) approximate fictitious mass.
        p["Ml"], # Large alpha (α→∞) approximate fictitious mass.
        p["M0a"], # Athermal asymptotic approximate fictitious mass (v0/w0).
        p["M0r"], # Athermal reduced mass of particle + fictitious particle system.
        p["R0"], # Athermal polaron radius (s. d. of a Gaussian wavefunction).
        p["Rs"], # Small alpha (α→0) approximate polaron radius.
        p["Rl"], # Large alpha (α→∞) approximate polaron radius.
        p["ΩFC"], # Large alpha (α → ∞)approximate Franck-Condon peak frequency.
        p["v"], # v params.
        p["w"], # w params.
        p["F"], # Energies.
        p["A"], # A parameter.
        p["B"], # B parameter.
        p["C"], # C parameter.
        p["κ"], # Spring constants.
        p["M"], # Fictitious masses.
        p["Ma"],  # Thermal asymptotic approximate fictitious mass (v/w).
        p["Mr"],  # Thermal reduced mass of particle + fictitious particle system.
        p["R"], # Polaron radii.
        p["μ"], # Mobilities.
        p["μFHIP"], # FHIP low-temperature DC mobiltiy approximation.
        p["μD"], # Kadanoff DC mobility, low-temperature approximation using Devreese2016.
        p["μK"], # Kadanoff DC mobility, low-temperature approximation.
        p["τ"], # relaxation time from Kadanoff Boltzmann transport equation.
        p["μH"], # Hellwarth DC mobility.
        p["μH0"], # Hellwarth DC mobility with b=0.
        p["χ"], # Memory function or self energy.
        p["z"], # Complex impedences.
        p["σ"], # Complex conductivities.
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
    
    # Extract material data from Material type.
    phonon_freqs = material.freqs

    # Generate polaron data from the arbitrary model constructor.
    p = polaron(material.α', Trange, Ωrange, phonon_freqs .* 2π |> u"THz2π", v_guesses, w_guesses, verbose=verbose)

    # Return material-specific, unitful Polaron type.
    return p
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

# Broadcast Polaron data.
function Base.show(io::IO, ::MIME"text/plain", x::Polaron)

    io_lim = IOContext(io, :limit => true, :compact => true)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                         Polaron Information:                          ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KPhonon frequencies             | ωeff = ", x.ω, " | ω = ", x.ω)
    println(io_lim, "\e[KFröhlich coupling              | αeff = ", x.αeff, " | α = ", x.α)
    println(io_lim, "\e[KSmall α→0 energy               | Fₛ = ", x.Fs)
    println(io_lim, "\e[KLarge α→∞ energy               | Fₗ = ", x.Fl)
    println(io_lim, "\e[KSmall α→0 fictitious mass      | Mₛ = ", x.Ms)
    println(io_lim, "\e[KLarge α→∞ fictitious mass      | Mₗ = ", x.Ml)
    println(io_lim, "\e[KSmall α→0 polaron radius       | Rₛ = ", x.Rs)
    println(io_lim, "\e[KLarge α→∞ polaron radius       | Rₗ = ", x.Rl)
    println(io_lim, "\e[KLarge α→∞ FC peak freq.        | ΩFC = ", x.ΩFC)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                     Zero Temperature Information:                     ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KVariational parameter          | v₀ = ", x.v0)
    println(io_lim, "\e[KVariational parameter          | w₀ = ", x.w0)
    println(io_lim, "\e[KEnergy                         | E₀ = ", x.F0)
    println(io_lim, "\e[KElectron energy                | A₀ = ", x.A0)
    println(io_lim, "\e[KInteraction energy             | B₀ = ", x.B0)
    println(io_lim, "\e[KTrial energy                   | C₀ = ", x.C0)
    println(io_lim, "\e[KFictitious spring constant     | κ₀ = ", x.κ0)
    println(io_lim, "\e[KFictitious mass                | M₀ = ", x.M0)
    println(io_lim, "\e[KFictitious mass (asymptotic)   | M₀ₐ = ", x.Ma)
    println(io_lim, "\e[KReduced mass                   | μ₀ᵣ = ", x.Mr)
    println(io_lim, "\e[KPolaron radius                 | R₀ = ", x.R0)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                    Finite Temperature Information:                    ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KTemperatures                   | T = ", x.T)
    println(io_lim, "\e[KReduced thermodynamic          | β = ", x.β)
    println(io_lim, "\e[KVariational parameter          | v = ", x.v)
    println(io_lim, "\e[KVariational parameter          | w = ", x.w)
    println(io_lim, "\e[KFree energy                    | F = ", x.F)
    println(io_lim, "\e[KElectron energy                | A = ", x.A)
    println(io_lim, "\e[KInteraction energy             | B = ", x.B)
    println(io_lim, "\e[KTrial energy                   | C = ", x.C)
    println(io_lim, "\e[KFictitious spring constant     | κ = ", x.κ)
    println(io_lim, "\e[KFictitious mass                | M = ", x.M)
    println(io_lim, "\e[KFictitious mass (asymptotic)   | Mₐ = ", x.Ma)
    println(io_lim, "\e[KReduced mass                   | μᵣ = ", x.Mr)
    println(io_lim, "\e[KPolaron radius                 | R = ", x.R)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                      DC Mobility Information:                         ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KFinite temperature mobility    | μ = ", x.μ)
    println(io_lim, "\e[KFHIP low-temp. mobility        | μFHIP = ", x.μFHIP)
    println(io_lim, "\e[KDevreese low-temp. mobility    | μDev = ", x.μD)
    println(io_lim, "\e[KKadanoff low-temp. mobility    | μKad = ", x.μK)
    println(io_lim, "\e[KHellwarth mobility             | μHel = ", x.μH)
    println(io_lim, "\e[KHellwarth mobility (b=0)       | μHel₀ = ", x.μH0)
    println(io_lim, "\e[KKadanoff relaxation time       | τ = ", x.τ)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                   Frequency Reponse Information:                      ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KElectric field frequency       | Ω = ", x.Ω)
    println(io_lim, "\e[KMemory function                | χ = ", x.χ)
    println(io_lim, "\e[KComplex impedance              | z = ", x.z)
    println(io_lim, "\e[KComplex conductivity           | σ = ", x.σ)

    println("\e[K-----------------------------------------------------------------------\n")
end

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
