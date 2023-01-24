# Polaron.jl

"""
    Polaron
"""
struct Polaron
    α       # Fröhlich coupling, unitless
    αeff    # Effective Fröhlich coupling summed for multiple modes, unitless
    T       # Temperature, K
    ω       # Phonon mode frequencies, 2π⋅THz
    β       # Reduced thermodynamic beta ħω₀/kBT, unitless
    Ω       # Electric field frequency, 2π⋅THz
    v0      # Variational parameter v, unitless
    w0      # Variational parameter w, unitless
    F0      # Polaron free energy, meV
    A0      # Bare electron free energy, meV
    B0      # ⟨S⟩ₜ interaction energy, meV
    C0      # ⟨Sₜ⟩ₜ free energy of trial system, meV
    v       # Variational parameter v, unitless
    w       # Variational parameter w, unitless
    F       # Polaron free energy, meV
    A       # Bare electron free energy, meV
    B       # ⟨S⟩ₜ interaction energy, meV
    C       # ⟨Sₜ⟩ₜ free energy of trial system, meV
    # Fw      # Weak coupling energy approximation, meV
    # Fs      # Strong coupling energy approximation, meV
    κ       # Fictitious spring constant, unitless (multiples of mₑω₀², electron band-mass and (2π⋅THz)²)
    M       # Fictitious mass, unitless (multiples of mₑ electron band-mass)
    R       # Schultz polaron radius, unitless
    z       # Complex impedence, V/A
    σ       # Complex conductivity, A/V
    μ       # DC mobility, cm²/Vs
    # μK      # Kadanoff DC mobility (cm²/Vs)
    # μH      # Hellwarth DC mobility (cm²/Vs)
    # τ       # relaxation time from Kadanoff Boltzmann transport equation

    function Polaron(x...)
        new(reduce_array.(x)...)
    end
end

"""
    polaron(αrange, Trange, Ωrange, ω, v_guesses, w_guesses; verbose=false)

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

            μ_FHIP = FHIP_mobility_lowT(v, w, α, ω, β)
            μ_Kadanoff = Kadanoff_mobility_lowT(v, w, α, ω, β)
            μ_Hellwarth = Hellwarth_mobility(v, w, α, ω, β)

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

# Broadcast Polaron data.
function Base.show(io::IO, ::MIME"text/plain", x::Polaron)
    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                         Polaron Information:                          ")
    println("\e[K-----------------------------------------------------------------------")

    println(IOContext(io, :limit => true, :compact => true), "\e[KPhonon frequencies         | ω = ", x.ω, " 2π THz")
    println(IOContext(io, :limit => true, :compact => true), "\e[KFröhlich coupling          | α = ", x.α, " | sum(α) = ", x.αeff)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                       Ground State Information:                       ")
    println("\e[K-----------------------------------------------------------------------")

    println(IOContext(io, :limit => true, :compact => true), "\e[KGS variational parameter   | v₀ = ", x.v0, " ω")
    println(IOContext(io, :limit => true, :compact => true), "\e[KGS variational parameter   | w₀ = ", x.w0, " ω")
    println(IOContext(io, :limit => true, :compact => true), "\e[KGS Energy                  | E₀ = ", x.F0, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KGS Electron energy         | A₀ = ", x.A0, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KGS Interaction energy      | B₀ = ", x.B0, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KGS Trial energy            | C₀ = ", x.C0, " meV")

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                    Finite Temperature Information:                    ")
    println("\e[K-----------------------------------------------------------------------")

    println(IOContext(io, :limit => true, :compact => true), "\e[KTemperatures               | T = ", x.T, " K")
    println(IOContext(io, :limit => true, :compact => true), "\e[KReduced thermodynamic      | β = ", x.β)
    println(IOContext(io, :limit => true, :compact => true), "\e[KVariational parameter      | v = ", x.v, " ω")
    println(IOContext(io, :limit => true, :compact => true), "\e[KVariational parameter      | w = ", x.w, " ω")
    println(IOContext(io, :limit => true, :compact => true), "\e[KFree energy                | F = ", x.F, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KElectron energy            | A = ", x.A, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KInteraction energy         | B = ", x.B, " meV")
    println(IOContext(io, :limit => true, :compact => true), "\e[KTrial energy               | C = ", x.C, " meV")

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                       Trial System Information:                       ")
    println("\e[K-----------------------------------------------------------------------")

    println(IOContext(io, :limit => true, :compact => true), "\e[KFictitious spring constant | κ = ", x.κ, " kg/s²")
    println(IOContext(io, :limit => true, :compact => true), "\e[KFictitious mass            | M = ", x.M, " mₑ")
    println(IOContext(io, :limit => true, :compact => true), "\e[KPolaron radius             | R = ", x.R, " rₚ")

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                     Linear Reponse Information:                       ")
    println("\e[K-----------------------------------------------------------------------")

    println(IOContext(io, :limit => true, :compact => true), "\e[KElectric field frequency   | Ω = ", x.Ω, " 2π THz")
    println(IOContext(io, :limit => true, :compact => true), "\e[KComplex impedance          | z = ", x.z, " V/A")
    println(IOContext(io, :limit => true, :compact => true), "\e[KComplex conductivity       | σ = ", x.σ, " A/V")
    println(IOContext(io, :limit => true, :compact => true), "\e[KMobility                   | μ = ", x.μ, " cm²/Vs")

    println("\e[K-----------------------------------------------------------------------")
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
