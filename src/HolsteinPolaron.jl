# Polaron.jl

"""
    Holstein(x...)

Type for storing the Holstein polaron information, `x...`.
"""
mutable struct Holstein
    α       # Holstein coupling.
    αeff    # Effective Holstein coupling summed for multiple modes.
    T       # Temperature.
    ω       # Phonon mode frequencies.
    ωeff    # Effective phonon mode frequency.
    β       # Reduced thermodynamic beta ħω₀/kBT.
    Ω       # Electric field frequency.
    d       # Number of spatial dimensions.       
    v0      # Athermal variational parameter v.
    w0      # Athermal variational parameter w.
    F0      # Polaron athermal total energy (enthalpy).
    K0      # Polaron athermal kinetic energy.
    P0      # Polaron athermal potential energy.
    κ0      # Fictitious spring constant.
    M0      # Athermal fictitious mass.
    M0a     # Athermal asymptotic approximate fictitious mass (v0/w0).
    M0r     # Athermal reduced mass of particle + fictitious particle system.
    R0      # Athermal polaron radius (s. d. of a Gaussian wavefunction).
    v       # Thermal variational parameter v.
    w       # Thermal variational parameter w.
    F       # Polaron free energy.
    K       # Polaron thermal kinetic energy.
    P       # Polaron thermal potential energy.
    κ       # Fictitious spring constant.
    M       # Thermal fictitious mass.
    Ma      # Thermal asymptotic approximate fictitious mass (v/w).
    Mr      # Thermal reduced mass of particle + fictitious particle system.
    R       # Schultz polaron radius.
    μ       # DC mobility.
    χ       # Memory function or polaron self energy.
    z       # Complex impedence.
    σ       # Complex conductivity.

    function Holstein(x...)
        new(reduce_array.(x)...)
    end
end

function holstein(αrange, Trange, Ωrange; ω=1, ωeff=1, mb=1, β0=1, v_guesses=4.0, w_guesses=3.9, dims=3, verbose=false)

    # v_guesses and w_guesses are initial values for v and w (including many v and w parameters).
    # These guesses are generally not needed unless instabilities are found in the minimisation and better initial values improve stability.

    # Get the length of any arrays etc.
    num_α = size(αrange, 1)
    num_T = length(Trange)
    num_Ω = length(Ωrange)
    num_ω = length(ω)

    ω = reduce_array(ω)
    
    # For multiple variational modes, ensure that the number of v and w parameters is the same.
    @assert length(v_guesses) == length(w_guesses) "v and w guesses must be the same length."
    num_vw = length(v_guesses)

    # Instantiate all the polaron data that will go into the Polaron type.
    p = Dict(
        "α"     => αrange, # Alphas.
        "αeff"  => sum(αrange, dims=2), # Alphas sums.
        "T"     => Trange, # Temperatures.
        "ω"     => ω, # Phonon frequencies.
        "ωeff"  => ωeff, # Effective Phonon frequency.
        "β"     => Matrix{Float64}(undef, num_T, num_ω), # Betas.
        "Ω"     => Ωrange, # Photon frequencies.
        "d"     => dims, # Number of spatial dimensions.
        "v0"    => Matrix{Float64}(undef, num_α, num_vw), # Athermal v params.
        "w0"    => Matrix{Float64}(undef, num_α, num_vw), # Athermal w params.
        "F0"    => Vector{Float64}(undef, num_α), # Athermal total energies.
        "K0"    => Vector{Float64}(undef, num_α), # Athermal kinetic energies.
        "P0"    => Vector{Float64}(undef, num_α), # Athermal potential energies.
        "κ0"    => Matrix{Float64}(undef, num_α, num_vw), # Fictitious spring constant.
        "M0"    => Matrix{Float64}(undef, num_α, num_vw), # Athermal fictitious mass.
        "M0a"   => Matrix{Float64}(undef, num_α, num_vw), # Athermal asymptotic approximate fictitious mass (v0/w0).
        "M0r"   => Matrix{Float64}(undef, num_α, num_vw), # Athermal reduced mass of particle + fictitious particle system.
        "R0"    => Matrix{Float64}(undef, num_α, num_vw), # Athermal polaron radius (s. d. of a Gaussian wavefunction).
        "v"     => Array{Float64,3}(undef, num_T, num_α, num_vw), # v params.
        "w"     => Array{Float64,3}(undef, num_T, num_α, num_vw), # w params.
        "F"     => Matrix{Float64}(undef, num_T, num_α), # Total energies.
        "K"     => Matrix{Float64}(undef, num_T, num_α), # Kinetic energies.
        "P"     => Matrix{Float64}(undef, num_T, num_α), # Potential energies.
        "κ"     => Array{Float64,3}(undef, num_T, num_α, num_vw), # Spring constants.
        "M"     => Array{Float64,3}(undef, num_T, num_α, num_vw), # Fictitious masses.
        "Ma"    => Array{Float64,3}(undef, num_T, num_α, num_vw),  # Thermal asymptotic approximate fictitious mass (v/w).
        "Mr"    => Array{Float64,3}(undef, num_T, num_α, num_vw),  # Thermal reduced mass of particle + fictitious particle system.
        "R"     => Array{Float64,3}(undef, num_T, num_α, num_vw), # Polaron radii.
        "μ"     => Matrix{Float64}(undef, num_T, num_α), # Mobilities.
        "χ"     => Array{ComplexF64,3}(undef, num_Ω, num_T, num_α), # 'Memory function' or self energy.
        "z"     => Array{ComplexF64,3}(undef, num_Ω, num_T, num_α), # Complex impedences.
        "σ"     => Array{ComplexF64,3}(undef, num_Ω, num_T, num_α), # Complex conductivities.
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

            # Different formatting for single vs multiple frequencies (limits an array to head and tail to limit prints).
            if num_ω == 1
                println(io, "\e[KPhonon frequencies             | ωeff = ", ωeff, " | ω = ", ω)
                println(io, "\e[KHolstein coupling              | αeff = ", αeff, " | α = ", join(α, ", ")...)
            else
                println(io, "\e[KPhonon frequencies             | ωeff = ", ωeff, " | ω = ", join(round.(first(ω, 2), digits=2), ", ")..., " ... ", join(round.(last(ω, 2), digits=2), ", ")...)
                println(io, "\e[KHolstein coupling              | αeff = ", αeff, " | α = ", join(round.(first(α, 2), digits=3), ", ")..., " ... ", join(round.(last(α, 2), digits=3), ", ")...)
            end
            println(io, "\e[KNumber of dimensions           | d = ", dims)
        end

        # Extract the ground-state, athermal polaron properties (energy (enthalpy) and variational parameters v and w).
        # w is also the frequency of oscillation of the SHM trial system composed of the bare particle and fictitous mass.
        # A, B, C are components of the total energy: A is the bare electron energy, B the electron-phonon interaction energy, C is the energy of the harmonic trial system.
        athermal_energy(v, w) = holstein_energy(v, w, α, ω; dims = dims)
        v_gs, w_gs, F_gs, K_gs, P_gs = vw_variation(athermal_energy, v_guesses, w_guesses)
 
        # Update the guesses to keep them close-ish to the true solutions during loops over alphas.
        v_guesses, w_guesses = v_gs, w_gs

        # Store the athermal data.
        p["v0"][j, :] .= v_gs
        p["w0"][j, :] .= w_gs
        p["F0"][j] = F_gs
        p["K0"][j] = K_gs
        p["P0"][j] = P_gs

        # Print athermal variational parameter and energy data.
        if verbose
            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[K                   Zero Temperature Information:                       ")
            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[KVariational parameter          | v0 = ", v_gs)
            println(io, "\e[KVariational parameter          | w0 = ", w_gs)
            println(io, "\e[KEnergy                         | F0 = ", F_gs)
            println(io, "\e[KKinetic energy                 | K0 = ", K_gs)
            println(io, "\e[KPotential energy               | P0 = ", P_gs)
            Tprocess = 1    # Counter for Trange.
            αprocess += 1   # Increment αrange iteration.
        end

        # Calculate and store fictitious spring constants. See after Eqn. (18), pg. 1007 of Feynman1962. Thermal
        κ_gs = (v_gs .^ 2 .- w_gs .^ 2)
        p["κ0"][j, :] .= κ_gs

        # Print athermal fictitious spring constant.
        if verbose
            println(io, "\e[KFictitious spring constant     | κ0 = ", κ_gs)
        end

        # Calculate and store fictitious masses. Athermal.
        M_gs = κ_gs ./ w_gs .^ 2
        p["M0"][j, :] .= M_gs

        # Print athermal fictitious mass.
        if verbose
            println(io, "\e[KFictitious mass                | M0 = ", M_gs)
        end

        # Approximate large coupling asymptotic limit for the polaron mass. Feynman1962. Athermal
        M_asymp_gs = v_gs ./ w_gs
        p["M0a"][j, :] .= M_asymp_gs

        # Print athermal asymptotic fictitious mass.
        if verbose
            println(io, "\e[KFictitious mass (asymptotic)   | M0a = ", M_asymp_gs)
        end

        # Reduced mass of particle and fictitious mass system. Before eqn. (2.4) in Schultz1959. Athermal.
        M_reduced_gs = (v_gs .^2 .- w_gs .^2) / v_gs .^ 2
        p["M0r"][j, :] .= M_reduced_gs

        # Print athermal reduced mass.
        if verbose
            println(io, "\e[KReduced mass                   | M0r = ", M_reduced_gs)
        end

        # Calculate and store polaron radii. Approximates the polaron wavefunction as a Gaussian and relates the size to the standard deviation. Eqn. (2.4) in Schultz1959. Athermal.
        R_gs = sqrt.(3 .* v_gs ./ (v_gs .^ 2 .- w_gs .^ 2) .^ 2)
        p["R0"][j, :] .= R_gs

        # Print athermal polaron radius.
        if verbose
            println(io, "\e[KPolaron radius                 | R0 = ", R_gs)
        end

        for i in eachindex(Trange)  # Temperatures loop.
            T = Trange[i]

            # Print temperature.
            if verbose
                println(io, "\e[K-----------------------------------------------------------------------")
                println(io, "\e[K         Finite Temperature Information: [$(Tprocess[]) / $(num_T) ($(round(Tprocess[] / (num_T) * 100, digits=1)) %)]")
                println(io, "\e[K-----------------------------------------------------------------------")
                println(io, "\e[KTemperatures                   | T = ", T)
            end

            # Calculate the reduced (unitless) thermodynamic betas for each phonon mode.
            β = 1 ./ T .* β0
            p["β"][i, :] .= β

            # Print thermodynamic betas.
            if verbose
                println(io, "\e[KReduced thermodynamic          | β = ", β)
            end

            # Calculate thermal polaron properties (energy (Gibbs free energy) and variational parameters v and w).
            # w is also the frequency of oscillation of the SHM trial system composed of the bare particle and fictitous mass.
            # A, B, C are components of the total energy: A is the bare electron energy, B the electron-phonon interaction energy, C is the energy of the harmonic trial system.
            thermal_energy(v, w) = holstein_energy(v, w, α, ω, β; dims = dims) 
            v, w, F, K, P = vw_variation(thermal_energy, v_guesses, w_guesses)

            # Update the guesses to keep them close-ish to the true solutions during loops over temperatures.
            v_guesses, w_guesses = v, w

            # Store thermal data.
            p["v"][i, j, :] .= v
            p["w"][i, j, :] .= w
            p["F"][i, j] = F
            p["K"][i, j] = K
            p["P"][i, j] = P

            # Print thermal data.
            if verbose
                println(io, "\e[KVariational parameter          | v = ", v)
                println(io, "\e[KVariational parameter          | w = ", w)
                println(io, "\e[KFree energy                    | F = ", F)
                println(io, "\e[KKinetic energy                 | K = ", K)
                println(io, "\e[KPotential energy               | P = ", P)
            end

            # Calculate and store fictitious spring constants. See after Eqn. (18), pg. 1007 of Feynman1962. Thermal
            κ = (v .^ 2 .- w .^ 2)
            p["κ"][i, j, :] .= κ

            # Print spring constants.
            if verbose
                println(io, "\e[KFictitious spring constant     | κ = ", κ)
            end

            # Calculate and store fictitious masses. Thermal.
            M = κ ./ w .^ 2
            p["M"][i, j, :] .= M

            # Print masses.
            if verbose
                println(io, "\e[KFictitious mass                | M = ", M)
            end

            # Approximate large coupling asymptotic limit for the polaron mass. Feynman1962. Thermal
            M_asymp = v ./ w 
            p["Ma"][i, j, :] .= M_asymp

            # Print asymptotic masses.
            if verbose
                println(io, "\e[KFictitious mass (asymptotic)   | Ma = ", M_asymp)
            end

            # Reduced mass of particle and fictitious mass system. Before eqn. (2.4) in Schultz1959. Athermal.
            M_reduced = (v .^2 .- w .^2) / v .^ 2
            p["Mr"][i, j, :] .= M_reduced

            # Print redcued masses.
            if verbose
                println(io, "\e[KReduced mass                   | Mr = ", M_reduced)
            end

            # Calculate and store polaron radii.
            R = sqrt.(3 .* v ./ (v .^ 2 .- w .^ 2) .^ 2)
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
            μ = holstein_mobility(v, w, α, ω, β; dims = dims) / mb
            p["μ"][i, j] = μ 

            # Print DC mobilities.
            if verbose
                println(io, "\e[KFinite temperature mobility    | μ = ", μ)
                Ωprocess = 1    # Counter for Ωrange.
                Tprocess += 1   # Increment Trange iterator.
            end

            for k in eachindex(Ωrange)  # E-field frequencies loop. 
                Ω = Ωrange[k] 

                # Print E-field frequency.
                if verbose
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K         Frequency Response Information: [$(Ωprocess[]) / $(num_Ω) ($(round(Ωprocess[] / (num_Ω) * 100, digits=1)) %)]")
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[KElectric field frequency       | Ω = ", Ω)
                end

                # Calculate and store polaron memory functions (akin to self energy).
                χ = holstein_memory_function(Ω, v, w, α, ω, β; dims = dims)
                p["χ"][k, i, j] = χ

                # Print memory function.
                if verbose
                    println(io, "\e[KMemory function                | χ = ", χ)
                end

                # Calculate and store polaron complex impedances.

                z = -(im * Ω + im * χ) .* mb
                p["z"][k, i, j] = z 

                # Print complex impedances.
                if verbose
                    println(io, "\e[KComplex impedance              | z = ", z)
                end

                # Calculate and store polaron complex conductivities.
                σ = 1 / z
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
            if verbose print(io, "\e[19F")end   # Move up 19 lines and erase.
        end 
        if verbose print(io, "\e[19F") end   # Move up 19 lines and erase. 
    end
    if verbose print("\e[?25h") end # Show cursor again.

    # Extract polaron data as Vector ready for Polaron type.
    polaron_data = [
        p["α"], # Alphas.
        p["αeff"], # Alphas sums.
        p["T"], # Temperatures.
        p["ω"], # Phonon frequencies.
        p["ωeff"], # Hellwarth eff phonon frequency.
        p["β"], # Betas.
        p["Ω"], # Photon frequencies.
        p["d"], # Number of spatial dimensions.
        p["v0"], # Athermal v params.
        p["w0"], # Athermal w params.
        p["F0"], # Athermal energies.
        p["K0"], # Athermal A parameter.
        p["P0"], # Athermal B parameter.
        p["κ0"], # Fictitious spring constant.
        p["M0"], # Athermal fictitious mass.
        p["M0a"], # Athermal asymptotic approximate fictitious mass (v0/w0).
        p["M0r"], # Athermal reduced mass of particle + fictitious particle system.
        p["R0"], # Athermal polaron radius (s. d. of a Gaussian wavefunction).
        p["v"], # v params.
        p["w"], # w params.
        p["F"], # Energies.
        p["K"], # A parameter.
        p["P"], # B parameter.
        p["κ"], # Spring constants.
        p["M"], # Fictitious masses.
        p["Ma"],  # Thermal asymptotic approximate fictitious mass (v/w).
        p["Mr"],  # Thermal reduced mass of particle + fictitious particle system.
        p["R"], # Polaron radii.
        p["μ"], # Mobilities.
        p["χ"], # Memory function or self energy.
        p["z"], # Complex impedences.
        p["σ"], # Complex conductivities.
    ]

    # Return Polaron type containing all generated data over any coupling strengths, temperatures and frequencies.
    return Holstein(polaron_data...)
end

"""
Single alpha parameter. holstein() expects alpha parameters to be in a Vector.
"""
holstein(α::Real, Trange, Ωrange; ω=1, ωeff=1, mb=1, β0=1, v_guesses=3.11, w_guesses=2.87, dims=3, verbose=false) = holstein([α], Trange, Ωrange; ω=ω, ωeff=ωeff, mb=mb, β0=β0, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, verbose=verbose)

"""
No frequency input => Phonon peak.
"""
holstein(αrange, Trange; ω=1, ωeff=1, mb=1, β0=1, v_guesses=3.11, w_guesses=2.87, dims=3, verbose=false) = holstein(αrange, Trange, ω; ω=ω, ωeff=ωeff, mb=mb, β0=β0, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, verbose=verbose)

"""
No temperature input => 300 K.
"""
holstein(αrange; ω=1, ωeff=1, mb=1, β0=1, v_guesses=3.11, w_guesses=2.87, dims=3, verbose=false) = holstein(αrange, 300, ω; ω=ω, ωeff=ωeff, mb=mb, β0=β0, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, verbose=verbose)

"""
No input => α = 1
"""
holstein(; ω=1, ωeff=1, mb=1, β0=1, v_guesses=3.11, w_guesses=2.87, dims=3, verbose=false) = holstein(1, 300, ω; ω=ω, ωeff=ωeff, mb=mb, β0=β0, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, verbose=verbose)

"""
    holstein(material::Material, TΩrange...; v_guesses=3.11, w_guesses=2.87, verbose=false)
Material specific constructors that use material specific parameters to parameterise the polaron.
Material data is inputted through the `Material` type.
Returns all data in either SI units or other common, suitable units otherwise.
"""
function holstein(material::Material, TΩrange...; v_guesses=3.11, w_guesses=2.87, dims=3, verbose=false)

    # Show material data.
    if verbose
        display(material)
    end
    
    # Extract material data from Material type.
    phonon_freqs = material.f
    phonon_eff_freq = material.feff
    mb = material.mb

    # Generate Holstein polaron data from the arbitrary model constructor.
    p = holstein(material.α', TΩrange...; ω=phonon_freqs, ωeff=phonon_eff_freq, mb=mb, β0=ħ/kB*1e12*2π, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, verbose=verbose)

    # Return material-specific, unitful Holstein type.
    return p
end

# Broadcast Polaron data.
function Base.show(io::IO, ::MIME"text/plain", x::Holstein)

    io_lim = IOContext(io, :limit => true, :compact => true)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                         Polaron Information:                          ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KPhonon frequencies             | ωeff = ", x.ωeff, " | ω = ", x.ω)
    println(io_lim, "\e[KHolstein coupling              | αeff = ", x.αeff, " | α = ", x.α)
    println(io_lim, "\e[KNumber of spatial dimensions   | d = ", x.d)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                     Zero Temperature Information:                     ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KVariational parameter          | v0 = ", x.v0)
    println(io_lim, "\e[KVariational parameter          | w0 = ", x.w0)
    println(io_lim, "\e[KTotal energy                   | F0 = ", x.F0)
    println(io_lim, "\e[KKinetic energy                 | K0 = ", x.K0)
    println(io_lim, "\e[KPotential energy               | P0 = ", x.P0)
    println(io_lim, "\e[KFictitious spring constant     | κ0 = ", x.κ0)
    println(io_lim, "\e[KFictitious mass                | M0 = ", x.M0)
    println(io_lim, "\e[KFictitious mass (asymptotic)   | M0a = ", x.M0a)
    println(io_lim, "\e[KReduced mass                   | M0r = ", x.M0r)
    println(io_lim, "\e[KPolaron radius                 | R0 = ", x.R0)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                    Finite Temperature Information:                    ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KTemperatures                   | T = ", x.T)
    println(io_lim, "\e[KReduced thermodynamic          | β = ", x.β)
    println(io_lim, "\e[KVariational parameter          | v = ", x.v)
    println(io_lim, "\e[KVariational parameter          | w = ", x.w)
    println(io_lim, "\e[KFree energy                    | F = ", x.F)
    println(io_lim, "\e[KKinetic energy                 | K = ", x.K)
    println(io_lim, "\e[KPotenital energy               | P = ", x.P)
    println(io_lim, "\e[KFictitious spring constant     | κ = ", x.κ)
    println(io_lim, "\e[KFictitious mass                | M = ", x.M)
    println(io_lim, "\e[KFictitious mass (asymptotic)   | Ma = ", x.Ma)
    println(io_lim, "\e[KReduced mass                   | Mr = ", x.Mr)
    println(io_lim, "\e[KPolaron radius                 | R = ", x.R)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                      DC Mobility Information:                         ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KFinite temperature mobility    | μ = ", x.μ)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                  Frequency Response Information:                      ")
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
function save_holstein_polaron(polaron::Holstein, prefix)

    println("Saving polaron data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "alpha", polaron.α,
        "alpha eff", polaron.αeff,
        "temperature", polaron.T,
        "phonon freq", polaron.ω,
        "phonon freq eff", polaron.ωeff,
        "beta", polaron.β,
        "E-field freq", polaron.Ω,
        "dims", polaron.d,
        "v0", polaron.v0,
        "w0", polaron.w0,
        "athermal energy", polaron.F0,
        "athermal K", polaron.K0,
        "athermal P", polaron.P0,
        "athermal spring", polaron.κ0,
        "athermal mass", polaron.M0,
        "athermal asympt mass", polaron.M0a,
        "athermal reduced mass", polaron.M0r,
        "athermal size", polaron.R0,
        "v", polaron.v,
        "w", polaron.w,
        "thermal energy", polaron.F,
        "thermal K", polaron.K,
        "thermal P", polaron.P,
        "thermal spring", polaron.κ,
        "thermal mass", polaron.M,
        "thermal asympt mass", polaron.Ma,
        "thermal reduced mass", polaron.Mr,
        "thermal size", polaron.R,
        "mobility", polaron.μ,
        "memory function", polaron.χ,
        "impedance", polaron.z,
        "conductivity", polaron.σ
    )

    println("... Polaron data saved.")
end

"""
    load_polaron(p::NewPolaron, prefix)

Loads data from file "polaron_file_path" into a NewPolaron type.
"""
function load_holstein_polaron(polaron_file_path)

    println("Loading polaron data from $polaron_file_path ...")

    data = JLD.load("$polaron_file_path")

    holstein_polaron = Holstein(
        data["alpha"],
        data["alpha eff"],
        data["temperature"],
        data["phonon freq"],
        data["phonon freq eff"],
        data["beta"],
        data["E-field freq"],
        data["dims"],
        data["v0"],
        data["w0"],
        data["athermal energy"],
        data["athermal K"],
        data["athermal P"],
        data["athermal spring"],
        data["athermal mass"],
        data["athermal asympt mass"],
        data["athermal reduced mass"],
        data["athermal size"],
        data["v"],
        data["w"],
        data["thermal energy"],
        data["thermal K"],
        data["thermal P"],
        data["thermal spring"],
        data["thermal mass"],
        data["thermal asympt mass"],
        data["thermal reduced mass"],
        data["thermal size"],
        data["mobility"],
        data["memory function"],
        data["impedance"],
        data["conductivity"]
    )
    println("... Polaron loaded.")

    return holstein_polaron
end
