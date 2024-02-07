# HolsteinPolaron.jl
# A do-it-all function for calculating properties of the Holstein polaron. Data is stored as a struct.

"""
    HolsteinPolaron(x...)

Type for storing the Holstein polaron information, `x...`.
"""
mutable struct HolsteinPolaron
    α       # Holstein coupling.
    αeff    # Effective Holstein coupling summed for multiple modes.
    mb      # Band mass.
    T       # Temperature.
    ω       # Phonon mode frequencies.
    ωeff    # Effective phonon mode frequency.
    γ       # Adiabaticity
    J       # Transfer integral.
    a       # Lattice constant.
    β       # Reduced thermodynamic beta ħω₀/kBT.
    Ω       # Electric field frequency.
    d       # Number of spatial dimensions.       
    v0      # Athermal variational parameter v.
    w0      # Athermal variational parameter w.
    F0      # Polaron athermal total energy (enthalpy).
    A0      # Bare electron enthalpy.
    B0      # ⟨S⟩ₜ interaction enthalpy.
    C0      # ⟨Sₜ⟩ₜ enthalpy of trial system.
    κ0      # Fictitious spring constant.
    M0      # Athermal fictitious mass.
    M0a     # Athermal asymptotic approximate fictitious mass (v0/w0).
    M0r     # Athermal reduced mass of particle + fictitious particle system.
    R0      # Athermal polaron radius (s. d. of a Gaussian wavefunction).
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
    χ       # Memory function or polaron self energy.
    z       # Complex impedence.
    σ       # Complex conductivity.

    function HolsteinPolaron(x...)
        new(reduce_array.(x)...)
    end
end

function holsteinpolaron(αrange, Trange, Ωrange; ω=1, ωeff=1, γ=1, mb=1, J=1, a=1, v_guesses=4.0, w_guesses=3.9, dims=3, kspace=false, verbose=false)

    # v_guesses and w_guesses are initial values for v and w (including many v and w parameters).
    # These guesses are generally not needed unless instabilities are found in the minimisation and better initial values improve stability.

    # Get the length of any arrays etc.
    num_α = size(αrange, 1)
    num_T = length(Trange)
    num_Ω = length(Ωrange)
    num_ω = length(ω)
    num_d = length(dims)

    ω = reduce_array(ω)
    
    # For multiple variational modes, ensure that the number of v and w parameters is the same.
    @assert length(v_guesses) == length(w_guesses) "v and w guesses must be the same length."
    num_vw = length(v_guesses)

    # Instantiate all the polaron data that will go into the Polaron type.
    p = Dict(
        "α"     => αrange, # Alphas.
        "αeff"  => sum(αrange, dims=2), # Alphas sums.
        "mb"    => mb, # Band mass
        "T"     => Trange, # Temperatures.
        "ω"     => ω, # Phonon frequencies.
        "ωeff"  => ωeff, # Effective Phonon frequency.
        "γ"     => γ, # Adiabaticity
        "J"     => J, # Transfer Integral.
        "a"     => a, # Lattice constant.
        "β"     => Vector{Float64}(undef, num_T), # Betas.
        "Ω"     => Ωrange, # Photon frequencies.
        "d"     => dims, # Number of spatial dimensions.
        "v0"    => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal v params.
        "w0"    => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal w params.
        "F0"    => Matrix{Float64}(undef, num_d, num_α), # Athermal total energies.
        "A0"    => Matrix{Float64}(undef, num_d, num_α), # Athermal A parameter.
        "B0"    => Matrix{Float64}(undef, num_d, num_α), # Athermal B parameter.
        "C0"    => Matrix{Float64}(undef, num_d, num_α), # Athermal C parameter.
        "κ0"    => Array{Float64,3}(undef, num_d, num_α, num_vw), # Fictitious spring constant.
        "M0"    => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal fictitious mass.
        "M0a"   => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal asymptotic approximate fictitious mass (v0/w0).
        "M0r"   => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal reduced mass of particle + fictitious particle system.
        "R0"    => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal polaron radius (s. d. of a Gaussian wavefunction).
        "v"     => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw), # v params.
        "w"     => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw), # w params.
        "F"     => Array{Float64,3}(undef, num_T, num_d, num_α), # Total energies.
        "A"     => Array{Float64,3}(undef, num_T, num_d, num_α), # A parameter.
        "B"     => Array{Float64,3}(undef, num_T, num_d, num_α), # B parameter.
        "C"     => Array{Float64,3}(undef, num_T, num_d, num_α), # C parameter.
        "κ"     => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw), # Spring constants.
        "M"     => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw), # Fictitious masses.
        "Ma"    => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw),  # Thermal asymptotic approximate fictitious mass (v/w).
        "Mr"    => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw),  # Thermal reduced mass of particle + fictitious particle system.
        "R"     => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw), # Polaron radii.
        "μ"     => Array{Float64,3}(undef, num_T, num_d, num_α), # Mobilities.
        "χ"     => Array{ComplexF64,4}(undef, num_Ω, num_T, num_d, num_α), # 'Memory function' or self energy.
        "z"     => Array{ComplexF64,4}(undef, num_Ω, num_T, num_d, num_α), # Complex impedences.
        "σ"     => Array{ComplexF64,4}(undef, num_Ω, num_T, num_d, num_α), # Complex conductivities.
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
            dprocess = 1    # Counter for drange iterations.

            # Hide cursor whilst printing live data.
            print(io, "\e[?25l") 

            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[K               Polaron Information: [$(αprocess[]) / $(num_α) ($(round(αprocess[] / (num_α) * 100, digits=1)) %)]")
            println(io, "\e[K-----------------------------------------------------------------------")

            αprocess += 1   # Increment αrange iteration.

            # Different formatting for single vs multiple frequencies (limits an array to head and tail to limit prints).
            if num_ω == 1
                println(io, "\e[KPhonon frequencies             | ωeff = ", ωeff, " | ω = ", ω)
                println(io, "\e[KHolstein coupling              | αeff = ", αeff, " | α = ", join(α, ", ")...)
            else
                println(io, "\e[KPhonon frequencies             | ωeff = ", ωeff, " | ω = ", join(round.(first(ω, 2), digits=2), ", ")..., " ... ", join(round.(last(ω, 2), digits=2), ", ")...)
                println(io, "\e[KHolstein coupling              | αeff = ", αeff, " | α = ", join(round.(first(α, 2), digits=3), ", ")..., " ... ", join(round.(last(α, 2), digits=3), ", ")...)
            end
        end

        for d in 1:num_d
            if verbose
                println(io, "\e[KNumber of dimensions", " [", dprocess, " / ", num_d, "]", "   | d = ", dims[d])
                dprocess += 1
            end

        # Extract the ground-state, athermal polaron properties (energy (enthalpy) and variational parameters v and w).
        # w is also the frequency of oscillation of the SHM trial system composed of the bare particle and fictitous mass.
        # A, B, C are components of the total energy: A is the bare electron energy, B the electron-phonon interaction energy, C is the energy of the harmonic trial system.
        athermal_energy(v, w) = !kspace ? holstein_energy(v, w, α, ω; dims = dims[d]) : holstein_energy_k_space(v, w, α, ω; dims = dims[d])
        v_gs, w_gs, F_gs, A_gs, B_gs, C_gs = vw_variation(athermal_energy, v_guesses, w_guesses)
 
        # Update the guesses to keep them close-ish to the true solutions during loops over alphas.
        v_guesses, w_guesses = v_gs, w_gs

        # Store the athermal data.
        p["v0"][d, j, :] .= v_gs ./ ωeff
        p["w0"][d, j, :] .= w_gs ./ ωeff
        p["F0"][d, j] = F_gs
        p["A0"][d, j] = A_gs
        p["B0"][d, j] = B_gs
        p["C0"][d, j] = C_gs

        # Print athermal variational parameter and energy data.
        if verbose
            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[K                   Zero Temperature Information:                       ")
            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[KVariational parameter          | v0 = ", v_gs ./ ωeff)
            println(io, "\e[KVariational parameter          | w0 = ", w_gs ./ ωeff)
            println(io, "\e[KEnergy                         | F0 = ", F_gs)
            println(io, "\e[KElectron energy                | A0 = ", A_gs)
            println(io, "\e[KInteraction energy             | B0 = ", B_gs)
            println(io, "\e[KTrial energy                   | C0 = ", C_gs)
            Tprocess = 1    # Counter for Trange.
        end

        # Calculate and store fictitious spring constants. See after Eqn. (18), pg. 1007 of Feynman1962. Thermal
        κ_gs = (v_gs .^ 2 .- w_gs .^ 2) 
        p["κ0"][d, j, :] .= κ_gs ./ ωeff .^2

        # Print athermal fictitious spring constant.
        if verbose
            println(io, "\e[KFictitious spring constant     | κ0 = ", κ_gs)
        end

        # Calculate and store fictitious masses. Athermal.
        M_gs = κ_gs ./ w_gs .^ 2
        p["M0"][d, j, :] .= M_gs

        # Print athermal fictitious mass.
        if verbose
            println(io, "\e[KFictitious mass                | M0 = ", M_gs)
        end

        # Approximate large coupling asymptotic limit for the polaron mass. Feynman1962. Athermal
        M_asymp_gs = v_gs ./ w_gs
        p["M0a"][d, j, :] .= M_asymp_gs

        # Print athermal asymptotic fictitious mass.
        if verbose
            println(io, "\e[KFictitious mass (asymptotic)   | M0a = ", M_asymp_gs)
        end

        # Reduced mass of particle and fictitious mass system. Before eqn. (2.4) in Schultz1959. Athermal.
        M_reduced_gs = (v_gs .^2 .- w_gs .^2) / v_gs .^ 2
        p["M0r"][d, j, :] .= M_reduced_gs

        # Print athermal reduced mass.
        if verbose
            println(io, "\e[KReduced mass                   | M0r = ", M_reduced_gs)
        end

        # Calculate and store polaron radii. Approximates the polaron wavefunction as a Gaussian and relates the size to the standard deviation. Eqn. (2.4) in Schultz1959. Athermal.
        R_gs = sqrt.(3 ./ 2 .* M_reduced_gs .* v_gs)
        p["R0"][d, j, :] .= R_gs ./ ωeff .^(1/2)

        # Print athermal polaron radius.
        if verbose
            println(io, "\e[KPolaron radius                 | R0 = ", R_gs)
        end

        for i in eachindex(Trange)  # Temperatures loop.
            T = Trange[i]

            if !iszero(T)

                # Print temperature.
                if verbose
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K         Finite Temperature Information: [$(Tprocess[]) / $(num_T) ($(round(Tprocess[] / (num_T) * 100, digits=1)) %)]")
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[KTemperatures                   | T = ", T)
                end

                # Calculate the reduced (unitless) thermodynamic betas for each phonon mode.
                β = 1 ./ T
                p["β"][i] = β

                # Print thermodynamic betas.
                if verbose
                    println(io, "\e[KReduced thermodynamic          | β = ", β)
                end

                # Calculate thermal polaron properties (energy (Gibbs free energy) and variational parameters v and w).
                # w is also the frequency of oscillation of the SHM trial system composed of the bare particle and fictitous mass.
                # A, B, C are components of the total energy: A is the bare electron energy, B the electron-phonon interaction energy, C is the energy of the harmonic trial system.
                thermal_energy(v, w) = !kspace ? holstein_energy(v, w, α, ω, β; dims = dims[d]) : holstein_energy_k_space(v, w, α, ω, β; dims = dims[d])
                v, w, F, A, B, C = vw_variation(thermal_energy, v_guesses, w_guesses)

                # Update the guesses to keep them close-ish to the true solutions during loops over temperatures.
                v_guesses, w_guesses = v, w

                # Store thermal data.
                p["v"][i, d, j, :] .= v ./ ωeff
                p["w"][i, d, j, :] .= w ./ ωeff
                p["F"][i, d, j] = F
                p["A"][i, d, j] = A
                p["B"][i, d, j] = B
                p["C"][i, d, j] = C

                # Print thermal data.
                if verbose
                    println(io, "\e[KVariational parameter          | v = ", v ./ ωeff)
                    println(io, "\e[KVariational parameter          | w = ", w ./ ωeff)
                    println(io, "\e[KFree energy                    | F = ", F)
                    println(io, "\e[KElectron energy                | A = ", A)
                    println(io, "\e[KInteraction energy             | B = ", B)
                    println(io, "\e[KTrial energy                   | C = ", C)
                end

                # Calculate and store fictitious spring constants. See after Eqn. (18), pg. 1007 of Feynman1962. Thermal
                κ = (v .^ 2 .- w .^ 2) 
                p["κ"][i, d, j, :] .= κ ./ ωeff .^2

                # Print spring constants.
                if verbose
                    println(io, "\e[KFictitious spring constant     | κ = ", κ)
                end

                # Calculate and store fictitious masses. Thermal.
                M = κ ./ w .^ 2
                p["M"][i, d, j, :] .= M

                # Print masses.
                if verbose
                    println(io, "\e[KFictitious mass                | M = ", M)
                end

                # Approximate large coupling asymptotic limit for the polaron mass. Feynman1962. Thermal
                M_asymp = v ./ w 
                p["Ma"][i, d, j, :] .= M_asymp

                # Print asymptotic masses.
                if verbose
                    println(io, "\e[KFictitious mass (asymptotic)   | Ma = ", M_asymp)
                end

                # Reduced mass of particle and fictitious mass system. Before eqn. (2.4) in Schultz1959. Athermal.
                M_reduced = (v .^2 .- w .^2) / v .^ 2 
                p["Mr"][i, d, j, :] .= M_reduced

                # Print redcued masses.
                if verbose
                    println(io, "\e[KReduced mass                   | Mr = ", M_reduced)
                end

                # Calculate and store polaron radii.
                R = sqrt.(3 ./ 2 .* v .* M_reduced) 
                p["R"][i, d, j, :] .= R ./ ωeff .^(1/2)

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
                μ = !kspace ? holstein_mobility(v, w, α, ω, β; dims = dims[d]) : holstein_mobility_k_space(v, w, α, ω, β; dims = dims[d])
                p["μ"][i, d, j] = μ 

                # Print DC mobilities.
                if verbose
                    println(io, "\e[KFinite temperature mobility    | μ = ", μ)
                end

            else
                # If zero temperature.
                v, w, β = v_gs, w_gs, Inf
                p["β"][i, :] .= β
                p["v"][i, d, j, :] .= v_gs ./ ωeff
                p["w"][i, d, j, :] .= w_gs ./ ωeff
                p["F"][i, d, j] = F_gs
                p["A"][i, d, j] = A_gs
                p["B"][i, d, j] = B_gs
                p["C"][i, d, j] = C_gs
                p["κ"][i, d, j, :] .= κ_gs
                p["M"][i, d, j, :] .= M_gs
                p["μ"][i, d, j] = Inf
                p["R"][i, d, j, :] .= R_gs
                p["Mr"][i, d, j, :] .= M_reduced_gs
                p["Ma"][i, d, j, :] .= M_asymp_gs
            end # End of zero temperature if statement.

            if verbose
                Ωprocess = 1    # Counter for Ωrange.
                Tprocess += 1   # Increment Trange iterator.
            end

            for k in eachindex(Ωrange)  # E-field frequencies loop. 
                Ω = Ωrange[k] 

                if !iszero(T) || !iszero(Ω)

                # Print E-field frequency.
                if verbose
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K         Frequency Response Information: [$(Ωprocess[]) / $(num_Ω) ($(round(Ωprocess[] / (num_Ω) * 100, digits=1)) %)]")
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[KElectric field frequency       | Ω = ", Ω)
                end

                # Calculate and store polaron memory functions (akin to self energy).
                χ = !kspace ? holstein_memory_function(Ω, v, w, α, ω, β; dims = dims[d]) : holstein_memory_function_k_space(Ω, v, w, α, ω, β; dims = dims[d])
                p["χ"][k, i, d, j] = χ

                # Print memory function.
                if verbose
                    println(io, "\e[KMemory function                | χ = ", χ)
                end

                # Calculate and store polaron complex impedances.

                z = -(im * Ω + im * χ) 
                p["z"][k, i, d, j] = z 

                # Print complex impedances.
                if verbose
                    println(io, "\e[KComplex impedance              | z = ", z)
                end

                # Calculate and store polaron complex conductivities.
                σ = 1 / z
                p["σ"][k, i, d, j] = σ 

                # Print complex conductivities and show total algorithm progress.
                if verbose
                    println(io, "\e[KComplex conductivity           | σ = ", σ)
                end

                else

                    # If zero frequency.
                    p["χ"][k, i, d, j] = Inf + 0 * im
                    p["z"][k, i, d, j] = 0 + im * Inf
                    p["σ"][k, i, d, j] = 0 + 0 * im

                end # End of zero frequency if statement.

                if verbose
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K[Total Progress: $(process[]) / $(num_α * num_T * num_Ω * num_d) ($(round(process[] / (num_α * num_T * num_Ω * num_d) * 100, digits=1)) %)]")
                    Ωprocess += 1   # Increment Ωrange iterator.
                    process += 1    # Increment total iterator.
                    print(io, "\e[2F")
                end

                if verbose && (!iszero(T) || !iszero(Ω)) print(io, "\e[7F") end
            end
            if verbose && !iszero(T) print(io, "\e[20F") end   # Move up 20 lines and erase.
        end 
        if verbose print(io, "\e[15F") end   # Move up 20 lines and erase. 
        end # End dimensions loop.

        if verbose print(io, "\e[5F") end
    end
    if verbose print("\e[?25h") end # Show cursor again.

    # Extract polaron data as Vector ready for Polaron type.
    polaron_data = [
        p["α"], # Alphas.
        p["αeff"], # Alphas sums.
        p["mb"], # Band mass.
        p["T"], # Temperatures.
        p["ω"], # Phonon frequencies.
        p["ωeff"], # Hellwarth eff phonon frequency.
        p["γ"], # Adiabaticity
        p["J"], # Transfer integral.
        p["a"], # Lattice constant.
        p["β"], # Inverse temperatures.
        p["Ω"], # Photon frequencies.
        p["d"], # Number of spatial dimensions.
        p["v0"], # Athermal v params.
        p["w0"], # Athermal w params.
        p["F0"], # Athermal energies.
        p["A0"], # Athermal A parameter.
        p["B0"], # Athermal B parameter.
        p["C0"], # Athermal C parameter.
        p["κ0"], # Fictitious spring constant.
        p["M0"], # Athermal fictitious mass.
        p["M0a"], # Athermal asymptotic approximate fictitious mass (v0/w0).
        p["M0r"], # Athermal reduced mass of particle + fictitious particle system.
        p["R0"], # Athermal polaron radius (s. d. of a Gaussian wavefunction).
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
        p["χ"], # Memory function or self energy.
        p["z"], # Complex impedences.
        p["σ"], # Complex conductivities.
    ]

    # Return Polaron type containing all generated data over any coupling strengths, temperatures and frequencies.
    return HolsteinPolaron(polaron_data...)
end

"""
Single alpha parameter. holstein() expects alpha parameters to be in a Vector.
"""
holsteinpolaron(α::Real, Trange, Ωrange; ω=1, ωeff=1, γ=1, mb=1, J=1, a=1, v_guesses=3.11, w_guesses=2.87, dims=3, kspace=false, verbose=false) = holsteinpolaron([α], Trange, Ωrange; ω=ω, ωeff=ωeff, γ=γ, mb=mb, J=J, a=a, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, kspace=kspace, verbose=verbose)

"""
No frequency input => Phonon peak.
"""
holsteinpolaron(αrange, Trange; ω=1, ωeff=1, γ=1, mb=1, J=1, a=1, v_guesses=3.11, w_guesses=2.87, dims=3, kspace=false, verbose=false) = holsteinpolaron(αrange, Trange, 0; ω=ω, ωeff=ωeff, γ=γ, mb=mb, J=J, a=a, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, kspace=kspace, verbose=verbose)

"""
No temperature input => 300 K.
"""
holsteinpolaron(αrange; ω=1, ωeff=1, γ=1, mb=1, J=1, a=1, v_guesses=3.11, w_guesses=2.87, dims=3, kspace=false, verbose=false) = holsteinpolaron(αrange, 0, 0; ω=ω, ωeff=ωeff, γ=γ, mb=mb, J=J, a=a, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, kspace=kspace, verbose=verbose)

"""
No input => α = 1
"""
holsteinpolaron(; ω=1, ωeff=1, γ=1, mb=1, J=1, a=1, v_guesses=3.11, w_guesses=2.87, dims=3, kspace=false, verbose=false) = holsteinpolaron(1, 0, 0; ω=ω, ωeff=ωeff, γ=γ, mb=mb, J=J, a=a, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, kspace=kspace, verbose=verbose)

"""
    holstein(material::Material, TΩrange...; v_guesses=3.11, w_guesses=2.87, verbose=false)
Material specific constructors that use material specific parameters to parameterise the polaron.
Material data is inputted through the `Material` type.
Returns all data in either SI units or other common, suitable units otherwise.
"""
function holsteinpolaron(material::Material, TΩrange...; v_guesses=3.11, w_guesses=2.87, kspace=false, verbose=false)

    # Show material data.
    if verbose
        display(material)
    end
    
    # Extract material data from Material type.
    ω = material.f
    ωeff = material.feff
    mb = material.mb
    J = material.J
    γ = material.γ
    a = material.a
    dims = material.z ./ 2

    TΩrange = length(TΩrange) == 1 ? TΩrange .* phustrip(1u"K") : TΩrange .* (phustrip(1u"K"),1)

    # Generate Holstein polaron data from the arbitrary model constructor.
    p = holsteinpolaron(material.α', TΩrange...; ω=ω, ωeff=ωeff, γ=γ, mb=mb, J=J, a=a, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, kspace=kspace, verbose=verbose)

    # Return material-specific, unitful Holstein type.
    return p
end

# Broadcast Polaron data.
function Base.show(io::IO, ::MIME"text/plain", x::HolsteinPolaron)

    io_lim = IOContext(io, :limit => true, :compact => true)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                         Polaron Information:                          ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KPhonon frequencies             | ωeff = ", x.ωeff, " | ω = ", x.ω)
    println(io_lim, "\e[KHolstein coupling              | αeff = ", x.αeff, " | α = ", x.α)
    println(io_lim, "\e[KAdiabaticity                   | γ = ", x.γ)
    println(io_lim, "\e[KNumber of spatial dimensions   | d = ", x.d)
    println(io_lim, "\e[KTransfer integral              | J = ", x.J)
    println(io_lim, "\e[KLattice constant               | a = ", x.a)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                     Zero Temperature Information:                     ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KVariational parameter          | v0 = ", x.v0)
    println(io_lim, "\e[KVariational parameter          | w0 = ", x.w0)
    println(io_lim, "\e[KTotal energy                   | F0 = ", x.F0)
    println(io_lim, "\e[KElectron energy                | A0 = ", x.A0)
    println(io_lim, "\e[KInteraction energy             | B0 = ", x.B0)
    println(io_lim, "\e[KTrial energy                   | C0 = ", x.C0)
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
    println(io_lim, "\e[KElectron energy                | A = ", x.A)
    println(io_lim, "\e[KInteraction energy             | B = ", x.B)
    println(io_lim, "\e[KTrial energy                   | C = ", x.C)
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
function save_holstein_polaron(polaron::HolsteinPolaron, prefix)

    println("Saving polaron data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "alpha", polaron.α,
        "alpha eff", polaron.αeff,
        "band mass", polaron.mb,
        "temperature", polaron.T,
        "phonon freq", polaron.ω,
        "phonon freq eff", polaron.ωeff,
        "adiabaticity", polaron.γ,
        "transfer integral", polaron.J,
        "lattice constant", polaron.a,
        "beta", polaron.β,
        "E-field freq", polaron.Ω,
        "dims", polaron.d,
        "v0", polaron.v0,
        "w0", polaron.w0,
        "athermal energy", polaron.F0,
        "athermal A", polaron.A0,
        "athermal B", polaron.B0,
        "athermal C", polaron.C0,
        "athermal spring", polaron.κ0,
        "athermal mass", polaron.M0,
        "athermal asympt mass", polaron.M0a,
        "athermal reduced mass", polaron.M0r,
        "athermal size", polaron.R0,
        "v", polaron.v,
        "w", polaron.w,
        "thermal energy", polaron.F,
        "thermal A", polaron.A,
        "thermal B", polaron.B,
        "thermal C", polaron.C,
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

    holstein_polaron = HolsteinPolaron(
        data["alpha"],
        data["alpha eff"],
        data["band mass"],
        data["temperature"],
        data["phonon freq"],
        data["phonon freq eff"],
        data["adiabaticity"],
        data["transfer integral"],
        data["lattice constant"],
        data["beta"],
        data["E-field freq"],
        data["dims"],
        data["v0"],
        data["w0"],
        data["athermal energy"],
        data["athermal A"],
        data["athermal B"],
        data["athermal C"],
        data["athermal spring"],
        data["athermal mass"],
        data["athermal asympt mass"],
        data["athermal reduced mass"],
        data["athermal size"],
        data["v"],
        data["w"],
        data["thermal energy"],
        data["thermal A"],
        data["thermal B"],
        data["thermal C"],
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
