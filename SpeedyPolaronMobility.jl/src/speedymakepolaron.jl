# speedymakepolaron.jl

"""
    speedymakepolaron(ϵ_optic, ϵ_static, phonon_freq, m_eff, T, Ωrange; v_guess = 5.4, w_guess = 3.4, volume = nothing, ir_activity = nothing, N_params = 1, verbose = false)

Solves the Feynman polaron problem variationally with finite temperature
Osaka energies. From the resulting v and w parameters, calculates polaron
structure (wave function size, etc.). 

Uses FHIP to predict a temperature- and frequency-dependent mobility, complex conductivity and impedance for the material system.

Can evaluate polaron properties for materials with multiple phonon branches using infrared activities and phonon branch frequencies.

Compared to `make_polaron`, this version multithreads.

# Arguments
- `ϵ_optic::Float64`: reduced optical dielectric constant.
- `ϵ_static::Float64`: reduced static dielectric constant.
- `phonon_freq::Vector{Float64}`: vector of characteristic dielectric phonon frequencies (THz).
- `m_eff::Float64`: bare-band effective-mass (mₑ).
- `T::Float64`: temperature value (K).
- `Ω::Float64`: electric field frequency vector or range (THz).
- `v_guess::Union{Float64, Vector{Float64}}`: initial guess for v params. Single value or vector of values for multiple parameters.
- `w_guess::Union{Float64, Vector{Float64}}`: initial guess for w params. Single value or vector of values for multiple parameters.
- `volume::Float64`: unit cell volume for the material (m³). `nothing` for one phonon mode.
- `ir_activity::Vector{Float64}`: vector of infrared activities. `nothing` for one phonon  mode.
- `N_params::Int`: number of variational parameters to minimise the polaron energy.
- 'verbose::Boolean': `true` to print progress. `false` to print nothing.

Returns a structure of type NewPolaron, containing arrays of useful
information.  Also prints a lot of information to the standard out - which
may be more useful if you're just inquiring as to a particular data point,
rather than plotting a temperature-dependent parameter.

As an example, to calculate the electron polaron in MAPI, at temperature 300.0 K and electric field frequency 2.25 THz, and inclusive of 15 optical phonon modes:

# Examples
```jldoctest
speedymakepolaron(
    4.5, 
    24.1, 
    [4.016471586720514, 3.887605410774121, 3.5313112232401513, 2.755392921480459, 2.4380741812443247, 2.2490917637719408, 2.079632190634424, 2.0336707697261187, 1.5673011873879714, 1.0188379384951798, 1.0022960504442775, 0.9970130778462072, 0.9201781906386209, 0.800604081794174, 0.5738689505255512], 
    0.12,
    300.0, 
    2.25, 
    volume = (6.29e-10)^3, 
    ir_activity = [0.08168931020200264, 0.006311654262282101, 0.05353548710183397, 0.021303020776321225, 0.23162784335484837, 0.2622203718355982, 0.23382298607799906, 0.0623239656843172, 0.0367465760261409, 0.0126328938653956, 0.006817361620021601, 0.0103757951973341, 0.01095811116040592, 0.0016830270365341532, 0.00646428491253749], 
    N_params = 1)
``` 
"""
function speedymakepolaron(ϵ_optic, ϵ_static, phonon_freq, m_eff, T::Float64, Ω::Float64; v_guess = 5.4, w_guess = 3.4, volume=nothing, ir_activity=nothing, N_params=1, verbose=false)

    # Number of phonon modes.
    N_modes = length(phonon_freq)

    # Convert THz phonon frequencies to 2π THz.
    ω = 2π .* phonon_freq

    if N_modes == 1
        # One phonon mode.

        # Calculate Frohlich alpha coupling parameters.
        α = [frohlichalpha(ϵ_optic, ϵ_static, phonon_freq, m_eff)]
    else
        # Multiple phonon modes

        # Calculate contribution to the ionic delectric constant for each phonon mode.
        ϵ_ionic = ϵ_ionic_mode.(phonon_freq, ir_activity, volume) 

        # Calculate contribution to Frohlich alpha for each phonon mode.
        α = multi_frohlichalpha.(ϵ_optic, ϵ_ionic, sum(ϵ_ionic), phonon_freq, m_eff)
    end

    # Calculate reduced thermodynamic betas for each phonon mode.
    # If temperature is zero, return Inf.
    betas = [T == 0.0 ? Inf64 : ħ * ω[j] / (kB * T) * 1e12 for j in eachindex(ω)]

    # Calculate variational parameters for each temperature from multiple phonon frequencies.
    params = T == 0.0 ? var_params(α; v=v_guess, w=v_guess, ω=ω, N=N_params) : var_params(α, betas; v=v_guess, w=w_guess, ω=ω, N=N_params)

    # Separate tuples of variational parameters into a list of 'v' and 'w' parameters for each temperature.
    v_params = params[1]
    w_params = params[2]

    # Calculate fictitious spring constants for each temperature.
    spring_constants = v_params .^ 2 .- w_params .^ 2

    # Calculate fictitious masses for each temperature.
    masses = spring_constants ./ w_params .^ 2

    # Calculate ground-state free energies for each temperature.
    energy = T == 0.0 ? multi_F(v_params, w_params, α; ω=ω) * 1000 * ħ / eV * 1e12 : multi_F(v_params, w_params, α, betas; ω=ω) * 1000 * ħ / eV * 1e12

    # Calculate complex impedances for each frequency and temperature. Returns a matrix.
    impedance = Ω == T == 0.0 ? 0.0 + 1im * 0.0 : polaron_complex_impedence(Ω, betas, α, v_params, w_params; ω=ω)

    # Calculate complex conductivities for each frequency and temperature. Returns a matrix.
    conductivity = Ω == T == 0.0 ? Inf64 + 1im * 0.0 : 1.0 / impedance

    # Calculates the dc mobility for each temperature.
    mobility = T == 0.0 ? Inf64 : polaron_mobility(betas, α, v_params, w_params; ω=ω) * eV / (1e12 * me * m_eff) * 100^2

    polaron = NewPolaron(α, T, betas, phonon_freq, hcat(v_params...), hcat(w_params...), hcat(spring_constants...), hcat(masses...), energy, Ω, impedance, conductivity, mobility)

    # Return Polaron mutable struct with evaluated data.
    return polaron
end

"""
    speedymakepolaron(ϵ_optic, ϵ_static, phonon_freq, m_eff, Trange, Ωrange; v_guess = 5.4, w_guess = 3.4, volume = nothing, ir_activity = nothing, N_params = 1, verbose = false)

Call the above `speedymakepolaron` function but for a range of temperatures `Trange::Union{Vector, StepRangeLen}`. Utilises the @tullio macro for multithreading and vectorisation.

# Arguments
- `ϵ_optic::Float64`: reduced optical dielectric constant.
- `ϵ_static::Float64`: reduced static dielectric constant.
- `phonon_freq::Vector{Float64}`: vector of characteristic dielectric phonon frequencies (THz).
- `m_eff::Float64`: bare-band effective-mass (mₑ).
- `Trange::Union{StepRangeLen, Vector}`: vector or range of temperature values (K).
- `Ωrange::Union{StepRangeLen, Vector}`: electric field frequency vector or range (THz).
- `v_guess::Union{Float64, Vector{Float64}}`: initial guess for v params. Single value or vector of values for multiple parameters.
- `w_guess::Union{Float64, Vector{Float64}}`: initial guess for w params. Single value or vector of values for multiple parameters.
- `volume::Float64`: unit cell volume for the material (m³). `nothing` for one phonon mode.
- `ir_activity::Vector{Float64}`: vector of infrared activities. `nothing` for one phonon  mode.
- `N_params::Int`: number of variational parameters to minimise the polaron energy.
- 'verbose::Boolean': `true` to print progress. `false` to print nothing.

Returns a structure of type NewPolaron, containing arrays of useful
information.  Also prints a lot of information to the standard out - which
may be more useful if you're just inquiring as to a particular data point,
rather than plotting a temperature-dependent parameter.

As an example, to calculate the electron polaron in MAPI, at temperatures 0:100:400 K and electric field frequencies 0.0:0.1:5.0 THz, and inclusive of 15 optical phonon modes:

# Examples
```jldoctest
speedymakepolaron(
    4.5, 
    24.1, 
    [4.016471586720514, 3.887605410774121, 3.5313112232401513, 2.755392921480459, 2.4380741812443247, 2.2490917637719408, 2.079632190634424, 2.0336707697261187, 1.5673011873879714, 1.0188379384951798, 1.0022960504442775, 0.9970130778462072, 0.9201781906386209, 0.800604081794174, 0.5738689505255512], 
    0.12,
    0.0:100.0:400.0, 
    0.0:0.1:5.0, 
    volume = (6.29e-10)^3, 
    ir_activity = [0.08168931020200264, 0.006311654262282101, 0.05353548710183397, 0.021303020776321225, 0.23162784335484837, 0.2622203718355982, 0.23382298607799906, 0.0623239656843172, 0.0367465760261409, 0.0126328938653956, 0.006817361620021601, 0.0103757951973341, 0.01095811116040592, 0.0016830270365341532, 0.00646428491253749], 
    N_params = 1)
``` 
"""
@noinline function speedymakepolaron(ϵ_optic, ϵ_static, phonon_freq, m_eff, Trange::Union{Vector, StepRangeLen}, Ωrange::Union{Vector, StepRangeLen}; volume=nothing, ir_activity=nothing, N_params=1, verbose=false)

    # Number of phonon modes.
    N_modes = length(phonon_freq)
    T_length = length(Trange)
    Ω_length = length(Ωrange)

    # Convert THz phonon frequencies to 2π THz.
    ω = 2π .* phonon_freq

    if N_modes == 1
        # One phonon mode.

        # Calculate Frohlich alpha coupling parameters.
        α = frohlichalpha(ϵ_optic, ϵ_static, phonon_freq, m_eff)
    else
        # Multiple phonon modes

        # Calculate contribution to the ionic delectric constant for each phonon mode.
        ϵ_ionic = ϵ_ionic_mode.(phonon_freq, ir_activity, volume) 

        # Calculate contribution to Frohlich alpha for each phonon mode.
        α = multi_frohlichalpha.(ϵ_optic, ϵ_ionic, sum(ϵ_ionic), phonon_freq, m_eff)
    end

    # Calculate reduced thermodynamic betas for each phonon mode.
    # If temperature is zero, return Inf.
    betas = Matrix{Float64}(undef, (T_length, N_modes))
    params = Vector{Any}(undef, T_length)
    v_params = Matrix{Float64}(undef, (T_length, N_params))
    w_params = Matrix{Float64}(undef, (T_length, N_params))
    energies = Vector{Float64}(undef, T_length)
    spring_constants = Matrix{Float64}(undef, (T_length, N_params))
    masses = Matrix{Float64}(undef, (T_length, N_params))
    impedances = Matrix{ComplexF64}(undef, (T_length, Ω_length))
    conductivities = Matrix{ComplexF64}(undef, (T_length, Ω_length))
    mobilities = Vector{Float64}(undef, T_length)

    if verbose
        process = Threads.Atomic{Int64}(1)
    end

    Threads.@threads for i in eachindex(Trange)

        if verbose
            println("[Progress : $(process[]) / $(T_length * Ω_length) ($(round(process[] / (T_length * Ω_length) * 100, digits=1)) %)] | Initialising... | Temperature T = $(Trange[i]) K | Thread #$(Threads.threadid())")
        end
        
        @fastmath @inbounds @simd for j in eachindex(phonon_freq)
            betas[i, j] = Trange[i] == 0.0 ? Inf64 : ħ * ω[j] / (kB * Trange[i]) * 1e12 
        end

        if !isassigned(params, i-1)
            v_guess = 5.6
            w_guess = 2.4
        else
            v_guess, w_guess, E_guess = params[i-1]
        end
        
        @fastmath @inbounds params[i] = Trange[i] == 0.0 ? var_params(α; v=v_guess, w=w_guess, ω=ω, N=N_params) : var_params(α, @view(betas[i, :]); v=v_guess, w=w_guess, ω=ω, N=N_params)

        @fastmath @inbounds @simd for j in 1:N_params
            v_params[i, j] = params[i][1][j]
            w_params[i, j] = params[i][2][j]
            spring_constants[i, j] = v_params[i, j] ^2 - w_params[i, j] ^2 
            masses[i, j] = spring_constants[i, j] / w_params[i, j]^2 
        end

        @fastmath @inbounds energies[i] = params[i][3] * 1000 * ħ / eV * 1e12

        @fastmath @inbounds @simd for j in eachindex(Ωrange)

            impedances[i, j] = Ωrange[j] == Trange[i] == 0.0 ? 0.0 + 1im * 0.0 : polaron_complex_impedence(Ωrange[j], @view(betas[i, :]), α, @view(v_params[i, :]), @view(w_params[i, :]); ω=ω) / eV^2 * 1e12 * me * m_eff * volume * 100^3
            
            conductivities[i, j] = Ωrange[j] == Trange[i] == 0.0 ? Inf64 + 1im * 0.0 : 1 / impedances[i, j]

            if verbose
                println("[Progress : $(process[]) / $(T_length * Ω_length) ($(round(process[] / (T_length * Ω_length) * 100, digits=1)) %)] | Temperature T = $(Trange[i]) K | Frequency Ω = $(Ωrange[j]) THz | Thread #$(Threads.threadid())")
                Threads.atomic_add!(process, 1)
            end
        end

        @fastmath @inbounds mobilities[i] = Trange[i] == 0.0 ? Inf64 : polaron_mobility(@view(betas[i, :]), α, @view(v_params[i, :]), @view(w_params[i, :]); ω=ω) * eV / (1e12 * me * m_eff) * 100^2
    end

    polaron = NewPolaron(α, Trange, betas, phonon_freq, v_params, w_params, spring_constants, masses, energies, Ωrange, impedances, conductivities, mobilities)

    # Return Polaron struct with evaluated data.
    return polaron
end

"""
    speedymakepolaron(α, T, Ωrange; ω = 1.0, v_guess = 5.4, w_guess = 3.4, verbose = false)

Same as first `speedymakepolaron` function but for a model system with specified alpha values rather than from material properties. Here we only have one phonon mode with a normalised frequency `ω = 1.0`.
"""
function speedymakepolaron(α::Float64, T::Float64, Ω::Float64; ω=1.0, v_guess = 5.4, w_guess = 3.4, N_params = 1, verbose=false)

    # Calculate reduced thermodynamic betas for each phonon mode.
    # If temperature is zero, return Inf.
    beta = T == 0.0 ? Inf64 : ω / T

    # Calculate variational parameters for each alpha parameter and temperature. Returns a Matrix of tuples.
    params = T == 0.0 ? var_params(α; v=v_guess, w=w_guess, ω=ω, N=N_params) : var_params(α, beta; v=v_guess, w=w_guess, ω=ω, N=N_params)

    # Separate tuples of variational parameters into a Matrices of 'v' and 'w' parameters for each alpha parameter and temperature.
    v_param = params[1]
    w_param = params[2]

    # Calculate fictitious spring constants for each alpha parameter and temperature. Returns a Matrix.
    spring_constant = v_param .^ 2 .- w_param .^ 2

    # Calculate fictitious masses for each alpha parameter and temperature. Returns a Matrix.
    mass = spring_constant ./ w_param .^ 2

    # Calculate ground-state free energies for each alpha parameter and temperature. Returns a Matrix.
    energy = params[3]

    # Calculate complex impedances for each alpha parameter, frequency and temperature. Returns a 3D Array.
    impedance =  Ω == T == 0.0 ? 0.0 + 1im * 0.0 : polaron_complex_impedence(Ω, beta, α, v_param, w_param; ω=ω) 

    # Calculate complex conductivities for each alpha parameter, frequency and temperature. Returns a 3D array.
    conductivity = Ω == T == 0.0 ? Inf64 + 1im * 0.0 : 1 / impedance

    # Calculates the dc mobility for each alpha parameter and each temperature.
    mobility = T == 0.0 ? Inf64 : polaron_mobility(beta, α, v_param, w_param; ω=ω)

    # Create Polaron type.
    polaron = NewPolaron(α, T, β, ω, v_param, w_param, spring_constant, mass, energy, Ω, impedance, conductivity, mobility)

    # Return Polaron struct with evaluated data.
    return polaron
end

"""
    speedymakepolaron(α, T, Ωrange; ω = 1.0, v_guess = 5.4, w_guess = 3.4, verbose = false)

Same as above `speedymakepolaron` function for a model system with specified alpha values but for a range or vector of temperatures.
"""
@noinline function speedymakepolaron(αrange::Union{Vector, StepRangeLen}, Trange::Union{Vector, StepRangeLen}, Ωrange::Union{Vector, StepRangeLen}; ω=1.0, N_params=1, verbose=false)

    # Number of phonon modes.
    T_length = length(Trange)
    Ω_length = length(Ωrange)
    α_length = length(αrange)

    # Calculate reduced thermodynamic betas for each phonon mode.
    # If temperature is zero, return Inf.
    betas = Vector{Float64}(undef, T_length)
    params = Matrix{Any}(undef, (α_length, T_length))
    v_params = Array{Float64, 3}(undef, (α_length, T_length, N_params))
    w_params = Array{Float64, 3}(undef, (α_length, T_length, N_params))
    energies = Matrix{Float64}(undef, (α_length, T_length))
    spring_constants = Array{Float64, 3}(undef, (α_length, T_length, N_params))
    masses = Array{Float64, 3}(undef, (α_length, T_length, N_params))
    impedances = Array{ComplexF64, 3}(undef, (α_length, T_length, Ω_length))
    conductivities = Array{ComplexF64, 3}(undef, (α_length, T_length, Ω_length))
    mobilities = Matrix{Float64}(undef, (α_length, T_length))

    if verbose
        process = Threads.Atomic{Int64}(1)
    end

    Threads.@threads for j in eachindex(Trange)

        @fastmath @inbounds betas[j] = Trange[j] == 0.0 ? Inf64 : ω / Trange[j]

        if !isassigned(params, j-1)
            v_guess = 5.6
            w_guess = 2.4
        else
            v_guess, w_guess, E_guess = params[j-1]
        end

        @fastmath @inbounds @simd for i in eachindex(αrange)

            # if verbose
            #     println("[Progress : $(process[]) / $(α_length * T_length * Ω_length) ($(round(process[] / (α_length * T_length * Ω_length) * 100, digits=1)) %)] | Initialising... | α = $(αrange[i]) | Temperature T = $(Trange[j]) K | Thread #$(Threads.threadid())")
            # end
        
            params[i, j] = Trange[j] == 0.0 ? var_params(αrange[i]; v=v_guess, w=w_guess, ω=ω, N=N_params) : var_params(αrange[i], betas[j]; v=v_guess, w=w_guess, ω=ω, N=N_params)

            for k in 1:N_params
                v_params[i, j, k] = params[i, j][1][k]
                w_params[i, j, k] = params[i, j][2][k]
                spring_constants[i, j, k] = v_params[i, j, k]^2 - w_params[i, j, k]^2 
                masses[i, j, k] = spring_constants[i, j, k] / w_params[i, j, k]^2 
            end

            energies[i, j] = params[i, j][3]

            @fastmath @inbounds @simd for k in eachindex(Ωrange)

                impedances[i, j, k] = Ωrange[k] == Trange[j] == 0.0 ? 0.0 + 1im * 0.0 : polaron_complex_impedence(Ωrange[k], betas[j], αrange[i], @view(v_params[i, j, :]), @view(w_params[i, j, :]); ω=ω)
                
                conductivities[i, j, k] = Ωrange[k] == Trange[j] == 0.0 ? Inf64 + 1im * 0.0 : 1 / impedances[i, j, k]

                if verbose
                    println("[Progress : $(process[]) / $(α_length * T_length * Ω_length) ($(round(process[] / (α_length * T_length * Ω_length) * 100, digits=1)) %)] | α = $(αrange[i]) | Temperature T = $(Trange[i]) K | Frequency Ω = $(Ωrange[k]) THz | Thread #$(Threads.threadid())")
                    Threads.atomic_add!(process, 1)
                end
            end

            mobilities[i, j] = Trange[j] == 0.0 ? Inf64 : polaron_mobility(betas[j], αrange[i], @view(v_params[i, j, :]), @view(w_params[i, j, :]); ω=ω)

        end
    end

    polaron = NewPolaron(αrange, Trange, betas, ω, v_params, w_params, spring_constants, masses, energies, Ωrange, impedances, conductivities, mobilities)

    # Return Polaron struct with evaluated data.
    return polaron
end