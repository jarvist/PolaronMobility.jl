# speedymakepolaron.jl

"""
    speedymakepolaron(ϵ_optic, ϵ_static, phonon_freq, m_eff, T, Ωrange; v_guess = 5.4, w_guess = 3.4, volume = nothing, ir_activity = nothing, N_params = 1, verbose = false)

Solves the Feynman polaron problem variationally with finite temperature
Osaka energies. From the resulting v and w parameters, calculates polaron
structure (wave function size, etc.). 

Uses FHIP to predict a temperature- and frequency-dependent mobility, complex conductivity and impedance for the material system.

Can evaluate polaron properties for materials with multiple phonon branches using infrared activities and phonon branch frequencies.

Compared to `make_polaron`, this version uses @tullio to multithread.

# Arguments
- `ϵ_optic::Float64`: reduced optical dielectric constant.
- `ϵ_static::Float64`: reduced static dielectric constant.
- `phonon_freq::Vector{Float64}`: vector of characteristic dielectric phonon frequencies (THz).
- `m_eff::Float64`: bare-band effective-mass (mₑ).
- `T::Float64`: temperature value (K).
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

As an example, to calculate the electron polaron in MAPI, at temperature 300.0 K and electric field frequencies 0.0:0.1:5.0 THz, and inclusive of 15 optical phonon modes:

# Examples
```jldoctest
speedymakepolaron(
    4.5, 
    24.1, 
    [4.016471586720514, 3.887605410774121, 3.5313112232401513, 2.755392921480459, 2.4380741812443247, 2.2490917637719408, 2.079632190634424, 2.0336707697261187, 1.5673011873879714, 1.0188379384951798, 1.0022960504442775, 0.9970130778462072, 0.9201781906386209, 0.800604081794174, 0.5738689505255512], 
    0.12,
    300.0, 
    0.0:0.1:5.0, 
    volume = (6.29e-10)^3, 
    ir_activity = [0.08168931020200264, 0.006311654262282101, 0.05353548710183397, 0.021303020776321225, 0.23162784335484837, 0.2622203718355982, 0.23382298607799906, 0.0623239656843172, 0.0367465760261409, 0.0126328938653956, 0.006817361620021601, 0.0103757951973341, 0.01095811116040592, 0.0016830270365341532, 0.00646428491253749], 
    N_params = 1)
``` 
"""
@noinline function speedymakepolaron(ϵ_optic, ϵ_static, phonon_freq, m_eff, T::Real, Ωrange; v_guess = 5.4, w_guess = 3.4, volume=nothing, ir_activity=nothing, N_params=1, verbose=false, processes=(1,1))

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
        ϵ_ionic = [ϵ_ionic_mode(i, j, volume) for (i, j) in zip(phonon_freq, ir_activity)]

        # Calculate contribution to Frohlich alpha for each phonon mode.
        α = [multi_frohlichalpha(ϵ_optic, i, sum(ϵ_ionic), j, m_eff) for (i, j) in zip(ϵ_ionic, phonon_freq)]
    end

    if verbose
        if Threads.threadid() == 1
            println("\e[0;0H\e[2J")
            print("Thread id: ", Threads.threadid(), " | Process: ", processes[1], " / ", processes[2], " ($(round.(processes[1] / processes[2] * 100, digits = 1)) %) | α = ", sum(α), " | T = ", T, "K")
        end
    end

    # Calculate reduced thermodynamic betas for each phonon mode.
    # If temperature is zero, return Inf.
    betas = [T == 0.0 ? Inf64 : ħ * ω[j] / (kB * T) * 1e12 for j in eachindex(ω)]

    # Calculate variational parameters for each temperature from multiple phonon frequencies.
    params = T == 0.0 ? var_params(α; v=5.6, w=3.4, ω=ω, N=N_params, verbose=verbose) : var_params(α, betas; v=v_guess, w=w_guess, ω=ω, N=N_params)

    # Separate tuples of variational parameters into a list of 'v' and 'w' parameters for each temperature.
    v_params = params[1]
    w_params = params[2]

    # Calculate fictitious spring constants for each temperature.
    spring_constants = v_params .^ 2 .- w_params .^ 2

    # Calculate fictitious masses for each temperature.
    masses = spring_constants ./ w_params .^ 2

    # Calculate ground-state free energies for each temperature.
    energies = T == 0.0 ? multi_F(v_params, w_params, α; ω=ω) * 1000 * ħ / eV * 1e12 : multi_F(v_params, w_params, α, betas; ω=ω) * 1000 * ħ / eV * 1e12

    # Calculate complex impedances for each frequency and temperature. Returns a matrix.
    @tullio impedances[i] := Ωrange[i] == T == 0.0 ? 0.0 + 1im * 0.0 : polaron_complex_impedence(Ωrange[i], betas, α, v_params, w_params; ω=ω) / eV^2 * 1e12 * me * m_eff * volume * 100^3
   
    # Calculate complex conductivities for each frequency and temperature. Returns a matrix.
    @tullio conductivities[i] := Ωrange[i] == T == 0.0 ? Inf64 + 1im * 0.0 : 1 / impedances[i]

    # Calculates the dc mobility for each temperature.
    mobilities = T == 0.0 ? Inf64 : polaron_mobility(betas, α, v_params, w_params; ω=ω) * eV / (1e12 * me * m_eff) * 100^2

    polaron = NewPolaron(α, [T], betas, phonon_freq, hcat(v_params...), hcat(w_params...), hcat(spring_constants...), hcat(masses...), energies, Ωrange, impedances, conductivities, mobilities)

    if verbose
        if Threads.threadid() == 1
            print("\e[0;0H\e[2J")
            println(polaron)
            print("\n")
        end
    end

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
@noinline function speedymakepolaron(ϵ_optic, ϵ_static, phonon_freq, m_eff, Trange::Union{Vector, StepRangeLen}, Ωrange; volume=nothing, ir_activity=nothing, N_params=1, verbose=false)

    # Initialise Vector for polarons at each temperature in Trange. 
    polarons = Vector{NewPolaron}(undef, length(Trange))

    # Fill polarons vector with Polaron types evaluated from `make_polaron` function. Uses @tullio macro for speedup.
    @tullio polarons[i] = speedymakepolaron(ϵ_optic, ϵ_static, phonon_freq, m_eff, Trange[i], Ωrange; volume=volume, ir_activity=ir_activity, N_params=N_params, verbose=verbose, processes=(1+sum([isassigned(polarons, j) for j in eachindex(polarons)]), length(Trange)))

    # Take Vector of Polaron types and combine them into a single Polaron type.
    polaron = combine_polarons(polarons)

    # Return Polaron struct with evaluated data.
    return polaron
end

"""
    speedymakepolaron(α, T, Ωrange; ω = 1.0, v_guess = 5.4, w_guess = 3.4, verbose = false)

Same as first `speedymakepolaron` function but for a model system with specified alpha values rather than from material properties. Here we only have one phonon mode with a normalised frequency `ω = 1.0`.
"""
@noinline function speedymakepolaron(α::Real, T::Real, Ωrange; ω=1.0, v_guess = 5.4, w_guess = 3.4, verbose=false, processes=(1, 1))

    if verbose
        if Threads.threadid() == 1
            println("\e[0;0H\e[2J")
            print("Thread id: ", Threads.threadid(), " | Process: ", processes[1], " / ", processes[2], " ($(round.(processes[1] / processes[2] * 100, digits = 1)) %) | α = ", α, " | T = ", T, "K")
        end
    end

    # Calculate reduced thermodynamic betas for each phonon mode.
    # If temperature is zero, return Inf.
    β = T == 0.0 ? Inf64 : ω / T

    # Calculate variational parameters for each alpha parameter and temperature. Returns a Matrix of tuples.
    params = T == 0.0 ? var_params(α; v=v_guess, w=w_guess, ω=ω) : var_params(α, β; v=v_guess, w=w_guess, ω=ω)

    # Separate tuples of variational parameters into a Matrices of 'v' and 'w' parameters for each alpha parameter and temperature.
    v_params = params[1]
    w_params = params[2]

    # Calculate fictitious spring constants for each alpha parameter and temperature. Returns a Matrix.
    spring_constants = v_params .^ 2 .- w_params .^ 2

    # Calculate fictitious masses for each alpha parameter and temperature. Returns a Matrix.
    masses = spring_constants ./ w_params .^ 2

    # Calculate ground-state free energies for each alpha parameter and temperature. Returns a Matrix.
    energies = T == 0.0 ? multi_F(v_params, w_params, α; ω=ω) : multi_F(v_params, w_params, α, β; ω=ω)

    # Calculate complex impedances for each alpha parameter, frequency and temperature. Returns a 3D Array.
    @tullio impedances[i] := Ωrange[i] == T == 0.0 ? 0.0 + 1im * 0.0 : polaron_complex_impedence(Ωrange[i], β, α, v_params, w_params; ω=ω) (i in eachindex(Ωrange))

    # Calculate complex conductivities for each alpha parameter, frequency and temperature. Returns a 3D array.
    @tullio conductivities[i] := Ωrange[i] == T == 0.0 ? Inf64 + 1im * 0.0 : 1 / impedances[i] (i in eachindex(Ωrange))

    # Calculates the dc mobility for each alpha parameter and each temperature.
    mobilities = T == 0.0 ? Inf64 : polaron_mobility(β, α, v_params, w_params; ω=ω)

    # Create Polaron type.
    polaron = NewPolaron(α, T, β, ω, hcat(v_params...), hcat(w_params...), hcat(spring_constants...), hcat(masses...), energies, Ωrange, impedances, conductivities, mobilities)

    if verbose
        if Threads.threadid() == 1
            print("\e[0;0H\e[2J")
            println(polaron)
            print("\n")
        end
    end

    # Return Polaron struct with evaluated data.
    return polaron
end

"""
    speedymakepolaron(α, T, Ωrange; ω = 1.0, v_guess = 5.4, w_guess = 3.4, verbose = false)

Same as above `speedymakepolaron` function for a model system with specified alpha values but for a range or vector of temperatures.
"""
@noinline function speedymakepolaron(αrange::Union{Vector, StepRangeLen}, Trange::Union{Vector, StepRangeLen}, Ωrange::Union{Vector, StepRangeLen}; ω=1.0, verbose=false)

    # Initialise Matrix for polarons at each α value and temperature T. 
    polarons = Matrix{NewPolaron}(undef, (length(Trange), length(αrange)))

    # Fill polarons Matrix with Polaron types evaluated from `make_polaron` function. Uses @tullio macro for speedup.
    @tullio polarons[j, i] = speedymakepolaron(αrange[i], Trange[j], Ωrange; ω=ω, verbose=verbose, processes=(1+sum([isassigned(polarons, j) for j in eachindex(polarons)]), length(Trange) * length(αrange)))  (j in eachindex(Trange), i in eachindex(αrange))

    # Take Matrix of Polaron types and combine them into a single Polaron type.
    polaron = combine_polarons(polarons)

    # Return Polaron struct with evaluated data.
    return polaron
end