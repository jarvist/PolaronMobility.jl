# MultipleBranches.jl
#  Extending the Feynman theory to multiple phonon branches

"""
frohlichPartial((f, ϵ_mode); ϵ_o,ϵ_s,meff)

Calculate a (partial) dielectric electron-phonon coupling element.
f - frequency of mode in THz
ϵ_mode - this mode's contribution to dielectric
ϵ_o - optical dielectric
ϵ_s - total static dielectric contribution
"""
function frohlichPartial((f, ϵ_mode); ϵ_o,ϵ_s,meff)
    ω=f*1E12*2π
    α=1/(4*π*ϵ_0) * ϵ_mode/(ϵ_o*ϵ_s) * (q^2/ħ) * sqrt(meff*me/(2*ω*ħ))
    return α
end

# deprecated signature wrapped via multiple dispatch
frohlichPartial(ϵ_o,ϵ_s,ϵ_mode,f,meff) = frohlichPartial((f, ϵ_mode), ϵ_o=ϵ_o,ϵ_s=ϵ_s,meff=meff)

rows(M::Matrix) = map(x->reshape(getindex(M, x, :), :, size(M)[2]), 1:size(M)[1])

"""
IRtoDielectric(IRmodes,volume)

From absolute value of IR activities of phonon modes, generate a per-mode
contribution to the low-frequency dielectric constant.

IRmodes are tuples f,S with Frequency in THz; InfraRed activity in e^2 amu^-1
"""
function IRtoDielectric(IRmodes,volume)
    ϵ=0.0 #* q^2 * THz^-2 * amu^-1 * m^-3
    for r in rows(IRmodes)
        f,S=r # frequency in THz; activity in e^2 amu^-1
        f=f * 1E12 #* THz
        ω=2π * f
        S=S * q^2 / amu
        ϵ_mode = S / ω^2 / volume
        ϵ_mode /= 3 # take isotropic average = divide by 3
        ϵ+=ϵ_mode
        println("Mode f= $f S= $S ϵ_mode = $(ϵ_mode/ϵ_0)")
    end
    println("Raw ionic dielectric contribution: $ϵ absolute $(ϵ/ϵ_0) relative")
    return ϵ/ϵ_0
end


"""
IRtoalpha(IR,volume)

Calculates contribution to dielectric constant for each polar phonon mode, and
thereby the Frohlich alpha contribution for this mode.
"""
function IRtoalpha(IR; volume, ϵ_o,ϵ_s,meff)
    ϵ=0.0 #* q^2 * THz^-2 * amu^-1 * m^-3
    α_sum=0.0
    for r in rows(IR)
        f,S=r # frequency in THz; activity in e^2 amu^-1
        f=f * 1E12 #* THz
        ω=2π * f
        S=S * q^2 / amu
        ϵ_mode = S / ω^2 / volume
        ϵ_mode /= 3 # take isotropic average = divide by 3
        ϵ_mode /= ϵ_0 # reduced dielectric constant
        ϵ+=ϵ_mode
        #println("Mode f= $f S= $S ϵ_mode = $(upreferred(ϵ_mode/u"ϵ0"))")
        α=frohlichPartial(ϵ_o,ϵ_s,ϵ_mode,f/1E12,meff)
        α_sum+=α
        if (α>0.1)
            println("Notable Mode f= $f    α_partial=$α")
        end
    end
    println("Sum alpha: $(α_sum)")
    return α_sum
end

"""
DieletricFromIRmode(IRmode)

Calculate dielectric from an individual mode.

IRmode is a tuple f,S with Frequency in THz; InfraRed activity in e^2 amu^-1
"""
function DielectricFromIRmode(IRmode; volume)
    f,S=IRmode # Assumes Frequency in THz; InfraRed activity in e^2 amu^-1
    f*=1E12 # convert to Hz from THz
    ω=2π * f
    S=S * q^2 / amu
    ϵ_mode = S / ω^2 / volume
    ϵ_mode /= 3 # take isotropic average = divide by 3
    ϵ_mode /= ϵ_0 # reduced dielectric
    return ϵ_mode
end

function Hellwarth1999mobilityRHS( (α, (v,w) ,f) , effectivemass, T)
    mb=effectivemass*MassElectron
    ω=f*1E12*2π
    βred=ħ*ω/(kB*T)

    R=(v^2-w^2)/(w^2*v) # inline, page 300 just after Eqn (2)
    b=R*βred/sinh(βred*v/2) # Feynman1962 version; page 1010, Eqn (47b)
    a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2))
    k(u,a,b,v) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u) # integrand in (2)
    K=quadgk(u->k(u,a,b,v),0,Inf)[1] # numerical quadrature integration of (2)

    #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
    RHS=α/(3*sqrt(π)) * βred^(5/2) / sinh(βred/2) * (v^3/w^3) * K
    μ=RHS^-1 * (q)/(ω*mb)

    return 1/μ
end

"""
----------------------------------------------------------------------
Multiple Branch Frohlich Alpha
----------------------------------------------------------------------

This section of the code is dedicated to determining the partial dielectric electron-phonon coupling parameter, decomposed into ionic dielectric contributions from each phonon mode of the matieral. 
"""

"""
ϵ_ionic_mode(phonon_mode_freq::Float64, ir_activity::Float64, volume::Float64)

    Calculate the ionic contribution to the dielectric function for a given phonon mode.
    phonon_mode_freq is the frequency of the mode in THz.

     - ir_activity is the infra-red activity of the mode in e^2 amu^-1.
     - volume is the volume of the unit cell of the material in m^3.
"""
function ϵ_ionic_mode(phonon_mode_freq, ir_activity, volume) # single ionic mode

    # Angular phonon frequency for the phonon mode (rad Hz)
    ω_j = 2π * phonon_mode_freq * 1e12 

    # Dielectric contribution from a single ionic phonon mode
    ϵ_mode = eV^2 * ir_activity / (3 * volume * ω_j^2 * amu)

    # Normalise ionic dielectric contribution with 1 / (4π ϵ_0) (NB: the 4π has been pre-cancelled)
    return ϵ_mode / ϵ_0 
end

"""
ϵ_total(freqs_and_ir_activity::Matrix{Float64}, volume::Float64)

    Calculate the total ionic contribution to the dielectric function from all phonon modes.

     - freqs_and_ir_activity is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e^2 amu^-1) in the second column.
     - volume is the volume of the unit cell of the material in m^3.
"""
function ϵ_total(freqs_and_ir_activity, volume) # total ionic contribution to dielectric

    # Extract phonon frequencies (THz)
    phonon_freqs = freqs_and_ir_activity[:, 1] 

    # Extra infra-red activities (e^2 amu^-1)
    ir_activity = freqs_and_ir_activity[:, 2]

    # Sum over all ionic contribution from each phonon mode
    total_ionic = 0.0
    for (f, r) in zip(phonon_freqs, ir_activity)
        total_ionic += ϵ_ionic_mode(f, r, volume) 
    end

    return total_ionic
end

"""
effective_freqs(freqs_and_ir_activity::Matrix{Float64}, num_var_params::Integer)

    Generates a matrix of effective phonon modes with frequencies and infra-red activities derived from a larger matrix using the Principal Component Analysis (PCA) method.

     - freqs_and_ir_activity: is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e^2 amu^-1) in the second column.
     - num_var_params: is the number of effective modes required (which needs to be less than the number of modes in freqs_and_ir_activity)

    *** POSSIBLY REDUNDANT ***
"""
function effective_freqs(freqs_and_ir_activity, num_var_params) # PCA Algorithm

    # Check that the number of effective modes is less than the number of actual phonon modes.
    if num_var_params >= size(freqs_and_ir_activity)[1]

        println("The number of effective phonon modes has to be less than the total number of phonon modes.")

    else

        # Centralise data by subtracting the columnwise mean
        standardized_matrix = freqs_and_ir_activity' .- mean(freqs_and_ir_activity', dims = 2) 

        # Calculate the covariance matrix S' * S. Matrix size is (n - 1) x (n - 1) for number of params (here n = 2)
        covariance_matrix = standardized_matrix' * standardized_matrix 

        # Extract eigenvectors of the covariance matrix
        eigenvectors = eigvecs(covariance_matrix) 

        # Project the original data along the covariance matrix eigenvectors and undo the centralisation
        reduced_matrix = standardized_matrix[:, 1:num_var_params] * eigenvectors[1:num_var_params, 1:num_var_params] * 
        eigenvectors[1:num_var_params, 1:num_var_params]' .+ mean(freqs_and_ir_activity', dims = 2)

        # Resultant matrix is positive definite and transposed.
        return abs.(reduced_matrix')
    end
end

"""
frohlich_α_j(ϵ_optic::Float64, ϵ_ionic::Float64, ϵ_total::Float64, phonon_mode_freq::Float64, m_eff::Float64)

    Calculates the partial dielectric electron-phonon coupling parameter for a given longitudinal optical phonon mode. This decomposes the original Frohlich alpha coupling parameter (defined for a single phonon branch) into contributions from multiple phonon branches.

     - ϵ_optic is the optical dielectric constant of the material.
     - ϵ_ionic is the ionic dielectric contribution from the phonon mode.
     - ϵ_total is the total ionic dielectric contribution from all phonon modes of the material.
     - phonon_mode_freq is the frequency of the phonon mode (THz).
     - m_eff is the band mass of the electron (in units of electron mass m_e)
"""
function frohlich_α_j(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff) 

    # The Rydberg energy unit
    Ry = eV^4 * me / (2 * ħ^2)

    # Angular phonon frequency for the phonon mode (rad Hz).
    ω = 2π * 1e12 * phonon_mode_freq 

    # The static dielectric constant. Calculated here instead of inputted so that ionic modes are properly normalised.
    ϵ_static = ϵ_total + ϵ_optic

    # The contribution to the electron-phonon parameter from the currrent phonon mode. 1 / (4π ϵ_0) is the dielectric normalisation.
    α_j  = (m_eff * Ry / (ħ * ω))^(1 / 2) * ϵ_ionic / (4π * ϵ_0) / (ϵ_optic * ϵ_static)

    return α_j
end

"""
----------------------------------------------------------------------
Multiple Branch Polaron Free Energy
----------------------------------------------------------------------

This section of the code is dedicated to calculating the polaron free energy, generalised from Osaka's expression to the case where multiple phonon modes are present in the material.
"""

"""
κ_i(i::Integer, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Calculates the spring-constant coupling the electron to the `ith' fictitious mass that approximates the exact electron-phonon interaction with a harmonic coupling to a massive fictitious particle. NB: Not to be confused with the number of physical phonon branches; many phonon branches could be approximated with one or two etc. fictitious masses for example. The number of fictitious mass does not necessarily need to match the number of phonon branches.

     - i enumerates the current fictitious mass.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function κ_i(i, v, w)
    κ = v[i]^2 - w[i]^2
    if length(v) > 1
        for j in 1:length(v)
            if j != i
                κ *= (v[j]^2 - w[i]^2) / (w[j]^2 - w[i]^2)
            end
        end
    end
    return κ
end

"""
h_i(i::Integer, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Calculates the normal-mode (the eigenmodes) frequency of the coupling between the electron and the `ith' fictitious mass that approximates the exact electron-phonon interaction with a harmonic coupling to a massive fictitious particle. NB: Not to be confused with the number of physical phonon branches; many phonon branches could be approximated with one or two etc. fictitious masses for example. The number of fictitious mass does not necessarily need to match the number of phonon branches.

     - i enumerates the current fictitious mass.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function h_i(i, v, w)
    h = v[i]^2 - w[i]^2
    if length(v) > 1
        for j in 1:length(v)
            if j != i
                h *= (w[j]^2 - v[i]^2) / (v[j]^2 - v[i]^2)
            end
        end
    end
    return h
end

"""
C_ij(i::Integer, j::Integer, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Calculates the element to the coupling matrix C_ij (a generalisation of Feynman's `C` coupling variational parameter) between the electron and the `ith' and `jth' fictitious masses that approximates the exact electron-phonon interaction with a harmonic coupling to a massive fictitious particle. NB: Not to be confused with the number of physical phonon branches; many phonon branches could be approximated with one or two etc. fictitious masses for example. The number of fictitious mass does not necessarily need to match the number of phonon branches.

     - i, j enumerate the current fictitious masses under focus (also the index of the element in the coupling matrix C)
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function C_ij(i, j, v, w)
    C = w[i] * κ_i(i, v, w) * h_i(j, v, w) / (4 * (v[j]^2 - w[i]^2))
    return C
end

"""
D_j(τ::Float64, β::Float64, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Calculates the recoil function that approximates the exact influence (recoil effects) of the phonon bath on the electron with the influence of the harmonicly coupled fictitious masses on the electron. It appears in the exponent of the intermediate scattering function.

     - τ is the imaginary time variable.
     - β is the reduced thermodynamic temperature ħ ω_j / (kB T) associated with the `jth` phonon mode.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function D_j(τ, β, v, w)
    D = τ * (1 - τ / β)
    for i in 1:length(v)
        if v[i] != w[i]
        D += (h_i(i, v, w) / v[i]^2) * (2 * sinh(v[i] * τ / 2) * sinh(v[i] * (β - τ) / 2) / (v[i] * sinh(v[i] * β / 2)) - τ * (1 - τ / β))
        end
    end
    return D
end

"""
B_j(α::Float64, β::Float64, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Generalisation of the B function from Equation 62c in Biaggio and Hellwarth []. This is the expected value of the exact action <S_j> taken w.r.t trial action, given for the `jth` phonon mode.

     - α is the partial dielectric electron-phonon coupling parameter for the `jth` phonon mode.
     - β is the reduced thermodynamic temperature ħ ω_j / (kB T) associated with the `jth` phonon mode.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function B_j(α, β, v, w)
    B_integrand(τ) = cosh(β / 2 - abs(τ)) / (sinh(β / 2) * sqrt(D_j(abs(τ), β, v, w)))
    B = α / √π * quadgk(τ -> B_integrand(τ), 0.0, β / 2)[1]
    return B
end

"""
C_j(β::Float64, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), n::Float64)

    Generalisation of the C function from Equation 62e in Biaggio and Hellwarth []. This is the expected value of the trial action <S_0> taken w.r.t trial action.

     - β is the reduced thermodynamic temperature ħ ω_j / (kB T) associated with the `jth` phonon mode.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - n is the total number of phonon modes.
"""
function C_j(β, v, w, n)

    s = 0.0

    # Sum over the contributions from each fictitious mass.
    for i in 1:length(v)
        for j in 1:length(v)
            s += C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2)  - 2 / (β * v[j]))
        end
    end

    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    return 3 * s / n
end

"""
A_j(β::Float64, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), n::Float64)

    Generalisation of the A function from Equation 62b in Biaggio and Hellwarth []. This is the Helmholtz free energy of the trial model.

     - β is the reduced thermodynamic temperature ħ ω_j / (kB T) associated with the `jth` phonon mode.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - n is the total number of phonon modes.
"""
function A_j(β, v, w, n)

    s = -log(2π * β) / 2

    # Sum over the contributions from each fictitious mass.
    for i in 1:length(v)
        if v[i] != w[i]
            s += log(v[i] / w[i]) - log(sinh(v[i] * β / 2) / sinh(w[i] * β / 2))
        end
    end

    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    3 / β * s / n
end

"""
multi_free_energy(v_params::Array{Float64}(undef, 1), w_params::Array{Float64}(undef, 1), T::Float64, ϵ_optic::Float64, m_eff::Float64, volume::Float64, freqs_and_ir_activity::Matrix{Float64})

    Calculates the Helmholtz free energy of the polaron for a material with multiple phonon branches. This generalises Osaka's free energy expression (below Equation (22) in []).

     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - T is the environment temperature in kelvin (K).
     - ϵ_optic is the optical dielectric constant of the material.
     - m_eff is the band mass of the electron (in units of electron mass m_e) .
     - volume is the volume of the unit cell of the material in m^3.
     - freqs_and_ir_activity is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e^2 amu^-1) in the second column.
"""
function multi_free_energy(v_params, w_params, T, ϵ_optic, m_eff, volume, freqs_and_ir_activity)

    # Speed up. Stops potential overflows.
    setprecision(BigFloat, 32) 

    # Extract phonon modes frequencies and ir activities.
    phonon_freqs = freqs_and_ir_activity[:, 1]
    ir_activity = freqs_and_ir_activity[:, 2]

    # Total number of phonon modes / branches.
    num_of_branches = length(phonon_freqs)
    
    # Total ionic dielectric contribution from all phonon modes (used as a normalisation).
    ϵ_tot = ϵ_total(freqs_and_ir_activity, volume)

    F = 0.0

    # Sum over the phonon modes.
	for j in 1:num_of_branches

        # Angular phonon frequency of the phonon mode (rad Hz).
		ω_j = 2π * 1e12 * phonon_freqs[j]

        # Reduced thermodynamic temperature of the phonon mode (unitless).
		β_j = ħ * ω_j / (kB * T) 

        # The ionic contribution to the dielectric function from the phonon mode.
		ϵ_ionic_j = ϵ_ionic_mode(phonon_freqs[j], ir_activity[j], volume)

        # The partial dielectric electron-phonon coupling parameter of the phonon mode.
		α_j = frohlich_α_j(ϵ_optic, ϵ_ionic_j, ϵ_tot, phonon_freqs[j], m_eff) 

        # Add contribution to the total free energy from the phonon mode.
		F += -(B_j(β_j, α_j, v_params, w_params) + C_j(β_j, v_params, w_params, num_of_branches) + A_j(β_j, v_params, w_params)) * ħ * ω_j
        
        # Prints out the frequency, reduced thermodynamic temperature, ionic dielectric and partial coupling for the phonon mode.
        # println("Free energy: Phonon freq = ", phonon_freqs[j], " | β = ", β_j, " | ϵ_ionic = ", ϵ_ionic_j, " | α_j = ", α_j)
    end
	
    # print out the total polaron free energy from all phonon modes.
    # println("Total free energy: ", F * ħ / eV * 1e3, " meV")

    # Free energy in units of meV
    return F / eV * 1e3 
end

"""
multi_variation(T::Float64, ϵ_optic::Float64, m_eff::Float64, volume::Float64, freqs_and_ir_activity::Matrix{Float64}; initial_vw::{Bool, Array{Float64}(undef, 1)}, N::Integer)

    Minimises the multiple phonon mode free energy function for a set of v_p and w_p variational parameters.
    The variational parameters follow the inequality: v_1 > w_1 > v_2 > w_2 > ... > v_N > w_N.

     - T is the environment temperature in kelvin (K).
     - ϵ_optic is the optical dielectric constant of the material.
     - m_eff is the band mass of the electron (in units of electron mass m_e) .
     - volume is the volume of the unit cell of the material in m^3.
     - freqs_and_ir_activity is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e^2 amu^-1) in the second column.
     - initial_vw determines if the function should start with a random initial set of variational parameters (Bool input) or a given set of variational parameter values (one dimensional array).
     - N specifies the number of variational parameter pairs, v_p and w_p, to use in minimising the free energy.
"""
function multi_variation(T, ϵ_optic, m_eff, volume, freqs_and_ir_activity; initial_vw = false, N = 1) # N number of v and w params

    # Speed up. Stops potential overflows.
    setprecision(BigFloat, 32) 

    # Use a random set of N initial v and w values.
    if initial_vw isa Bool

		# Intial guess for v and w parameters.
    	initial = sort(rand(2 * N)) .* 4.0 .+ 1.0 # initial guess around 4 and ≥ 1.
		
		# Limits of the optimisation.
		lower = repeat([0.1], 2 * N)
		upper = repeat([60.0], 2 * N)

    # Use the N given initial v and w values.
	else

		# Intial guess for v and w.
		initial = sort(vcat(initial_vw...))
		
		# Limits of the optimisation.
		lower = repeat([0.1], 2 * N)
		upper = repeat([60.0], 2 * N)
	end
	
    # Print out the initial v and w values.
	println("Initial guess: ", initial)

	# The multiple phonon mode free energy function to minimise.
	f(x) = multi_free_energy([x[2 * n] for n in 1:Int(N)], [x[2 * n - 1] for n in 1:Int(N)], T, ϵ_optic, m_eff, volume, freqs_and_ir_activity)

	# Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
	solution = Optim.optimize(
		Optim.OnceDifferentiable(f, initial; autodiff = :forward),
		lower,
		upper,
		initial,
		Fminbox(LBFGS()),
		# Optim.Options(time_limit = 20.0), # Set time limit for asymptotic convergence if needed.
	)

	# Extract the v and w parameters that minimised the free energy.
	var_params = Optim.minimizer(solution)

	# Separate the v and w parameters into one-dimensional arrays (vectors).
	v_params = [var_params[2 * n] for n in 1:N]
	w_params = [var_params[2 * n - 1] for n in 1:N]

	# Print the variational parameters that minimised the free energy.
	println("Variational parameters: ", var_params)

    # Return the variational parameters that minimised the free energy.
    return v_params, w_params
end

"""
----------------------------------------------------------------------
Multiple Branch Polaron Memory Function and Complex Conductivity
----------------------------------------------------------------------

This section of the code is dedicated to calculating the polaron memory function and complex conductivity, generalised from FHIP's expression to the case where multiple phonon modes are present in the material.
"""

"""
function multi_memory_function(Ω::Float64, β::Array{Float64}(undef, 1), α::Array{Float64}(undef, 1), v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), ω::Array{Float64}(undef, 1), m_eff::Float64)

    Calculate polaron complex memory function inclusive of multiple phonon branches j, each with angular frequency ω[j] (rad THz).

     - Ω is the frequency (THz) of applied electric field.
     - β is an array of reduced thermodynamic betas, one for each phonon frequency ω[j]. 
     - α is an array of decomposed Frohlich alphas, one for each phonon frequency ω[j]. 
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - m_eff is the is the conduction band mass of the particle (typically electron / hole, in units of electron mass m_e).
"""
function multi_memory_function(Ω, β, α, v, w, ω, m_eff)

    # FHIP1962, page 1009, eqn (36).
    S(t, β) = cos(t - 1im * β / 2) / sinh(β / 2) / D_j(-1im * t, β, v, w)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β, ν) = (1 - exp(1im * 2π * ν * t)) * imag(S(t, β))

    memory = 0.0

    # Sum over the phonon modes.
    for j in 1:length(ω) 
        
        # print out the current photon frequency and phonon mode frequency (THz).
        # println("Photon frequency = $ν, Phonon mode frequency = $(ω[j] / 2π)")

        # Add the contribution to the memory function from the `jth` phonon mode.
        memory += 2 * α[j] * ω[j]^2 * quadgk(t -> integrand(t, β[j], ν / ω[j]), 0.0, Inf)[1] / (3 * √π * 2π * ν)
    end

    # Print out the value of the memory function.
    # println("Memory function: ", memory)

    return memory
end

"""
function multi_conductivity(Ω::Float64, β::Array{Float64}(undef, 1), α::Array{Float64}(undef, 1), v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), ω::Array{Float64}(undef, 1), m_eff::Float64)

    Calculate polaron complex conductivity inclusive of multiple phonon branches j, each with angular frequency ω[j] (rad THz).

     - Ω is the frequency (THz) of applied electric field.
     - β is an array of reduced thermodynamic betas, one for each phonon frequency ω[j]. 
     - α is an array of decomposed Frohlich alphas, one for each phonon frequency ω[j]. 
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - m_eff is the is the conduction band mass of the particle (typically electron / hole, in units of electron mass m_e).
"""
function multi_conductivity(ν, β, α, v, w, ω, m_eff)

	# FHIP1962, page 1009, eqn (36).
    S(t, β) = cos(t - 1im * β / 2) / sinh(β / 2) / D_j(-1im * t, β, v, w)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β, ν) = (1 - exp(1im * 2π * ν * t)) * imag(S(t, β))

    impedence = 0.0

    # Sum over the phonon modes.
	for j in 1:length(ω)

        # print out the current photon frequency and phonon mode frequency (THz).
		# println("Photon frequency = $ν, Phonon mode frequency = $(ω[j] / 2π)")

        # Add the contribution to the complex impedence from the `jth` phonon mode.
		impedence += -1im * 2π * ν / length(ω) + 1im * 2 * α[j] * ω[j]^2 * quadgk(t -> integrand(t, β[j], ν / ω[j]), 0.0, Inf)[1] / (3 * √π * 2π * ν)
	end

    # Conductivity is the reciprocal of the impedence. 
	conductivity = 1 / impedence * eV * 100^2 / (m_eff * me * 1e12)

    # Print out the value of the complex conductivity.
    # println("Multiple complex conductivity: ", conductivity, " cm^2/Vs")

    return conductivity
end

