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

function Hellwarth1999mobilityRHS( (α, (v,w) ,f) ,effectivemass, T)
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

"Multiple branches frohlich α"

function ϵ_ionic_mode(phonon_mode_freq, ir_activity, volume) # single ionic mode
    ω_j = 2π * phonon_mode_freq * 1e12 # angular phonon freq in Hz
    ϵ_mode = eV^2 * ir_activity / (3 * volume * ω_j^2 * amu) # single ionic mode
    return ϵ_mode / ϵ_0 # normalise with 1 / (4π ϵ_0)
end

function ϵ_total(freqs_and_ir_activity, volume) # total ionic contribution to dielectric
    phonon_freqs = freqs_and_ir_activity[:, 1] 
    ir_activity = freqs_and_ir_activity[:, 2]
    result = 0.0
    for (f, r) in zip(phonon_freqs, ir_activity)
        result += ϵ_ionic_mode(f, r, volume) # sum over all ionic contributions
    end
    return result
end

function effective_freqs(freqs_and_ir_activity, num_var_params) #PCA Algorithm
    standardized_matrix = freqs_and_ir_activity' .- mean(freqs_and_ir_activity', dims = 2) # centralise data by subtracting columnwise mean
    covariance_matrix = standardized_matrix' * standardized_matrix # has 1 / (n - 1) for number of params n = 2
    eigenvectors = eigvecs(covariance_matrix) # eigenvectors to project data along
    reduced_matrix = # project data along eigenvectors and undo centralisation
    standardized_matrix[:, 1:num_var_params] * eigenvectors[1:num_var_params, 1:num_var_params] * 
    eigenvectors[1:num_var_params, 1:num_var_params]' .+ mean(freqs_and_ir_activity', dims = 2)
    return abs.(reduced_matrix')
end

function frohlich_α_j(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff) # Frohlich alpha decomposed into phonon branch contributions
    Ry = eV^4 * me / (2 * ħ^2) # Rydberg energy
    ω = 2π * 1e12 * phonon_mode_freq # angular phonon freq (Hz)
    ϵ_static = ϵ_total + ϵ_optic # static dielectric. Calculate here instead of input so that ionic modes properly normalised.
    return (m_eff * Ry / (ħ * ω))^(1 / 2) * ϵ_ionic / (4π * ϵ_0) / (ϵ_optic * ϵ_static) # 1 / (4π ϵ_0) dielectric normalisation
end

"Multiple Phonon Branches"

function κ_i(i, v, w) # fictitious spring constant, multiple variational params
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

function h_i(i, v, w) # some vector relating to harmonic eigenmodes
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

function C_ij(i, j, v, w) # generalised Feynman C variational parameter (inclusive of multiple v and w params)
    C = w[i] * κ_i(i, v, w) * h_i(j, v, w) / (4 * (v[j]^2 - w[i]^2))
    return C
end

function D_j(τ, β, v, w) # log of dynamic structure factor for polaron 
    D = τ * (1 - τ / β)
    for i in 1:length(v)
        if v[i] != w[i]
        D += (h_i(i, v, w) / v[i]^2) * (2 * sinh(v[i] * τ / 2) * sinh(v[i] * (β - τ) / 2) / (v[i] * sinh(v[i] * β / 2)) - τ * (1 - τ / β))
        end
    end
    return D
end

function multi_free_energy(v_params, w_params, T, ϵ_optic, m_eff, volume, freqs_and_ir_activity)

    setprecision(BigFloat, 32) # Speed up. Stops potential overflows.

    # Extract phonon frequencies and ir activities.
    phonon_freqs = freqs_and_ir_activity[:, 1]
    ir_activity = freqs_and_ir_activity[:, 2]

    num_of_branches = length(phonon_freqs)
    
    # total dielectric contribution from all phonon branches (used as a normalisation)
    ϵ_tot = ϵ_total(freqs_and_ir_activity, volume)

    # Generalisation of B i.e. Equation 62c in Hellwarth.
    B_j_integrand(τ, β, v, w) = cosh(β / 2 - abs(τ)) / (sinh(β / 2) * sqrt(D_j(abs(τ), β, v, w)))
    B_j(β, α, v, w) = α / √π * quadgk(τ -> B_j_integrand(τ, β, v, w), 0.0, β / 2)[1]

    # Generalisation of C i.e. Equation 62e in Hellwarth.
    function C_j(β, v, w)
        s = 0.0
        for i in 1:length(v)
            for j in 1:length(v)
                s += C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2)  - 2 / (β * v[j]))
            end
        end
        3 * s / num_of_branches
    end

    # Generalisation of A i.e. Equation 62b in Hellwarth.
    function A_j(β, v, w)
        s = -log(2π * β) / 2
        for i in 1:length(v)
            if v[i] != w[i]
                s += log(v[i] / w[i]) - log(sinh(v[i] * β / 2) / sinh(w[i] * β / 2))
            end
        end
        3 / β * s / num_of_branches
    end
	
	F = 0.0
	for j in 1:num_of_branches
		ω_j = 2π * 1e12 * phonon_freqs[j] # angular phonon freq im 2π Hz
		β_j = BigFloat(ħ * ω_j / (k_B * T)) # reduced thermodynamic beta
		ϵ_ionic_j = ϵ_ionic_mode(phonon_freqs[j], ir_activity[j], volume) # ionic dielectric contribution for current phonon branch
		α_j = frohlich_α_j(ϵ_optic, ϵ_ionic_j, ϵ_tot, phonon_freqs[j], m_eff) # decomposed alpha for current phonon branch
		F += -(B_j(β_j, α_j, v_params, w_params) + C_j(β_j, v_params, w_params) + A_j(β_j, v_params, w_params)) * ω_j # × ħω branch phonon energy
        println("Free energy:\n
                Phonon freq: ", phonon_freqs[j], "\n 
                β: ", β_j, "\n
                ϵ_ionic: ", ϵ_ionic_j, "\n
                α_j: ", α_j, "\n")
    end
		
    println("Total free energy: ", F * ħ / eV * 1e3, " meV")
    return F * ħ / eV * 1e3 # change to meV
end

"""
function multi_variation(T, ϵ_optic, m_eff, volume, freqs_and_ir_activity; N = 1)

Multiple branch variational theory.
"""
function multi_variation(T, ϵ_optic, m_eff, volume, freqs_and_ir_activity; initial_vw = false, N = 1) # N number of v and w params

    setprecision(BigFloat, 32) # Speed up. Stops potential overflows.

    if initial_vw isa Bool
		# Intial guess for v and w.
    	initial = sort(rand(2 * N)) .* 4.0 .+ 1.0 # initial guess around 4 and ≥ 1.
		
		# Limits of the optimisation.
		lower = repeat([0.1], 2 * N)
		upper = repeat([60.0], 2 * N)
	else
		# Intial guess for v and w.
		initial = sort(vcat(initial_vw...))
		
		# Limits of the optimisation.
		lower = repeat([0.1], 2 * N)
		upper = repeat([60.0], 2 * N)
	end
	
	println("Initial guess: ", initial)

	# Osaka Free Energy function to minimise.
	f(x) = multi_free_energy([x[2 * n] for n in 1:Int(N)], [x[2 * n - 1] for n in 1:Int(N)], T, ϵ_optic, m_eff, volume, freqs_and_ir_activity)

	# Use Optim to optimise the free energy function w.r.t v and w.
	solution = Optim.optimize(
		Optim.OnceDifferentiable(f, initial; autodiff = :forward),
		lower,
		upper,
		initial,
		Fminbox(BFGS()),
		Optim.Options(time_limit = 20.0),
	)

	# Get v and w params that minimised free energy.
	var_params = Optim.minimizer(solution)

	# Update matrices for v and w parameters.
	v_params = [var_params[2 * n] for n in 1:N]
	w_params = [var_params[2 * n - 1] for n in 1:N]

	# Show current v and w that minimise jth phonon branch.
	println("Variational parameters: ", var_params)

    # Return variational parameters that minimise the free energy.
    return v_params, w_params
end

"""
function multi_memory_function(ν, β::Array, α::Array, v, w, ω::Array, m_eff)

    Calculate polaron complex memory function inclusive of multiple phonon branches j, each with angular frequency ω[j] (rad THz).
    α is an array of decomposed Frohlich alphas, one for each phonon frequency ω[j].
    β is an array of reduced thermodynamic betas, one for each phonon frequency ω[j].
    ν is the frequency (THz) of applied electric field.
    v and w are variational parameters that minimise the polaron free energy for the system.
    m_eff is the conduction band mass of the particle (typically electron / hole)

Multiple phonon complex memory function.
"""
function multi_memory_function(ν, β, α, v, w, ω, m_eff)

    # FHIP1962, page 1009, eqn (36).
    S(t, β) = cos(t - 1im * β / 2) / sinh(β / 2) / D_j(-1im * t, β, v, w)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β, ν) = (1 - exp(1im * 2π * ν * t)) * imag(S(t, β))

    memory = 0.0

    for j in 1:length(f) # sum over phonon branches
        println("Photon frequency = $ν, Phonon mode frequency = $(ω[j] / 2π)")
        memory += 2 * α[j] * ω[j]^2 * quadgk(t -> integrand(t, β[j], ν / ω[j]), 0.0, Inf)[1] / (3 * √π * 2π * ν)
    end

    println("Memory function: ", memory)
    return memory
end

"""
function multi_conductivity(ν, β::Array, α::Array, v, w, ω::Array, m_eff)

    Calculate polaron complex conductivity inclusive of multiple phonon branches j, each with angular frequency ω[j] (rad THz).
    α is an array of decomposed Frohlich alphas, one for each phonon frequency ω[j].
    β is an array of reduced thermodynamic betas, one for each phonon frequency ω[j].
    ν is the frequency (THz) of applied electric field.
    v and w are variational parameters that minimise the polaron free energy for the system.
    m_eff is the conduction band mass of the particle (typically electron / hole)

Multiple phonon complex memory function.
"""
function multi_conductivity(ν, β, α, v, w, ω, m_eff)

	# FHIP1962, page 1009, eqn (36).
    S(t, β) = cos(t - 1im * β / 2) / sinh(β / 2) / D_j(-1im * t, β, v, w)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β, ν) = (1 - exp(1im * 2π * ν * t)) * imag(S(t, β))

    impedence = 0.0

	for j in 1:length(ω)
		println("Photon frequency = $ν, Phonon mode frequency = $(ω[j] / 2π)")
		impedence += -1im * 2π * ν / length(ω) + 1im * 2 * α[j] * ω[j]^2 * quadgk(t -> integrand(t, β[j], ν / ω[j]), 0.0, Inf)[1] / (3 * √π * 2π * ν)
	end

	conductivity = 1 / impedence * eV * 100^2 / (m_eff * m_e * 1e12)
    println("Multiple complex conductivity: ", conductivity, " cm^2/Vs")
end
