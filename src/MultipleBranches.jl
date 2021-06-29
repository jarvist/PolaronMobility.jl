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

"""
Brad's multiple phonon branch α, free energy and minimisation functions.
wip.
"""

MAPI= [
# 96.20813558773261 0.4996300522819191
# 93.13630357703363 1.7139631746083817
# 92.87834578121567 0.60108592692181
# 92.4847918585963 0.0058228799414729
# 92.26701437594754 0.100590086574602
# 89.43972834606603 0.006278895133832249
# 46.89209141511332 0.2460894564364346
# 46.420949316788 0.14174282581124137
# 44.0380222871706 0.1987196948553428
# 42.89702947649343 0.011159939465770681
# 42.67180170168193 0.02557751102757614
# 41.46971205834201 0.012555230726601503
# 37.08982543385215 0.00107488277468418
# 36.53555265689563 0.02126940080871224
# 30.20608114002676 0.009019481779712388
# 27.374810898415028 0.03994453721421388
# 26.363055017011728 0.05011922682554448
# 9.522966890022039 0.00075631870522737
4.016471586720514 0.08168931020200264
3.887605410774121 0.006311654262282101
3.5313112232401513 0.05353548710183397
2.755392921480459 0.021303020776321225
2.4380741812443247 0.23162784335484837
2.2490917637719408 0.2622203718355982
2.079632190634424 0.23382298607799906
2.0336707697261187 0.0623239656843172
1.5673011873879714 0.0367465760261409
1.0188379384951798 0.0126328938653956
1.0022960504442775 0.006817361620021601
0.9970130778462072 0.0103757951973341
0.9201781906386209 0.01095811116040592
0.800604081794174 0.0016830270365341532
0.5738689505255512 0.00646428491253749
# 0.022939578929507105 8.355742795827834e-06   # Acoustic modes!
# 0.04882611767873102 8.309858592685e-06
# 0.07575149723846182 2.778248540373041e-05
]

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

function D_ν(τ, β, v, w) # log of dynamic structure factor for polaron 
    D = τ * (1 - τ / β)
    for i in 1:length(v)
        if v[i] != w[i]
        D += (h_i(i, v, w) / v[i]^2) * (2 * sinh(v[i] * τ / 2) * sinh(v[i] * (β - τ) / 2) / (v[i] * sinh(v[i] * β / 2)) - τ * (1 - τ / β))
        end
    end
    return D
end

function multi_free_energy(v_params, w_params, T, ϵ_optic, m_eff, volume, freqs_and_ir_activity, phonon_branch)

    setprecision(BigFloat, 32) # Speed up. Stops potential overflows.

    # Extract phonon frequencies and ir activities.
    phonon_freqs = freqs_and_ir_activity[:, 1]
    ir_activity = freqs_and_ir_activity[:, 2]

    num_of_branches = length(phonon_freqs)
    j = phonon_branch # jth phonon branch
    
    # total dielectric contribution from all phonon branches (used as a normalisation)
    ϵ_tot = ϵ_total(freqs_and_ir_activity, volume)

    # Generalisation of B i.e. Equation 62c in Hellwarth.
    S_integrand(τ, β, v, w) = cosh(β / 2 - τ) / (sinh(β / 2) * sqrt(abs(D_ν(τ, β, v, w))))
    S(β, α, v, w) = α / √π * quadgk(τ -> S_integrand(τ, β, v, w), 0.0, β / 2)[1]

    # Generalisation of C i.e. Equation 62e in Hellwarth.
    function S_0(β, v, w)
        s = 0.0
        for i in 1:length(v)
            for j in 1:length(v)
                s += C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2)  - 2 / (β * v[j]))
            end
        end
        3 * s / num_of_branches 
    end

    # Generalisation of A i.e. Equation 62b in Hellwarth.
    function log_Z_0(β, v, w)
        s = -log(2π * β) / 2
        for i in 1:length(v)
            if v[i] != w[i]
                s += log(v[i] / w[i]) -log(sinh(v[i] * β / 2) / sinh(w[i] * β / 2))
            end
        end
        3 / β * s / num_of_branches
    end

    ω = 2π * 1e12 * phonon_freqs[j] # angular phonon freq im 2π Hz
    β = BigFloat(ħ * ω / (kB * T)) # reduced thermodynamic beta
    ϵ_ionic = ϵ_ionic_mode(phonon_freqs[j], ir_activity[j], volume) # ionic dielectric contribution for current phonon branch
    α = frohlich_α_j(ϵ_optic, ϵ_ionic, ϵ_tot, phonon_freqs[j], m_eff) # decomposed alpha for current phonon branch

    # Print out data.
    # println("α = $α, β = $β, f = $(phonon_freqs[j]), ")
    # print("F0 = $(log_Z_0(β, ω, v_params .* ω, w_params .* ω) / β), ")
    # print("<S0> = $(S_0(β, ω, v_params .* ω, w_params .* ω) / β), ")
    # print("<S> = $(S(β, α, ω, v_params .* ω, w_params .* ω) / β), ")
    # println("F = $(-log_Z_0(β, ω, v_params .* ω, w_params .* ω) / β - S(β, α, ω, v_params .* ω, w_params .* ω) / β + S_0(β, ω, v_params .* ω, w_params .* ω) / β)")
    
    # F = -(A + B + C) in Hellwarth.
    F = -(S(β, α, v_params, w_params) + S_0(β, v_params, w_params) + log_Z_0(β, v_params, w_params)) * ħ * ω # × ħω branch phonon energy
    
    return F / eV * 1e3 # change to meV
end

"""
function multi_variation(T, ϵ_optic, m_eff, volume, freqs_and_ir_activity; N = 1)

Multiple branch variational theory.
"""
function multi_variation(T, ϵ_optic, m_eff, volume, freqs_and_ir_activity; N = 1) # N number of v and w params

    setprecision(BigFloat, 32) # Speed up. Stops potential overflows.

    # Number of phonon branches.
    M = length(freqs_and_ir_activity[:, 1])

    # Initialise MxN matrices for v and w parameters. M is number of phonon branches. N is number of variational parameters (v & w) per branch.
    v_params = Matrix{Float64}(undef, M, N)
    w_params = Matrix{Float64}(undef, M, N)

    # Intial guess for v and w.
    initial = sort(rand(2 * N)) .* 4.0 .+ 1.0 # initial guess around 4 and ≥ 1.

    # Limits of the optimisation.
    lower = repeat([0.1], 2 * N) 
    upper = repeat([200.0], 2 * N)

    for j in 1:M # sum over phonon branches

        # Osaka Free Energy function to minimise.
        f(x) = multi_free_energy([x[2 * n] for n in 1:Int(N)], [x[2 * n - 1] for n in 1:Int(N)], T, ϵ_optic, m_eff, volume, freqs_and_ir_activity, j)

        # Use Optim to optimise the free energy function w.r.t v and w.
        solution = Optim.optimize(
            Optim.OnceDifferentiable(f, initial; autodiff = :forward),
            lower,
            upper,
            initial,
            Fminbox(BFGS()),
            # Optim.Options(time_limit = 10.0),
        )

        # Get v and w params that minimised free energy.
        var_params = Optim.minimizer(solution)

        # If v ≤ w quit as negative mass.
        # if any(sort([var_params[2 * n] for n in 1:Int(N)]) .<= sort([var_params[2 * n - 1] for n in 1:Int(N)]))
        #     return "v_i <= w_i"
        # end

        # Intialise next guess of v and w to be their values for the current phonon branch. (quicker convergence)
        initial = sort(var_params)

        # Update matrices for v and w parameters.
        v_params[j, :] .= [var_params[2 * n] for n in 1:N]
        w_params[j, :] .= [var_params[2 * n - 1] for n in 1:N]

        # Show current v and w that minimise jth phonon branch.
        println(var_params)
    end

    # Return variational parameters that minimise the free energy.
    return v_params, w_params
end


"""
function multi_susceptibility(Ω, β::Array, α::Array, v, w, f::Array)

    Calculate polaron complex susceptibility inclusive of multiple phonon branches j, each with frequency f[j] (THz).
    α is an array of decomposed Frohlich alphas, one for each phonon frequency f[j].
    β is an array of reduced thermodynamic betas, one for each phonon frequency f[j].
    Ω is the frequency (THz) of applied electric field.
    v and w are variational parameters that minimise the polaron free energy for the system.

Multiple branch complex susceptibility.
"""
function multi_impedence(Ω, β, α, v, w, f, m_eff)

    m_e = 9.10938188e-31; # Electron Mass (kg)
    eV = 1.602176487e-19; # Electron Volt (kg m^2 s^{-2})

    # FHIP1962, page 1009, eqn (36).
    S(t, β_j, v_j, w_j) = cos(t - 1im * β_j / 2) / sinh(β_j / 2) / D_ν(1im * t, β_j, v_j, w_j)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a). Scale Frequency Ω by phonon branch frequency f_j.
    integrand(t, β_j, v_j, w_j, Ω, ω_j) = (1 - exp(1im * Ω * t / ω_j)) * imag(S(t, β_j, v_j, w_j))

    impedence = 1im * Ω * 2π
    for j in 1:length(f) # sum over phonon branches
        impedence += 1im * (2 * α[j] * f[j] * quadgk(t -> integrand(t, β[j], v[j], w[j], Ω * 2π, 1), 0.0, Inf)[1] / (3 * √π * Ω))
        # Print data for current branch
        # println("$j, $(f[j])")
        # println("$Ω, $j, $(β[j]), $(α[j]), $(f[j]), $(v[j]), $(w[j]), $(impedence)")
    end
    println("$Ω")
    return impedence / (eV^2 / m_e / m_eff * 1e12 * 2π / 100^2)
end

"""
function multi_conductivity(susceptibility)

    Transforms complex susceptibility into the complex conductivity.
"""
function multi_conductivity(Ω, z)

    -1 ./ z
end
