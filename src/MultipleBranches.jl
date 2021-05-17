# MultipleBranches.jl
#  Extending the Feynman theory to multiple phonon branches

"""
frohlichPartial((f, ϵ_mode); ϵ_o,ϵ_s,meff)

Calculate a (partial) dielectric electron-phonon coupling element.
ϵ_o - optical dielectric
ϵ_s - total static dielectric contribution
ϵ_mode - this mode's contribution to dielectric
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
IRtoDielectric(IR,V)

From absolute value of IR activities of phonon modes, generate a per-mode
contribution to the low-frequencydielectric constant.
"""
function IRtoDielectric(IR,volume)
    ϵ=0.0 #* q^2 * THz^-2 * amu^-1 * m^-3
    for r in rows(IR)
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
    return f, α_sum
end

"""
DieletricFromIRmode(IRmode)

Calculate dielectric from an individual mode.
"""
function DielectricFromIRmode(IRmode; volume)
    f,S=IRmode # Assumes Frequency in THz; activity in e^2 amu^-1
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

# Physical constants
# const ħ = 1.05457162825e-34; # Reduced Planck's constant (kg m^2 s^{-1})
# const eV = 1.602176487e-19; # Electron Volt (kg m^2 s^{-2})
# const m_e = 9.10938188e-31; # Electron Mass (kg)
# const k_B = 1.3806504e-23; # Boltzmann's constant (kg m^2 K^{-1} s^{-2})
# const ϵ_0 = 8.854e-12 # Dielectric constant (C^2 N^{-1} m^{-2})
# const c = 2.99792458e8 # Speed of light (m s^{-1})
# const amu = 1.66053906660e-27 # Atomic Mass Unit (kg)

function ϵ_ionic_mode(phonon_mode_freq, ir_activity, volume)
    ω_ν = 2π * phonon_mode_freq * 1e12
    ϵ_mode = eV^2 * ir_activity / (3 * volume * ω_ν^2 * amu)
    return ϵ_mode / ϵ_0
end

function frohlich_α_ν(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff)
    Ry = eV^4 * m_e / (2 * ħ^2)
    ω = 2π * phonon_mode_freq * 1e12
    ϵ_static = ϵ_optic + ϵ_total
    (m_eff * Ry / (ħ * ω))^(1 / 2) * ϵ_ionic / (4π * ϵ_0 * ϵ_optic * ϵ_static)
end

MAPI= [
96.20813558773261 0.4996300522819191
93.13630357703363 1.7139631746083817
92.87834578121567 0.60108592692181
92.4847918585963 0.0058228799414729
92.26701437594754 0.100590086574602
89.43972834606603 0.006278895133832249
46.89209141511332 0.2460894564364346
46.420949316788 0.14174282581124137
44.0380222871706 0.1987196948553428
42.89702947649343 0.011159939465770681
42.67180170168193 0.02557751102757614
41.46971205834201 0.012555230726601503
37.08982543385215 0.00107488277468418
36.53555265689563 0.02126940080871224
30.20608114002676 0.009019481779712388
27.374810898415028 0.03994453721421388
26.363055017011728 0.05011922682554448
9.522966890022039 0.00075631870522737
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
0.022939578929507105 8.355742795827834e-06   # Acoustic modes!
0.04882611767873102 8.309858592685e-06
0.07575149723846182 2.778248540373041e-05
]

phonon_freqs = MAPI[:, 1]
ir_activity = MAPI[:, 2]
ir_total = sum(ir_activity)

function ϵ_total(phonon_freqs, ir_activity, volume)
    result = 0.0
    for (f, r) in zip(phonon_freqs, ir_activity)
        result += ϵ_ionic_mode(f, r, volume)
    end
    return result
end
"""
Check alpha parameter
"""
α = 0
s = 0
ϵ_t = ϵ_total(phonon_freqs, ir_activity, (6.29E-10)^3)
for (f, r) in eachrow(MAPI)
    global α
    global s
    ϵ_ionic = ϵ_ionic_mode(f, r, (6.29E-10)^3)
    α_ν = frohlich_α_ν(4.5, ϵ_ionic, ϵ_t, f, 0.12)
    # println("f: $f, ir: $r, α_ν: $α_ν, α: $α, ϵ: $ϵ_ionic")
    α += α_ν
    s += ϵ_ionic / ϵ_t
end
println("α = $(α), $s")

"""
Multiple Phonon Branch Free Energy
"""

using QuadGK

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

function C_ij(i, j, v, w)
    C = w[i] * κ_i(i, v, w) * h_i(j, v, w) / (4 * (v[j]^2 - w[i]^2))
    return C
end

function D_ν(τ, β, v, w)
    D = τ * (1 - τ / β)
    for i in 1:length(v)
        D += (h_i(i, v, w) / v[i]^2) * (2 * sinh(v[i] * τ / 2) * sinh(v[i] * (β - τ) / 2) / (v[i] * sinh(v[i] * β / 2)) - τ * (1 - τ / β))
    end
    return D
end

function multi_free_energy(v, w, T, ϵ_optic, ϵ_static, phonon_freqs, m_eff, ir_activity, volume)

    ϵ_tot = ϵ_total(phonon_freqs, ir_activity, volume)

    S_integrand(τ, β) = cosh(β / 2 - τ) / (sinh(β / 2) * sqrt(D_ν(τ, β, v, w)))
    S(β, α) = β * α / √π * QuadGK.quadgk(τ -> S_integrand(τ, β), 0, β / 2)[1]

    C(i, j) = C_ij(i, j, v, w)

    function S_0(β)
        s = 0.0
        for i in 1:length(v)
            for j in 1:length(v)
                s += C(i, j) / (v[j] * w[i]) * (coth(β * v[i] / 2) - 2 / (β * v[j]))
            end
        end
        -3 * s * β
    end

    function log_Z_0(β)
        s = 0.0
        for i in 1:length(v)
            s += log(v[i] * sinh(β * w[i] / 2) / (w[i] * sinh(β * v[i] / 2))) - log(2π * β) / 2
        end
        3 * s
    end

    F = 0.0
    for ν in 1:length(phonon_freqs)
        ω = 2π * 1e12 * phonon_freqs[ν]
        βred = ħ * ω / (k_B * T)
        β = ħ * 2π * 1e12 / (k_B * T) # Standard β to compare against.
        ϵ_ionic = ϵ_ionic_mode(phonon_freqs[ν], ir_activity[ν], volume)
        α = frohlich_α_ν(ϵ_optic, ϵ_ionic, ϵ_tot, phonon_freqs[ν], m_eff)
        # println("$α, $β, $(phonon_freqs[ν]), $(-log_Z_0(β) / β), $(-S(βred, α) / β), $(S_0(β) / β), $(-(log_Z_0(β) + S(βred, α) - S_0(β)) / β)")
        F += -(log_Z_0(β) + S(βred, α) - S_0(β)) / β * 2π * 1e12
    end
    return F * 1e3 * ħ / eV
end

# F = multi_free_energy([19.6], [16.3], 300, 4.5, 24.1, phonon_freqs, 0.12, ir_activity, (6.29E-10)^3)
# @show(F)
# 2.2490917637719408 0.2622203718355982

using Optim
using LineSearches

"""
Multiple Phonon Branch & Multiple variation parameter minimisation of free energy.
"""

function multi_variation(T, ϵ_optic, ϵ_static, phonon_freqs, m_eff, ir_activity, volume; N = 1)

    # Intial guess for v and w.
    initial = sort(rand(2 * N)) .* 10.0

    # Limits of the optimisation.
    lower = repeat([0.0], 2 * N)
    upper = repeat([Inf], 2 * N)

    # @show(initial, [initial[2 * n - 1] for n in 1:Int(N)], [initial[2 * n] for n in 1:Int(N)])

    # Osaka Free Energy function to minimise.
    f(x) = multi_free_energy([x[2 * n] for n in 1:Int(N)], [x[2 * n - 1] for n in 1:Int(N)], T, ϵ_optic, ϵ_static, phonon_freqs, m_eff, ir_activity, volume)

    err = false
    # Use Optim to optimise the free energy function w.r.t v and w.
    solution = try
        Optim.optimize(
            Optim.OnceDifferentiable(f, initial; autodiff = :forward),
            lower,
            upper,
            initial,
            Fminbox(BFGS(; linesearch = LineSearches.BackTracking())),
        )
    catch
        err = true
    end

    # # If optimisation does not converge or if v ≤ w, pick new random starting guesses for v and w between 1.0 and 11.0. Repeat until the optimisation converges with v > w.
    if err
        while err
            initial = sort(rand(2 * N)) .* 10.0
            err = false
            solution = try
                Optim.optimize(
                    Optim.OnceDifferentiable(f, initial; autodiff = :forward),
                    lower,
                    upper,
                    initial,
                    Fminbox(BFGS(; linesearch = LineSearches.BackTracking())),
                )
            catch
                err = true
            end
        end
    end

    # Return variational parameters that minimise the free energy.
    return Optim.minimizer(solution)
end

# N = 1
# r = sort(multi_variation(300, 4.5, 24.1, phonon_freqs, 0.12, ir_activity, (6.29E-10)^3; N = N))
# F = multi_free_energy([r[2 * n] for n in 1:N], [r[2 * n - 1] for n in 1:N], 300, 4.5, 24.1, phonon_freqs, 0.12, ir_activity, (6.29E-10)^3)
# println("$r, $F")
