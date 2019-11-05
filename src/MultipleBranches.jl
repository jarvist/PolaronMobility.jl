# MultipleBranches.jl
#  Extending the Feynman theory to multiple phonon branches

"""
frohlichPartial(ϵ_o,ϵ_s,ϵ_i,f,meff)

Calculate a (partial) dielectric electron-phonon coupling element. 
ϵ_o - optical dielectric
ϵ_s - total static dielectric contribution
ϵ_i - this mode contribution to dielectric
"""
function frohlichPartial(ϵ_o,ϵ_s,ϵ_i,f,meff)
    ω=f*2π
    α=1/(4*π*ɛ0) * ϵ_i/(ϵ_o*ϵ_s) * (q^2/ħ) * sqrt(meff*me/(2*ω*ħ))
    return α
end

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
        f=f #* THz
        ω=2π * f
        S=S * q^2 / amu
        ϵ_mode = S / ω^2 / volume
        ϵ_mode /= 3 # take isotropic average = divide by 3
        ϵ+=ϵ_mode
        println("Mode f= $f S= $S ϵ_mode = $(ϵ_mode/ɛ_0)")
    end
    println("Raw ionic dielectric contribution: $ϵ absolute $(ϵ/ɛ_0) relative")
    return ϵ
end


"""
IRtoalpha(IR,volume)

Calculates contribution to dielectric constant for each polar phonon mode, and
thereby the Frohlich alpha contribution for this mode.
"""
function IRtoalpha(IR,volume)
    ϵ=0.0 #* q^2 * THz^-2 * amu^-1 * m^-3
    α_sum=0.0
    for r in rows(IR)
        f,S=r # frequency in THz; activity in e^2 amu^-1
        f=f #* THz
        ω=2π * f
        S=S * q^2 / amu
        ϵ_mode = S / ω^2 / volume
        ϵ_mode /= 3 # take isotropic average = divide by 3
        ϵ+=ϵ_mode
        #println("Mode f= $f S= $S ϵ_mode = $(upreferred(ϵ_mode/u"ɛ0"))")
        α=frohlichPartial(ϵ_o,ϵ_s,ϵ_mode,f,meff)
        α_sum+=α
        if (α>0.1)
            println("Notable Mode f= $f    α_partial=$α")
        end
    end
    println("Raw ionic dielectric contribution: $ϵ absolute $(ϵ/ɛ_0) relative")
    return (ϵ/ɛ_0), α_sum
end

