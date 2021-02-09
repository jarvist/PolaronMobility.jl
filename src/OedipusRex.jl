# OedipusRex.jl
# - Polaron optical absorption
#     A WORK IN PROGRESS!

"""
    ReX(v=20,w=20)

Codes to implement 'Optical Absorption of Polarons in the
Feynman-Hellwarth-Iddings-Platzman Approximation',
https://doi.org/10.1103/PhysRevB.5.2367 ;

In particular we are trying to calculate Re[χ] , the real part of the polaron optical absorption.
# (13) - defn of Re[χ]

"""
function ReX(v=20,w=20)
    R=(v^2-w^2)/(w^2*v)

    Reχintegrand(u,Ω) =  (1-cos(Ω*u)*exp(im*u))/(R*(1-exp(im*v*u))-im*u)^(3/2)

    Ω=1
    #Reχ = quadgk(u->Reχintegrand(u,Ω), 0.0, Inf)
# OK, problematic as this is a complex (multi valued integration!)

    [  Reχintegrand(u,Ω) for u=0:20 ]
    # Oh yikes, and it explodes as u->0 takes it to ->0 in the denominator
end


"""
ℜχ_0(Ω::Float64, α::Float64, v::Float64, w::Float64, N::Int)

Calculate the real part of χ(Ω) in a zero temperature approximation (equation
(16) in Devreese's et al.) for a given frequency Ω. v and w are the variational
Polaron parameters that minimise the free energy, for the supplied
α Frohlich coupling. N is the upper limit of a sum that is analytically
infinite, however most of the information is encapuslated by N < 50.
"""
function ℜχ_0(Ω, α, v, w, N = 10) # Large N for greater accuracy.
    R = (v^2 - w^2) / (w^2 * v)
    integrand(x, n) =
        ((n + 1 / 2) * x^(n - 1 / 2) * exp(-R * x) - R * x^(n + 1 / 2) * exp(-R * x)) *
        log(abs((1 + n * v + x)^2 / (Ω^2 - (1 + n * v + x)^2))^(1 / 2))
    total_sum = 0.0
    for n = 0:Int(N)
        total_sum += -(1 / gamma(n + 3 / 2)) * QuadGK.quadgk(x -> integrand(x, n), 0, Inf)[1]
    end
    return total_sum
end

"""
----------------------------------------------------------------------
Polaron absorption coefficient Γ(Ω).
----------------------------------------------------------------------
"""

"""
optical_absorption(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64; N::Int)

Calculate the absorption coefficient Γ(Ω) for the polaron at finite
temperatures (equation (11b) in Devreese's et al.) for a given frequency Ω.
β is the thermodynamic beta. v and w are the variational Polaron
parameters that minimise the free energy, for the supplied α Frohlich
coupling. N is the upper limit of a sum within ℑχ(Ω).
"""
function optical_absorption(Ω, β, α, v, w; N=10)
    ((Ω^2 - w^2)^2 / (Ω^5 - Ω^3 * v^2)^2) * ℑχ(Ω, β, α, v, w, N)
end

"""
optical_absorption_zero(Ω::Float64, α::Float64, v::Float64, w::Float64, n::Float64)

Calculate the absorption coefficient Γ(Ω) for the polaron at zero-temperature
(equation (11a) in Devreese's et al.) for a given frequency Ω. v and w are the
variational Polaron parameters that minimise the free energy, for the supplied
α Frohlich coupling.
n is the index of refraction of the medium.
"""
function optical_absorption_zero(Ω, α, v, w, n)
    Reχ = ℜχ_0(Ω, α, v, w)
    Imχ = ℑχ_0(Ω, α, v, w)
    Ω * Imχ / ((c * ϵ_0 * n) * (Ω^4 - 2 * Ω^2 * Reχ + Reχ^2 + Imχ^2) )
end
