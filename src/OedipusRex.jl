# OedipusRex.jl
# - Polaron optical absorption
#     A WORK IN PROGRESS!

#using QuadGK
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

