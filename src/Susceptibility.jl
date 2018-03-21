# Susceptibility.jl
#  Frequency-dependent susceptibility (and thus mobility) of the Feynman polaron state
#    A WORK IN PROGRESS.

# Data structure to store the results
struct susceptibility
    nu
    ImX
    μ
end
Susceptibility()=susceptibility([],[],[])

"""
function ImX(nurange,v,w,βred,α,ω,mb)

Impedance in (47a) from Feynman1962, directly solving freq dep without taking
Hellwarth1999 limit of v->0 . 

Calculates a frequency dependent (over range of nu) susceptibility which can be linked back to mobility.

HERE BE DRAGONS! 
Not well tested or validated code; the central numeric integration is highly
oscillatory and becomes intractable for large nu.  
"""
function ImX(nurange,v,w,βred,α,ω,mb)
        @printf("\nAin't no Imaginary Impedance like a Feynman 1962 ImX...\n")
        println("ImX Call signature: v: $v w: $w βred: $βred α: $α ω: $ω mb: $mb")
    # Feynman, I love you - but using Nu, v; Omega, w in the same paper + formulas, for similar objects?!
    s=Susceptibility()
    for nu in nurange
        R=(v^2-w^2)/(w^2*v)     # inline, page 300 just after Eqn (2). Note this is wrong in some textbooks / 1990s PRB.
        b=R*βred/sinh(βred*v/2) # Feynman1962 version; page 1010, Eqn (47b)
        a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2))
        k(u,a,b,v,nu) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u)*cos(nu*u) # integrand with cos(vu) term, as (47a)

        @printf("Numerical integration of Feynman1962(42a): nu=%.2f ",nu)
        @time n=quadgk(u->k(u,a,b,v,nu),0,Inf,maxevals=10^9,reltol=0.1) # numerical quadrature integration of (2)
        K=n[1]
        err=n[2]
        @printf(" quadgk: K=%g err=%g\n",K,err)

        # Full 47a constructed here
        ImX= 2*α/(3*sqrt(π)) * βred^(3/2) * (sinh(βred*nu/2))/sinh(βred/2) * (v^3/w^3) * K

        μ=ImX^-1 * (q)/(ω*mb)

        @printf(" %.3f %g %g\n",nu,ImX,μ)

        append!(s.nu,nu)
        append!(s.ImX,ImX)
        append!(s.μ,μ)
    end

    @printf("\n\n")
    return(s)
end

