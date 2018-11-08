# FeynmanTheory.jl

"""
frohlichalpha(ε_Inf,ε_S,freq,m_eff)

    Calculates the Frohlich alpha parameter, for a given dielectric constant,
    frequency (f) of phonon in Hertz, and effective mass (in units of the
    bare electron mass).

    See Feynman 1955:
    http://dx.doi.org/10.1103/PhysRev.97.660

"""
function frohlichalpha(ε_optic,ε_static,freq,m_eff)
    ω=freq*2*pi #frequency to angular velocity
    # Note: we need to add a 4*pi factor to the permitivity of freespace.
    # This gives numeric agreement with literature values.  This is required as
    # the contemporary 1950s and 1960s literature implicitly used atomic units,
    # where the electric constant ^-1 has this factor baked in, k_e=1/(4πϵ_0).
    α= 0.5/(4*π*ε_0) *              # Units: m/F
       (1/ε_optic - 1/ε_static) *   # Units: none
       (q^2/(hbar*ω)) *             # Units: F
       sqrt(2*me*m_eff*ω/ħ)         # Units: 1/m
    return (α)
end

#####
# Athermal (Feynman 1955) model
# Set up equations for the polaron free energy, which we will variationally improve upon

# Integrand of (31) in Feynman I (Feynman 1955, Physical Review, "Slow electrons...")
fF(τ,v,w)=(w^2 * τ + (v^2-w^2)/v*(1-exp(-v*τ)))^-0.5 * exp(-τ)
# (31) in Feynman I
AF(v,w,α)=π^(-0.5) * α*v * quadgk(τ->fF(τ,v,w),0,Inf)[1]
# (33) in Feynman I
F(v,w,α)=(3/(4*v))*(v-w)^2-AF(v,w,α)

# Let's wrap the Feynman athermal variation approximation in a simple function
"""
    feynmanvw(α; v=7.0, w=6.0)

    Calculate v and w variational polaron parameters,
    for the supplied α Frohlich coupling.
    This version uses the original athermal action (Feynman 1955).
	Returns v,w.
"""
function feynmanvw(α; v=7.0, w=6.0) # v,w defaults
    initial=[v,w]
    # Main use of these bounds is stopping v or w going negative, at which you get a NaN error as you are evaluating log(-ve Real)
    lower=[0.1,0.1]
    upper=[1000.0,1000.0]

    myf(x) = F(x[1],x[2],α) # Wraps the function so just the two variational params are exposed

    # Now updated to use Optim.jl > 0.15.0 call signature (Julia >0.6 only)
    res=optimize(OnceDifferentiable(myf, initial; autodiff = :forward), 
                 lower, upper, initial, Fminbox( BFGS() ) )
    # specify Optim.jl optimizer. This is doing all the work.

    v,w=Optim.minimizer(res)

    return v,w
end



# Hellwarth et al. 1999 PRB - Part IV; T-dep of the Feynman variation parameter

# Originally a Friday afternoon of hacking to try and implement the T-dep electron-phonon coupling from the above PRB
# Which was unusually successful! And more or less reproduced Table III

# In Julia we have 'Multiple dispatch', so let's just construct the free
# energies (temperature-dependent) with the same name as before, but withthe
# thermodynamic beta where required

# Define Osaka's free-energies (Hellwarth1999 version) as Julia functions
# Equation numbers follow above Hellwarth et al. 1999 PRB
# 62b
A(v,w,β)=3/β*( log(v/w) - 1/2*log(2*π*β) - log(sinh(v*β/2)/sinh(w*β/2)))
# 62d
Y(x,v,β)=1/(1-exp(-v*β))*(1+exp(-v*β)-exp(-v*x)-exp(v*(x-β)))
# 62c integrand
f(x,v,w,β)=(exp(β-x)+exp(x))/(w^2*x*(1-x/β)+Y(x,v,β)*(v^2-w^2)/v)^(1/2)
# 62c
B(v,w,β,α) = α*v/(sqrt(π)*(exp(β)-1)) * quadgk(x->f(x,v,w,β),0,β/2)[1]
# 62e
C(v,w,β)=3/4*(v^2-w^2)/v * (coth(v*β/2)-2/(v*β))
# 62a
F(v,w,β,α)=-(A(v,w,β)+B(v,w,β,α)+C(v,w,β))

# Can now evaluate, e.g.
# F(v,w,β,α)=F(7.2,6.5,1.0,1.0)
# BUT - this is just the objective function! (The temperature-dependent free-energy.) Not the optimised parameters.
# Also there's a scary numeric integration (quadgk) buried within...

# Unbundled version of the finite temperature one
"""
    feynmanvw(α, βred; v=7.1, w=6.5, verbose::Bool=false)

    Calculate v and w variational polaron parameters,
    for the supplied α Frohlich coupling and βred reduced thermodynamic temperature.
    This uses the Osaka finite temperature action.
    Returns v,w.
"""
function feynmanvw(α, βred; v=7.1, w=6.5, verbose::Bool=false) # v,w defaults
    # Initial v,w to use
    initial=[v,w]
    # Main use of these bounds is stopping v or w going negative, at which you get a NaN error as you are evaluating log(-ve Real)
    lower=[0.1,0.1]
    upper=[100.0,100.0]

    myf(x) = F(x[1],x[2],βred,α) # Wraps the function so just the two variational params are exposed
    
    # Now updated to use Optim > 0.15.0 call signature (Julia >0.6 only)
    res=optimize(OnceDifferentiable(myf, initial; autodiff = :forward), 
                 lower, upper, initial, Fminbox( BFGS() ) )
    # specify Optim.jl optimizer. This is doing all the work.

    if Optim.converged(res) == false
        print("\tWARNING: Failed to converge to v,w soln? : ",Optim.converged(res) )
    end

    if verbose
        println()
        show(res)
    end

    v,w=Optim.minimizer(res)
    return v,w
end

