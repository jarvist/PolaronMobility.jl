# Forked from Julia Jupyter notebook 2019-12-PolaronMultipleBranches.ipynb

using PolaronMobility
using Test

using Plots
using Printf

const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const MassElectron = 9.10938188e-31;                          # kg

# ((freq THz)) ((IR Activity / e^2 amu^-1))
# These data from MAPbI3-Cubic_PeakTable.csv
# https://github.com/WMD-group/Phonons/tree/master/2015_MAPbI3/SimulatedSpectra
# Data published in Brivio2015 (PRB)
# https://doi.org/10.1103/PhysRevB.92.144308
MAPI= [
#96.20813558773261 0.4996300522819191
#93.13630357703363 1.7139631746083817
#92.87834578121567 0.60108592692181
#92.4847918585963 0.0058228799414729
#92.26701437594754 0.100590086574602
#89.43972834606603 0.006278895133832249
#46.89209141511332 0.2460894564364346
#46.420949316788 0.14174282581124137
#44.0380222871706 0.1987196948553428
#42.89702947649343 0.011159939465770681
#42.67180170168193 0.02557751102757614
#41.46971205834201 0.012555230726601503
#37.08982543385215 0.00107488277468418
#36.53555265689563 0.02126940080871224
#30.20608114002676 0.009019481779712388
#27.374810898415028 0.03994453721421388
#26.363055017011728 0.05011922682554448
#9.522966890022039 0.00075631870522737
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
#0.022939578929507105 8.355742795827834e-05   # Acoustic modes!
#0.04882611767873102 8.309858592685e-06
#0.07575149723846182 2.778248540373041e-05
]

vol=(6.29E-10)^3
ϵ_o=4.5
meff=0.12
mb=meff*MassElectron

ϵ_i=IRtoDielectric(MAPI,vol)
ϵ_s=sum(ϵ_i)+ϵ_o # total (static) dielectric = sum of ionic, and optical

IRtoalpha(MAPI, volume=vol, ϵ_o=ϵ_o, ϵ_s=ϵ_s, meff=meff)

ϵ_polar_modes=DielectricFromIRmode.(eachrow(MAPI), volume=vol)
f_dielectric=hcat( MAPI[:,1], ϵ_polar_modes)
alphas=frohlichPartial.(eachrow(f_dielectric), ϵ_o = ϵ_o, ϵ_s = ϵ_o+sum(ϵ_polar_modes), meff=meff)

# First idea: solve for each alpha_i as a separate variational problem, then
# aggregate solutions
function separate_variational_solutions(alphas, MAPI)
    mobilityproblem=hcat(alphas, feynmanvw.(alphas), MAPI[:,1])
    inverse_μ=Hellwarth1999mobilityRHS.(eachrow(mobilityproblem), meff, 300)
    μ=sum(inverse_μ)^-1
    @printf("\n\tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)

    # as published in FrostPolaronMobility2017 PRB
    polaronmobility(300, 4.5, 24.1, 2.25E12, 0.12)

    # re-use these data, directly from the polaron problem
    polaronmobility(300, ϵ_o, ϵ_s, 2.25E12, 0.12)

    for T in 10:10:500
        inverse_μ=Hellwarth1999mobilityRHS.(eachrow(mobilityproblem), meff, T)
        μ=sum(inverse_μ)^-1
        @printf("\nT %f \tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs",T,μ,μ*100^2)
    end

    normal=polaronmobility(10:10:500, ϵ_o, ϵ_s, 2.25E12, 0.12)

    hcat(normal.T, normal.Hμ)

    # OK, we have a problem I think in that we are optimising different things here
end

# Second idea: solve for variational problem with explicit set of interacting
# modes in the 'B' component of the free energy

# copy + paste from the code
using QuadGK
# Define Osaka's free-energies (Hellwarth1999 version) as Julia functions
# Equation numbers follow above Hellwarth et al. 1999 PRB
# 62b
A(v,w,β)=3/β*( log(v/w) - 1/2*log(2*π*β) - log(sinh(v*β/2)/sinh(w*β/2)))
# 62d
Y(x,v,β)=1/(1-exp(-v*β))*(1+exp(-v*β)-exp(-v*x)-exp(v*(x-β)))
# 62c integrand
#   Nb: Magic number 1e-10 adds stablity to optimisation; v,w never step -ve
f(x,v,w,β)=(exp(β-x)+exp(x))/sqrt(1e-10+ w^2*x*(1-x/β)+Y(x,v,β)*(v^2-w^2)/v)
# 62c
B(v,w,β,α) = α*v/(sqrt(π)*(exp(β)-1)) * quadgk(x->f(x,v,w,β),0,β/2)[1]
# 62e
C(v,w,β)=3/4*(v^2-w^2)/v * (coth(v*β/2)-2/(v*β))
# 62a
F(v,w,β,α)=-(A(v,w,β)+B(v,w,β,α)+C(v,w,β))


β(T,ω)=ħ*ω/(kB*T)

ω_hellwarth=2π*1E12*HellwarthBScheme(MAPI) # just low frequency branches, otherwise produces too high a value


F(v,w,ω_hellwarth,T,α_s,ω_s)=-(A(v,w,β(T,ω_hellwarth))+sum(B.(v,w,β.(T,ω_s),α_s))+C(v,w,β(T,ω_hellwarth)))

# OK, let's try this with a single value
@test F(7.6,6.5, ω_hellwarth, 300, 2.4, ω_hellwarth) ≈ F(7.6,6.5, β(300,ω_hellwarth), 2.4)
# And compare to the original route

# Split the strength, should still give same result
@test F(7.6,6.5, ω_hellwarth, 300, [1.2, 1.2], [ω_hellwarth, ω_hellwarth]) ≈ F(7.6,6.5, β(300,ω_hellwarth), 2.4)


# repack / rename variables from above
α_s=alphas
ω_s=2π*1E12*MAPI[:,1]

sum(α_s)

import PolaronMobility.feynmanvw
using Optim

function feynmanvw(T, ω_hellwarth, α_s, ω_s; v=7.1, w=6.5, verbose::Bool=true) # v,w defaults
    # Initial v,w to use
    initial=[v,w]
    # Main use of these bounds is stopping v or w going negative, at which you get a NaN error as you are evaluating log(-ve Real)
    lower=[0.1,0.1]
    upper=[100.0,100.0]
    
    
    
    myf(x) = F(x[1],x[2],ω_hellwarth, T, α_s, ω_s) 
    # Wraps the function so just the two variational params are exposed, so Optim can call it

    # Now updated to use Optim > 0.15.0 call signature (Julia >0.6 only)
    res=optimize(OnceDifferentiable(myf, initial; autodiff = :forward), 
                 lower, upper, initial, Fminbox( BFGS() )) 
#                ,Optim.Options(g_tol=1e-15, allow_f_increases=true))
    # specify Optim.jl optimizer. This is doing all the work.

    if Optim.converged(res) == false
        print("\tWARNING: Failed to converge to v,w soln? : ",Optim.converged(res) )
    end

    if verbose # pretty print Optim solution
        println()
        show(res)
    end

    v,w=Optim.minimizer(res)
    return v,w
end


function HellwarthMobility(v,w,βred, α, ω)        
# Hellwarth1999 - directly do contour integration in Feynman1962, for
        # finite temperature DC mobility
        # Hellwarth1999 Eqn (2) and (1) - These are going back to the general
        # (pre low-T limit) formulas in Feynman1962.  to evaluate these, you
        # need to do the explicit contour integration to get the polaron
        # self-energy
        R=(v^2-w^2)/(w^2*v) # inline, page 300 just after Eqn (2)

        #b=R*βred/sinh(b*βred*v/2) # This self-references b! What on Earth?
        # OK! I now understand that there is a typo in Hellwarth1999 and
        # Biaggio1997. They've introduced a spurious b on the R.H.S. compared to
        # the original, Feynman1962:
        b=R*βred/sinh(βred*v/2) # Feynman1962 version; page 1010, Eqn (47b)

        a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2))
        k(u,a,b,v) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u) # integrand in (2)
        K=quadgk(u->k(u,a,b,v),0,Inf)[1] # numerical quadrature integration of (2)

        #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS= α/(3*sqrt(π)) * βred^(5/2) / sinh(βred/2) * (v^3/w^3) * K
        μ=RHS^-1 * (q)/(ω*mb)
        @printf("\n\tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)
        μ
end

fout=open("MultipleBranches-MAPI.dat","w")
@printf(fout,"T ExplicitBranches-v w mob HellwarthBScheme-v w mob")

Trange=10:1:500
for T in Trange
    # With explicit modes
    v,w=feynmanvw(T, ω_hellwarth, α_s, ω_s)
    μ=HellwarthMobility(v,w, β(T, ω_hellwarth), 2.4, ω_hellwarth)
    @printf(fout,"%d    %f %f %f",T,v,w,10_000*μ)
    
    # With Hellwarth effective mode frequency
    v,w=feynmanvw(T, ω_hellwarth, 2.4, ω_hellwarth)
    μ=HellwarthMobility(v,w, β(T, ω_hellwarth), 2.4, ω_hellwarth)
    @printf(fout," %f %f %f\n",v,w,10_000*μ)


    println("T: $(T) v: $(v) w:$(w)") # Just for STDOUT so you can monitor what's going on...
end

close(fout)


