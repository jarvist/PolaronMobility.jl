# FeynmanAthermalAsymptoticComparison.jl - reproduce some 1950s numeric results

push!(LOAD_PATH,"../src/") # load module from local directory
using PolaronMobility 

using Plots

# Here we collect all the data we need for plotting
αrange=0.2:0.1:10 # unstable below α=0.2
vs=[ feynmanvw(α)[1]  for α in αrange] # takes a while; computes all vw
ws=[ feynmanvw(α)[2]  for α in αrange] 
# Urgh - can't seem to do this more elegantly; recalculating all the params just to strip the second component of the tuple

Es=map(F,vs,ws,αrange); # calculate polaron energy for these params

p=plot(xlim=αrange,ylim=(0,), 
    xlabel="α parameter", ylabel="Reduced polaron units", 
    title="Athermal (numeric, variational) v,w, and E for α" )

plot!(αrange,vs,label="v")
plot!(αrange,ws,label="w")
plot!(αrange,Es,label="E")
hline!([0.0],label="",color=:black)

# I like this plot!
# You can see how v=w=3 as α-->0.0
# w then tends down to 1.0 in a nice sigmoid form as α --> ~10+
# v increases monotonically

savefig("AthermalvwE.pdf")

# Feynman Stat Mech (1972), p. 240
smallαw(α)=3.0
P(α)=2/smallαw(α) * (sqrt(1+smallαw(α))-1) # Nb: error in Feynman Stat Mech: sqrt(1+W) , as in Feynman1955
smallαv(α)=3*(1 + 2*α*(1-P(α))/(3*smallαw(α))) 

c=0.5772 # Euler Mascheroni constant

largeαw(α)=1.0
largeαv(α)=(4*α^2/9π) - (4*(log(2) + c/2)-1)

p=plot(xlim=αrange,ylim=(0,15), xlabel="α parameter", ylabel="v,w, Reduced polaron units", title="Asymtotic limits for v,w" )

plot!(αrange,α->smallαw(α),label="w, smallαw",style=:dash)
plot!(αrange,α->smallαv(α),label="v, smallαv",style=:solid)

plot!(αrange,α->largeαw(α),label="w, largeαw",style=:dash)
plot!(αrange,α->largeαv(α),label="v, largeαv",style=:solid)

plot!(αrange,vs,label="v, numeric",w=3,style=:solid)
plot!(αrange,ws,label="w, numeric",w=3,style=:dash)

savefig("AthermalvwAsymptotic.pdf")

smallαE(α)=-α-α^2/81
FeynmanStatMechlargeαE(α)=-α^2/(3π) - 3/2 *(2*log(2)+c)-3/4 # As in Feynman Stat Mech; with 2π corrected to 3π ...

Feynman1955largeαE(α)=-α^2/3π - 3*log(2) # as in Feynman I, 1955

p=plot(xlim=αrange,ylim=(-Inf,0), xlabel="α parameter", ylabel="Energy, Reduced polaron units", title="Asymtotic limits for E(α)" )

plot!(αrange,α->smallαE(α),label="smallαE")
plot!(αrange,α->FeynmanStatMechlargeαE(α),label="FeynmanStatMechlargeαE")
plot!(αrange,α->Feynman1955largeαE(α),label="Feynman1955largeαE")
plot!(αrange,Es,label="E (numeric integration)")

savefig("AthermalAsymptoticEnergy.pdf")

