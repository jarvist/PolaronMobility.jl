push!(LOAD_PATH,"../src/") # load module from local directory
using PolaronMobility

# Ivan Biaggio emailed 7th Nov 2018, suggesting there may be an error in these
# codes, due to a typographic error in the equations (particularly 62c) in
# Hellwarth1999
# 
# As part of investigating this, these tests reproduce Table III in
# Hellwarth1999, to play with different ways of specifying the electron-phonon
# action.


"Print out F(alpha,beta) for a specific v,w; as a test"
function test_fns()
    @printf("\t\t")
    for α in 1:5
        @printf("α=%d\t\t",α)
    end
    @printf("\n")

    for β in 1:0.25:3.0
        v=w=4
        print("β: $β\t|\t")
        for α in 1:5
            @printf("%f\t",F(v,w,β,α))
        end
        println()
    end
end

# Our original definition of B, directly following the equation as written in Hellwarth1999
# 62c
println("Original B (62c, electron-phonon free energy):")
PolaronMobility.B(v,w,β,α) = α*v/(sqrt(π)*(exp(β)-1)) * quadgk(x->PolaronMobility.f(x,v,w,β),0,β/2)[1]
test_fns()  # OK - very primitive!
# Following private communication with Ivan Biaggio, and crosschecking against Osaka1959

println("Corrected B (62c, electron-phonon free energy, but corrected as Osaka.")
PolaronMobility.B(v,w,β,α) = exp(β)* α*v/(sqrt(π)*(exp(β)-1)) * quadgk(x->PolaronMobility.f(x,v,w,β),0,β/2)[1]
test_fns()

# Reproduce HellwarthTableIII, v and w solutions with finite Beta 
function HellwarthTableIII()
    @printf("\t\t")
    for α in 1:5
        @printf("α=%d\t\t",α)
    end
    @printf("\n")
    
    for β in 1:0.25:3.0
        print("β: $β\t|\t")
        for α in 1:5
            v,w = feynmanvw(α, β)
            @printf("%.2f %.2f\t",v,w)
        end
        println()
    end

    print("Athermal\n β=Inf\t|\t")
    for α in 1:5
        v,w = feynmanvw(α)
        @printf("%.2f %.2f\t",v,w)
    end
    println()
end

# Our original definition of B, directly following the equation as written in Hellwarth1999
# 62c
println("Original B (62c, electron-phonon free energy):")
PolaronMobility.B(v,w,β,α) = α*v/(sqrt(π)*(exp(β)-1)) * quadgk(x->PolaronMobility.f(x,v,w,β),0,β/2)[1]
HellwarthTableIII()  # OK - very primitive!
# Following private communication with Ivan Biaggio, and crosschecking against Osaka1959

println("Corrected B (62c, electron-phonon free energy, but corrected as Osaka.")
PolaronMobility.B(v,w,β,α) = exp(β)* α*v/(sqrt(π)*(exp(β)-1)) * quadgk(x->PolaronMobility.f(x,v,w,β),0,β/2)[1]
HellwarthTableIII()

end

# Taking it back to where it all started; directly doing the full integral in Osaka1959
OsakaIntegrand(t,v,w,β) = exp(-t) / sqrt( w^2*t*(1-t/β) + (v^2-w^2)/v  *(1-exp(-v*t)+(1-coth(v/2*β)*(cosh(v*t)-1))) )  
PolaronMobility.B(v,w,β,α) = 1/sqrt(pi) * α*v* exp(β)/(exp(β)-1) * quadgk(x->OsakaIntegrand(x,v,w,β),0,β)
# CURRENTLY BROKEN (!)
# sqrt(x) in the integrand is taken negative; so Julia throws a hissy fit - not expecting to find Imaginary values.
HellwarthTableIII()

