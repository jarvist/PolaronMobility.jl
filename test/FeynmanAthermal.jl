@testset "FeynmanAthermal tests" begin

println("FeynmanAthermal Tests - check athermal feynmanvw(α) -> v,w literature values.")

# Results from Feynman & Hibbs, Emended Edition, p 319.
# α v w E
Schultz=[
3.00 3.44 2.55 -3.1333
5.00 4.02 2.13 -5.4401
7.00 5.81 1.60 -8.1127 
9.00 9.85 1.28 -11.486
11.0 15.5 1.15 -15.710
]

rows(M::Matrix) = map(x->reshape(getindex(M, x, :), :, size(M)[2]), 1:size(M)[1])

for r in rows(Schultz)
    α,vSchultz,wSchultz,ESchultz=r # Unpack row of results
    v,w=feynmanvw(α) # performs the optimisation
    E=F(v,w,α) # energy at the optimised parameters
    
    println("α=$α v=$v w=$w E=$E  | Schultz: v=$vSchultz w=$wSchultz E=$ESchultz")

    @test v ≈ vSchultz atol=0.1 # Strangely these need more tolerance than the Energies
    @test w ≈ wSchultz atol=0.1
    @test E ≈ ESchultz atol=0.001
end

end

