println("Alpha-parameter, Cross check 'frohlichalpha()' fn vs. literature values.")

@testset "FrohlichAlpha" begin

α=frohlichalpha(2.3, 5.6, (4.9E13/(2*pi)), 1.0)
println("NaCl Frohlich paper, α=",α," should be ~about 5 (Feynman1955)")
@test α ≈ 5.0 atol=0.3


α=frohlichalpha(7.1, 10.4, 5.08E12, 0.095)
println("CdTe  α=",α," Stone 0.39 / Devreese 0.29")
@test α ≈ 0.3 atol=0.1


α=frohlichalpha(10.89, 12.9,  8.46E12, 0.063)
println("GaAs  α=",α," Devreese 0.068 ")
@test α ≈ 0.068 atol=0.01

end

