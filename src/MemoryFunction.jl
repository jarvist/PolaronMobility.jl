# MemoryFunction.jl
# References:
# [1] R. P. Feynman, R. W. Hellwarth, C. K. Iddings, and P. M. Platzman, Mobility of slow electrons in a polar crystal, PhysicalReview127, 1004 (1962) https://doi.org/10.1103/PhysRev.127.1004.
# [2] Devreese, J. De Sitter, and M. Goovaerts, “Optical absorption of polarons in thefeynman-hellwarth-iddings-platzman approximation,”Phys. Rev. B, vol. 5, pp. 2367–2381, Mar 1972.
# [3] F. Peeters and J. Devreese, “Theory of polaron mobility,” in Solid State Physics,pp. 81–133, Elsevier, 1984.

function frohlich_structure_factor(t, v, w, α, ω, β)
	
	coupling = 2 * α / 3 / √π * ω^2

	propagator = polaron_propagator(im * t, v, w, β * ω)
	
	integral = phonon_propagator(im * t / ω, ω, β) / propagator^(3/2)

	coupling * integral 
end

function frohlich_structure_factor(t, v, w, α, ω)
	
	coupling = 2 * α / 3 / √π * ω^2

	propagator = polaron_propagator(im * t, v, w)
	
	integral = phonon_propagator(im * t / ω, ω) / propagator^(3/2)

	coupling * integral 
end

function frohlich_memory_function(Ω, v, w, α, ω, β)
    structure_factor(t) = frohlich_structure_factor(t, v, w, α, ω, β)
    return general_memory_function(Ω / ω, structure_factor; limits = [0, 1e4])
end

function frohlich_memory_function(Ω, v, w, α, ω)
    structure_factor(t) = frohlich_structure_factor(t, v, w, α, ω)
    return general_memory_function(Ω / ω, structure_factor; limits = [0, 1e4])
end

frohlich_memory_function(Ω, v::Vector, w::Vector, α, ωβ...) = sum(frohlich_memory_function.(Ω, v, w, α, ωβ...))

frohlich_memory_function(Ω, v, w, α::Vector, ωβ...) = sum(frohlich_memory_function.(Ω, v, w, α, ωβ...))
