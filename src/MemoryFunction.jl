# MemoryFunction.jl
# References:
# [1] R. P. Feynman, R. W. Hellwarth, C. K. Iddings, and P. M. Platzman, Mobility of slow electrons in a polar crystal, PhysicalReview127, 1004 (1962) https://doi.org/10.1103/PhysRev.127.1004.
# [2] Devreese, J. De Sitter, and M. Goovaerts, “Optical absorption of polarons in thefeynman-hellwarth-iddings-platzman approximation,”Phys. Rev. B, vol. 5, pp. 2367–2381, Mar 1972.
# [3] F. Peeters and J. Devreese, “Theory of polaron mobility,” in Solid State Physics,pp. 81–133, Elsevier, 1984.

function frohlich_structure_factor(t, v, w, α, ω, β; dims = 3)
	if β == Inf
		return frohlich_structure_factor(t, v, w, α, ω; dims = dims)
	end
	coupling = frohlich_coupling(1, α, ω; dims = dims)
	propagator = polaron_propagator(im * t, v, w, β * ω) / 2
	integral = ball_surface(dims) / (2π)^dims * √π / 4 / propagator^(3/2)
	return coupling * integral * phonon_propagator(im * t / ω, ω, β) * 2 / dims / ω
end

function frohlich_structure_factor(t, v, w, α, ω; dims = 3)
	coupling = frohlich_coupling(1, α, ω; dims = dims)
	propagator = polaron_propagator(im * t, v, w) / 2
	integral = ball_surface(dims) / (2π)^dims * √π / 4 / propagator^(3/2)
	return coupling * integral * phonon_propagator(im * t / ω, ω) * 2 / dims / ω
end

function frohlich_memory_function(Ω, v, w, α, ω, β; dims = 3)
	if β == Inf
		return frohlich_memory_function(Ω, v, w, α, ω; dims = dims)
	end
    structure_factor(t) = frohlich_structure_factor(t, v, w, α, ω, β; dims = dims)
    return general_memory_function(Ω / ω, structure_factor)
end

function frohlich_memory_function(Ω, v, w, α, ω; dims = 3)
    structure_factor(t) = frohlich_structure_factor(t, v, w, α, ω; dims = dims)
    return general_memory_function(Ω / ω, structure_factor, limits = [0, 1e6])
end

frohlich_memory_function(Ω, v::Vector, w::Vector, α, ωβ...) = sum(frohlich_memory_function.(Ω, v, w, α, ωβ...))

frohlich_memory_function(Ω, v, w, α::Vector, ωβ...) = sum(frohlich_memory_function.(Ω, v, w, α, ωβ...))
