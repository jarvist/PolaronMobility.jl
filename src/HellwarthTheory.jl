# HellwarthTheory.jl

"""
    HellwarthBScheme(LO)

LO an array assumed to be of [freq ; absolute ir activity ]

Multiple phonon mode reduction to a single effective frequency.
Hellwarth et al. 1999 PRB, 'B scheme'; the athermal method.  Averaging procedure is
constructed by considering the average effect of the action of multiple branches.

Follows Eqn. (58) in this paper, assuming typo on LHS, should actually be W_e.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
function HellwarthBScheme(LO)
    println("Hellwarth B Scheme... (athermal)")
    H58 = sum((LO[:, 2] .^ 2) ./ LO[:, 1] .^ 2)
    println("Hellwarth (58) summation: ", H58)

    H59 = sum(LO[:, 2] .^ 2) # sum of total ir activity squarred
    println("Hellwarth (59) summation (total ir activity ^2): ", H59)
    println("Hellwarth (59) W_e (total ir activity ): ", sqrt(H59))

    omega = sqrt(H59 / H58)
    println("Hellwarth (61) Omega (freq): ", omega)

    return (omega)
end

# More complex scheme, involving thermodynamic Beta
# Hellwarth(50), RHS
"""
    HellwarthAScheme(phonon_modes; T = 295, convergence = 1e-6)

Multiple phonon mode reduction to a single effective frequency.
Temperature dependent, defaults to T = 295 K.

Solved iteratively by bisection until Δfreq < convergence.

Follows Hellwarth et al. 1999 PRB 'A' scheme, Eqn. (50) RHS.

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
function HellwarthAScheme(phonon_modes; T=295, convergence=1e-6)

    phonon_mode_freqs = phonon_modes[:, 1]
    ir_activities = phonon_modes[:, 2]

    condition(f) = coth(π * f * 1e12 * ħ / (kB * T)) / f
    -sum(ir_activities .*
         coth.(π .* phonon_mode_freqs .* 1e12 .* ħ ./ (kB * T)) ./ phonon_mode_freqs) / sum(ir_activities)

    # Solve by bisection
    minimum_frequency = minimum(phonon_mode_freqs)
    maximum_frequency = maximum(phonon_mode_freqs)
    middle_frequency = (maximum_frequency + minimum_frequency) / 2
    print("\n")

    while (maximum_frequency - minimum_frequency) / 2 > convergence

        if sign(condition(middle_frequency)) == sign(condition(minimum_frequency))
            minimum_frequency = middle_frequency
            middle_frequency = (maximum_frequency + minimum_frequency) / 2
        else
            maximum_frequency = middle_frequency
            middle_frequency = (maximum_frequency + minimum_frequency) / 2
        end
    end

    return middle_frequency
end
