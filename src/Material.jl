# Material.jl

"""
    mutable struct Material

    Structure to contain material-specific data. 

"""
mutable struct Material
    ϵo     # Optical dielectric constant
    ϵs     # Static dielectric constant
    ϵi       # Ionic dielectric contributions
    mb          # Effective band mass
    α           # Fröhlich coupling
    f       # Phonon frequencies
    feff    # Effective frequency
    ir          # Infrared activities
    V      # Unit cell volumes
    function Material(x...)
        new(reduce_array.(x)...)
    end
end

"""
    material(ϵ_optical, ϵ_static, m_eff, phonon_freq)

    Construct a 'struct Material' object from traditional Frohlich Hamiltonian parameters. 

    ϵ_ionic and α derived from these parmaters; other material properties filled with bogus values. 
"""
function material(ϵ_optical, ϵ_static, m_eff, phonon_freq)
    ϵ_ionic = ϵ_static - ϵ_optical
    α = frohlichalpha(ϵ_optical, ϵ_static, phonon_freq, m_eff)
    return Material(ϵ_optical, ϵ_static, ϵ_ionic, m_eff, α, phonon_freq, phonon_freq, 1, 1)
end

"""
    material(ϵ_optic, ϵ_static, m_eff, phonon_freqs, ir_activity, volume)

    Construct a 'struct Material' object from a set of infrared activities and phonon
    frequencies, for use with the 'multiple phonon branches' extension of the code. 

    effective_freq calculated with the Hellwarth 'B' scheme; ϵ_ionic and α calculated from
    the infrared activties.  
"""
function material(ϵ_optic, ϵ_static, m_eff, phonon_freqs, ir_activity, volume)
    effective_freq = HellwarthBScheme(hcat(phonon_freqs, ir_activity))
    ϵ_ionic = ϵ_ionic_mode.(phonon_freqs, ir_activity, volume)
    α = frohlichalpha.(ϵ_optic, ϵ_ionic, sum(ϵ_ionic), phonon_freqs, m_eff)
    return Material(ϵ_optic, ϵ_static, ϵ_ionic, m_eff, α, phonon_freqs, effective_freq, ir_activity, volume)
end

function Base.show(io::IO, ::MIME"text/plain", x::Material)
    io_limit = IOContext(io, :compact => true, :limit => true)
    println("\e[K------------------------------------------")
    println("\e[K           Material Information           ")
    println("\e[K------------------------------------------")
    println(io_limit, "\e[KOptic dielectric   | ϵo = ", x.ϵo)
    println(io_limit, "\e[KStatic dielectric  | ϵs = ", x.ϵs)
    println(io_limit, "\e[KIonic dielectric   | ϵi = ", x.ϵi)
    println(io_limit, "\e[KBand mass          | mb = ", x.mb)
    println(io_limit, "\e[KFröhlich coupling  | α = ", x.α)
    println(io_limit, "\e[KPhonon frequencies | f = ", x.f)
    println(io_limit, "\e[KEff Phonon freq    | feff = ", x.feff)
    println(io_limit, "\e[KIR activities      | ir = ", x.ir)
    println(io_limit, "\e[KUnit cell volume   | V = ", x.V)
    println("\e[K-------------------------------------------")
end

