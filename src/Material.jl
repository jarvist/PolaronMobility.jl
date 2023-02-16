# Material.jl

mutable struct Material
    optical     # Optical dielectric constant
    static      # Static dielectric constant
    ionic       # Ionic dielectric contributions
    mb          # Effective band mass
    α           # Fröhlich coupling
    freqs       # Phonon frequencies
    freq_eff    # Effective frequency
    ir          # Infrared activities
    volume      # Unit cell volumes
    function Material(x...)
        new(reduce_array.(x)...)
    end
end

function material(ϵ_optical, ϵ_static, m_eff, phonon_freq)
    ϵ_ionic = ϵ_static - ϵ_optical
    α = frohlichalpha(ϵ_optical, ϵ_static, phonon_freq, m_eff)
    return Material(ϵ_optical, ϵ_static, ϵ_ionic, m_eff, α, phonon_freq, phonon_freq, 1, 1)
end

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
    println(io_limit, "\e[KOptic dielectric   | ϵₒ = ", x.optical)
    println(io_limit, "\e[KStatic dielectric  | ϵₛ = ", x.static)
    println(io_limit, "\e[KIonic dielectric   | ϵᵢ = ", x.ionic)
    println(io_limit, "\e[KBand mass          | mb = ", x.mb)
    println(io_limit, "\e[KFröhlich coupling  | α = ", x.α)
    println(io_limit, "\e[KPhonon frequencies | f = ", x.freqs)
    println(io_limit, "\e[KIR activities      | IR = ", x.ir)
    println(io_limit, "\e[KUnit cell volume   | V₀ = ", x.volume)
    println("\e[K-------------------------------------------")
end

