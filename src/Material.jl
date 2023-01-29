# Material.jl

struct Material
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
    flush(io)
    println("\e[K------------------------------------------")
    println("\e[K           Material Information           ")
    println("\e[K------------------------------------------")
    println(IOContext(io, :compact => true, :limit => true), "\e[KOptic dielectric    | ϵ∞ = ", x.optical, " ")
    println(IOContext(io, :compact => true, :limit => true), "\e[KStatic dielectric   | ϵ0 = ", x.static, " ")
    println(IOContext(io, :compact => true, :limit => true), "\e[KIonic dielectric    | ϵᵢ = ", x.ionic, " ")
    println(IOContext(io, :compact => true, :limit => true), "\e[KBand mass           | mb = ", x.mb, " mₑ")
    println(IOContext(io, :compact => true, :limit => true), "\e[KFröhlich coupling   | α = ", x.α)
    println(IOContext(io, :compact => true, :limit => true), "\e[KPhonon frequencies  | f = ", x.freqs, " THz")
    println(IOContext(io, :compact => true, :limit => true), "\e[KEffective frequency | f = ", x.freq_eff, " THz")
    println(IOContext(io, :compact => true, :limit => true), "\e[KIR activities       | IR = ", x.ir, " ")
    println(IOContext(io, :compact => true, :limit => true), "\e[KUnit cell volume    | V₀ = ", x.volume, " m³")
    println("\e[K-------------------------------------------")
end
