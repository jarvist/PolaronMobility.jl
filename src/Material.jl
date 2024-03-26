# Material.jl
# Calculates inpit variables for the polaron theories from material data and saves the material as a Material type.

"""
    mutable struct Material

    Structure to contain material-specific data. 

"""
mutable struct Material
    # General polaron material properties
    α       # Unitless electron-phonon coupling (Unitless)
    mb      # Effective band mass (me)
    f       # Phonon frequencies (THz2π)
    feff    # Effective frequency (THz2π)

    # Frohlich specific properties for calculating electron-phonon coupling α
    ϵo      # Optical dielectric constant
    ϵs      # Static dielectric constant
    ϵi      # Ionic dielectric contributions

    # Frohlich specific properties for calculating multiple phonon electron-phonon coupling αⱼ
    ir      # Infrared activities
    V       # Unit cell volumes

    # Holstein specific properties
    γ       # Adiabaticity (Unitless)
    J       # Transfer integral (meV)
    a       # Lattice constant (Angstroms)
    g       # Electron-phonon coupling energy (Unitless)
    z       # Lattice dimensionality (Unitless)

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
    return Material(α, m_eff, phonon_freq, phonon_freq, ϵ_optical, ϵ_static, ϵ_ionic, 1, 1, 1, 1, 1, 1, 1)
end

"""
material(ϵ_optic, ϵ_static, m_eff, phonon_freqs, ir_activity, volume)

    Construct a 'struct Material' object, for a multiple phonon mode Frohlich polaron, from a set of infrared activities and phonon
    frequencies, for use with the 'multiple phonon branches' extension of the code. 

    effective_freq calculated with the Hellwarth 'B' scheme; ϵ_ionic and α calculated from
    the infrared activties.  
"""
function material(ϵ_optic, ϵ_static, m_eff, phonon_freqs, ir_activity, volume)
    effective_freq = HellwarthBScheme(hcat(phonon_freqs, ir_activity))
    ϵ_ionic = ϵ_ionic_mode.(phonon_freqs, ir_activity, volume)
    α = frohlichalpha.(ϵ_optic, ϵ_ionic, sum(ϵ_ionic), phonon_freqs, m_eff)
    return Material(α, m_eff, phonon_freqs, effective_freq, ϵ_optic, ϵ_static, ϵ_ionic, ir_activity, volume, 1, 1, 1, 1, 1)
end

function material(transfer_integral, coupling_energy, phonon_freq, lattice_constant, lattice_dimensionality)

    transfer_integral *= eV / 1000
    coupling_energy *= eV / 1000
    lattice_constant *= 1e-10

    phonon_energy = ħ * phonon_freq * 1e12 * 2π
    dimensionless_coupling = coupling_energy / phonon_energy
    adiabaticity = phonon_energy / transfer_integral
    α = dimensionless_coupling * 2 / lattice_dimensionality
    band_mass = ħ^2 / transfer_integral / lattice_constant^2 / 2 / me
    
    return Material(α, band_mass, phonon_freq, phonon_freq, 1, 1, 1, 1, lattice_constant^(lattice_dimensionality/2), adiabaticity, transfer_integral / eV * 1000, lattice_constant * 1e10, dimensionless_coupling, lattice_dimensionality)
end

function Base.show(io::IO, ::MIME"text/plain", x::Material)
    io_limit = IOContext(io, :compact => true, :limit => true)
    println("\e[K------------------------------------------")
    println("\e[K           Material Information           ")
    println("\e[K------------------------------------------")
    println(io_limit, "\e[KUnitless Coupling  | α = ", x.α)
    println(io_limit, "\e[KBand mass          | mb = ", x.mb)
    println(io_limit, "\e[KPhonon frequencies | f = ", x.f)
    println(io_limit, "\e[KEff Phonon freq    | feff = ", x.feff)
    println(io_limit, "\e[KOptic dielectric   | ϵo = ", x.ϵo)
    println(io_limit, "\e[KStatic dielectric  | ϵs = ", x.ϵs)
    println(io_limit, "\e[KIonic dielectric   | ϵi = ", x.ϵi)
    println(io_limit, "\e[KIR activities      | ir = ", x.ir)
    println(io_limit, "\e[KUnit cell volume   | V = ", x.V)
    println(io_limit, "\e[KAdiabaticity       | γ = ", x.γ)
    println(io_limit, "\e[KTransfer integral  | J = ", x.J)
    println(io_limit, "\e[KLattice constant   | a = ", x.a)
    println(io_limit, "\e[KEl-ph coupling     | g = ", x.g)
    println(io_limit, "\e[KLattice Dimensions | z = ", x.z)

    println("\e[K-------------------------------------------")
end

function save_material(material::Material, prefix)

    println("Saving material data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "alpha", material.α,
        "band mass", material.mb,
        "phonon freq", material.f,
        "phonon freq eff", material.feff,
        "optical dielectric", material.ϵo,
        "static dielectric", material.ϵs,
        "ionic dielectric", material.ϵi,
        "infrared activity", material.ir,
        "unitcell volume", material.V,
        "adiabaticity", material.γ,
        "transfer integral", material.J,
        "lattice constant", material.a,
        "dimensionless coupling", material.g,
        "lattice dimension", material.z
    )

    println("... Material data saved.")
end

function load_material(material_file_path)

    println("Loading material data from $material_file_path ...")

    data = JLD.load("$material_file_path")

    material = Material(
        data["alpha"],
        data["band mass"],
        data["phonon freq"],
        data["phonon freq eff"],
        data["optical dielectric"],
        data["static dielectric"],
        data["ionic dielectric"],
        data["infrared activity"],
        data["unitcell volume"],
        data["adiabaticity"],
        data["transfer integral"],
        data["lattice constant"],
        data["dimensionless coupling"],
        data["lattice dimension"]
    )
    println("... Material loaded.")

    return material
end
