Perhaps the easiest way to understand the code is to see how it can be used for
science. 
As an example system, we are going to look at some of the basic polaron
properties of methylammonium lead-iodide perovskite.

## Loading the Module

However you run Julia (whether at the REPL; within a Jupyter notebook; or in
a standalone script/program), you will first need to load the module. 
As the Module is not yet registered, you will first need to tell Julia where
to find the codes.

```julia
push!(LOAD_PATH,"../src/") # load module from local directory
using PolaronMobility 
```

## α/alpha parameter 

The Frohlich electron-phonon coupling parameter can be characterised by
a dimensionless coupling, alpha (`α`). 
This gives the long-range ('non analytic') contribution from electrodynamic
coupling into infrared active phonon modes.

As an example, let us calculate α for CdTe, and compare it to literature
values.

The call signature is (ϵ-optical, ϵ-static, phonon-frequency (Hz),
effective-mass (in mass-of-electron units)).

```
α=frohlichalpha(7.1, 10.4, 5.08E12, 0.095)
println("CdTe  α=",α," Stone 0.39 / Devreese 0.29")
@test α ≈ 0.3 atol=0.1
```

We get a value of `0.351`. 

## Single temperature phonon properties

Let us calculate the room-temperature (300 K) character of the electron-polaron
in methylammonium lead iodide perovskite (MAPI). 
The parameters we use are as in Frost2017PRB.

The call signature is (Temperature range, ϵ-optical, ϵ-static, phonon-frequency
(Hz), effective-mass (in mass-of-electron units)). 
For electrons in MAPI, these are ϵ=4.5/24.1, f=2.25 THz, me=0.12 electron
masses. 


```
MAPIe=polaronmobility(300, 4.5, 24.1, 2.25E12, 0.12)
```

This will think for a bit (as Julia just-in-time compiles the required
functions), and then spits out a considerable amount of information to
`STDOUT`. 

```
Polaron mobility for system ε_Inf=4.5, ε_S=24.1, freq=2.25e12, 
                 effectivemass=0.12; with Trange 300 ...
Polaron mobility input parameters: ε_Inf=4.500000 ε_S=24.100000 freq=2.25e+12 α=2.393991 
Derived params in SI: ω =1.41372e+13 mb=1.09313e-31 
T: 300.000000 β: 2.41e+20 βred: 0.36 ħω  = 9.31 meV		Converged? : true
 VariationalParams v= 19.86  w= 16.96	||	 M=0.371407 k=106.835753	
 Polaron frequency (SI) v=  4.5e+13 Hz 	 w=  3.8e+13 Hz	
 Feynman1955(46,47): meSmallAlpha(α)= 0.542 meLargeAlpha(α)= 0.067
 Feynman1962: Approximate ~ Large alpha limit, v/w = 1.17  =~approx~= alpha^2 = 5.73 
 POLARON SIZE (rf), following Schultz1959. (s.d. of Gaussian polaron ψ )
	 Schultz1959(2.4): rf= 0.528075 (int units) = 2.68001e-09 m [SI]
	 Schultz1959(2.5a) with 0.44: Feynman α→0 expansion: rfa= 1.68761 (int units) = 2.96691e-09 m [SI]
	 Schultz1959(2.5a) with 4/9 re-derivation: Feynman α→0 expansion: rfa= 1.67915 (int units) = 2.95204e-09 m [SI]
	 Schultz1959(2.5b): Feynman α→∞ expansion: rf= 9.00127 (int units) = 1.58247e-08 m [SI]
	 Schultz1959(5.7-5.8): fixed-e: phononfreq= 16.9603 (int units) = 3.81607e+13 [SI, Hz] = 157.82 [meV]
	 Schultz1959(5.7-5.8): reducd mass: phononfreq= 19.8617 (int units) = 4.46888e+13 [SI, Hz] = 184.818 [meV]
	 Schultz1959: electronfreq= 10.3361 (int units) = 2.32563e+13 [SI, Hz] = 96.1804 [meV]
	 Schultz1959: combinedfreq= 8.82623 (int units) = 1.9859e+13 [SI, Hz] = 82.1303 [meV]
 Devreese1972: (Large Alpha) Franck-Condon frequency = 0.81
 Polaron Free Energy: A= -6.448815 B= 7.355626 C= 2.911977 F= -3.818788	 = -35.534786 meV
Polaron Mobility theories:
	μ(FHIP)= 0.082049 m^2/Vs 	= 820.49 cm^2/Vs
	μ(Kadanoff,via Devreese2016)= 0.019689 m^2/Vs 	= 196.89 cm^2/Vs
	 Eqm. Phonon. pop. Nbar: 2.308150 
	exp(Bred): 1.433247 exp(-Bred): 0.697716 exp(Bred)-1: 0.433247
	μ(Kadanoff1963 [Eqn. 25]) = 0.019689 m^2/Vs 	 = 196.89 cm^2/Vs
	Gamma0 = 5.4282e+13 rad/s = 8.63925e+12 /s 
	Tau=1/Gamma0 = 1.15751e-13 = 0.115751 ps
	Energy Loss = 1.288e-08 J/s = 80.3904 meV/ps
	μ(Hellwarth1999)= 0.013642 m^2/Vs 	= 136.42 cm^2/Vs
	μ(Hellwarth1999,b=0)= 0.013663 m^2/Vs 	= 136.63 cm^2/Vs
	Error due to b=0; 0.001539
```

The output is a little ad-hoc, and specific values are perhaps best understood
with comparison to the code, and to the references to the original papers!

Initially the polaron calculation is made; this involves varying `v` and `w` to
minimise the miss-match between the trial (analytically solvable) polaron
Hamiltonian energy, and the true temperature-dependent free-energy (as
specified by Osaka). 
The method uses automatic differentiation to get gradients for the optimisation
procedure. 

'Textbook' expressions that predict polaron character and mobilities make
assumptions about `v` and `w` (usually that either `v` is small, or `v=w`).

Values that can be directly derived from these `v` and `w` variational
parameters are then displayed.  Essentially we are just using `Julia` as
a glorified scientific calculator at this point, but with the units checked. 
Perhaps most of interest are the Schultz polaron 'size', various resonant
frequencies, and the polaron free energies. 

The polaron theories are constructed in reduced units. Generally this means
that energy is in units of ħω, and frequencies in a unit of ω (of the input
phonon frequency). For convenience, these are re-printed in SI or more standard
units. 

Beyond `Polaron Mobility theories:`, the code enters its final phase and uses
the `v` and `w` parameters specifying the polaron to calculate a charge
carrier mobility. 
The asymtotic 'FHIP' mobility (low T) is calculated, this can be most easily
related to textbook expressions that directly infer a mobility from an `α`
parameter. It lacks optical phonon emission, and so shows pathological high
temperature (kT > ħω) behaviour. 
The Kadanoff mobility (see the original paper) improves on this by assuming
a Boltzmann process (independent scattering events). 
From this theory we can also get an average scattering time, which we relate to
the time-scale of the polaron interacting with the phonon cloud, and so to the
rate of polaron cooling. 
Finally the Hellwarth1999 scheme is used, which goes back to the original 1962
FHIP paper, and directly carries out the contour integral for the polaron
impedence function. We improve on this slightly by explicitly calculating with
`b`, though the approximation `b=0` makes very little difference for any so-far
tested materials. 

## Hellwarth's multi-mode scheme

The above examples are slightly back-to-front - in that we've specified
a single mode frequency, as if the material were a simple tetrahedral
semiconductor with only one infrared active mode. 
(The Linear Optical 'LO' phonon mode.)

In order to use these theories with more complex (many atoms in a unit cell)
materials of technological relevance, we must first reduce all of these
Infrared-active phonon responses to a single effective one. 

For this we will use the averaging scheme described in Hellwarth1999. 
Currently only the B scheme (athermal) is correctly implemented; a partial
A scheme implementation is present.

Let's test it against the Hellwarth1999 literature data. 
The argument to the function is a table of frequencies (cm^-1) and infrared
activities (unit does not matter, as long as it is consistent). 

```julia
# Hellwarth et al. PRB 1999 Table II - BiSiO frequencies and activities
HellwarthII = [
    106.23 8.86
    160.51 9.50
    180.33 20.85
    206.69 10.05
    252.76 27.00
    369.64 61.78
    501.71 52.87
    553.60 86.18
    585.36 75.41
    607.29 98.15
    834.53 89.36
]

println("Attempting to reproduce Hellwarth et al.'s data.")
println("\nB scheme: (athermal)")
HellwarthBScheme(HellwarthII)
println("    ... should agree with values given in Hellwarth(60) W_e=196.9 cm^-1 and Hellwarth(61) Ω_e=500 cm^-1")
```

The output agrees to within three significant figures with the literatures values;
```
Hellwarth (58) summation: 0.15505835776181887
Hellwarth (59) summation (total ir activity ^2): 38777.7725
Hellwarth (59) W_e (total ir activity ): 196.92072643579192
Hellwarth (61) Omega (freq): 500.08501275972833
```

## Temperature-dependent behaviour

Getting temperature-dependent behaviour is a matter of sending a temperature
range to the polaronmobility funtion.

```julia
MAPIe=polaronmobility(10:10:1000, 4.5, 24.1, 2.25E12, 0.12)
savepolaron("MAPI-electron",MAPIe)
```

`savepolaron` saves a column-delimited text file for post-production plotting
(with gnuplot) or similar.

```julia
using PlotPolaron
plotpolaron("MAPI-electron", MAPIe)
```
The convenience function `plotpolaron` generates (and saves) a number of
`Plots.jl` figures of the temperature dependent behaviour.
It has been separated off into its own submodule (`PlotPolaron`), so that the
`Plots.jl` dependency does not slow down general loading of the package.

## Further examples

More complete examples are provided in the examples/ and HalidePerovskites/
folders of the main PolaronMobility.jl repository.

