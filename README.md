# PolaronMobility-FeynmanKadanoffOsakaHellwarth

[![Build Status](https://travis-ci.org/jarvist/PolaronMobility.jl.svg?branch=master)](https://travis-ci.org/jarvist/PolaronMobility.jl)
[![Coverage Status](https://coveralls.io/repos/jarvist/PolaronMobility.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jarvist/PolaronMobility.jl?branch=master)
[![codecov.io](http://codecov.io/github/jarvist/PolaronMobility.jl/coverage.svg?branch=master)](http://codecov.io/github/jarvist/PolaronMobility.jl?branch=master)



These codes calculate the temperature-dependent polaron mobility for
a material.
We have parameters for various metal-halide Perovskite.

More extensive documentation is available
[![here](https://img.shields.io/badge/docs-latest-blue.svg)](https://jarvist.github.io/PolaronMobility.jl/),
which is perhaps easiest to read and understand alongside the pre-print:
[ArXiv:1704.05404](https://arxiv.org/abs/1704.05404) .

The scientific inputs are the dielectric constants, a characteristic phonon
frequency, and the bare-electron band effective-mass.

From this you can solve a temperature-dependent polaron model (variationally
optimising temperature-dependent free-energies for the coupled electron-phonon
system), and from this optimised parameters calculate derivative quantities
such as the polaron mobility, and polaron features (effective masses, energy
loss rates, oscillation etc.)

May your phonons drag in a manner truly sublime.

![MAPI Polaron mobility, plotted vs expt data](mobility-calculated-experimental.png)

## Research outputs

Polaron mobilities, three different ways
![Polaron mobilities, three different ways](mobility-calculated.png)

Effective mass of phonon cloud
![Effective mass of phonon cloud](mass.png)

Spring constant for coupling to phonon cloud
![Spring constant for coupling to phonon cloud](spring.png)

Variational (temperature-dependent free-energy) parameters for the coupled system
![Variational (temperature-dependent free-energy) parameters for the coupled system](variational.png)


