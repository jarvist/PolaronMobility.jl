These codes calculate the temperature-dependent polaron mobility for
a material.
We have parameters for various metal-halide Perovskite.

It is perhaps easiest to read and understand this documentation alongside the paper:
[Frost2017PRB](https://doi.org/10.1103/PhysRevB.96.195202)
([ArXiv:1704.05404](https://arxiv.org/abs/1704.05404) ).

The required inputs are the dielectric constants (ϵ-static and ϵ-optic)
, a characteristic phonon frequency (ω), and the bare-electron band
effective-mass (me).  These values can be relatively easily calculated in the
ab-initio electronic structure package of your choosing, or measured directly.

From these four values, the code solves a temperature-dependent polaron model.
This is done by variationally optimising the temperature-dependent
free-energies for the coupled
electron-phonon system.
These optimised parameters describe the polaron with the infinite quantum field
of lattice vibrations 'integrated through', and replaced with a phonon-drag
term.
From this the polaron features such as effective-mass, size of the
wavefunction, frequency of energy oscillation etc. can be calculated.

This polaron state can then be used as an input to further models for polaron
mobility.

The codes are designed to produce a set of temperature-dependent mobilities and
other data, for direct incorporation into a scientific publication.

May your phonons drag in a manner truly sublime.

