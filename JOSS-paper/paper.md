---
title: 'PolaronMobility.jl: Implementation of the Feynman variational polaron model'
tags:
  - polaron
  - path-integral
  - mobility
authors:
 - name: Jarvist Moore Frost 
   orcid: 0000-0003-1938-4430 
   affiliation: 1, 2
affiliations:
 - name: Department of Chemistry, University of Bath, UK.
   index: 1
 - name: Department of Materials, Imperial College London, UK.
   index: 2
date: 28 October 2017
bibliography: paper.bib
---

# Summary

An additional electron (i.e. a charge carrier) in a material applies
a polarisation field. 
In a polar material (one where the ions of the lattice have different charge),
this results in a large coupling into the lattice motion. 
Small distortions in the lattice can be described in terms of a harmonic
restoring force, and thereby a harmonic quasi-particle of vibration termed
a phonon. 
Electron-phonon coupling provides a scattering mechanism by which momentum
(kinetic-energy) of an electron is dissipated. 
These scattering processes limit charge-carrier mobility in materials. 
Being able to predict charge-carrier mobility helps the computational design of
new technological materials. 
Many new materials of potential utility and interest are polar. 
These include oxides used as battery anodes and transparent conductors for
electronic displays; and chalcogenides and halides indicated for light emission
(displays, lighting) and absorption (photovoltaic solar cells).

In a polar material the dielectric electron-phonon coupling dominates the
electronic scattering. 
Unusually, this scattering process can be modelled without any empirical
parameters, and so the temperature-dependent absolute-mobility of a polar
material can be calculated. 

This package is an implementation in the Julia programming language of solving
the Feynman [@Feynman1955] variational path-integral solution to the Frohlich
[@Frohlich1952] Hamiltonian specifying the polaron problem. 
The Frohlich Hamiltonian is very simple. 
Electrons are treated at the effective-mass (quadratic dispersion relationship)
level, and the vibrational response of the material as a single
effective-frequency harmonic mode. 
Finite-temperature (free) energies of Osaka [@Osaka1959] are used, as tabulated
in the modern presentation of Hellwarth et al. [@Hellwarth1999]. 
Hellwarth et al. also provides a rigorous method to calculate an effective
frequency. 

The physical system is specified by four parameters. 
1) bare-band effective-mass; 2) high-frequency and 3) zero-frequency
dielectric constants, and 4) an effective dielectric phonon frequency. 
These are most easily calculated by electronic structure calculations.
Components 3 and 4 are dependent on the lattice response, and can be derived
from calculation of the infrared (dielectric) properties of the harmonic
phonons. 

The Feynman polaron model integrates through the infinite quantum field of
these (as specified) lattice vibrations. 
The method is variational, and consists of an optimisation of the
finite-temperature (free) energies. 

Having solved for the finite temperature polaron state,
the codes can then calculate various parameters of interest for device physics. 
Most notably polaron mobilities in the original FHIP [@Feynman1962] asymtotic
limit, the Kadanoff [@Kadanoff1963] Boltzmann formulation and the most recent
Hellwarth et al. [@Hellwarth1999] explicit contour integral forms. 
The size and nature of the polaron state is also described, most of which was
previously investigated by Schultz [@Schultz1959].

These codes were developed for and enabled a recent publication by Frost
[@FrostPolaronMobility2017],  which provided calculated temperature-dependent
mobilities and polaron configurations for the halide perovskite family of
semiconductors.

In providing robust codes to calculate the polaron state, this work enables
calculation of further parameters such as the nature of polaron scattering,
frequency-dependent mobility and polaron optical absorption. 

# References
