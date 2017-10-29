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

This package is an implementation in the Julia programming language of solving
the Feynman @Feynman1955 variational path-integral solution to the Frohlich
@Frohlich1952 Hamiltonian specifying the polaron problem. 
Finite-temperature (free) energies of Osaka @Osaka1959 are used, in the modern
presentation of Hellwarth et al.@Hellwarth1999. 

The codes are designed to have the physical system specified from ab-initio
calculations. 
The inputs are a bare-band effective-mass; high-frequency and zero-frequency
dielectric constants, and an effective dielectric phonon frequency. 
These lattice-dependent characteristics can all be derived from an ab-initio
calculation of the infrared (dielectric) properties of the harmonic phonons. 

Having solved for this ab-initio specified, finite temperature, polaron state,
the codes can then calculate various parameters of interest for device physics. 
Most notably polaron mobilities in the original FHIP @Feynman1962 asymtotic
limit, the Kadanoff @Kadanoff1963 Boltzmann formulation and the most recent
Hellwarth et al. @Hellwarth1999 explicit contour integral forms. 
The size and nature of the polaron state is also described, most of which was
previously investigated by Schultz @Schultz1959.

These codes were developed for and enabled a recent publication by Frost
@FrostPolaronMobility2017, which provided calculated temperature-dependent
mobilities and polaron configurations in the halide perovskite family of
semiconductors.

In providing robust codes to calculate the polaron state, this work enables
calculation of further parameters such as the nature of polaron scattering,
frequency-dependent mobility and polaron optical absorption.

((Spiel from original website below.))

- A summary describing the high-level functionality and purpose of the software
for a diverse, non-specialist audience
- A clear statement of need that illustrates the purpose of the software
- A list of key references including a link to the software archive
- Mentions (if applicable) of any ongoing research projects using the software
or recent scholarly publications enabled by it

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

Figures can be included like this: ![Fidgit deposited in figshare.](figshare_article.png)

# References
