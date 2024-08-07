---
title: 'TSE: A triple stellar evolution code'
tags:
  - Python
  - Fortran
  - astronomy
  - dynamics
  - few-body dynamics
  - stars
authors:
  - name: Jakob Stegmann
    orcid: 0000-0003-2340-8140
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Max Planck Research Fellow, Max Planck Institute for Astrophysics, Germany
   index: 1
date: 16 July 2024
bibliography: paper.bib

---

# Summary

Most massive stars are found in hierarchical triples or higher multiplicity systems
in which a close inner binary is orbited by one or more distant companions [@Moe:2017].
A distant companion may significantly alter the evolution of the system by causing
large-amplitude oscillations of the inner binary eccentricity (so-called von Zeipel-
Kozai-Lidov oscillations) [@Naoz:2016]. These oscillations determine how close the 
inner binary stars pass each other at periapsis and are therefore essential for
our understanding of interactions between massive binary stars [@Sana:2012] as well as 
their compact remnants. Understanding the role of tertiary companions for the evolution 
of massive stars requires an efficient numerical tool to simulate the complex interplay 
between the dynamics of hierachical triples and their stellar evolution.

# Statement of need

`TSE` is a Python code for hierarchical triple stellar dynamics. At its core, 
the secular equations of motions for the orbital elements of the inner and outer
orbits [@Liu:2015] and the spin vectors of the inner binary stars [@Hamers:2021] 
are evolved using the IVP solver implemented in `scipy.integrate`. Optionally, it
follows the trajectory of the triple barycentre throughout the Milky Way and includes 
the perturbative effect of Galactic tides on the orbital elements of the inner and 
outer binary [@Bub:2020; @Stegmann:2024]. The evolution of stellar properties such as 
masses and radii are followed by using the stellar evolution code `MOBSE` 
[@Giacobbo:2018; @Giacobbo:2019; @Giacobbo:2023] which builds upon the binary stellar 
evolution code `BSE` written in Fortran [@Hurley:2000; @Hurley:2002]. While integrating 
the equations of motion `TSE` constantly checks for events such as Roche-lobe overflow, 
supernova explosions, stellar or compact object mergers, orbital disruptions and other 
user-defined custom events. In either case, `TSE` either models the impact of the event 
on the triple, e.g., by simulating the result of Roche-lobe overflow using `MOBSE`,
or terminates the evolution and stores the final outcome. 

`TSE` was designed to study the impact of a tertiary companion on massive stellar
evolution, and has been already used in a number of scientific publications 
[@Stegmann:2022a; @Stegmann:2022b]. While previous work mostly focused on the long-term 
evolution of black hole triples towards gravitational-wave sources (e.g., 
@Rodriguez:2018), already including the dynamical effect of a tertiary companion during
the lifetime of the progenitors has been largely neglected. Thus, `TSE` can be 
used to investigate the role of tertiary companions for a range of massive star 
phenomena, such as X-ray binaries, stellar mergers, and gravitational-wave
sources. There are only few other codes published which combine the dynamical evolution
of triples with that of stars in a comparable way than `TSE` (`TrES` 
[@Toonen:2016; @Toonen:2023] and `MSE` [@Hamers:2021]). Comparing the outcomes of the 
codes will be an essential part in order to understand how different model assumptions
determine the evolution of stellar triples.

# Acknowledgements

I acknowledge contributions from Fabio Antonini and Jordan Barber and support from 
Michela Mapelli and Max Moe during the genesis of this project. This work was largely
carried out at Cardiff University and was funded in parts by the Netherlands Organisation 
for Scientific Research (NWO), as part of the Vidi research program BinWaves (project 
number 639.042.728, PI: de Mink).

# References