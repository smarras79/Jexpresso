---
title: 'Jexpresso: a fast, modular and general solver of partial differential equations on CPUs and GPUs using Julia'
tags:
  - julia
  - pdes
  - partial differential equations
  - spectral elements
  - finite differences
  - finite volumes
authors:
  - name: Simone Marras
    orcid: 0000-0002-7498-049X
    affiliation: "1, 2"
  - name: Yassine Tissaoui
    orcid: 0000-0000-0000-0000
    affiliation: 3
  - name: Hang Wang
     orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: Department of Mechanical and Industrial Engineering, New Jersey Institute of Technology, Newark (NJ), USA
   index: 1
 - name: Center for Applied Mathematics and Statistics, New Jersey Institute of Technology, Newark (NJ), USA
   index: 2
 - name: Department of Mathematics, University of Wisconsit, Madison (WI), USA
   index: 3
date: 19 Mau 2024
bibliography: paper.bib
---

# Summary

 (see [@tissaoui2024])

![](code.pdf)

Other FE packages like FEniCS also achieve such compact user interfaces, but in contrast to Gridap, they are based on a sophisticated compiler of variational forms [@Kirby2006], which generates, compiles and links a specialized C++ back-end for the problem at hand. One of the limitations of this approach is that the form compiler is a rigid system that is not designed to be extended by average users.

Jexpresso is an open-source software hosted at Github and distributed under an MIT license.

# References
