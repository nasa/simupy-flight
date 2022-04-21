---
title: 'SimuPy Flight Vehicle Toolkit'
tags:
  - Python
  - flight vehicle
  - simulation
authors:
  - name: Benjamin W. L. Margolis
    orcid: 0000-0001-5602-1888
    affiliation: 1
  - name: Kenneth R. Lyons
    orcid: 0000-0002-9143-8459
    affiliation: 1
affiliations:
 - name:  NASA Ames Research Center, Systems Analysis Office
   index: 1
date: 06 September 2020
bibliography: references.bib
---

# Summary

Vehicle flight simulation is an important part of the innovation of aerospace vehicle technology. Built on the SimuPy modeling and simulation framework [@simupy], the SimuPy Flight Vehicle Toolkit provides a modular framework for rapid implementation of simulations for novel flight vehicle concepts, such as hypersonic re-entry vehicles or urban air mobility vehicles. The open source repository of the source code includes implementations for the sixteen atmospheric test cases defined by the NASA Engineering Safety Center (NESC), which serve as validation of the simulation framework and examples of its usage.

# Statement of Need

The NESC has identified and addressed the need to verify flight vehicle simulations through their work on the flight simulation test cases [@nesc1]. In that work, the NESC established flight vehicle simulation test cases to compare and validate a suite of simulation tools, several from within NASA and one external, open-source tool. Implementations of the NESC test cases via the SimuPy Flight Vehicle Toolkit's API help verify its correctness and demonstrate its effectiveness in succinctly constructing flight vehicle simulations.

One author has used a precursor to this software package to simulate control system performance for a novel mechancially deployed hypersonic entry vehicle [@d2019developing; @margolis2019control; @margolis2019iac; @okolo2020scitech; @margolis2020fuzzy; @d2021scitech; @margolis2021scitech].

# Description

The SimuPy Flight Vehicle toolkit provides a modular programming interface for specifying vehicle and planet characteristics, such as aerodynamic coefficients and a gravity model, depicted in \autoref{fig:diagram}. Equations of motion are composed in blocks using the SimuPy library [@simupy], an open source alternative to Simulink. These blocks can be formed into a standalone block diagram to simulate free behavior of the vehicle, or they can be incorporated into a more complex model. Implementations of the NESC test cases provided with the toolkit demonstrate usage for increasingly complex models, from a free-falling sphere in Earth atmosphere to a maneuvered F-16.

![Simupy Flight Vehicle Toolkit architecture\label{fig:diagram}](simupy_flight_diagram.pdf)

The SimuPy Flight Vehicle Toolkit leverages open source scientific computing tools to implement an efficient simulation framework for flight vehicles in Python. Differential equations are solved using SciPy's wrappers for standard Fortran implementations [@scipy]. Equations of motion for the inertial state of a rigid-body model of the vehicle representing the position, orientation, and their corresponding rates for integration are developed using the SymPy symbolic library [@sympy] and implemented using code generation. Kinematics equations are implemented through symbolic definition and code generation. Open-source scientific libraries are leveraged where possible, such as solving the inverse geodesy problem [@pyerfa] and implementing a standard atmosphere model [@fluids]. The library also provides a parser for the American Institute of Aeronautics and Astronautics's (AIAA) simulation description mark-up language standard [@daveml] using code generation. Aerodynamic data table interpolation is implemented using ndsplines [@ndsplines].

# References
