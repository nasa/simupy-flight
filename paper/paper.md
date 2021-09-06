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

Vehicle flight simulation is an important part of the innovation of aerospace vehicle technology. The NASA Engineering Safety Center (NESC) has identified and addressed the need to verify flight vehicle simulations through their work on the flight simulation test cases [@nesc1]. The SimuPy Flight Vehicle Toolkit provides a modular framework for the rapid implementation of simulations for novel flight vehicle concepts including hypersonic re-entry vehicles and urban air mobility vehicles. The open source repository of the source code includes implementations for the sixteen atmospheric test cases defined by the NESC. These implementations serve as validation of the simulation framework and provides example code for its usage.

This software library leverages open source scientific computing tools to implement an efficient simulation framework for flight vehicles in Python. Equations of motion are composed in blocks using the SimuPy library [@simupy], an open source Python alternative to Simulink. The resulting differential equations are solved using SciPy’s wrappers for standard Fortran implementations [@scipy]. Equations of motion for the inertial state of a rigid-body model of the vehicle representing the position, orientation, and their corresponding rates for integration are developed using the SymPy symbolic library [@sympy] and implemented using code generation. Kinematics equations are implemented through symbolic definition and code generation. Open-source scientific libraries are leveraged where possible, such as solving the inverse geodesy problem [@pyerfa]. The library also provides a parser for the American Institute of Aeronautics and Astronautics's (AIAA) simulation description mark-up language standard [@daveml] using code generation. Aerodynamic data table interpolation is implemented using ndsplines [@ndsplines].


# References
