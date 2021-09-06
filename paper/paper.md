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

Vehicle flight simulation is an important part of the innovation of aerospace vehicle technology. The NASA Engineering Safety Center (NESC) has identified and addressed the need to verify flight vehicle simulations through their work on the `“Six degree-of-freedom (6-DOF) Flight Simulation Test Cases.” <https://nescacademy.nasa.gov/flightsim/>`_ The author was prompted to develop tools that would allow for the rapid implementation of simulations for novel vehicle concepts including hypersonic re-entry vehicles and urban air mobility vehicles.

This software library leverages open source scientific computing tools to implement an efficient simulation framework for flight vehicles in Python. Equations of motion are composed in blocks using the SimuPy library [@margolis2017simupy], an open source Python alternative to Simulink, solves the resulting differential equations using SciPy’s [@scipy] wrappers for standard Fortran implementations. Dynamics equations of the inertial state variables for the position, orientation, and their corresponding rates for integration are developed using the SymPy symbolic library and implemented using code generation. Kinematics equations are implemented through symbolic definition and code generation as well as leveraging other open source software that implements useful functions, such as the solutions to the inverse geodesy problem.

Also provides a parser for the American Institute of Aeronautics and Astronautics's (AIAA) simulation description mark-up language standard [@daveml] which usess SymPy [] to generate code. Aerodynamic data table interpolation is implemented with ndsplines [].

# References
