SimuPy Flight Vehicle Toolkit
=============================

Vehicle flight simulation is an important part of the innovation of aerospace vehicle technology. The NASA Engineering Safety Center (NESC) has identified and addressed the need to verify flight vehicle simulations through their work on the `“Six degree-of-freedom (6-DOF) Flight Simulation Test Cases.” <https://nescacademy.nasa.gov/flightsim/>`_ The author was prompted to develop tools that would allow for the rapid implementation of simulations for novel vehicle concepts including hypersonic re-entry vehicles and urban air mobility vehicles.

This software library leverages open source scientific computing tools to implement an efficient simulation framework for flight vehicles in Python. Equations of motion are composed in blocks using the SimuPy library, an open source Python alternative to Simulink, solves the resulting differential equations using SciPy’s wrappers for standard Fortran implementations. Dynamics equations of the inertial state variables for the position, orientation, and their corresponding rates for integration are developed using the SymPy symbolic library and implemented using code generation. Kinematics equations are implemented through symbolic definition and code generation as well as leveraging other open source software that implements useful functions, such as the solutions to the inverse geodesy problem.

Library Design
--------------

.. |SimuPyAPI| replace:: SimuPy's ``DynamicalSystem`` interface
.. _SimuPyAPI: https://simupy.readthedocs.io/en/latest/api/api.html>

The SimuPy Flight Vehicle Toolkit API is designed to prioritize explicit and modular models, to facilitate adherance to best practices from the scientific Python community, and to follow SimuPy design principles. The library provides a ``Planet`` and ``Vehicle`` that implement |SimuPyAPI|_. This API separation was designed to reflect the modeling distinction between the two bodies being simuated. This separation could have been reflected in other ways; indeed within the ``Planet`` and ``Vehicle`` classes there is a namespace separation between model components. However, it is often convenient to break SimuPy models into distinct systems when intermediate signals may be used for multiple purposes. In particular, the separation allows easy access to the vehicle's sensed acceleration which is an important measurement for flight vehicles, particularly entry vehicles. Due to the disciplined modeling approach set out by SimuPy where the output of systems with state are not dependent on the input, a vehicle whose input changed acceleration could not be modeled without separating the ``Vehicle``. In the future, convenience functions will provided for common ``BlockDiagram`` manipulations like connecting guidance, navigation, or control models to reduce the learning curve for new users.

The ``Vehicle`` class provides the state-less dynamics (angular and translational acceleration) of the flight vehicle. It is constructed by

.. code:: python

    vehicle = Vehicle(
        m, I_xx, I_yy, I_zz, I_xy, I_yz, I_xz, x_com, y_com, z_com, # intertia
        base_aero_coeffs, x_mrc, y_mrc, z_mrc, S_A, a_l, b_l, c_l, d_l, # aero
        input_aero_coeffs, input_force_moment, # extra callbacks for control modeling
        input_aero_coeffs_idx, input_force_moment_idx, # routing for control callbacks
        dim_extra_input, # total number of extra (control) inputs
    )

The ``Planet`` class provides the kinematic equations of motion used to drive the position,  orientaiton, and their corresponding rates based on the dynamics modeled by the ``Vehicle``. The ``Plenet`` is constructed by 

.. code:: python

    planet = Planet(gravity, winds, atmosphere, planetodetics)

For more details, see the docstrings (via ``help(Vehicle)`` and ``help(Planet)`` respectively) and NESC Test Case examples.

NESC Test Cases
---------------

..
    TODO: update this section -- retain reference to NESC data and explain
    the local copy; possibly provide the argument information

A number of the NESC Atmospheric test cases have been implemented to verify the implementation and derivation of the equations of motion. These are located in the `nesc_test_cases` directory. To run, simply execute any of `nesc_case##.py` files or the `run_nesc_cases.py` which will iterate through test cases that have been implemented. These scripts will attempt to load the NESC reference results from the parent directory and plot the results along with the results from the SimuPy implemntation. To include the NESC results in the comparison plots, download the `Atmospheric trajectory data <https://nescacademy.nasa.gov/src/flightsim/Datasets/Atmospheric_checkcases.zip>`_ and unzip the `Atmospheric_checkcases` directory to the root `simupy_flight` directory. You can place the `Atmospheric_checkcases` directory in different location by changing the `data_relative_path` variable in the `nesc_testcase_helper.py` script.

License
-------

This software is released under the `NASA Open Source Agreement Version 1.3 <https://github.com/nasa/simupy-flight/raw/master/license.pdf>`_.


Notices
-------

Copyright © 2021 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  All Rights Reserved.

Disclaimers
-----------

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
