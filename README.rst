SimuPy Flight Vehicle Toolkit
=============================

.. image:: https://github.com/nasa/simupy-flight/actions/workflows/docs.yml/badge.svg
   :target: https://nasa.github.io/simupy-flight
.. image:: https://joss.theoj.org/papers/10.21105/joss.04299/status.svg
   :target: https://doi.org/10.21105/joss.04299

Vehicle flight simulation is an important part of the innovation of aerospace vehicle technology. The NASA Engineering Safety Center (NESC) has identified and addressed the need to verify flight vehicle simulations through their work on the `"Six degree-of-freedom (6-DOF) Flight Simulation Test Cases." <https://nescacademy.nasa.gov/flightsim/>`_ The author was prompted to develop tools that would allow for the rapid implementation of simulations for novel vehicle concepts including hypersonic re-entry vehicles and urban air mobility vehicles.

This software library leverages open source scientific computing tools to implement an efficient simulation framework for flight vehicles in Python. Equations of motion are composed in blocks using the SimuPy library, an open source Python alternative to Simulink, solves the resulting differential equations using SciPy’s wrappers for standard Fortran implementations. Dynamics equations of the inertial state variables for the position, orientation, and their corresponding rates for integration are developed using the SymPy symbolic library and implemented using code generation. Kinematics equations are implemented through symbolic definition and code generation as well as leveraging other open source software that implements useful functions, such as the solutions to the inverse geodesy problem.

Library Design
--------------

.. |SimuPyAPI| replace:: SimuPy's ``DynamicalSystem`` interface
.. _SimuPyAPI: https://simupy.readthedocs.io/en/latest/api/api.html

The SimuPy Flight Vehicle Toolkit API is designed to prioritize explicit and modular models, to facilitate adherence to best practices from the scientific Python community, and to follow SimuPy design principles. The library provides a ``Planet`` and ``Vehicle`` that implement |SimuPyAPI|_. This API separation was designed to reflect the modeling distinction between the two bodies being simulated. This separation could have been reflected in other ways; indeed within the ``Planet`` and ``Vehicle`` classes there is a namespace separation between model components. However, it is often convenient to break SimuPy models into distinct systems when intermediate signals may be used for multiple purposes. In particular, the separation allows easy access to the vehicle's sensed acceleration which is an important measurement for flight vehicles, particularly entry vehicles. Due to the disciplined modeling approach set out by SimuPy where the output of systems with state are not dependent on the input, a vehicle whose input changed acceleration could not be modeled without separating the ``Vehicle``. In the future, convenience functions will provided for common ``BlockDiagram`` manipulations like connecting guidance, navigation, or control models to reduce the learning curve for new users.

The ``Vehicle`` class provides the state-less dynamics (angular and translational acceleration) of the flight vehicle. It is constructed by

.. code:: python

    vehicle = Vehicle(
        m, I_xx, I_yy, I_zz, I_xy, I_yz, I_xz, x_com, y_com, z_com, # inertia
        base_aero_coeffs, x_mrc, y_mrc, z_mrc, S_A, a_l, b_l, c_l, d_l, # aero
        input_aero_coeffs, input_force_moment, # extra callbacks for control modeling
        input_aero_coeffs_idx, input_force_moment_idx, # routing for control callbacks
        dim_extra_input, # total number of extra (control) inputs
    )

The ``Planet`` class provides the kinematic equations of motion used to drive the position,  orientation, and their corresponding rates based on the dynamics modeled by the ``Vehicle``. The ``Planet`` is constructed by

.. code:: python

    planet = Planet(gravity, winds, atmosphere, planetodetics)

For more details, see the docstrings (via ``help(Vehicle)`` and ``help(Planet)`` respectively) and NESC Test Case examples.

Installation
------------

To install, clone the repository and install with pip

.. code:: bash

   $ git clone https://github.com/nasa/simupy-flight.git
   $ cd simupy-flight/
   $ pip install .

NESC Test Cases
---------------

A number of the `NESC Atmospheric test cases <https://nescacademy.nasa.gov/flightsim>`_ have been implemented to verify the implementation and derivation of the equations of motion. These are located in the ``nesc_test_cases`` directory. The implementations and their outputs can also be viewed in `the example gallery <https://nasa.github.io/simupy-flight/nesc_test_cases/index.html>`_ in the documentation. Additional information about each test case can be found in `the appendices of the NESC report <https://ntrs.nasa.gov/citations/20150001263>`_.

Running the examples requires a few extra dependencies. They can be be installed with:

.. code:: bash

    $ pip install .[test]

To run an example, simply execute any of the ``nesc_case##.py`` files or the ``run_nesc_cases.py`` script, which will iterate through test cases that have been implemented. These scripts load the NESC reference results from the ``NESC_data/`` directory and plot them along with the results from the SimuPy Flight implementation. The reference results are included in the SimuPy Flight repository, but they can be obtained directly from the NESC `here <https://nescacademy.nasa.gov/src/flightsim/Datasets/Atmospheric_checkcases.zip>`_.

.. note::

    The generated code performs a divide by zero if the velocity is zero,
    generating ``RuntimeWarning``\s. However, this condition is checked and
    handled correctly.

The SimuPy Flight results from running all NESC test cases are also included in the repository. By default, running any or all of the tests cases will perform a regression test against this data and report the result(s).

To re-generate the regression data, pass the ``--write-regression-data`` flag::

    $ python nesc_test_cases/run_nesc_cases.py --write-regression-data

Use ``-h`` or ``--help`` to see additional options that can be passed to the test case scripts.

Every case is annotated with at least a basic description adapted from the NESC reports. Cases 1-3 have moderate annotations to highlight basic API usage and modeling approaches. `Case 11 <https://nasa.github.io/simupy-flight/nesc_test_cases/nesc_case11.html#sphx-glr-nesc-test-cases-nesc-case11-py>`_, which demonstrates the trimming and straight and level flight of an F-16 model, is thoroughly annotated to illustrate how this simulation framework can be used for a sophisticated simulation. The F-16 vehicle model itself is also thoroughly annotated because it highlights how the ``Vehicle`` API can be adapted to alternate modeling approaches like the one used for the F-16 model implementation provided by the NESC.

DaveML Parsing
--------------

The American Institute of Aeronautics and Astronautics (AIAA) has developed a XML exchange format for aircraft simulation flight dynamics models called the `Dynamic Aerospace Vehicle Exchange Markup Language (DAVE-ML) <https://daveml.org/>`_. The ``parse_daveml`` submodule implements a parser that can be used to generate python code from valid DaveML. To use it, call the ``ProcessDaveML`` with a filename to the DaveML file. A python file will be created in the working directory with the same base-name as the DaveML file (replacing the extension, if any, with ``.py``). This feature was used to generate the vehicle models for the NESC test cases using the ``nesc_test_cases/process_NESC_DaveML.py`` script.

The DaveML specification includes elements for check-case data sets to assist in verification and debugging. The parser adds a function to each generated file, called ``run_checks``, which is executed when the file is run as a script. The NESC-provided F16 models include such data sets so they can be checked by running the generated files themselves. For example, to check the F16 aerodynamics model::

    $ python nesc_test_cases/F16_aero.py
    All checks for F16_aero passed.

Contributing
------------

Please feel free to share any thoughts or opinions about the design and
implementation of this software by `opening an issue on GitHub
<https://github.com/nasa/simupy-flight/issues/new>`_. Constructive feedback is
welcomed and appreciated.

Bug fix pull requests are always welcome. For feature additions, breaking
changes, etc. check if there is an open issue discussing the change and
reference it in the pull request. If there isn't one, it is recommended to open
one with your rationale for the change before spending significant time
preparing the pull request.

Ideally, new/changed functionality should come with tests and documentation. If
you are new to contributing, it is perfectly fine to open a work-in-progress
pull request and have it iteratively reviewed.

For pull requests to be accepted, all contributors must have a contributor's agreement on file with NASA. We will provide contributors with additional information during the review process.

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
