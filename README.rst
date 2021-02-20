SimuPy Flight Vehicle Toolbox
=============================

Vehicle flight simulation is an important part of the innovation of aerospace vehicle technology. The NASA Engineering Safety Center (NESC) has identified and addressed the need to verify flight vehicle simulations through their work on the `“Six degree-of-freedom (6-DOF) Flight Simulation Test Cases.” <https://nescacademy.nasa.gov/flightsim/>`_ The author was prompted to develop tools that would allow for the rapid implementation of simulations for novel vehicle concepts including hypersonic re-entry vehicles and urban air mobility vehicles.

This software library leverages open source scientific computing tools to implement an efficient simulation framework for flight vehicles in Python. Equations of motion are composed in blocks using the SimuPy library, an open source Python alternative to Simulink, and integrated using SciPy’s wrappers for standard Fortran implementations of ordinary differential equation solvers. Dynamics equations of the inertial state variables for the position, orientation, and their corresponding rates for integration are developed using the SymPy symbolic library and implemented using code generation. Kinematics equations are implemented through symbolic definition and code generation as well as leveraging other open source software that implements useful functions, such as the solutions to the inverse geodesy problem.


NESC Test Cases
---------------

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