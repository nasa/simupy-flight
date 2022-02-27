"""
=========================================================
Case 13.4: Lateral offset manuever of a subsonic aircraft
=========================================================

==============  ===============
Verifies        Multidimensional table look-up
Gravitation     J2
Geodesy         WGS-84 rotating
Atmosphere      US 1976 STD
Winds           still air
Vehicle         F-16 with simple auto-pilot
Notes           Initially straight & level. t=20s, 2000-ft lateral course offset.
==============  ===============

For the manuevering examples, the BlockDiagram from case 11 is modified to replace the
controller with the auto-pilot configuration and generate the appropriate command
signals.
"""


from simupy import systems
import numpy as np
from scipy import interpolate

from nesc_testcase_helper import plot_nesc_comparisons, plot_F16_controls, benchmark
from nesc_case11 import int_opts, get_controller_function, BD, spec_ic_args, opt_ctrl, dim_feedback, planet, baseChiCmdBlock, latOffsetStateEquation, latOffsetOutputEquation


def latOffsetOutputEquationShift(t, x):
    latOffset = latOffsetOutputEquation(t, x)
    if t < 20:
        return latOffset
    else:
        return latOffset - 2000


latOffsetBlock = systems.DynamicalSystem(state_equation_function=latOffsetStateEquation, output_equation_function=latOffsetOutputEquationShift, dim_state=1, dim_input=3, dim_output=1)

BD.systems[-2] = latOffsetBlock
BD.systems[2] = systems.SystemFromCallable(get_controller_function(*opt_ctrl, sasOn=True, apOn=True), dim_feedback + 4, 4)
BD.connect(baseChiCmdBlock, latOffsetBlock, inputs=[0])
BD.connect(planet, latOffsetBlock, outputs=[planet.V_N_idx, planet.V_E_idx], inputs=[1,2])

with benchmark() as b:
    res = BD.simulate(60, integrator_options=int_opts)

plot_nesc_comparisons(res, '13p4')
plot_F16_controls(res, '13p4')
