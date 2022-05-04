"""
===================================================
Case 13.3: Course change of a subsonic aircraft
===================================================

==============  ===============
Verifies        Multidimensional table look-up
Gravitation     J2
Geodesy         WGS-84 rotating
Atmosphere      US 1976 STD
Winds           still air
Vehicle         F-16 with simple auto-pilot
Notes           Initially straight & level. t=15s, command 15 degree right heading change
==============  ===============

For the manuevering examples, the BlockDiagram from case 11 is modified to replace the
controller with the auto-pilot configuration and generate the appropriate command
signals.
"""

from simupy import systems
import numpy as np
from scipy import interpolate

from nesc_testcase_helper import plot_nesc_comparisons, plot_F16_controls, benchmark
from nesc_case11 import (
    int_opts,
    F16ControllerBlock,
    BD,
    spec_ic_args,
    opt_ctrl,
    dim_feedback,
)

baseChiCmdBlock = systems.SystemFromCallable(
    interpolate.make_interp_spline(
        [0, 15],
        [spec_ic_args["psi"] * 180 / np.pi, spec_ic_args["psi"] * 180 / np.pi + 15.0],
        k=0,
    ),
    0,
    1,
)

BD.systems[-1] = baseChiCmdBlock
BD.systems[2] = F16ControllerBlock(*opt_ctrl, sasOn=True, apOn=True, event_t=15.)

with benchmark() as b:
    res = BD.simulate(30, integrator_options=int_opts)

# %%

plot_nesc_comparisons(res, "13p3")
plot_F16_controls(res, "13p3")
