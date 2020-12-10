from simupy import systems
import os
import numpy as np
from scipy import interpolate
from nesc_testcase_helper import plot_nesc_comparisons, plot_F16_controls, data_relative_path
from nesc_case11 import int_opts, get_controller_function, BD, spec_ic_args, opt_ctrl, dim_feedback


glob_path = os.path.join(data_relative_path, 'Atmospheric_checkcases', 'Atmos_13p3_SubsonicHeadingChangeF16', 'Atmos_13p3_sim_*.csv')

baseChiCmdBlock = systems.SystemFromCallable(interpolate.make_interp_spline([0, 15], [spec_ic_args['psi']*180/np.pi, spec_ic_args['psi']*180/np.pi+15.0], k=0), 0, 1)
    
BD.systems[-1] = baseChiCmdBlock
BD.systems[2] = systems.SystemFromCallable(get_controller_function(*opt_ctrl, sasOn=True, apOn=True), dim_feedback + 4, 4)
res = BD.simulate(30, integrator_options=int_opts)

plot_nesc_comparisons(res, glob_path, '13p3')
plot_F16_controls(res, '13p3')

