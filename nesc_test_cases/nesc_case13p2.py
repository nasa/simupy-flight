from simupy import systems
import os
from scipy import interpolate
from nesc_testcase_helper import plot_nesc_comparisons, plot_F16_controls, data_relative_path
from nesc_case11 import int_opts, get_controller_function, BD, opt_ctrl, dim_feedback, trimmed_KEAS


glob_path = os.path.join(data_relative_path, 'Atmospheric_checkcases', 'Atmos_13p2_SubsonicAirspeedChangeF16', 'Atmos_13p2_sim_*.csv')

keasCmdBlock = systems.SystemFromCallable(interpolate.make_interp_spline([0, 5], [trimmed_KEAS, trimmed_KEAS-5.0], k=0), 0, 1)
    
BD.systems[-4] = keasCmdBlock
BD.systems[2] = systems.SystemFromCallable(get_controller_function(*opt_ctrl, sasOn=True, apOn=True), dim_feedback + 4, 4)
res = BD.simulate(20, integrator_options=int_opts)

plot_nesc_comparisons(res, glob_path, '13p2')
plot_F16_controls(res, '13p2')

